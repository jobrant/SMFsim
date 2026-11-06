#!/usr/bin/env Rscript

# calibrate_bio_variation_2v2.R
#
# Which bio_region_sd makes SMFsim's null look like the REAL null? Runs the
# same 2-vs-2 splits of one group's replicates (default M1) on the real data
# and on parametric simulations of it at several bio_region_sd values, and
# counts false-positive DMRs for raw and SMFnorm (1.0 / 1.0) in each.
#
# Why a dedicated script: the main null (validate_null_calibration.R) is 4 vs 4
# -- every source replicate is drawn once per group -- so its false-positive
# counts cannot be compared with the real 2 vs 2 null (validate_real_null.R).
# Here each real library is simulated ONCE per bio_region_sd, from its own
# coverage and its own efficiency (mean rate relative to the highest), and the
# splits reuse those simulated libraries exactly as the real splits reuse the
# real ones. Real and simulated arms then differ only in where the counts
# come from.
#
# Read the result per split. For M1, N2+N3 vs N4+N5 is a sequencing-run split
# (data/M-Series_batches.tsv): it carries a locus-level run effect SMFsim does
# not model, so expect the real count to exceed every simulated one there.
# Calibrate on the other two splits.
#
# Usage, from the manuscript root (~15 min for chr1):
#   Rscript SMFsim/inst/scripts/calibrate_bio_variation_2v2.R
#
# Environment variables (all optional):
#   BIO2V2_GROUP   group in the sample sheet           (default M1)
#   BIO2V2_CHR     chromosome regex                    (default '^(chr)?1$';
#                  '' = all standard chromosomes, several hours)
#   BIO2V2_SD      comma-separated bio_region_sd grid  (default 0,0.1,0.25,0.5)
#   BIO2V2_OUT     output directory  (default results/bio_variation_2v2_<chr>)
#   SMFSIM_OVERWRITE  1 to write into a non-empty output directory
#
# Other settings (noise, dispersion, efficiency model, metilene) come from
# parse_args(), i.e. the same defaults as the main simulations.
#
# Outputs: bio_2v2_fp.csv (one row per source x split x arm) and
# replicate_efficiency.csv. metilene_input/ is deleted after each run.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
# Internal package functions (visible under load_all; fall back to the
# namespace when SMFsim is installed rather than loaded from source).
.internal <- function(f) if (exists(f)) get(f) else getFromNamespace(f, "SMFsim")
.n_tests_of   <- .internal(".metilene_n_tests")
.draw_rep     <- .internal(".draw_parametric_rep")
.noise_of     <- .internal(".noise_from_config")
.pooled       <- .internal(".pooled_rate")

# --- Settings --------------------------------------------------------------
group       <- Sys.getenv("BIO2V2_GROUP", "M1")
chr_pattern <- Sys.getenv("BIO2V2_CHR", "^(chr)?1$")
sd_grid     <- as.numeric(strsplit(Sys.getenv("BIO2V2_SD", "0,0.1,0.25,0.5"), ",")[[1]])
if (anyNA(sd_grid) || any(sd_grid < 0)) {
    stop("BIO2V2_SD must be comma-separated non-negative numbers.", call. = FALSE)
}
out <- Sys.getenv("BIO2V2_OUT", file.path(
    "results", paste0("bio_variation_2v2",
                      if (nzchar(chr_pattern)) paste0("_", gsub("[^A-Za-z0-9]", "", chr_pattern)))))

cfg <- parse_args(character(0))
cfg$wt_group_id <- group
if (nzchar(chr_pattern)) cfg$chr_pattern <- chr_pattern

message(sprintf("bio_region_sd 2v2 calibration: group %s, %s, sd grid %s -> %s",
                group, if (nzchar(chr_pattern)) chr_pattern else "all standard chromosomes",
                paste(sd_grid, collapse = "/"), out))
message(sprintf("  dispersion_s %s, fragment_rho %s / %s bp, bio_region_bp %s, efficiency %s",
                cfg$dispersion_s, cfg$fragment_rho, cfg$fragment_bp,
                cfg$bio_region_bp, cfg$efficiency_model))

if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 &&
    !identical(Sys.getenv("SMFSIM_OVERWRITE"), "1")) {
    stop("Output directory is not empty: ", normalizePath(out, mustWork = FALSE),
         "\nArchive it, set BIO2V2_OUT, or re-run with SMFSIM_OVERWRITE=1.",
         call. = FALSE)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# --- Data -------------------------------------------------------------------
reps <- prepare_wt_replicates(cfg)
k <- length(reps)
if (k < 4) {
    stop(sprintf("Group %s has %d replicates; a 2 vs 2 null needs at least 4.", group, k),
         call. = FALSE)
}
nm <- names(reps)

# Efficiency of each library relative to the best one, from its mean rate. The
# true rate is the pooled rate divided back by the mean efficiency, so the
# simulated libraries reproduce the real mean rates.
rate_mean <- vapply(reps, function(dt) mean(dt$rate), numeric(1))
eff <- rate_mean / max(rate_mean)
p_true <- pmin(1, .pooled(reps) / mean(eff))
fwrite(data.table(replicate = nm, mean_rate = rate_mean, efficiency = eff),
       file.path(out, "replicate_efficiency.csv"))
for (i in seq_len(k)) message(sprintf("  %-14s mean rate %.4f  efficiency %.3f",
                                      nm[i], rate_mean[i], eff[i]))

half <- k %/% 2
splits <- combn(k, half, simplify = FALSE)
if (k %% 2 == 0) splits <- Filter(function(a) 1L %in% a, splits)

arms <- list(
    raw     = function(pg) run_raw(pg)$data,
    SMFnorm = function(pg) run_SMFnorm(pg, within_alpha = 1.0, between_alpha = 1.0,
                                       min_coverage = cfg$min_coverage %||% 10,
                                       rate_between_groups = TRUE)$data
)

# Sources: the real libraries, then one simulated set per bio_region_sd.
sources <- c(list(real = reps), setNames(vector("list", length(sd_grid)),
                                         sprintf("sim_sd%s", sd_grid)))

res <- list()
for (si in seq_along(sources)) {
    src <- names(sources)[si]
    libs <- sources[[si]]
    if (is.null(libs)) {
        bsd <- sd_grid[si - 1]
        noise <- .noise_of(modifyList(cfg, list(bio_region_sd = bsd)))
        message("\n", strrep("=", 64), "\nSIMULATING ", src, "\n", strrep("=", 64))
        libs <- lapply(seq_len(k), function(j) {
            .draw_rep(reps[[j]], p_true, eff[j], seed_i = 7000L + 100L * si + j,
                      dispersion_s = cfg$dispersion_s, noise = noise,
                      efficiency_model = cfg$efficiency_model %||% "multiplicative")
        })
        names(libs) <- nm
    }

    for (s in seq_along(splits)) {
        ia <- splits[[s]]
        ib <- setdiff(seq_len(k), ia)
        label <- paste(paste(sub("^allc_", "", nm[ia]), collapse = "+"), "vs",
                       paste(sub("^allc_", "", nm[ib]), collapse = "+"))
        pg <- list(PseudoA = setNames(lapply(libs[ia], copy), paste0("PseudoA_", nm[ia])),
                   PseudoB = setNames(lapply(libs[ib], copy), paste0("PseudoB_", nm[ib])),
                   params = list())

        for (a in names(arms)) {
            dir <- file.path(out, src, paste0("split", s), a)
            dmrs <- suppressMessages(call_dmrs_metilene(
                arms[[a]](pg), dir, cfg$metilene_path,
                min_cpg = cfg$metilene_min_cpg, min_diff = cfg$metilene_min_diff,
                metilene_max_dist = cfg$metilene_max_dist,
                q_cutoff = cfg$metilene_qval, q_source = cfg$metilene_q_source,
                p_column = cfg$metilene_p_column, min_effect = cfg$metilene_min_effect,
                mtc = cfg$metilene_mtc))
            # An empty table from a crashed metilene would read as a perfect null.
            n_tests <- .n_tests_of(file.path(dir, "metilene_log.txt"))
            if (is.na(n_tests)) {
                stop("metilene did not complete for ", src, " / split ", s, " / ", a,
                     " (no 'Number of Tests' line). See ", dir, call. = FALSE)
            }
            unlink(file.path(dir, "metilene_input"), recursive = TRUE)

            res[[length(res) + 1]] <- data.table(
                source = src, bio_region_sd = if (src == "real") NA_real_ else sd_grid[si - 1],
                split = s, split_label = label, arm = a, FP = nrow(dmrs), n_tests = n_tests)
            message(sprintf("  %-10s split %d (%s) %-8s %6d FP DMRs", src, s, label, a,
                            nrow(dmrs)))
        }
        rm(pg); invisible(gc())
    }
    fwrite(rbindlist(res), file.path(out, "bio_2v2_fp.csv"))   # checkpoint per source
    if (src != "real") rm(libs)
}

# --- Summary ----------------------------------------------------------------
res <- rbindlist(res)
wide <- dcast(res, split + split_label + arm ~ source, value.var = "FP")
setcolorder(wide, c("split", "split_label", "arm", names(sources)))
message("\n", strrep("=", 64),
        "\nFALSE-POSITIVE DMRs, real vs simulated at each bio_region_sd (2 vs 2)\n",
        strrep("=", 64))
print(wide)
message("\nCalibrate on splits that are not sequencing-run splits; see the header.")
message("Results: ", normalizePath(file.path(out, "bio_2v2_fp.csv")))
