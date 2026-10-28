#!/usr/bin/env Rscript

# validate_real_null.R
#
# Specificity on REAL data: split one condition's replicates into two halves
# (e.g. M1 2 vs 2), so the two "groups" are the same biology and every DMR
# called is a false positive. Runs every distinct split, for raw, downsampled
# and several SMFnorm settings.
#
# This is the primary specificity evidence. SMFsim's simulated null cannot
# stand in for it: simulated replicates share one true rate, sites are drawn
# independently, and the efficiency artifact is a pure multiplier that
# SMFnorm's per-bin offset alignment inverts by construction.
#
# Usage, from the manuscript root (batch: inst/slurm/real_null.sbatch):
#   Rscript SMFsim/inst/scripts/validate_real_null.R
#
# Environment variables (all optional):
#   REAL_NULL_GROUP   group in the sample sheet to split      (default M1)
#   REAL_NULL_CHR     chromosome regex, e.g. '^(chr)?1$'      (default: all
#                     standard chromosomes; chr1 takes ~10 min, the genome
#                     several hours)
#   REAL_NULL_DATA    allc directory, e.g. data/allc_methylseq4.2_clip5
#                     (default: parse_args() data_dir, i.e. data/allc)
#   REAL_NULL_OUT     output directory  (default results/real_null_<group>[_chr])
#   SMFSIM_OVERWRITE  1 to write into a non-empty output directory
#
# Outputs under the output directory:
#   replicate_mean_rates.csv   per-replicate mean GpC rate (efficiency proxy)
#   real_null_fp.csv           one row per split x arm: FP DMRs, metilene test
#                              count, each half's mean rate and the gap
#   split<i>/<arm>/            metilene input, output and log for each run
#
# Genome-wide result (2026-10-27; 12 min, 12.6 GB), FP DMRs; replicate mean
# rates N2 .140, N3 .231, N4 .143, N5 .169:
#                     raw  downsampled  w1.0+btw1.0  w0.5+btw0.95  w1.0 btw-off
#   N2+N4 vs N3+N5  11654     2432            6          366          18599
#   N2+N5 vs N3+N4      7        0            5          482           1531
#   N2+N3 vs N4+N5   2460      757         3305         6386           7859
# In the last split SMFnorm's calls overlap raw's only 44% and are 78% "B
# higher" (raw: 59% "A higher"): new, directional false positives, not the
# same difference detected more sensitively. Cause: M1 behaves as two batches,
# {N2, N3} and {N4, N5}, sharing locus-level deviations that a rate-bin
# correction cannot remove, so that split is a batch split rather than an
# efficiency split (diagnose_replicate_pairing.R). Finer bins do not help
# (diagnose_bin_resolution.R).

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
.n_tests_of <- if (exists(".metilene_n_tests")) .metilene_n_tests else
    getFromNamespace(".metilene_n_tests", "SMFsim")

# --- Settings --------------------------------------------------------------
group       <- Sys.getenv("REAL_NULL_GROUP", "M1")
chr_pattern <- Sys.getenv("REAL_NULL_CHR", "")
out <- Sys.getenv("REAL_NULL_OUT", file.path(
    "results", paste0("real_null_", group,
                      if (nzchar(chr_pattern)) paste0("_", gsub("[^A-Za-z0-9]", "", chr_pattern)))))

cfg <- parse_args(character(0))      # metilene, coverage and q-value settings
cfg$wt_group_id <- group
if (nzchar(chr_pattern)) cfg$chr_pattern <- chr_pattern
if (nzchar(Sys.getenv("REAL_NULL_DATA"))) cfg$data_dir <- Sys.getenv("REAL_NULL_DATA")
message("allc directory: ", cfg$data_dir)

# Arms. "btw" = rate_between_groups with that between_alpha; "off" = within-group
# normalization only. alpha only scales site deviations from each bin's mean;
# the offset removal that corrects efficiency happens at any alpha.
arms <- list(
    raw                      = list(type = "raw"),
    downsampled              = list(type = "downsampled"),
    `SMFnorm w1.0 btw1.0`    = list(type = "SMFnorm", wa = 1.0, ba = 1.0,  rbg = TRUE),
    `SMFnorm w0.5 btw1.0`    = list(type = "SMFnorm", wa = 0.5, ba = 1.0,  rbg = TRUE),
    `SMFnorm w0.5 btw0.95`   = list(type = "SMFnorm", wa = 0.5, ba = 0.95, rbg = TRUE),
    `SMFnorm w1.0 btw-off`   = list(type = "SMFnorm", wa = 1.0, ba = 1.0,  rbg = FALSE)
)

message(sprintf("Real-data null: group %s, %s -> %s", group,
                if (nzchar(chr_pattern)) chr_pattern else "all standard chromosomes", out))

# --- Overwrite guard --------------------------------------------------------
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 &&
    !identical(Sys.getenv("SMFSIM_OVERWRITE"), "1")) {
    stop("Output directory is not empty: ", normalizePath(out, mustWork = FALSE),
         "\nArchive it, set REAL_NULL_OUT, or re-run with SMFSIM_OVERWRITE=1.",
         call. = FALSE)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# --- Data -------------------------------------------------------------------
reps <- prepare_wt_replicates(cfg)
k <- length(reps)
if (k < 4) {
    stop(sprintf("Group %s has %d replicates; a split null needs at least 4 (2 vs 2).",
                 group, k), call. = FALSE)
}
nm <- names(reps)

rate_mean <- vapply(reps, function(dt) mean(dt$rate), numeric(1))
rates_dt <- data.table(replicate = nm, mean_rate = rate_mean, n_sites = nrow(reps[[1]]))
fwrite(rates_dt, file.path(out, "replicate_mean_rates.csv"))
message("Per-replicate mean rate (efficiency proxy):")
for (i in seq_len(k)) message(sprintf("  %-14s %.4f", nm[i], rate_mean[i]))

# Every distinct split into two halves. With equal halves (even k), keeping
# only those that contain replicate 1 drops each split's mirror image; with odd
# k the halves differ in size, so every combination is distinct.
half <- k %/% 2
splits <- combn(k, half, simplify = FALSE)
if (k %% 2 == 0) splits <- Filter(function(a) 1L %in% a, splits)
message(sprintf("%d replicates -> %d distinct %d-vs-%d splits", k, length(splits),
                half, k - half))

arm_data <- function(arm, pg) {
    switch(arm$type,
           raw         = run_raw(pg)$data,
           downsampled = run_downsampled(pg)$data,
           SMFnorm     = run_SMFnorm(pg, within_alpha = arm$wa, between_alpha = arm$ba,
                                     min_coverage = cfg$min_coverage %||% 10,
                                     rate_between_groups = arm$rbg)$data)
}

# --- Run --------------------------------------------------------------------
res <- list()
for (s in seq_along(splits)) {
    ia <- splits[[s]]
    ib <- setdiff(seq_len(k), ia)
    label <- paste(paste(sub("^allc_", "", nm[ia]), collapse = "+"), "vs",
                   paste(sub("^allc_", "", nm[ib]), collapse = "+"))
    message("\n", strrep("=", 64), "\nSPLIT ", s, ": ", label, "\n", strrep("=", 64))

    pg <- list(PseudoA = setNames(lapply(reps[ia], copy), paste0("PseudoA_", nm[ia])),
               PseudoB = setNames(lapply(reps[ib], copy), paste0("PseudoB_", nm[ib])),
               params = list())

    for (a in names(arms)) {
        message("\n--- ", a, " ---")
        t0 <- Sys.time()
        dat <- arm_data(arms[[a]], pg)
        dir <- file.path(out, paste0("split", s), make.names(a))
        dmrs <- call_dmrs_metilene(
            dat, dir, cfg$metilene_path,
            min_cpg = cfg$metilene_min_cpg, min_diff = cfg$metilene_min_diff,
            metilene_max_dist = cfg$metilene_max_dist,
            q_cutoff = cfg$metilene_qval, q_source = cfg$metilene_q_source,
            p_column = cfg$metilene_p_column, min_effect = cfg$metilene_min_effect,
            mtc = cfg$metilene_mtc)

        # call_dmrs_metilene() returns an empty table when metilene crashes,
        # which would read as a perfect null. Refuse to record it.
        n_tests <- .n_tests_of(file.path(dir, "metilene_log.txt"))
        if (is.na(n_tests)) {
            stop("metilene did not complete for split ", s, " / ", a,
                 " (no 'Number of Tests' line). See ", dir, call. = FALSE)
        }

        res[[length(res) + 1]] <- data.table(
            split = s, split_label = label, arm = a, FP = nrow(dmrs),
            n_tests = n_tests,
            mean_rate_A = mean(rate_mean[ia]), mean_rate_B = mean(rate_mean[ib]),
            rate_gap = abs(mean(rate_mean[ia]) - mean(rate_mean[ib])),
            minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
        message(sprintf("  %s: %d false-positive DMRs (%d tests)", a, nrow(dmrs), n_tests))
        rm(dat); invisible(gc())
    }
    fwrite(rbindlist(res), file.path(out, "real_null_fp.csv"))   # checkpoint per split
}

# --- Summary ----------------------------------------------------------------
res <- rbindlist(res)
fwrite(res, file.path(out, "real_null_fp.csv"))

wide <- dcast(res, split + split_label + rate_gap ~ arm, value.var = "FP")
setcolorder(wide, c("split", "split_label", "rate_gap", intersect(names(arms), names(wide))))
setorder(wide, -rate_gap)
message("\n", strrep("=", 64), "\nFALSE-POSITIVE DMRs (every call is false)\n",
        "rate_gap = difference in the two halves' mean rate (efficiency imbalance)\n",
        strrep("=", 64))
print(wide)
message("\nResults: ", normalizePath(file.path(out, "real_null_fp.csv")))
