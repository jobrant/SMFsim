#!/usr/bin/env Rscript

# diagnose_bin_resolution.R
#
# Does SMFnorm's rate-bin resolution explain the false positives in the
# N2+N3 vs N4+N5 split of the M1 real-data null (validate_real_null.R)?
#
# Hypothesis: SMFnorm corrects each sample with ONE additive shift per rate bin,
# but labeling efficiency acts multiplicatively. Bins hold equal numbers of
# sites and most GpC rates are low, so the high-rate bins span a wide range of
# rates. Inside a wide bin a single shift over-corrects sites at one end and
# under-corrects the other, worst for the most extreme replicate (N3). If so,
# finer bins should remove those false positives without harming the splits
# that already sit near zero.
#
# Result (2026-10-27, chr1): hypothesis NOT supported. N2+N3 split FP DMRs at
# 50 / 200 / 1000 bins: 379 / 342 / 400 (other splits 0-4 throughout). Only 23
# of 379 sit above rate 0.5, where bins are widest; they spread over rates
# 0.1-0.5 and stay ~78% "N4+N5 higher" at every setting. The cause is a
# locus-level batch effect pairing {N2, N3} vs {N4, N5}; see
# diagnose_replicate_pairing.R.
#
# Usage, from the manuscript root (genome-wide takes a few minutes per setting):
#   Rscript SMFsim/inst/scripts/diagnose_bin_resolution.R
#
# Environment variables (optional):
#   BINRES_GROUP    group to split                      (default M1)
#   BINRES_CHR      chromosome regex                    (default: all standard)
#   BINRES_BINS     comma-separated max_quantiles values (default 50,200,1000)
#   BINRES_OUT      output directory (default results/diagnostics/bin_resolution_<group>)
#
# Outputs under the output directory:
#   bin_resolution_fp.csv   split x max_quantiles: FP DMRs, bins realised, widest bin
#   fp_by_bin_width.csv     where the false positives sit, by pooled-rate level
#   split<i>/mq<n>/         metilene input, output and log

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
.ns <- function(f) if (exists(f)) get(f) else getFromNamespace(f, "SMFsim")
.n_tests_of   <- .ns(".metilene_n_tests")
.filter_cov   <- .ns(".filter_pseudo_groups_by_min_coverage")
.format_smf   <- .ns("format_for_smfnorm")

group <- Sys.getenv("BINRES_GROUP", "M1")
chr_pattern <- Sys.getenv("BINRES_CHR", "")
bins <- as.integer(strsplit(Sys.getenv("BINRES_BINS", "50,200,1000"), ",")[[1]])
out <- Sys.getenv("BINRES_OUT", file.path("results", "diagnostics",
                                          paste0("bin_resolution_", group)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cfg <- parse_args(character(0))
cfg$wt_group_id <- group
if (nzchar(chr_pattern)) cfg$chr_pattern <- chr_pattern

reps <- prepare_wt_replicates(cfg)
k <- length(reps)
if (k < 4) stop("Need at least 4 replicates for a 2-vs-2 split.", call. = FALSE)
nm <- names(reps)
splits <- Filter(function(a) 1L %in% a, combn(k, k %/% 2, simplify = FALSE))

# Pooled raw average rate per site (the same for every split), as intervals so
# false-positive DMRs can be mapped back to the sites they cover.
avg <- rowMeans(do.call(cbind, lapply(reps, `[[`, "rate")))
site_iv <- data.table(chr = as.character(reps[[1]]$chr), start = reps[[1]]$pos,
                      end = reps[[1]]$pos, avg = avg)
setkey(site_iv, chr, start, end)

# SMFnorm with within = between = 1.0 (no shrinkage) at a given bin cap,
# otherwise exactly as run_SMFnorm() calls it.
smfnorm_mq <- function(pg, mq) {
    pg <- .filter_cov(pg, cfg$min_coverage)
    split_data <- SMFnorm::split_by_groups(.format_smf(pg))
    SMFnorm::normalize_methylation_data(
        data_list = split_data, do_coverage_norm = TRUE, normalize_rates = TRUE,
        coverage_between_groups = FALSE, rate_within_groups = TRUE,
        rate_between_groups = TRUE, within_alpha = 1.0, between_alpha = 1.0,
        min_coverage = cfg$min_coverage, max_quantiles = mq)
}

# Bins the between-group step actually uses: equal-count quantiles of the
# pooled average rate, de-duplicated (ties at low rates collapse edges).
realised_breaks <- function(avg, mq) {
    nq <- as.integer(max(5, min(length(avg) / 2000, mq)))
    unique(quantile(avg, probs = seq(0, 1, length.out = nq + 1), names = FALSE))
}

res <- list()
where <- list()
for (s in seq_along(splits)) {
    ia <- splits[[s]]; ib <- setdiff(seq_len(k), ia)
    label <- paste(paste(sub("^allc_", "", nm[ia]), collapse = "+"), "vs",
                   paste(sub("^allc_", "", nm[ib]), collapse = "+"))
    pg <- list(PseudoA = setNames(lapply(reps[ia], copy), paste0("PseudoA_", nm[ia])),
               PseudoB = setNames(lapply(reps[ib], copy), paste0("PseudoB_", nm[ib])),
               params = list())

    for (mq in bins) {
        message(sprintf("\n=== split %d (%s), max_quantiles = %d ===", s, label, mq))
        br <- realised_breaks(avg, mq)
        dat <- suppressMessages(smfnorm_mq(pg, mq))
        dir <- file.path(out, paste0("split", s), paste0("mq", mq))
        dmrs <- call_dmrs_metilene(
            dat, dir, cfg$metilene_path,
            min_cpg = cfg$metilene_min_cpg, min_diff = cfg$metilene_min_diff,
            metilene_max_dist = cfg$metilene_max_dist,
            q_cutoff = cfg$metilene_qval, q_source = cfg$metilene_q_source,
            p_column = cfg$metilene_p_column, min_effect = cfg$metilene_min_effect,
            mtc = cfg$metilene_mtc)
        n_tests <- .n_tests_of(file.path(dir, "metilene_log.txt"))
        if (is.na(n_tests)) {
            stop("metilene did not complete: ", dir, call. = FALSE)
        }

        widths <- diff(br)
        res[[length(res) + 1]] <- data.table(
            split = s, split_label = label, max_quantiles = mq,
            bins_realised = length(br) - 1L, widest_bin = max(widths),
            widest_bin_range = sprintf("%.3f-%.3f", br[which.max(widths)],
                                       br[which.max(widths) + 1]),
            FP = nrow(dmrs), n_tests = n_tests)
        message(sprintf("  %d bins realised (widest %.3f), %d false-positive DMRs",
                        length(br) - 1L, max(widths), nrow(dmrs)))

        # Where do the false positives sit? Pooled raw rate over each call's
        # sites (close to the scale SMFnorm bins on), mapped to its bin width.
        if (nrow(dmrs)) {
            d <- copy(dmrs)[, `:=`(chr = as.character(chr), dmr = .I)]
            setkey(d, chr, start, end)
            ov <- foverlaps(site_iv, d[, .(chr, start, end, dmr)], nomatch = 0L)
            lvl <- ov[, .(pooled_rate = mean(avg)), keyby = dmr]
            d <- lvl[d, on = "dmr"]
            b <- findInterval(d$pooled_rate, br, rightmost.closed = TRUE, all.inside = TRUE)
            where[[length(where) + 1]] <- data.table(
                split = s, max_quantiles = mq, mean_diff = d$mean_diff,
                pooled_rate = d$pooled_rate, bin_width = widths[b])
        }
        rm(dat); invisible(gc())
    }
    fwrite(rbindlist(res), file.path(out, "bin_resolution_fp.csv"))   # checkpoint
}

res <- rbindlist(res)
fwrite(res, file.path(out, "bin_resolution_fp.csv"))
if (length(where)) fwrite(rbindlist(where), file.path(out, "fp_by_bin_width.csv"))

message("\n", strrep("=", 64), "\nFALSE-POSITIVE DMRs by bin resolution ",
        "(SMFnorm, within = between = 1.0)\n", strrep("=", 64))
print(dcast(res, split + split_label ~ max_quantiles, value.var = "FP"))
message("\nBins actually realised (requested -> realised, widest bin):")
print(unique(res[, .(max_quantiles, bins_realised, widest_bin, widest_bin_range)]))
