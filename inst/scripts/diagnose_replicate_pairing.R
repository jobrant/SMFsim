#!/usr/bin/env Rscript

# diagnose_replicate_pairing.R
#
# Do a group's replicates share LOCUS-LEVEL deviations in subsets -- a batch
# effect beyond global labeling efficiency, which no rate-bin correction
# (SMFnorm included) can remove?
#
#   Rscript SMFsim/inst/scripts/diagnose_replicate_pairing.R [group] [chr_pattern]
#   # defaults: M1, chromosome 1. Run from the manuscript root. ~2 minutes.
#
# Each replicate's per-rate-bin offset is removed (SMFnorm's within-group step,
# alpha 1, all replicates in one set). Rates are averaged in 1 kb windows with
# >= 10 GpC sites, each replicate's deviation from the window mean is taken,
# and deviations are correlated between replicates. Deviations sum to zero, so
# the no-structure baseline is -1/(k-1) (-0.33 for 4 replicates); pairs above
# it share locus-level effects.
#
# Finding (2026-10-27, M1, chr1, 17,103 windows), after offset alignment:
#   N2-N3 -0.06   N4-N5 +0.17   cross pairs -0.40 to -0.63
# M1 behaves as two batches, {N2, N3} and {N4, N5}. In the real-data null
# (validate_real_null.R) the N2+N3 vs N4+N5 split is therefore a batch split:
# raw calls 2,460 DMRs there against 7 in the balanced split, and SMFnorm,
# which corrects efficiency only, calls 3,305. Finer rate bins do not help
# (diagnose_bin_resolution.R).

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)

args <- commandArgs(TRUE)
group <- if (length(args) >= 1) args[1] else "M1"
chr_pattern <- if (length(args) >= 2) args[2] else "^(chr)?1$"

cfg <- parse_args(character(0))
cfg$wt_group_id <- group
cfg$chr_pattern <- chr_pattern
reps <- suppressMessages(prepare_wt_replicates(cfg))
k <- length(reps)
if (k < 3) stop("Need at least 3 replicates to look for shared structure.", call. = FALSE)
nm <- sub("^allc_", "", names(reps))

# Within-group offset alignment over all replicates at once. run_SMFnorm()
# needs two groups, so group B is a throwaway copy of two replicates.
pg <- list(PseudoA = setNames(reps, paste0("PseudoA_", names(reps))),
           PseudoB = setNames(lapply(reps[1:2], copy), c("PseudoB_x1", "PseudoB_x2")),
           params = list())
aligned <- suppressMessages(run_SMFnorm(
    pg, within_alpha = 1, between_alpha = 1, min_coverage = 0,
    rate_between_groups = FALSE))$data$PseudoA

window_dev_cor <- function(R) {
    win <- paste(reps[[1]]$chr, reps[[1]]$pos %/% 1000)
    sums <- rowsum(R, win)
    n_in <- as.vector(rowsum(rep(1, nrow(R)), win))
    W <- (sums / n_in)[n_in >= 10, , drop = FALSE]
    C <- cor(W - rowMeans(W))
    dimnames(C) <- list(nm, nm)
    list(C = C, n = nrow(W))
}

raw <- do.call(cbind, lapply(reps, `[[`, "rate"))
raw <- sweep(raw, 2, colMeans(raw), "/") * mean(raw)     # remove global level only
out <- list(`raw (global mean removed)` = window_dev_cor(raw),
            `after per-bin offset alignment` =
                window_dev_cor(do.call(cbind, lapply(aligned, `[[`, "rate"))))

for (lbl in names(out)) {
    cat(sprintf("\n%s, %d 1-kb windows (no-structure baseline %.2f):\n",
                lbl, out[[lbl]]$n, -1 / (k - 1)))
    print(round(out[[lbl]]$C, 3))
}
