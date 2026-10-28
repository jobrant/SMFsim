#!/usr/bin/env Rscript

# diagnose_noise_correlation.R
#
# Is replicate-to-replicate noise in real data independent from site to site,
# as SMFsim's simulator assumes, or spatially correlated?
#
#   Rscript SMFsim/inst/scripts/diagnose_noise_correlation.R [group] [chr_pattern]
#   # defaults: M1, chromosome 1. Run from the manuscript root. ~5 minutes.
#
# For each replicate it standardises the residual from the pooled site rate by
# its binomial SD, removes the replicate's global offset (efficiency), and
# correlates residuals of neighbouring sites in distance bins. A random-pair
# control should be ~0.
#
# Finding (2026-10-27, M1, chr1):
#   1-50 bp apart:     0.327   <- shared sequenced fragments; SMFsim has none
#   51-150 bp:         0.006
#   151-300 bp:        0.031
#   1-5 kb:            0.003   <- no region-level biological replicate variation
#   random pairs:      0.001
# M1 replicates are cell-line replicates; tissue or patient samples would be
# expected to show region-level variation that M1 does not.

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
if (length(reps) < 2) stop("Need at least 2 replicates in group ", group, call. = FALSE)

R <- do.call(cbind, lapply(reps, `[[`, "rate"))
C <- do.call(cbind, lapply(reps, `[[`, "cov"))
pos <- reps[[1]]$pos
n <- nrow(R)
k_rep <- ncol(R)

m <- rowSums(R * C) / rowSums(C)                    # pooled site rate
Z <- (R - m) / sqrt(pmax(m * (1 - m), 1e-4) / C)    # binomial-standardised residuals
Z <- sweep(Z, 2, colMeans(Z))                       # drop per-replicate global offset

cat(sprintf("%s, %s: %d sites shared by all %d replicates\n", group, chr_pattern, n, k_rep))
cat(sprintf("residual variance / binomial expectation: %.2f  (about (k-1)/k = %.2f for pure binomial noise)\n",
            mean(Z^2), (k_rep - 1) / k_rep))

d <- diff(pos)
for (rng in list(c(1, 50), c(51, 150), c(151, 300), c(1000, 5000), c(20000, 1e6))) {
    k <- which(d >= rng[1] & d <= rng[2])
    r <- mean(sapply(seq_len(k_rep), function(j) cor(Z[k, j], Z[k + 1, j])))
    cat(sprintf("neighbour-site residual correlation, %6.0f-%7.0f bp apart (n=%7d): %.3f\n",
                rng[1], rng[2], length(k), r))
}

k <- which(d <= 150)
set.seed(1)
sh <- sample(n)
cat(sprintf("control, random site pairs: %.3f\n",
            mean(sapply(seq_len(k_rep), function(j) cor(Z[k, j], Z[sh[k], j])))))
