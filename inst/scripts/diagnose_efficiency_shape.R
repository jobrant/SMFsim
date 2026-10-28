#!/usr/bin/env Rscript

# diagnose_efficiency_shape.R
#
# What shape do real labeling-efficiency differences have across rate levels?
# For every pair of replicates in a group, sites are binned into deciles of the
# group's mean rate, and the pair's difference in each decile is expressed as
#   ratio            rate_hi / rate_lo           constant if MULTIPLICATIVE
#   difference       rate_hi - rate_lo           constant if ADDITIVE
#   logit_difference qlogis(hi) - qlogis(lo)     constant if a LOGIT SHIFT
# plus the logit shift SMFsim's matched-mean rule assigns to the pair's
# efficiency ratio (see .apply_efficiency()), for comparison.
#
#   Rscript SMFsim/inst/scripts/diagnose_efficiency_shape.R [group] [chr_pattern]
#   # defaults: M1, chromosome 1. Run from the manuscript root. ~2 minutes.
#
# Result (2026-10-28, M1-M4, chr1), coefficient of variation across deciles
# 4-10 for pairs with a real efficiency difference (ratio vs logit difference):
#   M1 N3/N2 0.080 vs 0.021 (logit better)   M1 N3/N4 0.069 vs 0.090
#   M1 N3/N5 0.074 vs 0.114                  M1 N5/N2 0.018 vs 0.185
#   M1 N5/N4 0.024 vs 0.296                  M2 0.068 vs 0.108
#   M4 0.041 vs 0.166
# A constant ratio (multiplicative) fits 6 of 7; additive fits worst in all.
# Matched-mean logit shift vs observed: 0.72/0.69, 0.41/0.41, 0.34/0.33.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
.shift <- if (exists(".efficiency_shift")) .efficiency_shift else
    getFromNamespace(".efficiency_shift", "SMFsim")

args <- commandArgs(TRUE)
group <- if (length(args) >= 1) args[1] else "M1"
chr_pattern <- if (length(args) >= 2) args[2] else "^(chr)?1$"

cfg <- parse_args(character(0))
cfg$wt_group_id <- group
cfg$chr_pattern <- chr_pattern
reps <- suppressMessages(prepare_wt_replicates(cfg))
nm <- sub("^allc_", "", names(reps))
R <- do.call(cbind, lapply(reps, `[[`, "rate"))
colnames(R) <- nm
m <- rowMeans(R)
decile <- cut(m, unique(quantile(m, seq(0, 1, 0.1))), include.lowest = TRUE,
              labels = FALSE)

cat(sprintf("%s, %s: %d replicates, %d shared sites (%s)\n\n", group, chr_pattern,
            ncol(R), nrow(R), cfg$data_dir))

pairs <- combn(ncol(R), 2, simplify = FALSE)
summary_rows <- list()
for (pr in pairs) {
    means <- colMeans(R[, pr])
    hi <- pr[which.max(means)]
    lo <- pr[which.min(means)]
    tab <- data.table(decile = decile, hi = R[, hi], lo = R[, lo])[
        , .(pooled = mean(m[.I]), hi = mean(hi), lo = mean(lo)), keyby = decile]
    tab[, `:=`(ratio = hi / lo, difference = hi - lo,
               logit_difference = qlogis(hi) - qlogis(lo))]
    eff <- mean(R[, lo]) / mean(R[, hi])
    rule <- -.shift(R[, hi], eff)
    cat(sprintf("== %s (higher) vs %s: efficiency ratio %.3f, matched-mean logit shift %.2f\n",
                nm[hi], nm[lo], eff, rule))
    print(tab[, lapply(.SD, function(x) round(x, 3))], row.names = FALSE)
    upper <- tab[decile >= 4]
    summary_rows[[length(summary_rows) + 1]] <- data.table(
        pair = paste(nm[hi], "vs", nm[lo]), efficiency_ratio = round(eff, 3),
        cv_ratio = round(sd(upper$ratio) / mean(upper$ratio), 3),
        cv_difference = round(sd(upper$difference) / mean(upper$difference), 3),
        cv_logit_difference = round(sd(upper$logit_difference) /
                                    mean(upper$logit_difference), 3),
        logit_shift_observed = round(mean(upper$logit_difference), 2),
        logit_shift_rule = round(rule, 2))
    cat("\n")
}

cat("Which scale is the difference constant on? (coefficient of variation\n",
    "across deciles 4-10; smaller = more constant)\n", sep = "")
print(rbindlist(summary_rows), row.names = FALSE)
