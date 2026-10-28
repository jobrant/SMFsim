#!/usr/bin/env Rscript

# diagnose_alpha_semantics.R
#
# What do SMFnorm's within_alpha and between_alpha actually do? A controlled
# synthetic test (no real data, runs in about a minute):
#
#   Rscript SMFsim/inst/scripts/diagnose_alpha_semantics.R      # from the manuscript root
#
# Group B is simulated ~0.30 less efficient than group A (an aligned bias, like
# aligned_strong), and 1% of sites carry a true +0.10 difference in group B.
# For each SMFnorm setting it reports:
#   null A-B gap    mean group difference at non-spiked sites (the artifact)
#   within-grp SD   replicate SD at non-spiked sites (pure binomial noise here)
#   spike diff      mean B-A difference at spiked sites (the real signal)
#
# Finding (2026-10-27): in each rate bin SMFnorm sets every sample's bin mean to
# the common bin mean -- that offset removal IS the efficiency correction and
# happens fully at ANY alpha -- then multiplies each site's deviation from the
# sample's bin mean by alpha. So:
#   * between_alpha = 1.0 is NOT the identity / "no between-group correction":
#     it removes the whole group gap (0.069 -> 0.001) and recovers the spike.
#   * within_alpha < 1 corrects nothing; it only shrinks the noise (SD halves at
#     0.5), which makes downstream tests overconfident.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)

set.seed(7)
n <- 60000
cov <- 30L
p <- rbeta(n, 1.2, 4)                          # true accessibility
spike <- seq_len(n) <= 600                      # 1% of sites: +0.10 in group B
pB <- pmin(p + 0.10 * spike, 1)

mk <- function(prob, effs, prefix) {
    setNames(lapply(seq_along(effs), function(i) {
        mc <- rbinom(n, cov, pmin(prob * effs[i], 1))
        data.table(chr = "chr1", pos = seq_len(n), strand = "+", site = seq_len(n),
                   mc = mc, cov = cov, rate = mc / cov)
    }), paste0(prefix, "_", seq_along(effs)))
}
pg <- list(PseudoA = mk(p,  c(0.92, 0.88, 0.90, 0.86), "PseudoA"),
           PseudoB = mk(pB, c(0.62, 0.58, 0.60, 0.56), "PseudoB"),
           params = list())

summ <- function(d, label) {
    R <- function(g) do.call(cbind, lapply(d[[g]], `[[`, "rate"))
    A <- R("PseudoA"); B <- R("PseudoB"); nul <- !spike
    cat(sprintf("%-36s null A-B gap %+.4f | within-grp SD (null) %.4f | spike diff %+.4f\n",
                label,
                mean(rowMeans(A[nul, ]) - rowMeans(B[nul, ])),
                mean(c(apply(A[nul, ], 1, sd), apply(B[nul, ], 1, sd))),
                mean(rowMeans(B[spike, ]) - rowMeans(A[spike, ]))))
}
run <- function(wa, ba, rbg) {
    suppressMessages(run_SMFnorm(pg, within_alpha = wa, between_alpha = ba,
                                 min_coverage = 0, rate_between_groups = rbg))$data
}

summ(pg, "raw")
summ(run(1.0, 1.0,  FALSE), "within 1.0, between OFF")
summ(run(0.5, 1.0,  FALSE), "within 0.5, between OFF")
summ(run(1.0, 1.0,  TRUE),  "within 1.0, between ON  alpha 1.0")
summ(run(1.0, 0.95, TRUE),  "within 1.0, between ON  alpha 0.95")
summ(run(0.5, 0.95, TRUE),  "within 0.5, between ON  alpha 0.95")
