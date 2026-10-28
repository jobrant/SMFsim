#!/usr/bin/env Rscript

# diagnose_bio_variation_null.R
#
# Does SMFnorm's within_alpha shrinkage produce false positives once replicates
# vary biologically? Synthetic null (no real data), raw vs SMFnorm at
# within_alpha 1.0 / 0.5 / 0.3, at three levels of region-level biological
# variation between replicates.
#
#   Rscript SMFsim/inst/scripts/diagnose_bio_variation_null.R [out_dir]
#   # default out_dir: results/diagnostics/bio_variation_null. ~10 minutes.
#
# Each replicate gets its own random effect per 25-site region (SD on the logit
# scale), on top of binomial sampling. Both groups draw from the same
# distribution, so every DMR is a false positive.
#
# Finding (2026-10-27), FP DMRs:
#   region SD 0    : raw 0, w1.0 0,  w0.5 0,   w0.3 0     <- what SMFsim simulates
#   region SD 0.25 : raw 0, w1.0 0,  w0.5 0,   w0.3 0
#   region SD 0.5  : raw 7, w1.0 19, w0.5 110, w0.3 177
#
# Sites come in 40-site clusters separated by 2 kb gaps. Without gaps larger
# than metilene's -M (300 bp) metilene tries to segment the whole chromosome
# as one block and is killed for running out of memory.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
# Internal package function (visible under load_all; fall back to the namespace
# when SMFsim is installed rather than loaded from source).
.n_tests_of <- if (exists(".metilene_n_tests")) .metilene_n_tests else
    getFromNamespace(".metilene_n_tests", "SMFsim")

args <- commandArgs(TRUE)
out <- if (length(args) >= 1) args[1] else "results/diagnostics/bio_variation_null"
metilene <- parse_args(character(0))$metilene_path

n <- 200000
cov_mean <- 30
block <- 25                                        # region = 25 sites, ~750 bp
set.seed(11)
pos <- cumsum(sample(15:45, n, TRUE) +
              ifelse(seq_len(n) %% 40 == 1, 2000L, 0L))   # ~30 bp within clusters, 2 kb gaps
lp  <- qlogis(rbeta(n, 1.2, 4))                    # true accessibility (logit)
blk <- (seq_len(n) - 1) %/% block
cov <- pmax(5L, rpois(n, cov_mean))

mkrep <- function(eff, bio_sd, seed) {
    set.seed(seed)
    b  <- rnorm(max(blk) + 1, 0, bio_sd)[blk + 1]  # replicate-specific REGION effect
    pr <- plogis(lp + b) * eff
    mc <- rbinom(n, cov, pr)
    data.table(chr = "chr1", pos = pos, strand = "+", site = seq_len(n),
               mc = mc, cov = cov, rate = mc / cov)
}

fp <- function(bio_sd) {
    effA <- c(0.95, 0.70, 0.85, 0.75)
    effB <- c(0.90, 0.65, 0.95, 0.75)
    pg <- list(
        PseudoA = setNames(lapply(1:4, function(i) mkrep(effA[i], bio_sd, 100 + i)),
                           paste0("PseudoA_", 1:4)),
        PseudoB = setNames(lapply(1:4, function(i) mkrep(effB[i], bio_sd, 200 + i)),
                           paste0("PseudoB_", 1:4)),
        params = list())
    arms <- list(raw = pg[c("PseudoA", "PseudoB")])
    for (wa in c(1.0, 0.5, 0.3)) {
        arms[[sprintf("SMFnorm_w%.1f", wa)]] <- suppressMessages(run_SMFnorm(
            pg, within_alpha = wa, between_alpha = 1, min_coverage = 0,
            rate_between_groups = FALSE))$data
    }
    vapply(names(arms), function(a) {
        dir <- file.path(out, bio_sd, a)
        d <- suppressMessages(call_dmrs_metilene(
            arms[[a]], dir, metilene, min_cpg = 10, min_diff = 0,
            metilene_max_dist = 300, q_cutoff = 0.05, q_source = "BH",
            p_column = "2dks", min_effect = 0.1, mtc = 2))
        # call_dmrs_metilene returns an empty table if metilene crashes; refuse
        # to report that as zero false positives.
        if (is.na(.n_tests_of(file.path(dir, "metilene_log.txt")))) {
            stop("metilene did not complete for ", a, " (bio_sd ", bio_sd,
                 "); see ", dir, call. = FALSE)
        }
        nrow(d)
    }, integer(1))
}

res <- rbindlist(lapply(c(0, 0.25, 0.5), function(s) {
    as.data.table(c(list(region_bio_sd_logit = s), as.list(fp(s))))
}))
print(res)
