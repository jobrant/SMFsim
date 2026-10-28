#!/usr/bin/env Rscript

# calibrate_replicate_noise.R
#
# Calibrate the simulator's replicate noise to a real group of replicates:
#   dispersion_s  - Beta precision (overall overdispersion)
#   fragment_rho  - shared-fragment noise (correlation between nearby sites)
#
#   Rscript SMFsim/inst/scripts/calibrate_replicate_noise.R [group] [chr_pattern]
#   # defaults: M1, chromosome 1. Run from the manuscript root. A few minutes.
#
# Two statistics are computed identically on the real replicates and on
# replicates simulated from the same sites (pooled true rate, real coverage,
# each replicate's real relative efficiency), using .draw_parametric_rep():
#   var_ratio  mean squared binomial-standardised residual, after removing
#              each replicate's efficiency per rate bin with SMFnorm's own
#              within-set step (alpha = 1), so efficiency differences of any
#              shape do not count as noise. Pure binomial noise gives about
#              (k-1)/k; overdispersion raises it. Set by dispersion_s.
#
# Why per-bin: in real M1 the N3/N2 rate ratio is ~1.4 at low rates, ~1.8 in
# the middle and ~1.4 at high rates, so efficiency is neither a constant
# offset nor a single multiplier. Removing it either way leaves efficiency in
# the residuals (var_ratio 1.01 vs 10.3 for the same data).
#   cor_50     correlation of those residuals between neighbouring sites
#              1-50 bp apart. Set by fragment_rho (and fragment_bp); the
#              per-site distribution, and so var_ratio, is unaffected because
#              the copula keeps each count Binomial given its latent rate.
# dispersion_s is chosen first with fragment_rho = 0, then fragment_rho for
# each fragment_bp at that dispersion. cor_150 (51-150 bp) checks the decay:
# real data show ~0.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
.ns <- function(f) if (exists(f)) get(f) else getFromNamespace(f, "SMFsim")
.draw   <- .ns(".draw_parametric_rep")
.pooled <- .ns(".pooled_rate")

args <- commandArgs(TRUE)
group <- if (length(args) >= 1) args[1] else "M1"
chr_pattern <- if (length(args) >= 2) args[2] else "^(chr)?1$"

cfg <- parse_args(character(0))
cfg$wt_group_id <- group
cfg$chr_pattern <- chr_pattern
reps <- suppressMessages(prepare_wt_replicates(cfg))
k <- length(reps)
message(sprintf("%s, %s: %d replicates, %d shared sites (data: %s)",
                group, chr_pattern, k, nrow(reps[[1]]), cfg$data_dir))

align <- getFromNamespace(".normalize_methylation_within_set", "SMFnorm")

noise_stats <- function(rs) {
    rs <- suppressMessages(align(lapply(rs, copy), alpha = 1, diagnostics = FALSE))
    R <- do.call(cbind, lapply(rs, `[[`, "rate"))
    C <- do.call(cbind, lapply(rs, `[[`, "cov"))
    m <- rowSums(R * C) / rowSums(C)
    m <- pmin(pmax(m, 1e-4), 1 - 1e-4)
    Z <- (R - m) / sqrt(m * (1 - m) / C)
    d <- diff(rs[[1]]$pos)
    same_chr <- rs[[1]]$chr[-1] == rs[[1]]$chr[-nrow(R)]
    nb <- function(lo, hi) {
        i <- which(same_chr & d >= lo & d <= hi)
        mean(vapply(seq_len(ncol(Z)), function(j) cor(Z[i, j], Z[i + 1, j]), numeric(1)))
    }
    c(var_ratio = mean(Z^2), cor_50 = nb(1, 50), cor_150 = nb(51, 150))
}

real <- noise_stats(reps)
message(sprintf("REAL       var_ratio %.3f  cor_50 %.3f  cor_150 %.3f",
                real["var_ratio"], real["cor_50"], real["cor_150"]))

pooled <- .pooled(reps)
rel_eff <- vapply(reps, function(dt) mean(dt$rate), numeric(1))
rel_eff <- rel_eff / max(rel_eff)

simulate <- function(s, rho, seed = 1, bp = cfg$fragment_bp) {
    noise <- list(fragment_bp = bp, fragment_rho = rho,
                  bio_region_bp = 1000L, bio_region_sd = 0)
    sim <- lapply(seq_len(k), function(i) {
        .draw(reps[[i]], pooled, rel_eff[i], seed + i, s,
              if (rho > 0) noise else NULL)
    })
    noise_stats(sim)
}

# Step 1: dispersion_s, fragment noise off.
s_grid <- c(15, 26, 50, 100, 200, 400, 1000, 1e5)
s_tab <- rbindlist(lapply(s_grid, function(s) {
    st <- simulate(s, 0)
    message(sprintf("  s = %5.1f            var_ratio %.3f", s, st["var_ratio"]))
    data.table(s = s, var_ratio = st["var_ratio"])
}))
# var_ratio falls as s rises; interpolate on log(s).
s_best <- exp(approx(s_tab$var_ratio, log(s_tab$s), xout = real["var_ratio"],
                     rule = 2)$y)
if (real["var_ratio"] <= min(s_tab$var_ratio)) {
    message("   real data are at or below the pure-binomial floor: no Beta ",
            "overdispersion needed")
    s_best <- Inf
}
message(sprintf("-> dispersion_s ~ %s", format(round(s_best, 1))))

# Step 2: fragment_rho at that dispersion, for each candidate block width.
rho_grid <- c(0, 0.2, 0.4, 0.6, 0.8, 0.95)
fits <- rbindlist(lapply(c(50L, 100L), function(bp) {
    rho_tab <- rbindlist(lapply(rho_grid, function(r) {
        st <- simulate(s_best, r, bp = bp)
        message(sprintf("  bp = %3d  rho = %.2f  var_ratio %.3f  cor_50 %.3f  cor_150 %.3f",
                        bp, r, st["var_ratio"], st["cor_50"], st["cor_150"]))
        data.table(rho = r, cor_50 = st["cor_50"], cor_150 = st["cor_150"])
    }))
    rho_best <- approx(rho_tab$cor_50, rho_tab$rho, xout = real["cor_50"], rule = 2)$y
    reached <- max(rho_tab$cor_50) >= real["cor_50"]
    final <- simulate(s_best, rho_best, seed = 99, bp = bp)
    message(sprintf("-> bp = %d: fragment_rho ~ %.2f%s | check (new seed): var_ratio %.3f vs %.3f, cor_50 %.3f vs %.3f, cor_150 %.3f vs %.3f",
                    bp, rho_best, if (reached) "" else " (target NOT reached)",
                    final["var_ratio"], real["var_ratio"], final["cor_50"], real["cor_50"],
                    final["cor_150"], real["cor_150"]))
    data.table(fragment_bp = bp, fragment_rho = rho_best, reached = reached,
               var_ratio = final["var_ratio"], cor_50 = final["cor_50"],
               cor_150 = final["cor_150"])
}))
fits[, error := abs(cor_50 - real["cor_50"]) + abs(cor_150 - real["cor_150"])]
best <- fits[which.min(error)]
cat(sprintf("\nCALIBRATION %s: dispersion_s = %s, fragment_bp = %d, fragment_rho = %.2f\n",
            group, format(round(s_best, 1)), best$fragment_bp, best$fragment_rho))
cat(sprintf("  real      var_ratio %.3f  cor_50 %.3f  cor_150 %.3f\n",
            real["var_ratio"], real["cor_50"], real["cor_150"]))
cat(sprintf("  simulated var_ratio %.3f  cor_50 %.3f  cor_150 %.3f\n",
            best$var_ratio, best$cor_50, best$cor_150))
