#!/usr/bin/env Rscript

# diagnose_variance.R
#
# Diagnostic for the SMFsim noise model. Answers one question:
#   "How far below real data does the simulated noise floor sit, and what
#    beta-binomial dispersion would I need to close the gap?"
#
# It measures three things and writes a report + CSVs + plots:
#
#   1. REAL overdispersion. For the source replicates (e.g. PrEC), compares the
#      observed across-replicate variance of the methylation rate to the
#      variance a pure Binomial would predict. The ratio (VIF) is 1 under
#      Binomial and > 1 under real beta-binomial overdispersion. Back-solves a
#      calibrated precision parameter s = alpha + beta for a Beta prior on the
#      per-site rate, which is what the simulator would need to match reality.
#
#   2. BETWEEN-GROUP deficit. The null simulation clones the SAME source reps
#      into both pseudo-groups, so the A-vs-B contrast carries no INDEPENDENT
#      biological variation. This compares the simulated between-group per-site
#      rate difference against the difference between two INDEPENDENT real
#      same-condition replicates. If the simulated contrast is much tighter,
#      that is why metilene never fires on raw data in the null.
#
#   3. WITHIN-GROUP spread. Real across-replicate SD vs simulated pseudo-group
#      SD, per efficiency scenario. Clone mode inherits the real within-group
#      spread, so this isolates what the efficiency knob adds on top.
#
# Usage:
#   devtools::load_all()          # makes internal sim functions available
#   source("inst/scripts/diagnose_variance.R")
#   # edit the config block below to point at your real data first

suppressWarnings(suppressMessages({
    library(data.table)
    library(ggplot2)
}))

# --- Config (mirror your normal run) -------------------------------------
config <- list(
    data_dir     = "../m-series-data/allc",
    sample_sheet = "../m-series-data/sample_sheet.csv",
    output_dir   = "results/simulation",
    wt_group_id  = "M1",
    scenarios    = c("mild", "moderate", "severe"),
    seed         = 42,
    min_coverage = 5
)

# Keep only standard chromosomes. Non-standard contigs (e.g. *_random, Un_*)
# break SMFnorm's numeric site id (.chr_to_int -> NA) and the shared-site
# intersection, leaving replicates with mismatched row counts.
STD_CHR_PATTERN <- "^(chr)?([0-9]{1,2}|X|Y|M|MT)$"

diag_dir <- file.path(config$output_dir, "m-series-diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# Internal sim functions must be loaded (devtools::load_all()).
for (fn in c("prepare_wt_replicates", "create_pseudo_groups",
             "get_efficiency_scenarios")) {
    if (!exists(fn, mode = "function")) {
        stop("`", fn, "` not found. Run devtools::load_all() before sourcing.")
    }
}

# Number of rate strata for binned, coverage-aware moment estimates.
N_BINS <- 20


# --- Replicate filtering / alignment -------------------------------------

# Drop non-standard contigs so every replicate covers the same chromosome set.
filter_standard_chr <- function(reps, pattern = STD_CHR_PATTERN) {
    lapply(reps, function(dt) dt[grepl(pattern, chr)])
}

# Inner-join all replicates onto their common sites, returning equal-length,
# identically ordered data.tables. Keys on chr+pos (coordinate based) to avoid
# the numeric uid path that NA-coerces on non-standard contigs. Does not trust
# find_shared_sites() to have aligned rows.
align_replicates <- function(reps, key = c("chr", "pos")) {
    # allc should hold one row per cytosine; drop accidental duplicate keys.
    reps <- lapply(reps, function(dt) {
        if (anyDuplicated(dt, by = key)) dt <- unique(dt, by = key)
        dt
    })
    common <- Reduce(function(a, b) merge(a, b, by = key),
                     lapply(reps, function(dt) dt[, ..key]))
    aligned <- lapply(reps, function(dt) merge(common, dt, by = key, sort = TRUE))
    lens <- vapply(aligned, nrow, integer(1))
    if (length(unique(lens)) != 1L) {
        stop("Alignment failed; replicate lengths differ: ",
             paste(lens, collapse = ", "))
    }
    aligned
}

# Fit an efficiency vector to n replicates (scenarios are length-3; a source
# group may have a different rep count, e.g. M1 has 4).
.fit_eff <- function(eff, n) if (length(eff) == n) eff else rep_len(eff, n)


# --- Matrix helpers ------------------------------------------------------

# Stack a list of aligned replicate data.tables into a site x replicate matrix
# for one column. Rows must already be aligned (see align_replicates()).
.rep_matrix <- function(reps, col) {
    m <- vapply(reps, function(dt) as.numeric(dt[[col]]),
                FUN.VALUE = numeric(nrow(reps[[1]])))
    if (is.null(dim(m))) m <- matrix(m, ncol = length(reps))
    m
}

# Row-wise sample variance (denominator R - 1) of a site x replicate matrix.
.row_var <- function(m) {
    R <- ncol(m)
    if (R < 2) return(rep(NA_real_, nrow(m)))
    mu <- rowMeans(m)
    rowSums((m - mu)^2) / (R - 1)
}


# --- 1. Real overdispersion ----------------------------------------------

# For each site: pooled rate p, observed across-rep variance of the rate, and
# the Binomial expectation p(1-p) * mean(1/cov). Aggregated per rate bin with a
# coverage-aware method-of-moments dispersion (rho) and precision (s = 1/rho-1).
#
# Model: rate_i ~ mean p, beta-binomial intra-class correlation rho, so
#   Var(rate_i) = p(1-p) * [1/n_i + (1 - 1/n_i) * rho].
# Solving the bin-level moment equation for rho:
#   rho = (mean[s2] - mean[p(1-p)/n]) / mean[p(1-p)(1 - 1/n)].
measure_overdispersion <- function(reps, n_bins = N_BINS) {
    R <- length(reps)
    if (R < 2) stop("Need >= 2 replicates to estimate overdispersion.")

    mc  <- .rep_matrix(reps, "mc")
    cov <- .rep_matrix(reps, "cov")
    rate <- mc / cov

    pooled_p <- rowSums(mc) / rowSums(cov)
    obs_var  <- .row_var(rate)
    mean_inv_cov <- rowMeans(1 / cov)
    pq <- pooled_p * (1 - pooled_p)
    binom_var <- pq * mean_inv_cov

    site <- data.table(
        pooled_p = pooled_p,
        obs_var = obs_var,
        binom_var = binom_var,
        pq = pq,
        mean_inv_cov = mean_inv_cov,
        mean_cov = rowMeans(cov),
        bin = cut(pooled_p, breaks = seq(0, 1, length.out = n_bins + 1),
                  include.lowest = TRUE)
    )
    # Drop invariant sites (p = 0 or 1: no variance to explain).
    site <- site[is.finite(obs_var) & pq > 0]

    by_bin <- site[, {
        m_s2  <- mean(obs_var)
        m_bin <- mean(binom_var)                       # mean p(1-p)/n
        m_pq_within <- mean(pq * (1 - mean_inv_cov))   # mean p(1-p)(1 - 1/n)
        rho <- (m_s2 - m_bin) / m_pq_within
        rho <- max(min(rho, 0.999), 0)
        .(
            n_sites   = .N,
            mid_rate  = mean(pooled_p),
            mean_cov  = mean(mean_cov),
            obs_var   = m_s2,
            binom_var = m_bin,
            VIF       = m_s2 / m_bin,                  # variance inflation factor
            rho       = rho,
            precision_s = if (rho > 0) (1 - rho) / rho else Inf
        )
    }, by = bin][order(mid_rate)]

    # Overall coverage-weighted moment estimate (weight by site count).
    overall_rho <- site[, {
        m_s2  <- mean(obs_var)
        m_bin <- mean(binom_var)
        m_pq_within <- mean(pq * (1 - mean_inv_cov))
        max(min((m_s2 - m_bin) / m_pq_within, 0.999), 0)
    }]
    overall <- list(
        n_sites = nrow(site),
        median_VIF = median(site$obs_var / site$binom_var, na.rm = TRUE),
        rho = overall_rho,
        precision_s = if (overall_rho > 0) (1 - overall_rho) / overall_rho else Inf
    )

    list(by_bin = by_bin, overall = overall, site = site)
}


# --- 2. Between-group deficit --------------------------------------------

# All unique pairwise |rate_i - rate_j| among real replicates: how far apart two
# INDEPENDENT same-condition samples sit. This is the realistic null contrast.
real_independent_diffs <- function(reps) {
    rate <- .rep_matrix(reps, "rate")
    R <- ncol(rate)
    pairs <- combn(R, 2)
    diffs <- lapply(seq_len(ncol(pairs)), function(k) {
        i <- pairs[1, k]; j <- pairs[2, k]
        data.table(pair = paste0(i, "v", j),
                   abs_diff = abs(rate[, i] - rate[, j]))
    })
    rbindlist(diffs)
}

# Simulated between-group per-site contrast from clone-mode pseudo-groups:
# |mean(rate_A) - mean(rate_B)| plus the within-group SDs that metilene sees.
simulated_group_contrast <- function(pseudo) {
    rate_A <- .rep_matrix(pseudo$PseudoA, "rate")
    rate_B <- .rep_matrix(pseudo$PseudoB, "rate")
    data.table(
        meanA   = rowMeans(rate_A),
        meanB   = rowMeans(rate_B),
        abs_diff = abs(rowMeans(rate_A) - rowMeans(rate_B)),
        sd_A    = sqrt(.row_var(rate_A)),
        sd_B    = sqrt(.row_var(rate_B))
    )
}


# --- Run -----------------------------------------------------------------

message("Loading real replicates (", config$wt_group_id, ")...")
wt_reps <- prepare_wt_replicates(config)

# Defensive: filter to standard chromosomes and re-align rows ourselves, since
# upstream shared-site filtering can leave ragged replicates on data with
# non-standard contigs.
n_before <- vapply(wt_reps, nrow, integer(1))
wt_reps <- filter_standard_chr(wt_reps)
wt_reps <- align_replicates(wt_reps)
R <- length(wt_reps)
message(sprintf("Replicate rows before align: %s",
                paste(n_before, collapse = ", ")))
message(sprintf("Loaded %d replicates, %d aligned shared sites each.",
                R, nrow(wt_reps[[1]])))

# 1. Real overdispersion
message("\n[1/3] Measuring real overdispersion...")
od <- measure_overdispersion(wt_reps)
fwrite(od$by_bin, file.path(diag_dir, "overdispersion_by_rate_bin.csv"))

# 2/3. Per-scenario simulated contrast + within-group spread
message("[2/3] Simulating null pseudo-groups per scenario...")
scenarios <- get_efficiency_scenarios()[config$scenarios]

real_diffs <- real_independent_diffs(wt_reps)
real_within_sd <- sqrt(.row_var(.rep_matrix(wt_reps, "rate")))

contrast_summ <- list()
diff_density  <- list(real_diffs[, .(source = "real_independent_reps",
                                     abs_diff)])
# Most-severe scenario in the requested set: used for the ECDF plot so we can
# compare clone vs parametric without overcrowding the figure.
plot_scenario <- tail(names(scenarios), 1)

s_label <- function(s) if (is.finite(s)) as.character(s) else "Inf"

# (mode, dispersion_s) configurations to compare: clone vs parametric at the
# binomial baseline (Inf) and the calibrated precision (26 ~ M-series).
sweep_configs <- list(
    list(mode = "clone",      s = Inf),
    list(mode = "clone",      s = 26),
    list(mode = "parametric", s = Inf),
    list(mode = "parametric", s = 26)
)

for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    for (cfg in sweep_configs) {
        pseudo <- create_pseudo_groups(
            replicates = wt_reps,
            efficiency_A = .fit_eff(sc$efficiency_A, R),
            efficiency_B = .fit_eff(sc$efficiency_B, R),
            mode = cfg$mode,
            seed = config$seed,
            dispersion_s = cfg$s
        )
        gc_dt <- simulated_group_contrast(pseudo)
        tag <- sprintf("%s_s%s", cfg$mode, s_label(cfg$s))
        key <- paste0(sc_name, "_", tag)

        contrast_summ[[key]] <- data.table(
            scenario = sc_name,
            mode = cfg$mode,
            dispersion_s = cfg$s,
            median_between_diff_sim = median(gc_dt$abs_diff),
            q90_between_diff_sim     = quantile(gc_dt$abs_diff, 0.90),
            median_within_sd_sim     = median(c(gc_dt$sd_A, gc_dt$sd_B), na.rm = TRUE),
            median_within_sd_real    = median(real_within_sd, na.rm = TRUE)
        )
        if (sc_name == plot_scenario) {
            diff_density[[key]] <- gc_dt[, .(source = tag, abs_diff)]
        }
    }
}
contrast_dt <- rbindlist(contrast_summ)
fwrite(contrast_dt, file.path(diag_dir, "between_group_contrast.csv"))


# --- Report --------------------------------------------------------------

cat("\n", strrep("=", 70), "\n", sep = "")
cat("SMFsim VARIANCE DIAGNOSTIC\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("1. REAL OVERDISPERSION (", config$wt_group_id, ", ", R, " reps)\n", sep = "")
cat("   A pure Binomial simulator produces VIF = 1 (rho = 0, precision = Inf).\n")
cat(sprintf("   Real median VIF   : %.2f   (1.0 = no overdispersion)\n",
            od$overall$median_VIF))
cat(sprintf("   Real dispersion   : rho = %.4f\n", od$overall$rho))
cat(sprintf("   Calibrated Beta precision s = alpha+beta = %.1f\n",
            od$overall$precision_s))
cat("   -> To match reality, draw per-site rates from Beta(mu*s, (1-mu)*s)\n")
cat("      BEFORE Binomial sampling, with s of this order. Smaller s = noisier.\n")
cat("   Per-rate-bin detail: diagnostics/overdispersion_by_rate_bin.csv\n\n")

cat("2. CLONE vs PARAMETRIC (within-group SD + between-group contrast)\n")
cat(sprintf("   Real within-group SD (median): %.4f  <- target for sim within-SD\n",
            median(real_within_sd, na.rm = TRUE)))
cat(sprintf("   Real single-rep |Δrate| (median): %.4f  (scale ref only; not a\n",
            median(real_diffs$abs_diff)))
cat("                                            target for a diff-of-means)\n")
cat(sprintf("     %-9s %-11s %5s  %10s  %10s\n",
            "scenario", "mode", "s", "within_SD", "between"))
for (i in seq_len(nrow(contrast_dt))) {
    r <- contrast_dt[i]
    cat(sprintf("     %-9s %-11s %5s  %10.4f  %10.4f\n",
                r$scenario, r$mode, s_label(r$dispersion_s),
                r$median_within_sd_sim, r$median_between_diff_sim))
}
cat("\n   Goal: parametric@s=26 matches the real within-SD (~0.083) AND keeps\n")
cat("   a healthy between-group contrast. Clone@s=26 overshoots within-SD\n")
cat("   (double-counts) while its between-group is driven only by the\n")
cat("   independent dispersion, since shared biology cancels.\n")
cat("\n", strrep("=", 70), "\n", sep = "")


# --- Plots ---------------------------------------------------------------

message("Writing plots to ", diag_dir, " ...")

# VIF vs rate bin: how overdispersion varies with methylation level.
p1 <- ggplot(od$by_bin, aes(mid_rate, VIF)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line() + geom_point(aes(size = n_sites)) +
    labs(title = "Real overdispersion vs methylation rate",
         subtitle = "Dashed line = Binomial floor (current simulator)",
         x = "Pooled methylation rate", y = "Variance inflation factor (VIF)",
         size = "n sites") +
    theme_minimal()
ggsave(file.path(diag_dir, "fig_vif_by_rate.png"), p1,
       width = 7, height = 4.5, dpi = 150)

# Distribution of per-site rate differences: real independent vs simulated.
diff_dt <- rbindlist(diff_density)
p2 <- ggplot(diff_dt, aes(abs_diff, colour = source)) +
    stat_ecdf() +
    coord_cartesian(xlim = c(0, quantile(diff_dt$abs_diff, 0.99))) +
    labs(title = "Per-site rate difference: real vs simulated null contrast",
         subtitle = sprintf(
             "Clone vs parametric (scenario '%s'); closer to real = better",
             plot_scenario),
         x = "|rate difference|", y = "ECDF", colour = NULL) +
    theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(diag_dir, "fig_between_group_diff_ecdf.png"), p2,
       width = 7, height = 5, dpi = 150)

message("Done. See ", diag_dir)
