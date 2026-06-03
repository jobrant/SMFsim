#!/usr/bin/env Rscript

# run_alpha_sweep.R
#
# Sweep SMFnorm's `between_alpha` on the systematic-bias scenarios to map the
# sensitivity/precision tradeoff. With rate_between_groups = TRUE, between_alpha
# sets how hard SMFnorm shrinks the between-group correction:
#   - high alpha (e.g. 0.5): aggressive -> removes the artifact but also most of
#     the spike signal (high precision, low sensitivity).
#   - low alpha  (e.g. 0.1): gentle -> retains more spike signal but leaves more
#     residual artifact (higher sensitivity, more false positives).
# The goal is to find the alpha that best trades sensitivity (spike-in) against
# specificity (null false positives).
#
# Only SMFnorm is run (the other methods don't depend on alpha); compare against
# the raw/downsampled baselines from your main bias run. Uses the corrected
# pipeline (parametric mode, s=26, metilene-q, region-level recall).
#
# Usage (from the package root):
#   Rscript inst/scripts/run_alpha_sweep.R

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".")
} else {
    library(SMFsim)
}
suppressWarnings(suppressMessages(library(data.table)))

# === EDIT THESE FOR YOUR ENVIRONMENT =====================================
config <- parse_args()
config$data_dir      <- "../m-series-data/allc"
config$sample_sheet  <- "../m-series-data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "results/alpha_sweep"

# between_alpha values to test. Start coarse; extend toward 0.9 if the optimum
# looks to be at the high end.
ALPHA_GRID <- seq(0.1, 0.5, 0.1)
# =========================================================================

# Fixed settings — match the main bias run.
config$sim_mode            <- "parametric"
config$dispersion_s        <- 26
config$rate_between_groups <- TRUE
config$within_alpha        <- 0.3            # held fixed; only between_alpha sweeps
config$metilene_q_source   <- "metilene"
config$methods             <- c("SMFnorm")   # only alpha-dependent method
config$scenarios           <- c("aligned_moderate", "aligned_strong",
                                "imbalanced_moderate", "imbalanced_strong")
config$effect_sizes        <- c(0.10, 0.15, 0.20, 0.30)
config$seed                <- 42

dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Load source replicates once.
wt_reps <- prepare_wt_replicates(config)

null_all <- list()
spk_all  <- list()

for (a in ALPHA_GRID) {
    message("\n", strrep("#", 64))
    message(sprintf("# between_alpha = %.2f", a))
    message(strrep("#", 64))

    cfg <- config
    cfg$between_alpha <- a
    cfg$output_dir    <- file.path(base_output, sprintf("alpha_%.2f", a))
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

    null <- run_null_simulation(wt_reps, cfg)
    null[, between_alpha := a]
    null_all[[length(null_all) + 1]] <- null

    spk <- run_spikein_simulation(wt_reps, cfg)
    spk[, between_alpha := a]
    spk_all[[length(spk_all) + 1]] <- spk
}

null_dt <- rbindlist(null_all, fill = TRUE)
spk_dt  <- rbindlist(spk_all, fill = TRUE)
fwrite(null_dt, file.path(base_output, "alpha_sweep_null.csv"))
fwrite(spk_dt,  file.path(base_output, "alpha_sweep_spikein.csv"))

message("\nDONE. Combined results:")
message("  null    -> ", file.path(base_output, "alpha_sweep_null.csv"))
message("  spikein -> ", file.path(base_output, "alpha_sweep_spikein.csv"))

# Console summary: sensitivity / precision / F1 by alpha.
print(spk_dt[order(efficiency_scenario, effect_size, between_alpha),
             .(efficiency_scenario, effect_size, between_alpha,
               sens = sensitivity, prec = precision, F1)])

# Optional diagnostic plot: sensitivity & precision vs between_alpha.
tryCatch({
    suppressWarnings(suppressMessages(library(ggplot2)))
    pdt <- melt(spk_dt,
                id.vars = c("efficiency_scenario", "effect_size", "between_alpha"),
                measure.vars = c("sensitivity", "precision"),
                variable.name = "metric", value.name = "value")
    p <- ggplot(pdt, aes(between_alpha, value, colour = metric)) +
        geom_line() + geom_point() +
        facet_grid(effect_size ~ efficiency_scenario) +
        labs(title = "SMFnorm sensitivity/precision vs between_alpha",
             x = "between_alpha", y = NULL, colour = NULL) +
        theme_minimal()
    ggsave(file.path(base_output, "alpha_sweep_tradeoff.png"), p,
           width = 11, height = 8, dpi = 150)
    message("  plot    -> ", file.path(base_output, "alpha_sweep_tradeoff.png"))
}, error = function(e) message("  [plot skipped] ", conditionMessage(e)))
