#!/usr/bin/env Rscript

# run_method_comparison_bias.R
#
# Compare normalization methods on the systematic between-group bias scenarios.
# This script is intended to test the methods in the regime where SMFnorm is
# expected to provide benefit.
#
# Usage (from the package root):
#   Rscript inst/scripts/run_method_comparison_bias.R

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
  devtools::load_all(".")
} else {
  library(SMFsim)
}

suppressWarnings(suppressMessages(library(data.table)))

config <- list(
  data_dir            = "data/allc",
  sample_sheet        = "data/sample_sheet.csv",
  output_dir          = "results/method_comparison_bias",
  metilene_path       = "/apps/metilene/0.2.8/metilene",
  scenarios           = c("aligned_moderate", "aligned_strong",
                            "imbalanced_moderate", "imbalanced_strong"),
  effect_sizes        = c(0.10, 0.15, 0.20, 0.30),
  n_spikein_regions   = 300,
  region_width_bp     = 2000,
  max_median_distance = 200,
  min_sites           = 10,
  wt_group_id         = "M1",
  seed                = 42,
  methods             = c("raw", "downsampled", "SMFnorm", "ComBatMet"),
  metilene_max_dist   = 300,
  metilene_min_cpg    = 10,
  metilene_min_diff   = 0.1,
  metilene_qval       = 0.05,
  metilene_q_source   = "metilene",
  within_alpha        = 0.3,
  between_alpha       = 0.8,
  min_coverage        = 5,
  dispersion_s        = Inf,
  sim_mode            = "parametric",
  rate_between_groups = TRUE
)

message(sprintf("Running method comparison on bias scenarios with between_alpha = %.1f and within_alpha = %.1f",
                config$between_alpha, config$within_alpha))
message(sprintf("Scenarios: %s", paste(config$scenarios, collapse = ", ")))

# Prepare WT replicates once.
wt_reps <- prepare_wt_replicates(config)

# Run null simulation (false positive assessment).
null_results <- run_null_simulation(wt_reps = wt_reps, config = config)

# Run spike-in simulation (sensitivity assessment).
spikein_results <- run_spikein_simulation(wt_reps = wt_reps, config = config)

# Save primary results.
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
results_rds <- file.path(config$output_dir, "method_comparison_bias_results.rds")
saveRDS(list(null = null_results, spikein = spikein_results, config = config), results_rds)

message(sprintf("Results written to %s", config$output_dir))

# Generate figures.
fig_dir <- file.path(config$output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tryCatch({
  generate_all_figures(results_rds, fig_dir)
  generate_additional_figures(
    file.path(config$output_dir, "spikein_simulation_results.csv"),
    fig_dir
  )
  message(sprintf("Figures written to %s", fig_dir))
}, error = function(e) {
  message("[figures skipped] ", conditionMessage(e))
})

# Print a quick summary table.
for (method in unique(spikein_results$method)) {
  dt <- spikein_results[method == !!method]
  tp <- sum(dt$TP)
  fp <- sum(dt$FP)
  fn <- sum(dt$FN)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  message(sprintf("%s: TP=%d, FP=%d, FN=%d, sens=%.4f, prec=%.4f",
                  method, tp, fp, fn, sens, prec))
}

message("Null false positives by method:")
for (method in unique(null_results$method)) {
  fp <- sum(null_results[method == !!method]$FP)
  message(sprintf("  %s: FP=%d", method, fp))
}
