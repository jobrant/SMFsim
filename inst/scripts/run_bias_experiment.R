#!/usr/bin/env Rscript

# run_bias_experiment.R
#
# Non-interactive driver for the between-group efficiency-bias benchmark.
# Designed for batch execution on an HPC node, e.g. from the package root:
#
#   Rscript inst/scripts/run_bias_experiment.R
#
# It runs three logical blocks and writes, under <base_output>/<block>/:
#   - null_simulation_results.csv        (FP counts by scenario x method)
#   - spikein_simulation_results.csv     (TP/FP/FN + metrics)
#   - <block>_results.rds                (list(null=, spikein=, config=))
# plus an all_blocks.rds at the top level.
#
# The blocks:
#   control_within  -- matched-mean within-group variation (specificity control)
#   bias_parametric -- systematic between-group efficiency bias (the artifact)
#   bias_clone      -- optional: same bias under the old clone mode, for contrast
#
# Calibration (s = 26, parametric mode) comes from inst/scripts/diagnose_variance.R.

# --- Load package --------------------------------------------------------
if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all("SMFsim/")
} else {
    library(SMFsim)
}

# === EDIT THESE FOR YOUR ENVIRONMENT =====================================
config <- parse_args()
config$data_dir      <- "data/allc"
config$sample_sheet  <- "data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "results/bias_experiment"

# Optional extra runs / outputs (off by default to keep the job lean).
RUN_CLONE_BASELINE <- FALSE   # rerun bias scenarios under clone mode for contrast
MAKE_FIGURES       <- TRUE   # generate manuscript figures after each block
# =========================================================================

# --- Simulation settings (shared across blocks) --------------------------
config$sim_mode          <- "parametric"            # independent group construction
config$dispersion_s      <- 26                      # calibrated to M-series
config$effect_sizes      <- c(0.10, 0.15, 0.20, 0.30)
config$seed              <- 42
config$standard_chr_only <- TRUE
# DMR-calling thresholds and method list keep parse_args() defaults; override
# here if needed, e.g.:
# config$methods         <- c("raw", "downsampled", "SMFnorm", "ComBatMet")
# config$metilene_min_cpg <- 10

dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# --- Load and prepare source replicates once -----------------------------
wt_reps <- prepare_wt_replicates(config)

# --- Block runner --------------------------------------------------------
# Each block gets its own output_dir so per-run summary CSVs do not collide.
run_block <- function(block, scenarios, rate_between,
                      sim_mode = config$sim_mode, do_spikein = TRUE) {
    cfg <- config
    cfg$scenarios           <- scenarios
    cfg$rate_between_groups <- rate_between
    cfg$sim_mode            <- sim_mode
    cfg$output_dir          <- file.path(base_output, block)
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

    message("\n", strrep("#", 64))
    message(sprintf("# BLOCK: %s  (mode=%s, rate_between_groups=%s)",
                    block, sim_mode, rate_between))
    message(sprintf("#   scenarios: %s", paste(scenarios, collapse = ", ")))
    message(strrep("#", 64))

    null    <- run_null_simulation(wt_reps, cfg)
    spikein <- if (do_spikein) run_spikein_simulation(wt_reps, cfg) else NULL

    rds <- file.path(cfg$output_dir, paste0(block, "_results.rds"))
    saveRDS(list(null = null, spikein = spikein, config = cfg), rds)

    if (MAKE_FIGURES && do_spikein) {
        fig_dir <- file.path(cfg$output_dir, "figures")
        tryCatch({
            generate_all_figures(rds, fig_dir)
            generate_additional_figures(
                file.path(cfg$output_dir, "spikein_simulation_results.csv"),
                fig_dir)
        }, error = function(e)
            message("  [figures skipped] ", conditionMessage(e)))
    }

    list(null = null, spikein = spikein)
}

# --- Blocks --------------------------------------------------------------
results <- list()

# (1) Specificity control: matched-mean within-group variation. Parametric so
#     the null is correctly calibrated; nothing systematic to correct, so
#     rate_between_groups = FALSE. Expect few/no DMRs for every method.
results$control <- run_block(
    "control_within",
    scenarios    = c("mild", "moderate", "severe"),
    rate_between = FALSE)

# (2) The artifact: systematic between-group efficiency bias. rate_between_groups
#     = TRUE so SMFnorm can correct the between-group shift. Expect RAW to call
#     false positives (null) and biased DMRs (spike-in); normalization to reduce
#     them. The imbalanced_* scenarios (overlapping per-sample efficiencies)
#     stress label-based correction (ComBatMet) the most.
results$bias <- run_block(
    "bias_parametric",
    scenarios    = c("aligned_strong", "aligned_moderate",
                     "imbalanced_strong", "imbalanced_moderate"),
    rate_between = TRUE)

# (3) Optional contrast: the same bias scenarios under the old CLONE mode, to
#     show how the mis-specified (too-tight) null understates raw's FPs.
if (RUN_CLONE_BASELINE) {
    results$bias_clone <- run_block(
        "bias_clone",
        scenarios    = c("aligned_strong", "imbalanced_strong"),
        rate_between = TRUE,
        sim_mode     = "clone")
}

saveRDS(results, file.path(base_output, "all_blocks.rds"))
message("\nDONE. Results under: ",
        normalizePath(base_output, mustWork = FALSE))
