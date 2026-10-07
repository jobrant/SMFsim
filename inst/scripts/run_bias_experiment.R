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
if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
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
# NOTE: do not set SMFSIM_ALLOW_METHOD_FAILURE here. ComBatMet's refusal on the
# aligned_* scenarios (perfect batch/group confounding) is now signalled as
# INAPPLICABLE and skipped automatically -- see .method_inapplicable() -- so the
# override is no longer needed to get past it. Setting it globally would also
# let a genuine SMFnorm failure through silently, which is how the 2026-09-02
# run finished with SMFnorm absent from every block.
# =========================================================================

# --- Simulation settings (shared across blocks) --------------------------
config$sim_mode          <- "parametric"            # independent group construction
config$dispersion_s      <- 26                      # calibrated to M-series
config$effect_sizes      <- c(0.10, 0.15, 0.20, 0.30)
config$seed              <- 42
config$standard_chr_only <- TRUE
config$metilene_min_cpg  <- 10
config$methods           <- c("raw", "downsampled", "SMFnorm", "ComBatMet")

# --- Which blocks to run ------------------------------------------------
# All of them by default. To rerun a subset -- e.g. after changing only the
# between-group scenarios -- name them in SMFSIM_BLOCKS:
#
#   SMFSIM_BLOCKS=bias_parametric Rscript SMFsim/inst/scripts/run_bias_experiment.R
#
# Each block writes only under <base_output>/<block>/, so a subset rerun leaves
# the other blocks' results untouched, and all_blocks.rds is rebuilt at the end
# from whatever block results are on disk.
known_blocks <- c("control_within", "bias_parametric", "bias_clone")
default_blocks <- c("control_within", "bias_parametric",
                    if (RUN_CLONE_BASELINE) "bias_clone")

blocks_env <- trimws(Sys.getenv("SMFSIM_BLOCKS"))
run_blocks <- if (nzchar(blocks_env)) {
    trimws(strsplit(blocks_env, ",")[[1]])
} else {
    default_blocks
}
unknown <- setdiff(run_blocks, known_blocks)
if (length(unknown)) {
    stop("SMFSIM_BLOCKS names unknown block(s): ", paste(unknown, collapse = ", "),
         ". Valid: ", paste(known_blocks, collapse = ", "), call. = FALSE)
}
message("Blocks to run: ", paste(run_blocks, collapse = ", "))

# --- Guard: never silently overwrite a previous run ----------------------
# Every block writes under a fixed path, so re-running lands on top of whatever
# is already there and destroys it in place. That is how an earlier run was
# lost -- not by anyone deleting it. Abort instead, and require the overwrite
# to be deliberate. Checked per block, so rerunning one block only requires
# archiving that block's directory.
#
# To overwrite on purpose:  SMFSIM_OVERWRITE=1 Rscript inst/scripts/run_bias_experiment.R
overwrite_ok <- identical(Sys.getenv("SMFSIM_OVERWRITE"), "1")

occupied <- Filter(function(b) {
    d <- file.path(base_output, b)
    dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE)) > 0
}, run_blocks)

if (length(occupied)) {
    paths <- normalizePath(file.path(base_output, occupied), mustWork = FALSE)
    if (!overwrite_ok) {
        stop(sprintf(
            paste0("Output for these blocks already exists and is not empty:\n  %s\n",
                   "Refusing to overwrite a previous run. Archive each under a ",
                   "dated name (e.g. mv %s %s_$(date +%%F)_superseded), point ",
                   "`base_output` at a new path, or re-run with ",
                   "SMFSIM_OVERWRITE=1 to overwrite deliberately."),
            paste(paths, collapse = "\n  "), paths[1], paths[1]),
            call. = FALSE)
    }
    message("SMFSIM_OVERWRITE=1 set -- overwriting existing results in:\n  ",
            paste(paths, collapse = "\n  "))
}

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

# (1) Specificity control: matched-mean within-group variation. Parametric so
#     the null is correctly calibrated; nothing systematic to correct, so
#     rate_between_groups = FALSE. Expect few/no DMRs for every method.
if ("control_within" %in% run_blocks) {
    run_block("control_within",
              scenarios    = c("mild", "moderate", "severe"),
              rate_between = FALSE)
}

# (2) The artifact: systematic between-group efficiency bias. rate_between_groups
#     = TRUE so SMFnorm can correct the between-group shift. Expect RAW to call
#     false positives (null) and biased DMRs (spike-in); normalization to reduce
#     them. ComBatMet is inapplicable in aligned_* (batch == group) and runs in
#     imbalanced_*, whose explicit efficiencies make the batch cross the groups.
if ("bias_parametric" %in% run_blocks) {
    run_block("bias_parametric",
              scenarios    = c("aligned_strong", "aligned_moderate",
                               "imbalanced_strong", "imbalanced_moderate"),
              rate_between = TRUE)
}

# (3) Optional contrast: the same bias scenarios under the old CLONE mode, to
#     show how the mis-specified (too-tight) null understates raw's FPs.
if ("bias_clone" %in% run_blocks) {
    run_block("bias_clone",
              scenarios    = c("aligned_strong", "imbalanced_strong"),
              rate_between = TRUE,
              sim_mode     = "clone")
}

# --- Combined results ----------------------------------------------------
# Rebuilt from every block result on disk, not just the blocks run this time,
# so a subset rerun still leaves a complete all_blocks.rds.
block_keys <- c(control_within = "control", bias_parametric = "bias",
                bias_clone = "bias_clone")
results <- list()
for (b in names(block_keys)) {
    rds <- file.path(base_output, b, paste0(b, "_results.rds"))
    if (file.exists(rds)) {
        r <- readRDS(rds)
        results[[block_keys[[b]]]] <- list(null = r$null, spikein = r$spikein)
    }
}
saveRDS(results, file.path(base_output, "all_blocks.rds"))
message("\nDONE (ran: ", paste(run_blocks, collapse = ", "), "). Results under: ",
        normalizePath(base_output, mustWork = FALSE))
