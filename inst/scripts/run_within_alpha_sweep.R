#!/usr/bin/env Rscript

# run_within_alpha_sweep.R
#
# Sweep SMFnorm's `within_alpha` on the within-group variation scenarios to map the
# sensitivity/precision tradeoff. With rate_between_groups = FALSE, within_alpha
# sets how hard SMFnorm shrinks the within-group rate variation:
#   - low alpha  (e.g. 0.1): strongest correction -> groups pulled toward pooled
#     mean, reducing both signal and noise.
#   - high alpha (e.g. 0.9): weakest correction -> more original signal retained
#     but also more within-group noise/variation.
# The goal is to find the alpha that best trades sensitivity (spike-in) against
# specificity (null false positives) -- the highest alpha that keeps null FP low.
#
# Only SMFnorm is run (the other methods don't depend on within_alpha); uses the
# within-group scenarios (mild, moderate, severe) with rate_between_groups = FALSE.
# The between_alpha is held fixed at 0.8 (from prior tuning).
#
# Usage (from the package root):
#   Rscript inst/scripts/run_within_alpha_sweep.R

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".")
} else {
    library(SMFsim)
}
suppressWarnings(suppressMessages(library(data.table)))

# === EDIT THESE FOR YOUR ENVIRONMENT =====================================
config <- parse_args()
config$data_dir      <- "data/allc"
config$sample_sheet  <- "data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "results/within_alpha_sweep"

# within_alpha values to test.
ALPHA_GRID <- seq(0.1, 0.9, 0.1)

# FAST_TUNE: run on a subset of regions so the whole sweep finishes quickly.
# Use it to find the rough optimum, then re-confirm the chosen alpha on the
# full set with FAST_TUNE <- FALSE.
FAST_TUNE <- FALSE
# =========================================================================

# Fixed settings — within-group scenarios only, no between-group artifact.
config$sim_mode            <- "parametric"
# s = 26 is the M-series calibration (diagnose_variance.R). This was previously
# Inf (pure binomial), which gives a within-group SD of ~0.064 against the real
# ~0.083 - i.e. the within-group correction was being tuned on data ~20% less
# variable than the data it is meant to correct. run_alpha_sweep.R already used
# 26; this script was the odd one out.
config$dispersion_s        <- 26
config$rate_between_groups <- FALSE
config$between_alpha       <- 0.9            # held fixed; only within_alpha sweeps
config$methods             <- c("SMFnorm")   # only within_alpha-dependent method
config$scenarios           <- c("mild", "moderate", "severe")
config$effect_sizes        <- c(0.10, 0.15, 0.20, 0.30)
config$seed                <- 42

# Fast-tuning: restrict to chr1 and fewer spike-in regions (~10x faster).
if (FAST_TUNE) {
    config$chr_pattern       <- "^(chr)?1$"   # chr1 only
    config$n_spikein_regions <- 100
    base_output              <- "results/within_alpha_sweep_fast"
    message(sprintf("FAST_TUNE on: chr1 only, %d regions, output -> %s",
                    config$n_spikein_regions, base_output))
}

dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Load source replicates once.
wt_reps <- prepare_wt_replicates(config)

# SLURM array support: if SLURM_ARRAY_TASK_ID is set, this process runs ONLY
# that one alpha (one alpha per array task).
task_id    <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
array_mode <- nzchar(task_id)
if (array_mode) {
    idx <- suppressWarnings(as.integer(task_id))
    if (is.na(idx) || idx < 1L || idx > length(ALPHA_GRID)) {
        stop(sprintf("SLURM_ARRAY_TASK_ID=%s out of range 1..%d",
                     task_id, length(ALPHA_GRID)))
    }
    ALPHA_GRID <- ALPHA_GRID[idx]
    message(sprintf("Array mode: task %s -> within_alpha = %.2f",
                    task_id, ALPHA_GRID))
}

null_all <- list()
spk_all  <- list()
null_csv <- file.path(base_output, "within_alpha_sweep_null.csv")
spk_csv  <- file.path(base_output, "within_alpha_sweep_spikein.csv")

for (a in ALPHA_GRID) {
    message("\n", strrep("#", 64))
    message(sprintf("# within_alpha = %.2f", a))
    message(strrep("#", 64))

    cfg <- config
    cfg$within_alpha <- a
    cfg$output_dir   <- file.path(base_output, sprintf("alpha_%.2f", a))
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

    a_null <- file.path(cfg$output_dir, "null_simulation_results.csv")
    a_spk  <- file.path(cfg$output_dir, "spikein_simulation_results.csv")

    # Resume: skip an alpha whose per-alpha CSVs already exist (completed).
    if (file.exists(a_null) && file.exists(a_spk)) {
        message(sprintf("  [resume] alpha %.2f already done; reading cached", a))
        null <- fread(a_null)
        spk  <- fread(a_spk)
    } else {
        null <- run_null_simulation(wt_reps, cfg)
        spk  <- run_spikein_simulation(wt_reps, cfg)
    }
    null[, within_alpha := a]
    spk[, within_alpha := a]
    null_all[[length(null_all) + 1]] <- null
    spk_all[[length(spk_all) + 1]]  <- spk

    # Checkpoint the combined CSVs after EVERY alpha.
    if (!array_mode) {
        fwrite(rbindlist(null_all, fill = TRUE), null_csv)
        fwrite(rbindlist(spk_all,  fill = TRUE), spk_csv)
        message(sprintf("  [checkpoint] combined CSVs updated through alpha %.2f", a))
    }
}

if (array_mode) {
    message(sprintf(
        "\nArray task done (alpha %.2f). Run once WITHOUT SLURM_ARRAY_TASK_ID ",
        ALPHA_GRID))
    message("to stitch the combined CSVs + plot from all alpha_*/ dirs.")
    quit(save = "no")
}

null_dt <- rbindlist(null_all, fill = TRUE)
spk_dt  <- rbindlist(spk_all,  fill = TRUE)

message("\nDONE. Combined results:")
message("  null    -> ", null_csv)
message("  spikein -> ", spk_csv)

# Console summary: sensitivity / precision / F1 by alpha.
print(spk_dt[order(efficiency_scenario, effect_size, within_alpha),
             .(efficiency_scenario, effect_size, within_alpha,
               sens = sensitivity, prec = precision, F1)])

# Optional diagnostic plot: sensitivity & precision vs within_alpha.
tryCatch({
    suppressWarnings(suppressMessages(library(ggplot2)))
    pdt <- melt(spk_dt,
                id.vars = c("efficiency_scenario", "effect_size", "within_alpha"),
                measure.vars = c("sensitivity", "precision"),
                variable.name = "metric", value.name = "value")
    p <- ggplot(pdt, aes(within_alpha, value, colour = metric)) +
        geom_line() + geom_point() +
        facet_grid(effect_size ~ efficiency_scenario) +
        labs(title = "SMFnorm sensitivity/precision vs within_alpha",
             x = "within_alpha", y = NULL, colour = NULL) +
        theme_minimal()
    ggsave(file.path(base_output, "within_alpha_sweep_tradeoff.png"), p,
           width = 11, height = 8, dpi = 150)
    message("  plot    -> ", file.path(base_output, "within_alpha_sweep_tradeoff.png"))
}, error = function(e) message("  [plot skipped] ", conditionMessage(e)))
