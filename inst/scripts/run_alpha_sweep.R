#!/usr/bin/env Rscript

# run_alpha_sweep.R
#
# Sweep SMFnorm's `between_alpha` on the systematic-bias scenarios to map the
# sensitivity/precision tradeoff. With rate_between_groups = TRUE, between_alpha
# sets how hard SMFnorm shrinks the between-group correction (verified
# empirically: LOWER alpha = STRONGER correction):
#   - low alpha  (e.g. 0.1): strongest correction -> groups pulled to equal,
#     removing the artifact AND the spike signal (no detection at the extreme).
#   - high alpha (e.g. 0.9): weakest correction -> more spike signal retained
#     (higher sensitivity) but more residual artifact (more false positives).
# The goal is to find the alpha that best trades sensitivity (spike-in) against
# specificity (null false positives) -- the highest alpha that keeps null FP low.
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
config$data_dir      <- "../m-series-data/data/allc"
config$sample_sheet  <- "../m-series-data/data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "results/alpha_sweep"

# between_alpha values to test. Start coarse; extend toward 0.9 if the optimum
# looks to be at the high end.
ALPHA_GRID <- seq(0.5, 0.9, 0.1)

# FAST_TUNE: run on chr1 only with fewer spike-in regions so the whole sweep
# finishes in ~1-2 h instead of ~8 h/alpha. Use it to find the rough optimum,
# then re-confirm the chosen alpha on the full genome with FAST_TUNE <- FALSE.
FAST_TUNE <- F
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

# Fast-tuning: restrict to chr1 and fewer spike-in regions (~10x faster). The
# optimal between_alpha transfers to the full genome; confirm it there after.
# Separate output dir so fast and full results never collide.
if (FAST_TUNE) {
    config$chr_pattern       <- "^(chr)?1$"   # chr1 only
    config$n_spikein_regions <- 100
    base_output              <- "results/alpha_sweep_fast"
    message(sprintf("FAST_TUNE on: chr1 only, %d regions, output -> %s",
                    config$n_spikein_regions, base_output))
}

dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Load source replicates once.
wt_reps <- prepare_wt_replicates(config)

# SLURM array support: if SLURM_ARRAY_TASK_ID is set, this process runs ONLY
# that one alpha (one alpha per array task). The combined CSVs/plot are then
# produced by a final plain run (no array var), which resume-skips the finished
# alphas and stitches them together.
task_id    <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
array_mode <- nzchar(task_id)
if (array_mode) {
    idx <- suppressWarnings(as.integer(task_id))
    if (is.na(idx) || idx < 1L || idx > length(ALPHA_GRID)) {
        stop(sprintf("SLURM_ARRAY_TASK_ID=%s out of range 1..%d",
                     task_id, length(ALPHA_GRID)))
    }
    ALPHA_GRID <- ALPHA_GRID[idx]
    message(sprintf("Array mode: task %s -> between_alpha = %.2f",
                    task_id, ALPHA_GRID))
}

null_all <- list()
spk_all  <- list()
null_csv <- file.path(base_output, "alpha_sweep_null.csv")
spk_csv  <- file.path(base_output, "alpha_sweep_spikein.csv")

for (a in ALPHA_GRID) {
    message("\n", strrep("#", 64))
    message(sprintf("# between_alpha = %.2f", a))
    message(strrep("#", 64))

    cfg <- config
    cfg$between_alpha <- a
    cfg$output_dir    <- file.path(base_output, sprintf("alpha_%.2f", a))
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
    null[, between_alpha := a]
    spk[, between_alpha := a]
    null_all[[length(null_all) + 1]] <- null
    spk_all[[length(spk_all) + 1]]  <- spk

    # Checkpoint the combined CSVs after EVERY alpha, so a killed job still
    # leaves usable partial results and the next run only does what's missing.
    # Skipped in array mode (each task does one alpha) to avoid write races.
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
