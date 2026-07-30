#!/usr/bin/env Rscript

# run_alpha_sweep.R
#
# Sweep SMFnorm's `within_alpha` AND `between_alpha` together, as a 2-D grid, on
# the systematic-bias scenarios.
#
# WHY A GRID. The earlier design swept each alpha with the other held fixed
# (between at within = 0.3; within at between = 0.8). Two problems with that:
#   - each 1-D slice is only valid near the other alpha's true optimum, and the
#     held-fixed values have since changed;
#   - the old within sweep ran with rate_between_groups = FALSE, so within_alpha
#     was never tuned while between-correction was active -- which is exactly
#     the regime the headline bias experiment uses.
# This grid tunes both in the regime they are actually used in, and reports
# whether the optimum is SEPARABLE (best within_alpha the same at every
# between_alpha). If it is, you can go back to quoting two 1-D sweeps.
#
# Direction of the alphas (verified empirically): LOWER alpha = STRONGER
# correction. between_alpha -> 0.1 pulls the groups to equal, removing the
# artifact AND the spike signal; -> 0.9 retains more spike signal but leaves
# more residual artifact. The goal is the highest alpha that keeps null FP low.
#
# Only SMFnorm is run (no other method depends on alpha); compare against the
# raw/downsampled baselines from the main bias run.
#
# Q-VALUE SETTINGS ARE INHERITED FROM parse_args() ON PURPOSE. The sweep must
# use exactly the same significance path as run_bias_experiment.R, so neither
# script pins them -- parse_args() is the single source of truth. Pinning
# metilene_q_source here previously silently held the sweep on the old
# (Bonferroni / wrong-p-value) path.
#
# Usage (from the package root):
#   Rscript inst/scripts/run_alpha_sweep.R
# SLURM array (one CELL per task, 1..nrow(GRID)):
#   sbatch --array=1-15 ... ; then one plain run to stitch + plot.

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

# The 2-D grid. Keep it coarse: the point is to locate the optimum and test
# separability, not to resolve it to two decimal places.
WITHIN_GRID  <- seq(0.1, 0.9, 0.2)   # 0.1 0.3 0.5 0.7 0.9
BETWEEN_GRID <- c(0.5, 0.7, 0.9)

# FAST_TUNE: chr1 only with fewer spike-in regions (~10x faster per cell).
# Strongly recommended for a 2-D grid -- 15 cells genome-wide is a very long
# job. Find the optimum here, then confirm the single chosen pair with
# FAST_TUNE <- FALSE and a 1x1 grid.
FAST_TUNE <- TRUE
# =========================================================================

# Fixed settings — match the main bias run.
config$sim_mode            <- "parametric"
config$dispersion_s        <- 26
config$rate_between_groups <- TRUE
config$methods             <- c("SMFnorm")   # only alpha-dependent method
config$scenarios           <- c("aligned_moderate", "aligned_strong",
                                "imbalanced_moderate", "imbalanced_strong")
config$effect_sizes        <- c(0.10, 0.15, 0.20, 0.30)
config$seed                <- 42

GRID <- CJ(within_alpha = WITHIN_GRID, between_alpha = BETWEEN_GRID)
n_cells <- nrow(GRID)

if (FAST_TUNE) {
    config$chr_pattern       <- "^(chr)?1$"   # chr1 only
    config$n_spikein_regions <- 100
    base_output              <- "results/alpha_sweep_fast"
    message(sprintf("FAST_TUNE on: chr1 only, %d regions, output -> %s",
                    config$n_spikein_regions, base_output))
} else if (n_cells > 4) {
    message(sprintf(
        "!! FAST_TUNE is off with %d grid cells. Genome-wide this is roughly ",
        n_cells))
    message("   8 h per cell. Consider FAST_TUNE <- TRUE, or shrink the grid.")
}

message(sprintf("Grid: %d within x %d between = %d cells",
                length(WITHIN_GRID), length(BETWEEN_GRID), n_cells))

# --- Guard: never silently overwrite a previous sweep ---------------------
# Same rationale as run_bias_experiment.R: base_output is a fixed path, so a
# re-run lands on top of whatever is there. Per-cell resume (below) is the
# intended way to continue an interrupted sweep; this guard only fires when the
# directory holds something that is not a resumable sweep.
overwrite_ok <- identical(Sys.getenv("SMFSIM_OVERWRITE"), "1")
if (dir.exists(base_output) && !overwrite_ok) {
    stray <- setdiff(list.files(base_output),
                     grep("^(w[0-9.]+_b[0-9.]+|alpha_sweep_.*)$",
                          list.files(base_output), value = TRUE))
    if (length(stray) > 0) {
        stop(sprintf(paste0(
            "Output directory holds files that are not from this sweep:\n  %s\n",
            "  (%s)\nArchive it, point `base_output` elsewhere, or re-run with ",
            "SMFSIM_OVERWRITE=1."),
            normalizePath(base_output, mustWork = FALSE),
            paste(utils::head(stray, 5), collapse = ", ")), call. = FALSE)
    }
}
dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Load source replicates once.
wt_reps <- prepare_wt_replicates(config)

# SLURM array support: one CELL per task, indexing rows of GRID.
task_id    <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
array_mode <- nzchar(task_id)
if (array_mode) {
    idx <- suppressWarnings(as.integer(task_id))
    if (is.na(idx) || idx < 1L || idx > n_cells) {
        stop(sprintf("SLURM_ARRAY_TASK_ID=%s out of range 1..%d",
                     task_id, n_cells))
    }
    GRID <- GRID[idx]
    message(sprintf("Array mode: task %s -> within=%.2f between=%.2f",
                    task_id, GRID$within_alpha, GRID$between_alpha))
}

null_all <- list()
spk_all  <- list()
null_csv <- file.path(base_output, "alpha_sweep_null.csv")
spk_csv  <- file.path(base_output, "alpha_sweep_spikein.csv")

for (i in seq_len(nrow(GRID))) {
    w <- GRID$within_alpha[i]
    b <- GRID$between_alpha[i]

    message("\n", strrep("#", 64))
    message(sprintf("# cell %d/%d: within_alpha = %.2f, between_alpha = %.2f",
                    i, nrow(GRID), w, b))
    message(strrep("#", 64))

    cfg <- config
    cfg$within_alpha  <- w
    cfg$between_alpha <- b
    cfg$output_dir    <- file.path(base_output, sprintf("w%.2f_b%.2f", w, b))
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

    c_null <- file.path(cfg$output_dir, "null_simulation_results.csv")
    c_spk  <- file.path(cfg$output_dir, "spikein_simulation_results.csv")

    # Resume: skip a cell whose per-cell CSVs already exist (completed).
    if (file.exists(c_null) && file.exists(c_spk)) {
        message(sprintf("  [resume] cell w%.2f_b%.2f already done; reading cached",
                        w, b))
        null <- fread(c_null)
        spk  <- fread(c_spk)
    } else {
        null <- run_null_simulation(wt_reps, cfg)
        spk  <- run_spikein_simulation(wt_reps, cfg)
    }
    null[, `:=`(within_alpha = w, between_alpha = b)]
    spk[,  `:=`(within_alpha = w, between_alpha = b)]
    null_all[[length(null_all) + 1]] <- null
    spk_all[[length(spk_all) + 1]]  <- spk

    # Checkpoint after EVERY cell, so a killed job leaves usable partials.
    # Skipped in array mode (one cell per task) to avoid write races.
    if (!array_mode) {
        fwrite(rbindlist(null_all, fill = TRUE), null_csv)
        fwrite(rbindlist(spk_all,  fill = TRUE), spk_csv)
        message(sprintf("  [checkpoint] combined CSVs updated through cell %d/%d",
                        i, nrow(GRID)))
    }
}

if (array_mode) {
    message(sprintf(
        "\nArray task done (within=%.2f between=%.2f). Run once WITHOUT ",
        GRID$within_alpha, GRID$between_alpha))
    message("SLURM_ARRAY_TASK_ID to stitch the combined CSVs + plots.")
    quit(save = "no")
}

null_dt <- rbindlist(null_all, fill = TRUE)
spk_dt  <- rbindlist(spk_all,  fill = TRUE)

message("\nDONE. Combined results:")
message("  null    -> ", null_csv)
message("  spikein -> ", spk_csv)

# --- Is the optimum separable? -------------------------------------------
# The question the grid exists to answer: does the best within_alpha depend on
# between_alpha? If not, the two can be tuned independently and reported as two
# 1-D sweeps. Specificity comes first -- a cell is only eligible if its null
# false-positive count is acceptable.
NULL_FP_TOL <- 0

message("\n=== Separability of the two alphas ===")
tryCatch({
    fp <- null_dt[, .(null_FP = sum(FP, na.rm = TRUE)),
                  by = .(within_alpha, between_alpha)]
    perf <- spk_dt[, .(mean_F1 = mean(F1, na.rm = TRUE),
                       mean_sens = mean(sensitivity, na.rm = TRUE),
                       min_prec = min(precision, na.rm = TRUE)),
                   by = .(within_alpha, between_alpha)]
    cells <- merge(perf, fp, by = c("within_alpha", "between_alpha"), all = TRUE)
    cells[is.na(null_FP), null_FP := 0L]
    setorder(cells, between_alpha, within_alpha)

    message("\nPer-cell summary (null_FP first -- specificity gates eligibility):")
    print(cells)

    eligible <- cells[null_FP <= NULL_FP_TOL]
    if (nrow(eligible) == 0) {
        message("\n!! No cell met null_FP <= ", NULL_FP_TOL,
                ". Relax NULL_FP_TOL or inspect the null directly.")
    } else {
        best_per_b <- eligible[, .SD[which.max(mean_F1)], by = between_alpha]
        message("\nBest within_alpha at each between_alpha:")
        print(best_per_b[, .(between_alpha, best_within = within_alpha,
                             mean_F1, mean_sens, null_FP)])

        if (length(unique(best_per_b$within_alpha)) == 1L) {
            message(sprintf(
              "\n==> SEPARABLE: best within_alpha = %.2f at every between_alpha.",
              best_per_b$within_alpha[1]))
            message("    The two can be tuned and reported independently.")
        } else {
            message("\n==> NOT SEPARABLE: the best within_alpha changes with ",
                    "between_alpha.")
            message("    Report the joint optimum; do not quote two 1-D sweeps.")
        }

        overall <- eligible[which.max(mean_F1)]
        message(sprintf(
            "\nJoint optimum: within = %.2f, between = %.2f  (mean F1 %.3f, ",
            overall$within_alpha, overall$between_alpha, overall$mean_F1))
        message(sprintf("  mean sensitivity %.3f, min precision %.3f, null FP %d)",
                        overall$mean_sens, overall$min_prec, overall$null_FP))
    }
}, error = function(e) message("  [separability check skipped] ",
                               conditionMessage(e)))

# --- Diagnostic plots ----------------------------------------------------
tryCatch({
    suppressWarnings(suppressMessages(library(ggplot2)))

    # F1 across the grid, per scenario x effect size.
    p1 <- ggplot(spk_dt, aes(factor(between_alpha), factor(within_alpha),
                             fill = F1)) +
        geom_tile(colour = "white") +
        facet_grid(effect_size ~ efficiency_scenario) +
        scale_fill_gradient2(low = "#E41A1C", mid = "#FFFFBF", high = "#4DAF4A",
                             midpoint = 0.5, limits = c(0, 1),
                             na.value = "grey90") +
        labs(title = "SMFnorm F1 across the alpha grid",
             x = "between_alpha", y = "within_alpha") +
        theme_minimal()
    ggsave(file.path(base_output, "alpha_grid_F1.png"), p1,
           width = 11, height = 8, dpi = 150)

    # Null false positives across the grid -- the specificity gate.
    fp_dt <- null_dt[, .(null_FP = sum(FP, na.rm = TRUE)),
                     by = .(within_alpha, between_alpha, scenario)]
    p2 <- ggplot(fp_dt, aes(factor(between_alpha), factor(within_alpha),
                            fill = null_FP)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = null_FP), size = 3) +
        facet_wrap(~ scenario) +
        scale_fill_gradient(low = "#FFFFBF", high = "#E41A1C") +
        labs(title = "Null false positives across the alpha grid",
             x = "between_alpha", y = "within_alpha") +
        theme_minimal()
    ggsave(file.path(base_output, "alpha_grid_null_FP.png"), p2,
           width = 10, height = 7, dpi = 150)

    message("  plots   -> ", file.path(base_output, "alpha_grid_*.png"))
}, error = function(e) message("  [plots skipped] ", conditionMessage(e)))
