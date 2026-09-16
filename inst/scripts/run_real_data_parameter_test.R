#!/usr/bin/env Rscript

# run_real_data_parameter_test.R
#
# Test all four methods on a real two-group comparison and compare SMFnorm with
# within-group correction only vs both within + between-group correction.
#
# Usage (defaults to M1 vs M2, between_alpha 0.8):
#   Rscript inst/scripts/run_real_data_parameter_test.R
#
# Parameterized for any pair (positional groups, optional --key value):
#   Rscript inst/scripts/run_real_data_parameter_test.R M3 M4
#   Rscript inst/scripts/run_real_data_parameter_test.R M3 M4 --between_alpha 0.9
#   Rscript inst/scripts/run_real_data_parameter_test.R --group_A M3 --group_B M4 \
#       --between_alpha 0.9 --output_dir results/real_data_M3_M4

# Locate the package root whether this is launched from the package directory
# itself or from one level up (e.g. the parent "simulation-testing" folder).
pkg_dir <- if (file.exists("DESCRIPTION")) "." else "SMFsim"
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  devtools::load_all(pkg_dir)
} else {
  library(SMFsim)
}
suppressWarnings(suppressMessages(library(data.table)))

# Parse command-line arguments: first two bare args are group_A/group_B; any
# --key value pairs override (group_A, group_B, between_alpha, output_dir).
cli <- commandArgs(trailingOnly = TRUE)
positional <- character(0)
opts <- list()
i <- 1
while (i <= length(cli)) {
  if (startsWith(cli[i], "--")) {
    opts[[sub("^--", "", cli[i])]] <- cli[i + 1]
    i <- i + 2
  } else {
    positional <- c(positional, cli[i])
    i <- i + 1
  }
}

group_A <- opts$group_A %||% (if (length(positional) >= 1) positional[1] else "M1")
group_B <- opts$group_B %||% (if (length(positional) >= 2) positional[2] else "M2")
between_alpha <- as.numeric(opts$between_alpha %||% "0.8")
output_dir <- opts$output_dir %||% sprintf("results/real_data_%s_%s", group_A, group_B)

config <- list(
  data_dir        = "data/allc",
  sample_sheet    = "data/sample_sheet.csv",
  output_dir      = output_dir,
  metilene_path   = "/apps/metilene/0.2.8/metilene",
  metilene_min_cpg = 10,
  metilene_min_diff = 0.1,
  min_coverage    = 10,
  seed            = 42
)

message(sprintf("Real-data comparison: %s vs %s  (between_alpha = %.2f, output = %s)",
                group_A, group_B, between_alpha, output_dir))

# min_coverage is applied inside load_real_data() -> SMFnorm::load_data() at
# read time (per file, before find_shared_sites), so no manual post-filtering
# is needed here.
real_data <- load_real_data(config, group_A, group_B)

sample_rows <- vapply(c(real_data$group_A, real_data$group_B), nrow, integer(1))
if (length(unique(sample_rows)) != 1L) {
  stop("Shared-site filtering failed: sample row counts differ after find_shared_sites()")
}
message(sprintf("min_coverage = %d, shared sites after loading: %d",
                config$min_coverage, sample_rows[1]))

# Two correction regimes to compare. within-group correction is always on; the
# modes differ only in whether between-group correction is added on top.
modes <- list(
  within_only = list(
    within_alpha = 0.3,
    between_alpha = between_alpha,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE
  ),
  both = list(
    within_alpha = 0.3,
    between_alpha = between_alpha,
    rate_within_groups = TRUE,
    rate_between_groups = TRUE
  )
)

results_summary <- list()

for (mode_name in names(modes)) {
  mode_dir <- file.path(config$output_dir, mode_name)
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n", strrep("=", 70))
  message(sprintf("Running real-data comparison: %s", mode_name))
  message(strrep("=", 70))

  SMFnorm_params <- list(
    within_alpha = modes[[mode_name]]$within_alpha,
    between_alpha = modes[[mode_name]]$between_alpha,
    min_coverage = config$min_coverage,
    rate_within_groups = modes[[mode_name]]$rate_within_groups,
    rate_between_groups = modes[[mode_name]]$rate_between_groups
  )

  method_results <- run_all_methods_real(
    real_data = real_data,
    methods = c("raw", "downsampled", "SMFnorm", "ComBatMet"),
    SMFnorm_params = SMFnorm_params
  )

  dmr_list <- list()
  for (method_name in names(method_results)) {
    message(sprintf("\nCalling DMRs for method: %s", method_name))
    dmr_dir <- file.path(mode_dir, "dmrs", method_name)
    dmrs <- tryCatch({
      call_dmrs_real(
        split_data = method_results[[method_name]]$data,
        out_dir = dmr_dir,
        metilene_path = config$metilene_path,
        min_cpg = config$metilene_min_cpg,
        min_diff = config$metilene_min_diff
      )
    }, error = function(e) {
      warning(sprintf("DMR calling failed for %s: %s", method_name, e$message))
      data.table(chr = character(), start = integer(), end = integer(),
                 q_value = numeric(), mean_diff = numeric(), n_sites = integer())
    })
    dmr_list[[method_name]] <- dmrs
  }

  summary_dt <- data.table(
    method = names(dmr_list),
    n_dmrs = sapply(dmr_list, nrow),
    median_effect = sapply(dmr_list, function(dt) if (nrow(dt) > 0) round(median(abs(dt$mean_diff)), 4) else NA_real_),
    mean_effect = sapply(dmr_list, function(dt) if (nrow(dt) > 0) round(mean(abs(dt$mean_diff)), 4) else NA_real_)
  )
  summary_dt[, mode := mode_name]
  fwrite(summary_dt, file.path(mode_dir, "dmr_summary.csv"))

  overlap_dt <- compute_dmr_overlap(dmr_list)
  fwrite(overlap_dt, file.path(mode_dir, "dmr_overlap.csv"))

  classification <- classify_method_dmrs(dmr_list, reference_method = "raw")

  fig_dir <- file.path(mode_dir, "figures")
  plot_dmr_counts(dmr_list, fig_dir)
  plot_effect_size_distribution(dmr_list, fig_dir)
  plot_overlap_heatmap(overlap_dt, fig_dir)
  plot_lost_dmr_effects(dmr_list, classification, fig_dir)

  saveRDS(list(
    dmr_list = dmr_list,
    summary = summary_dt,
    overlap = overlap_dt,
    classification = classification,
    config = list(config = config, mode = mode_name, SMFnorm_params = SMFnorm_params)
  ), file.path(mode_dir, "real_data_results.rds"))

  results_summary[[mode_name]] <- summary_dt
}

combined_summary <- rbindlist(results_summary)
fwrite(combined_summary, file.path(config$output_dir, "dmr_summary_by_mode.csv"))

message(sprintf("\nCompleted real-data parameter test for %s vs %s", group_A, group_B))
message("Results written to: ", normalizePath(config$output_dir, winslash = "/"))
message("\nNote: 'both' adds between-group correction on top of within-group, and is")
message("only appropriate if there is a true systematic technical bias aligned with the group labels.")
