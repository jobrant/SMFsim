#!/usr/bin/env Rscript
# =============================================================================
# 04_run_simulation.R
# Main orchestration script for the full simulation benchmark
#
# Usage:
#   Rscript 04_run_simulation.R \
#     --data_dir path/to/allc_files \
#     --sample_sheet path/to/sample_sheet.csv \
#     --output_dir results/simulation \
#     --metilene_path metilene \
#     --scenarios mild,moderate,severe \
#     --seed 42
# =============================================================================

library(data.table)
library(SMFnorm)

# Source simulation modules
source("SMFsim/01_simulate_efficiency.R")
source("SMFsim/02_run_methods.R")
source("SMFsim/03_evaluate.R")

# Configuration -----------------------------------------------------------

#' Parse command line arguments or use defaults
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults — adjust these to your local setup
  config <- list(
    data_dir      = "data/allc",
    sample_sheet  = "data/sample_sheet.csv",
    output_dir    = "results/simulation",
    metilene_path = "/apps/metilene/0.2.8/metilene",
    scenarios     = c("mild", "moderate", "severe"),
    effect_sizes  = c(0.10, 0.15, 0.20, 0.30),
    n_spikein_regions = 300,
    region_width_bp   = 500,
    wt_group_id       = "PrEC",   # normal prostate epithelial cells (3 reps)
    seed              = 42,
    methods       = c("raw", "downsampled", "SMFsim", "ComBatMet"),
    metilene_max_dist = 500,
    metilene_min_cpg  = 5,
    metilene_min_diff = 0.1,
    within_alpha = 0.03, 
    between_alpha = 0.5, 
    min_coverage = 5
  )

  # Override from command line if provided
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key %in% c("scenarios", "methods")) {
      config[[key]] <- strsplit(args[i + 1], ",")[[1]]
    } else if (key %in% c("seed", "n_spikein_regions", "region_width_bp",
                           "metilene_min_cpg", "metilene_max_dist")) {
      config[[key]] <- as.integer(args[i + 1])
    } else if (key %in% c("metilene_min_diff")) {
      config[[key]] <- as.numeric(args[i + 1])
    } else if (key %in% c("effect_sizes")) {
      config[[key]] <- as.numeric(strsplit(args[i + 1], ",")[[1]])
    } else if (key %in% names(config)) {
      config[[key]] <- args[i + 1]
    }
    i <- i + 2
  }

  return(config)
}


# Data preparation --------------------------------------------------------

#' Load WT replicates and prepare for simulation
#'
#' @param config Configuration list from parse_args().
#' @return Named list of WT replicate data.tables (post shared-site filtering).
prepare_wt_replicates <- function(config) {
  message("\n", strrep("=", 60))
  message("STEP 0: Loading and preparing WT replicates")
  message(strrep("=", 60))

  # Load all data
  all_samples <- load_data(
    dir_path = config$data_dir,
    sample_sheet = config$sample_sheet,
    groups = config$wt_group_id
  )

  message(sprintf("Loaded %d WT samples", length(all_samples)))

  # Find shared sites across WT replicates
  all_samples <- find_shared_sites(all_samples, filter = TRUE)

  message(sprintf("After shared-site filtering: %d sites per sample",
                  nrow(all_samples[[1]])))

  # Print summary
  # summarize_meth(all_samples)

  return(all_samples)
}


# Layer 1: Null simulation ------------------------------------------------

#' Run the null simulation (false positive assessment)
#'
#' @param wt_reps Named list of WT replicate data.tables.
#' @param config Configuration list.
#' @return data.table of null simulation results (FP counts by scenario × method).
run_null_simulation <- function(wt_reps, config) {
  message("\n", strrep("=", 60))
  message("LAYER 1: NULL SIMULATION (False Positive Assessment)")
  message(strrep("=", 60))

  scenarios <- get_efficiency_scenarios()
  scenarios <- scenarios[config$scenarios]

  null_results <- list()

  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    message(sprintf("\n--- Scenario: %s ---", sc$label))

    # Create pseudo-groups with efficiency distortion (clone mode: 3 vs 3)
    pseudo <- create_pseudo_groups(
      replicates = wt_reps,
      efficiency_A = sc$efficiency_A,
      efficiency_B = sc$efficiency_B,
      mode = "clone",
      seed = config$seed
    )

    # Run all normalization methods
    method_results <- run_all_methods(
      pseudo_groups = pseudo,
      methods = config$methods, 
      SMFnorm_params = list(
        within_alpha = config$within_alpha,
        between_alpha = config$between_alpha,
        min_coverage = config$min_coverage %||% 5
      )
    )

    # Call DMRs and classify for each method
    null_results[[sc_name]] <- list()

    for (m_name in names(method_results)) {
      m_res <- method_results[[m_name]]

      dmr_dir <- file.path(config$output_dir, "null", sc_name, m_name)

      dmrs <- tryCatch({
        call_dmrs_metilene(
          split_data = m_res$data,
          out_dir = dmr_dir,
          metilene_path = config$metilene_path,
          min_cpg = config$metilene_min_cpg,
          min_diff = config$metilene_min_diff, 
          metilene_max_dist = config$metilene_mix_dist
        )
      }, error = function(e) {
        warning(sprintf("DMR calling failed for %s/%s: %s", sc_name, m_name, e$message))
        data.table(chr = character(), start = integer(), end = integer(),
                   q_value = numeric(), mean_diff = numeric(), n_sites = integer())
      })

      # In null simulation, all DMRs are FP
      null_results[[sc_name]][[m_name]] <- classify_dmrs(
        called_dmrs = dmrs,
        truth_regions = NULL
      )

      message(sprintf("  %s: %d false positive DMRs", m_name, nrow(dmrs)))
    }
  }

  # Compile results
  null_summary <- compile_null_results(null_results)
  out_path <- file.path(config$output_dir, "null_simulation_results.csv")
  fwrite(null_summary, out_path)
  message(sprintf("\nNull simulation results saved to: %s", out_path))

  return(null_summary)
}


# Layer 2: Spike-in simulation --------------------------------------------

#' Run the spike-in simulation (sensitivity assessment)
#'
#' @param wt_reps Named list of WT replicate data.tables.
#' @param config Configuration list.
#' @return data.table of spike-in results by scenario × effect_size × method.
run_spikein_simulation <- function(wt_reps, config) {
  message("\n", strrep("=", 60))
  message("LAYER 2: SPIKE-IN SIMULATION (Sensitivity Assessment)")
  message(strrep("=", 60))

  scenarios <- get_efficiency_scenarios()
  scenarios <- scenarios[config$scenarios]

  # Select spike-in regions once (same regions for all scenarios)
  ref_dt <- wt_reps[[1]]
  spikein_regions <- select_dense_spikein_regions(
    dt = ref_dt,
    n_regions = 300,
    region_width_bp = 2000,
    min_sites = 10,
    max_median_spacing = 200,
    seed = config$seed + 1000
  )
  
  message("Regions class: ", class(spikein_regions)[1])
  message("Regions nrow: ", nrow(spikein_regions))
  message("Regions names: ", paste(names(spikein_regions), collapse = ", "))

  # Save regions for reference
  regions_path <- file.path(config$output_dir, "spikein_regions.csv")
  dir.create(dirname(regions_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(spikein_regions, regions_path)

  all_results <- list()

  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    message(sprintf("\n--- Scenario: %s ---", sc$label))

    for (effect in config$effect_sizes) {
      result_key <- paste(sc_name, effect, sep = "_effect")
      message(sprintf("\n  Effect size: %.2f", effect))

      # Step 1: Build pseudo-groups using clone mode
      # PseudoA: all reps with efficiency_A, no spike-in (original biology)
      # PseudoB: all reps with spike-in DMRs injected, then efficiency_B applied
      n_reps <- length(wt_reps)
      rep_names <- names(wt_reps)

      # Inject spike-in DMRs into all reps for group B
      reps_B_injected <- lapply(seq_len(n_reps), function(i) {
        injection <- inject_spikein_dmrs(
          dt = wt_reps[[i]],
          regions = spikein_regions,
          effect_size = effect,
          direction = "both",
          seed = config$seed + 2000 + i
        )
        injection$data
      })
      names(reps_B_injected) <- rep_names

      # Step 2: Apply efficiency distortion to both groups
      set.seed(config$seed + 3000)
      eff_A <- runif(n_reps, sc$efficiency_A[1], sc$efficiency_A[2])
      eff_B <- runif(n_reps, sc$efficiency_B[1], sc$efficiency_B[2])

      distorted_A <- mapply(function(dt, eff, i) {
        simulate_efficiency(dt, eff, seed = config$seed + 4000 + i)
      }, wt_reps, eff_A, seq_along(wt_reps), SIMPLIFY = FALSE)
      names(distorted_A) <- paste0("PseudoA_", rep_names)

      distorted_B <- mapply(function(dt, eff, i) {
        simulate_efficiency(dt, eff, seed = config$seed + 5000 + i)
      }, reps_B_injected, eff_B, seq_along(reps_B_injected), SIMPLIFY = FALSE)
      names(distorted_B) <- paste0("PseudoB_", rep_names)

      pseudo <- list(
        PseudoA = distorted_A,
        PseudoB = distorted_B,
        params = list(
          efficiency_A = eff_A,
          efficiency_B = eff_B,
          effect_size = effect,
          n_spikein_regions = nrow(spikein_regions)
        )
      )

      # Step 3: Run all normalization methods
      method_results <- run_all_methods(
        pseudo_groups = pseudo,
        methods = config$methods, 
        SMFnorm_params = list(
          within_alpha = config$within_alpha,
          between_alpha = config$between_alpha,
          min_coverage = config$min_coverage %||% 5
        )
      )

      # Step 4: Call DMRs and evaluate
      all_results[[result_key]] <- list()

      for (m_name in names(method_results)) {
        m_res <- method_results[[m_name]]

        dmr_dir <- file.path(config$output_dir, "spikein", sc_name,
                             paste0("effect_", effect), m_name)

        dmrs <- tryCatch({
          call_dmrs_metilene(
            split_data = m_res$data,
            out_dir = dmr_dir,
            metilene_path = config$metilene_path,
            min_cpg = config$metilene_min_cpg,
            min_diff = config$metilene_min_diff, 
            metilene_max_dist = config$metilene_max_dist
          )
        }, error = function(e) {
          warning(sprintf("DMR calling failed for %s/%s/effect%.2f: %s",
                          sc_name, m_name, effect, e$message))
          data.table(chr = character(), start = integer(), end = integer(),
                     q_value = numeric(), mean_diff = numeric(), n_sites = integer())
        })

        eval_result <- classify_dmrs(
          called_dmrs = dmrs,
          truth_regions = spikein_regions
        )

        all_results[[result_key]][[m_name]] <- eval_result
        metrics <- compute_metrics(eval_result$summary)

        message(sprintf("    %s: TP=%d FP=%d FN=%d | sens=%.3f prec=%.3f F1=%.3f",
                        m_name, metrics$TP, metrics$FP, metrics$FN,
                        metrics$sensitivity, metrics$precision, metrics$F1))
      }
    }
  }

  # Compile results
  spikein_summary <- compile_results(all_results)

  # Parse scenario and effect size from the result key
  spikein_summary[, c("efficiency_scenario", "effect_info") :=
                    tstrsplit(scenario, "_effect", fixed = TRUE)]
  spikein_summary[, effect_size := as.numeric(effect_info)]
  spikein_summary[, effect_info := NULL]

  out_path <- file.path(config$output_dir, "spikein_simulation_results.csv")
  fwrite(spikein_summary, out_path)
  message(sprintf("\nSpike-in simulation results saved to: %s", out_path))

  return(spikein_summary)
}


# Main entry point --------------------------------------------------------

main <- function() {
  config <- parse_args()

  message("Simulation Configuration:")
  message(sprintf("  Data directory:   %s", config$data_dir))
  message(sprintf("  Sample sheet:     %s", config$sample_sheet))
  message(sprintf("  Output directory: %s", config$output_dir))
  message(sprintf("  WT group:         %s", config$wt_group_id))
  message(sprintf("  Scenarios:        %s", paste(config$scenarios, collapse = ", ")))
  message(sprintf("  Methods:          %s", paste(config$methods, collapse = ", ")))
  message(sprintf("  Effect sizes:     %s", paste(config$effect_sizes, collapse = ", ")))
  message(sprintf("  Seed:             %d", config$seed))

  # Create output directory
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Save config for reproducibility
  saveRDS(config, file.path(config$output_dir, "config.rds"))

  # Load data
  wt_reps <- prepare_wt_replicates(config)

  # Run simulations
  null_results   <- run_null_simulation(wt_reps, config)
  spikein_results <- run_spikein_simulation(wt_reps, config)

  # Save combined results
  saveRDS(list(
    null = null_results,
    spikein = spikein_results,
    config = config
  ), file.path(config$output_dir, "all_results.rds"))

  message("\n", strrep("=", 60))
  message("SIMULATION COMPLETE")
  message(strrep("=", 60))

  message("\nNull simulation summary:")
  print(null_results)

  message("\nSpike-in simulation summary:")
  print(spikein_results)
}


# Run if called as a script
if (!interactive()) {
  main()
}
