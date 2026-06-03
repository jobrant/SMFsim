#!/usr/bin/env Rscript

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
#
# Updated run_spikein_simulation for 04_run_simulation.R
#
# Key changes:
# 1. Uses per-replicate efficiency vectors from get_efficiency_scenarios()
# 2. rate_between_groups = FALSE (within-group normalization only)
# 3. Dense region selection for spike-in placement


# Configuration -----------------------------------------------------------

#' Parse command line arguments or use defaults
#' @export
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
    region_width_bp   = 2000,
    max_median_spacing = 200,
    min_sites = 10, 
    wt_group_id       = "PrEC",   # normal prostate epithelial cells (3 reps)
    seed              = 42,
    methods       = c("raw", "downsampled", "SMFnorm", "ComBatMet"),
    metilene_max_dist = 300,
    metilene_min_cpg  = 10,
    metilene_min_diff = 0.1,
    metilene_qval     = 0.05,   # q-value cutoff for DMR significance
    metilene_q_source = "metilene",  # "metilene" (col-4 q) or "BH" (needs min_diff=0)
    within_alpha = 0.3,
    between_alpha = 0.5,
    min_coverage = 5,
    dispersion_s = Inf,    # Beta precision for overdispersion; Inf = pure binomial
    sim_mode = "clone",    # pseudo-group construction: "clone" or "parametric"
    rate_between_groups = FALSE  # SMFnorm: correct between-group rate differences
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
    } else if (key %in% c("metilene_min_diff", "dispersion_s", "metilene_qval")) {
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
#' Non-standard contigs (e.g. `*_random`, `Un_*`) break the numeric site id used
#' by `find_shared_sites()` and can leave replicates with mismatched row counts.
#' By default this filters to standard chromosomes before shared-site filtering
#' (override the pattern with `config$chr_pattern`, or disable with
#' `config$standard_chr_only = FALSE`), then verifies all replicates are aligned.
#'
#' @param config Configuration list from parse_args().
#' @return Named list of WT replicate data.tables (aligned, post shared-site
#'   filtering).
#' @export
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

  # Filter to standard chromosomes (autosomes + X/Y/M, with/without "chr")
  if (config$standard_chr_only %||% TRUE) {
    pattern <- config$chr_pattern %||% "^(chr)?([0-9]{1,2}|X|Y|M|MT)$"
    n_before <- nrow(all_samples[[1]])
    all_samples <- lapply(all_samples, function(dt) dt[grepl(pattern, chr)])
    message(sprintf("Standard-chromosome filter: %d -> %d rows in sample 1",
                    n_before, nrow(all_samples[[1]])))
  }

  # Find shared sites across WT replicates
  all_samples <- find_shared_sites(all_samples, filter = TRUE)

  # Verify alignment: downstream code (.pooled_rate, run_downsampled) assumes
  # all replicates share the same sites in the same row order.
  lens <- vapply(all_samples, nrow, integer(1))
  if (length(unique(lens)) != 1L) {
    stop(sprintf(
      "Replicates are not aligned after find_shared_sites (rows: %s). Check for non-standard contigs or duplicate sites.",
      paste(lens, collapse = ", ")))
  }

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
#' @export
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
      mode = config$sim_mode %||% "clone",
      seed = config$seed,
      dispersion_s = config$dispersion_s
    )

    # Run all normalization methods
    method_results <- run_all_methods(
      pseudo_groups = pseudo,
      methods = config$methods,
      SMFnorm_params = list(
        within_alpha = config$within_alpha,
        between_alpha = config$between_alpha,
        min_coverage = config$min_coverage %||% 5,
        rate_between_groups = config$rate_between_groups %||% FALSE
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
          metilene_max_dist = config$metilene_max_dist
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
#' @export
run_spikein_simulation <- function(wt_reps, config) {
  message("\n", strrep("=", 60))
  message("LAYER 2: SPIKE-IN SIMULATION (Sensitivity Assessment)")
  message(strrep("=", 60))
  
  scenarios <- get_efficiency_scenarios()
  scenarios <- scenarios[config$scenarios]

  sim_mode <- config$sim_mode %||% "clone"
  message(sprintf("Spike-in construction mode: %s", sim_mode))

  # Select spike-in regions once (same regions for all scenarios)
  ref_dt <- wt_reps[[1]]
  spikein_regions <- select_dense_spikein_regions(
    dt = ref_dt,
    n_regions = config$n_spikein_regions %||% 300,
    region_width_bp = config$region_width_bp %||% 2000,
    min_sites = config$min_sites %||% 10,
    max_median_spacing = config$max_median_spacing %||% 200,
    seed = config$seed + 1000
  )
  
  # Fix data.table internal reference
  spikein_regions <- spikein_regions[]
  
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
      
      n_reps <- length(wt_reps)
      rep_names <- names(wt_reps)
      
      # Resolve per-replicate efficiencies. Seed first so range-based scenarios
      # (which sample via runif) give the same efficiencies across effect sizes
      # and are reproducible across runs.
      set.seed(config$seed)
      eff_A <- .resolve_efficiencies(sc$efficiency_A, n_reps, "A")
      eff_B <- .resolve_efficiencies(sc$efficiency_B, n_reps, "B")
      
      if (sim_mode == "parametric") {
        # Parametric: inject the effect into the pooled TRUE rate for group B,
        # then draw both groups independently from their respective rates.
        pooled_p <- .pooled_rate(wt_reps)
        inj <- inject_spikein_parametric(
          ref_dt = wt_reps[[1]],
          p = pooled_p,
          regions = spikein_regions,
          effect_size = effect,
          direction = "auto",
          seed = config$seed + 2000
        )

        distorted_A <- mapply(function(dt, eff, i) {
          .draw_parametric_rep(dt, pooled_p, eff,
                               config$seed + 4000 + i, config$dispersion_s)
        }, wt_reps, eff_A, seq_along(wt_reps), SIMPLIFY = FALSE)
        names(distorted_A) <- paste0("PseudoA_", rep_names)

        distorted_B <- mapply(function(dt, eff, i) {
          .draw_parametric_rep(dt, inj$p, eff,
                               config$seed + 5000 + i, config$dispersion_s)
        }, wt_reps, eff_B, seq_along(wt_reps), SIMPLIFY = FALSE)
        names(distorted_B) <- paste0("PseudoB_", rep_names)

      } else {
        # Clone: inject spike-ins into realized group-B counts BEFORE efficiency,
        # then thin/resample both groups from the same source replicates.
        reps_B_injected <- lapply(seq_len(n_reps), function(i) {
          injection <- inject_spikein_dmrs(
            dt = wt_reps[[i]],
            regions = spikein_regions,
            effect_size = effect,
            direction = "auto",
            seed = config$seed + 2000 + i
          )
          injection$data
        })
        names(reps_B_injected) <- rep_names

        distorted_A <- mapply(function(dt, eff, i) {
          simulate_efficiency(dt, eff, seed = config$seed + 4000 + i,
                              dispersion_s = config$dispersion_s)
        }, wt_reps, eff_A, seq_along(wt_reps), SIMPLIFY = FALSE)
        names(distorted_A) <- paste0("PseudoA_", rep_names)

        distorted_B <- mapply(function(dt, eff, i) {
          simulate_efficiency(dt, eff, seed = config$seed + 5000 + i,
                              dispersion_s = config$dispersion_s)
        }, reps_B_injected, eff_B, seq_along(reps_B_injected), SIMPLIFY = FALSE)
        names(distorted_B) <- paste0("PseudoB_", rep_names)
      }
      
      pseudo <- list(
        PseudoA = distorted_A,
        PseudoB = distorted_B,
        params = list(
          efficiency_A = eff_A,
          efficiency_B = eff_B,
          effect_size = effect,
          dispersion_s = config$dispersion_s,
          sim_mode = sim_mode,
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
          min_coverage = config$min_coverage,
          rate_between_groups = config$rate_between_groups %||% FALSE
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
            metilene_max_dist = config$metilene_max_dist,
            q_cutoff = config$metilene_qval %||% 0.05,
            q_source = config$metilene_q_source %||% "metilene"
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
        
        message(sprintf(
          "    %s: calls TP=%d FP=%d | regions det=%d/%d | sens=%.3f prec=%.3f F1=%.3f",
          m_name, metrics$TP, metrics$FP,
          metrics$regions_detected, metrics$n_regions,
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

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x # now in R >= 4.4


# Main entry point --------------------------------------------------------

#' Run the full simulation benchmark
#' 
#' @export 
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
