setwd("/orange/paefron/jobrant/NOME-SEQ/mapitnorm_testing/simulation")

devtools::install_github("jobrant/SMFnorm")
devtools::install_github("jobrant/SMFsim")


x <- c("ggplot2", "RColorBrewer", "data.table", "sva", "GenomicRanges", 
       "TxDb.Hsapiens.UCSC.hg38.knownGene", "IRanges", "S4Vectors", 
       "SMFnorm", "SMFsim")
lapply(x, require, character.only = TRUE)

#rm(list = ls())

# Set up config
config <- list(
  data_dir            = "data/allc",
  sample_sheet        = "data/sample_sheet.csv",
  output_dir          = "results/simulation",
  metilene_path       = "/apps/metilene/0.2.8/metilene",
  scenarios           = c("mild", "moderate", "severe"),
  effect_sizes        = c(0.10, 0.15, 0.20, 0.30),
  n_spikein_regions   = 300,
  region_width_bp     = 2000, 
  max_median_distance = 200, 
  min_sites           = 10,
  wt_group_id         = "PrEC",   # normal prostate epithelial cells (3 reps)
  seed                = 42,
  methods             = c("raw", "downsampled", "SMFnorm", "ComBatMet"),
  metilene_max_dist   = 300,
  metilene_min_cpg    = 10,
  metilene_min_diff   = 0.1,
  within_alpha        = 0.3, 
  between_alpha       = 0.5, 
  min_coverage        = 5
)


# Load PrEC replicates
wt_reps <- prepare_wt_replicates(config = config)

# Run null simulation (false positive assessment)
null_results <- run_null_simulation(wt_reps = wt_reps, config = config)

# Run Spike-in simulation (sensitivity assessment)
spikein_results <- run_spikein_simulation(wt_reps = wt_reps, config = config)

# Save primary results
saveRDS(list(null = null_results, spikein = spikein_results, config = config),
        file.path(config$output_dir, "primary_results.rds"))


# Alpha sweep (SMFnorm only, spike-in only)
config$methods <- c("SMFnorm")
alphas <- seq(0.1, 0.9, by = 0.1)
alpha_results <- list()

for (a in alphas) {
  message(sprintf("\n=== Alpha = %.1f ===", a))
  config$within_alpha <- a
  result <- run_spikein_simulation(wt_reps, config)
  alpha_results[[as.character(a)]] <- result
}

# save alpha sweep results object
saveRDS(alpha_results, file.path(config$output_dir, "alpha_sweep_results.rds"))

# summary table of alpha sweep results
for (a in names(alpha_results)) {
  dt <- alpha_results[[a]]
  row <- dt[efficiency_scenario == "severe" & effect_size == 0.3 & method == "SMFnorm"]
  cat(sprintf("alpha=%-3s  TP=%3d  FP=%3d  sensitivity=%.3f  precision=%s  F1=%s\n",
              a, row$TP, row$FP, row$sensitivity,
              ifelse(is.na(row$precision), "NA", sprintf("%.3f", row$precision)),
              ifelse(is.na(row$F1), "NA", sprintf("%.3f", row$F1))))
}


# Generate figures
generate_all_figures(results_file = "results/simulation/primary_results.rds", out_dir = "results/figures")

generate_additional_figures("results/simulation/spikein_simulation_results.csv", "results/figures/")


# Real world data test ----------------------------------------------------
# make sure to reset alpha to 0.3 if running after sweep
config$within_alpha = 0.3
config$metilene_min_cpg = 5

results <- run_real_data_comparison(
  config,
  group_A = "PrEC",
  group_B = "PC",
  out_dir = "results/real_data"
)


