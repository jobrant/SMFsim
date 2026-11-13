#!/usr/bin/env Rscript

# summarize_seeds.R
#
# Combine per-seed runs of run_bias_experiment.R into across-seed summaries.
#
# Expects  <base>/seed_<n>/<block>/{null_simulation_results.csv,
#                                   spikein_simulation_results.csv,
#                                   spikein_by_density.csv}
# as written by the seeds array job (inst/slurm/bias_experiment_seeds), and
# writes, per block,
#   <base>/all_seeds/<block>/<table>_by_seed.csv        every seed, stacked
#   <base>/all_seeds/<block>/<table>_across_seeds.csv   mean, sd, min, max, n_seeds
#
# Usage, from the manuscript root:
#   Rscript SMFsim/inst/scripts/summarize_seeds.R [base_dir]
#
# base_dir defaults to SMFSIM_SEEDS_OUTPUT or results/bias_experiment_seeds.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else
    Sys.getenv("SMFSIM_SEEDS_OUTPUT", "results/bias_experiment_seeds")
if (!dir.exists(base)) stop("No such directory: ", base, call. = FALSE)

seed_dirs <- list.dirs(base, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("^seed_-?[0-9]+$", basename(seed_dirs))]
if (!length(seed_dirs)) stop("No seed_<n> directories under ", base, call. = FALSE)

blocks <- unique(unlist(lapply(seed_dirs, function(d) {
    b <- list.dirs(d, recursive = FALSE, full.names = FALSE)
    b[file.exists(file.path(d, b, "null_simulation_results.csv"))]
})))
if (!length(blocks)) stop("No block results found under ", base, call. = FALSE)

for (b in blocks) {
    message("\n== Block: ", b)
    res <- combine_seed_dirs(base, block = b)
    for (nm in names(res$across_seeds)) {
        message("\n-- ", nm, " (mean across seeds)")
        print(res$across_seeds[[nm]])
    }
}
