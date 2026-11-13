# 08_seeds.R
# Repeat the simulation over several seeds and summarise across them.
#
# Every result used to come from one seed (42), so a reported count was a
# single realization with no measure of its run-to-run variation. A seed
# changes the replicate noise, the efficiencies of range-based scenarios, the
# spike-in regions and which group carries each effect.


#' Summarise result columns across seeds
#'
#' @param dt data.table with a `seed` column.
#' @param by Columns defining one result (e.g. scenario and method).
#' @param value_cols Numeric columns to summarise.
#' @return data.table with one row per `by` group: `n_seeds`, and for each
#'   value column its mean (under the original column name) plus `<col>_sd`,
#'   `<col>_min` and `<col>_max`. `n_seeds` counts the seeds with a row for that
#'   group, so a method that was inapplicable under some seeds shows fewer.
#' @keywords internal
summarize_across_seeds <- function(dt, by, value_cols) {
  stopifnot("seed" %in% names(dt), all(by %in% names(dt)))
  value_cols <- intersect(value_cols, names(dt))
  dt[, {
    out <- list(n_seeds = uniqueN(seed))
    for (v in value_cols) {
      x <- as.numeric(get(v))
      ok <- x[!is.na(x)]
      out[[v]] <- if (length(ok)) mean(ok) else NA_real_
      out[[paste0(v, "_sd")]]  <- if (length(ok) > 1L) stats::sd(ok) else NA_real_
      out[[paste0(v, "_min")]] <- if (length(ok)) min(ok) else NA_real_
      out[[paste0(v, "_max")]] <- if (length(ok)) max(ok) else NA_real_
    }
    out
  }, by = by]
}


#' Combine per-seed result directories
#'
#' Reads the result tables written by [run_null_simulation()] and
#' [run_spikein_simulation()] from `<base_dir>/seed_<n>/` (or
#' `<base_dir>/seed_<n>/<block>/`), stacks them with a `seed` column and
#' summarises across seeds. Works for seeds run in one R session by
#' [run_simulation_seeds()] and for seeds run as separate jobs.
#'
#' Written to `<base_dir>/` (or `<base_dir>/all_seeds/<block>/`): for each of
#' `null`, `spikein` and `spikein_by_density` that is present, a
#' `<name>_by_seed.csv` and a `<name>_across_seeds.csv`.
#'
#' @param base_dir Directory holding the `seed_<n>` directories.
#' @param block Optional sub-directory inside each seed directory.
#' @return Invisibly, a list with elements `by_seed` and `across_seeds`, each a
#'   list of data.tables named `null`, `spikein`, `spikein_by_density`.
#' @export
combine_seed_dirs <- function(base_dir, block = NULL) {
  seed_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed_-?[0-9]+$", basename(seed_dirs))]
  if (!length(seed_dirs)) {
    stop("combine_seed_dirs: no seed_<n> directories under ", base_dir,
         call. = FALSE)
  }
  seeds <- as.integer(sub("^seed_", "", basename(seed_dirs)))
  seed_dirs <- seed_dirs[order(seeds)]
  seeds <- sort(seeds)
  if (!is.null(block)) seed_dirs <- file.path(seed_dirs, block)
  out_dir <- if (is.null(block)) base_dir else file.path(base_dir, "all_seeds", block)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  spec <- list(
    null = list(file = "null_simulation_results.csv",
                by = c("scenario", "method"),
                values = c("FP", "domain_calls", "domain_retained")),
    spikein = list(file = "spikein_simulation_results.csv",
                   by = c("efficiency_scenario", "effect_size", "method"),
                   values = c("sensitivity", "precision", "FDR", "F1", "TP", "FP",
                              "regions_detected", "total_called",
                              "domain_calls", "domain_retained")),
    spikein_by_density = list(file = "spikein_by_density.csv",
                              by = c("efficiency_scenario", "effect_size",
                                     "method", "density_stratum"),
                              values = c("sensitivity", "regions_detected"))
  )

  by_seed <- list()
  across <- list()
  for (nm in names(spec)) {
    paths <- file.path(seed_dirs, spec[[nm]]$file)
    have <- file.exists(paths)
    if (!any(have)) next
    if (!all(have)) {
      warning("combine_seed_dirs: ", spec[[nm]]$file, " is missing for seed(s) ",
              paste(seeds[!have], collapse = ", "),
              "; those seeds are left out of the ", nm, " summary.", call. = FALSE)
    }
    stacked <- rbindlist(Map(function(p, s) {
      x <- fread(p)
      x[, seed := s]
      x
    }, paths[have], seeds[have]), fill = TRUE)
    setcolorder(stacked, "seed")
    by_seed[[nm]] <- stacked
    across[[nm]] <- summarize_across_seeds(stacked, spec[[nm]]$by, spec[[nm]]$values)
    fwrite(stacked, file.path(out_dir, paste0(nm, "_by_seed.csv")))
    fwrite(across[[nm]], file.path(out_dir, paste0(nm, "_across_seeds.csv")))
  }
  if (!length(by_seed)) {
    stop("combine_seed_dirs: no result tables found under ",
         paste(seed_dirs, collapse = ", "), call. = FALSE)
  }

  message(sprintf("Combined %d seed(s) (%s) into %s", length(seeds),
                  paste(seeds, collapse = ", "), out_dir))
  invisible(list(by_seed = by_seed, across_seeds = across))
}


#' Run the null and spike-in simulations over several seeds
#'
#' Runs each seed into `<config$output_dir>/seed_<n>/` and then calls
#' [combine_seed_dirs()]. Seeds run one after another; for genome-wide runs it
#' is usually better to run one seed per job (see the seeds array job script
#' in `inst/slurm/`) and combine afterwards.
#'
#' @param wt_reps Named list of source replicate data.tables.
#' @param config Configuration list from [parse_args()]. `config$seed` is
#'   replaced by each element of `seeds` in turn.
#' @param seeds Integer vector of seeds (no duplicates).
#' @param layers Which simulations to run: "null", "spikein" or both.
#' @return Invisibly, the list returned by [combine_seed_dirs()].
#' @export
run_simulation_seeds <- function(wt_reps, config, seeds,
                                 layers = c("null", "spikein")) {
  layers <- match.arg(layers, several.ok = TRUE)
  seeds <- as.integer(seeds)
  if (!length(seeds) || anyNA(seeds) || anyDuplicated(seeds)) {
    stop("run_simulation_seeds: `seeds` must be distinct whole numbers.",
         call. = FALSE)
  }
  base_dir <- config$output_dir

  for (s in seeds) {
    message("\n", strrep("#", 60), sprintf("\n# SEED %d\n", s), strrep("#", 60))
    cfg <- config
    cfg$seed <- s
    cfg$output_dir <- file.path(base_dir, sprintf("seed_%d", s))
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(cfg, file.path(cfg$output_dir, "config.rds"))
    if ("null" %in% layers) run_null_simulation(wt_reps, cfg)
    if ("spikein" %in% layers) run_spikein_simulation(wt_reps, cfg)
  }

  combine_seed_dirs(base_dir)
}
