#!/usr/bin/env Rscript

# 02_run_methods.R
# Apply normalization methods to simulated data and call DMRs


# Format helpers ----------------------------------------------------------

#' Convert pseudo-group data to the format expected by SMFnorm
#'
#' Build a nested list with sample_metadata attribute, like load_data() output.
#'
#' @param pseudo_groups List from create_pseudo_groups() with PseudoA and PseudoB.
#' @return List structured as load_data() output (flat sample list with metadata).
format_for_smfnorm <-
  function(pseudo_groups) {
  all_samples <- c(pseudo_groups$PseudoA, pseudo_groups$PseudoB)

  # Build sample metadata
  sample_ids <- names(all_samples)
  group_ids <- ifelse(grepl("^PseudoA", sample_ids), "PseudoA", "PseudoB")

  metadata <- data.table(
    group_id = group_ids,
    replicate = seq_along(sample_ids),
    sample_id = sample_ids,
    file_name = paste0(sample_ids, ".tsv.gz")  # placeholder
  )

  attr(all_samples, "sample_metadata") <- metadata
  return(all_samples)
}


#' Filter pseudo-group samples to shared sites meeting the minimum coverage.
#'
#' SMFnorm requires all input samples to share the same site set. If we only
#' apply `min_coverage` inside SMFnorm, samples with uneven coverage may be
#' filtered differently and trigger unequal row counts.
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param min_coverage Minimum per-sample coverage threshold.
#' @return Pseudo-group list with all samples filtered to the same shared sites.
.filter_pseudo_groups_by_min_coverage <- function(pseudo_groups, min_coverage) {
  if (is.null(min_coverage) || min_coverage <= 0) {
    return(pseudo_groups)
  }

  all_samples <- c(pseudo_groups$PseudoA, pseudo_groups$PseudoB)
  keep <- Reduce(`&`, lapply(all_samples, function(dt) dt$cov >= min_coverage))

  if (!any(keep)) {
    stop(sprintf("No sites remain after applying min_coverage = %d", min_coverage))
  }

  filtered <- lapply(all_samples, function(dt) dt[keep])
  names(filtered) <- names(all_samples)

  n_A <- length(pseudo_groups$PseudoA)
  list(
    PseudoA = filtered[seq_len(n_A)],
    PseudoB = filtered[seq(n_A + 1, length(filtered))],
    params = pseudo_groups$params
  )
}


#' Export samples to allc-format TSV files for external tools (metilene, ComBatMet)
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param out_dir Directory to write files to.
#' @return data.table mapping sample_id to file_path and group.
export_allc_files <- function(pseudo_groups, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  all_samples <- c(pseudo_groups$PseudoA, pseudo_groups$PseudoB)
  file_map <- data.table(
    sample_id = character(),
    group = character(),
    file_path = character()
  )

  for (nm in names(all_samples)) {
    dt <- all_samples[[nm]]
    group <- ifelse(grepl("^PseudoA", nm), "PseudoA", "PseudoB")

    # Write in allc format: chr, pos, strand, site, mc, cov
    out_path <- file.path(out_dir, paste0(nm, ".tsv.gz"))
    fwrite(dt[, .(chr, pos, strand, site, mc, cov)],
           out_path, sep = "\t", col.names = FALSE)

    file_map <- rbind(file_map, data.table(
      sample_id = nm, group = group, file_path = out_path
    ))
  }

  message(sprintf("Exported %d files to %s", nrow(file_map), out_dir))
  return(file_map)
}


# Method: Raw (no normalization) ------------------------------------------

#' Run analysis on raw (unnormalized) data
#'
#' Simply reformats data for DMR calling without any normalization.
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @return Named list with method = "raw" and the split data.
run_raw <- function(pseudo_groups) {
  message("Method: Raw (no normalization)")

  split_data <- list(
    PseudoA = pseudo_groups$PseudoA,
    PseudoB = pseudo_groups$PseudoB
  )

  return(list(method = "raw", data = split_data))
}


# Method: Downsampling ----------------------------------------------------

#' Downsample coverage to the minimum across all samples
#'
#' For each site, downsample mc and cov to the minimum total coverage observed
#' across all samples at that site. Uses hypergeometric sampling.
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @return Named list with method = "downsampled" and normalized split data.
run_downsampled <- function(pseudo_groups) {
  message("Method: Downsampling to minimum coverage")

  all_samples <- c(pseudo_groups$PseudoA, pseudo_groups$PseudoB)

  # Find minimum coverage at each site across all samples
  # Assume all samples have the same sites in the same order (post find_shared_sites)
  n_sites <- nrow(all_samples[[1]])
  min_cov <- rep(Inf, n_sites)

  for (s in all_samples) {
    min_cov <- pmin(min_cov, s$cov)
  }
  min_cov <- as.integer(pmax(0, min_cov))

  # Downsample each sample
  downsampled <- lapply(all_samples, function(dt) {
    out <- copy(dt)

    # For each site, sample min_cov reads from cov total reads,
    # where mc of them are methylated → hypergeometric
    new_mc <- rhyper(
      nn = n_sites,
      m  = pmax(0L, as.integer(out$mc)),        # white balls (methylated)
      n  = pmax(0L, as.integer(out$cov - out$mc)),  # black balls (unmethylated)
      k  = min_cov                               # draws
    )

    out[, `:=`(
      mc  = new_mc,
      cov = min_cov,
      rate = ifelse(min_cov > 0, new_mc / min_cov, 0)
    )]

    return(out)
  })

  names(downsampled) <- names(all_samples)
  n_A <- length(pseudo_groups$PseudoA)

  split_data <- list(
    PseudoA = downsampled[seq_len(n_A)],
    PseudoB = downsampled[seq(n_A + 1, length(downsampled))]
  )

  message(sprintf("Downsampled to median min coverage: %d", median(min_cov)))
  return(list(method = "downsampled", data = split_data))
}


# Method: SMFnorm -------------------------------------------------------

#' Run SMFnorm normalization pipeline
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param within_alpha Alpha for within-group rate normalization.
#' @param between_alpha Alpha for between-group rate normalization.
#' @param min_coverage Minimum coverage threshold.
#' @param rate_between_groups Logical: correct between-group rate differences.
#'   Default FALSE (within-group only). Set TRUE to correct a systematic
#'   between-group efficiency artifact (e.g. the `aligned_*` bias scenarios);
#'   leave FALSE for the matched-mean within-group scenarios.
#' @return Named list with method = "SMFnorm" and normalized split data.
run_SMFnorm <- function(pseudo_groups,
                        within_alpha = 0.3,
                        between_alpha = 0.9,
                        min_coverage = 10,
                        rate_between_groups = FALSE) {
  message(sprintf("Method: SMFnorm (rate_between_groups = %s)",
                  rate_between_groups))

  # Ensure all samples share the same site set at the requested coverage.
  pseudo_groups <- .filter_pseudo_groups_by_min_coverage(pseudo_groups,
                                                         min_coverage)
  formatted <- format_for_smfnorm(pseudo_groups)

  # Split by groups
  split_data <- SMFnorm::split_by_groups(formatted)

  # Run full normalization
  normalized <- SMFnorm::normalize_methylation_data(
    data_list = split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    coverage_between_groups = FALSE,
    rate_within_groups = TRUE,
    rate_between_groups = rate_between_groups,
    within_alpha = within_alpha,
    between_alpha = between_alpha,
    min_coverage = min_coverage
  )

  return(list(method = "SMFnorm", data = normalized))
}


# Method: ComBatMet (via external call) -----------------------------------

#' Derive a ComBat batch label from per-replicate enzyme efficiency
#'
#' ComBat requires a discrete batch variable, but the artifact SMFsim simulates
#' is a continuous per-replicate efficiency. Splitting replicates at the median
#' efficiency gives the batch label an analyst could plausibly construct from
#' QC (high- vs low-efficiency samples).
#'
#' The split is deliberately scenario-dependent, and that is the point:
#' \itemize{
#'   \item `imbalanced_*` scenarios draw the two groups' efficiencies from
#'     OVERLAPPING ranges, so the split crosses the group boundary and ComBat
#'     has a genuine, separable artifact to remove.
#'   \item `aligned_*` scenarios use non-overlapping ranges, so the split
#'     reproduces the group labels exactly. Batch is then perfectly confounded
#'     with group and no batch-correction method can separate artifact from
#'     biology.
#' }
#'
#' @param params `params` element of a [create_pseudo_groups()] result.
#' @param n_A,n_B Replicate counts for the two pseudo-groups.
#' @return Factor of batch labels, PseudoA replicates first then PseudoB.
#' @keywords internal
.efficiency_batch <- function(params, n_A, n_B) {
  eff <- c(params$efficiency_A, params$efficiency_B)
  if (length(eff) != n_A + n_B) {
    stop("run_combatmet: params$efficiency_A/B hold ", length(eff),
         " values but there are ", n_A + n_B, " replicates. Pass an explicit ",
         "`batch` instead.", call. = FALSE)
  }
  factor(ifelse(eff > stats::median(eff), "high_eff", "low_eff"),
         levels = c("low_eff", "high_eff"))
}


#' Run ComBatMet normalization
#'
#' Wraps [ComBatMet::ComBat_met()], the published ComBat-met beta-regression
#' batch correction, on the per-site rate matrix.
#'
#' The biological group is passed as `group` with `full_mod = TRUE` so it enters
#' the model and is PROTECTED; only batch effects are removed. Batch defaults to
#' an efficiency stratum from [.efficiency_batch()].
#'
#' Note: an earlier version of this function hand-rolled a "mean-only"
#' adjustment that subtracted each group's own mean and added the grand mean.
#' That forces both group means to the grand mean at every site, making the
#' group difference identically zero, so no DMR could ever be called. It was
#' also not ComBat-met. Both problems are fixed here.
#'
#' @param pseudo_groups List from [create_pseudo_groups()].
#' @param batch Optional explicit batch factor, one entry per replicate ordered
#'   PseudoA then PseudoB. Defaults to the median-efficiency split.
#' @return Named list with method = "ComBatMet" and normalized split data.
run_combatmet <- function(pseudo_groups, batch = NULL) {
  if (!requireNamespace("ComBatMet", quietly = TRUE)) {
    stop("run_combatmet: package 'ComBatMet' is not installed. Install it with ",
         "remotes::install_github('JmWangBio/ComBatMet'), or drop 'ComBatMet' ",
         "from config$methods.", call. = FALSE)
  }

  n_A <- length(pseudo_groups$PseudoA)
  n_B <- length(pseudo_groups$PseudoB)

  all_samples <- c(lapply(pseudo_groups$PseudoA, data.table::copy),
                   lapply(pseudo_groups$PseudoB, data.table::copy))
  names(all_samples) <- c(names(pseudo_groups$PseudoA),
                          names(pseudo_groups$PseudoB))

  group <- factor(rep(c("PseudoA", "PseudoB"), c(n_A, n_B)),
                  levels = c("PseudoA", "PseudoB"))
  if (is.null(batch)) {
    batch <- .efficiency_batch(pseudo_groups$params, n_A, n_B)
  }
  batch <- as.factor(batch)

  message("Method: ComBatMet (ComBat_met, group protected)")
  message("  batch: ",
          paste(sprintf("%s=%s", names(all_samples), as.character(batch)),
                collapse = ", "))

  # ComBat cannot separate a batch effect from biology when every batch level
  # sits entirely inside one group. Fail loudly rather than silently returning
  # group-centred (zero-difference) data, which is what the old implementation
  # did for every scenario.
  confounded <- all(rowSums(table(batch, group) > 0) == 1)
  if (nlevels(batch) < 2 || confounded) {
    stop("run_combatmet: batch is perfectly confounded with group (each batch ",
         "level falls entirely within one group), so ComBat cannot separate ",
         "the efficiency artifact from biological signal. This is expected for ",
         "the aligned_* scenarios, whose efficiency ranges do not overlap - ",
         "report the method as inapplicable there, not as zero sensitivity.",
         call. = FALSE)
  }

  rate_matrix <- do.call(cbind,
                         lapply(all_samples, function(dt) as.numeric(dt$rate)))
  colnames(rate_matrix) <- names(all_samples)

  # ComBat_met swaps in pseudo_beta for exact 0/1 but errors on anything outside
  # (0, 1); keep inputs strictly inside the beta support.
  eps <- 1e-6
  rate_matrix[rate_matrix < eps] <- eps
  rate_matrix[rate_matrix > 1 - eps] <- 1 - eps

  adjusted <- ComBatMet::ComBat_met(
    vmat = rate_matrix,
    dtype = "b-value",
    batch = batch,
    group = group,
    full_mod = TRUE,
    ncores = 1
  )

  if (!identical(dim(adjusted), dim(rate_matrix))) {
    stop("run_combatmet: ComBat_met returned a ", nrow(adjusted), " x ",
         ncol(adjusted), " matrix, expected ", nrow(rate_matrix), " x ",
         ncol(rate_matrix), ".", call. = FALSE)
  }

  result_samples <- lapply(seq_along(all_samples), function(i) {
    out <- all_samples[[i]]
    out[, rate := as.numeric(adjusted[, i])]
    out[, mc := as.integer(round(cov * rate))]
    out
  })
  names(result_samples) <- names(all_samples)

  split_data <- list(
    PseudoA = result_samples[seq_len(n_A)],
    PseudoB = result_samples[n_A + seq_len(n_B)]
  )

  return(list(method = "ComBatMet", data = split_data))
}

# DMR calling with metilene -----------------------------------------------

#' Prepare bedGraph files for metilene input
#'
#' Metilene expects per-group bedGraph files with one column per replicate.
#'
#' @param split_data Named list with PseudoA and PseudoB, each a list of data.tables.
#' @param out_dir Output directory for bedGraph files.
#' @return List with file paths for group A and B bedGraphs.
prepare_metilene_input <- function(split_data, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  .write_group_bedgraph <- function(group_data, group_name) {
    # Merge all replicates on chr + pos
    merged <- Reduce(function(x, y) {
      merge(x, y, by = c("chr", "pos"), suffixes = c("", paste0("_", ncol(x))))
    }, lapply(seq_along(group_data), function(i) {
      dt <- group_data[[i]][, .(chr, pos, rate)]
      setnames(dt, "rate", paste0("rate_", names(group_data)[i]))
      dt
    }))

    # Sort by chr and pos
    setorder(merged, chr, pos)

    # Write as bedGraph: chr, start (0-based), end, rate_1, rate_2, ...
    out <- copy(merged)
    out[, start := pos - 1L]
    out[, end := pos]

    # Reorder columns
    rate_cols <- grep("^rate_", names(out), value = TRUE)
    out <- out[, c("chr", "start", "end", rate_cols), with = FALSE]

    out_path <- file.path(out_dir, paste0(group_name, "_metilene.bedGraph"))
    fwrite(out, out_path, sep = "\t", col.names = FALSE)

    return(out_path)
  }

  path_A <- .write_group_bedgraph(split_data$PseudoA, "PseudoA")
  path_B <- .write_group_bedgraph(split_data$PseudoB, "PseudoB")

  return(list(groupA = path_A, groupB = path_B))
}


#' Call DMRs using metilene
#'
#' DMRs are kept if their q-value is below `q_cutoff`. Two q-values are stored:
#' \describe{
#'   \item{`q_metilene` (default)}{metilene's column-4 q-value, corrected over
#'     ALL tested segments using the method chosen by `mtc`: Bonferroni
#'     (`mtc = 1`, metilene's own default) or Benjamini-Hochberg / FDR
#'     (`mtc = 2`). Because the denominator covers every tested segment rather
#'     than only the emitted ones, this is the statistically valid choice for
#'     the standard `min_diff > 0` workflow. Validated against the null (gives
#'     ~0 false positives where it should).}
#'   \item{`q_bh`}{Benjamini-Hochberg recomputed in R from metilene's
#'     Mann-Whitney U p-value (column 7), corrected over the FULL number of
#'     tests metilene reports (`Number of Tests:` on stderr, read by
#'     [.metilene_n_tests()]), not over the emitted rows.}
#' }
#'
#' The denominator is the whole problem here. metilene emits far fewer segments
#' than it tests - on one comparison, 985949 tests against 32162 emitted rows at
#' `-d 0`, and only 568 rows at `-d 0.1`. Running `p.adjust()` over the emitted
#' rows therefore under-corrects by roughly 30x at `-d 0` and 1700x at
#' `-d 0.1`; the latter produced 819 false positives in a true null where the
#' column-4 q gave 0.
#'
#' To compute your own FDR in R rather than relying on metilene's correction:
#' ```
#' call_dmrs_metilene(..., min_diff = 0, q_source = "BH", min_effect = 0.1)
#' ```
#' `min_diff` is metilene's `-d` (which segments get emitted), while
#' `min_effect` is an R-side `|mean_diff|` filter applied AFTER q-value
#' filtering, so significance and effect size stay separable.
#'
#' Assumption worth stating in any write-up: supplying `p.adjust()` a subset of
#' p-values with `n =` the full test count treats the withheld tests as
#' non-significant. That is appropriate when the withheld tests are less
#' significant than the emitted ones, which is why `min_diff = 0` matters - it
#' makes the emitted set as complete as metilene will report. Validate against
#' the null (a true null should yield ~0 false positives) rather than assuming.
#'
#' @param split_data Named list with PseudoA and PseudoB.
#' @param out_dir Output directory.
#' @param metilene_path Path to metilene binary (default: assumes in PATH).
#' @param min_cpg Minimum number of CpGs/sites in a DMR (default: 5).
#' @param min_diff Minimum mean difference to report (default: 0.1).
#' @param metilene_max_dist Maximum distance between CpGs within a DMR.
#' @param q_cutoff Significance threshold on the chosen q-value (default 0.05).
#'   Set to 1 to keep all candidate DMRs.
#' @param q_source Which q-value to filter on: "metilene" (default, column-4 q)
#'   or "BH" (R-recomputed; only valid with `min_diff = 0`).
#' @param p_column Which metilene p-value the R-side BH is computed from:
#'   "2dks" (column 8, 2D Kolmogorov-Smirnov; the default) or "mwu" (column 7,
#'   Mann-Whitney U). metilene's manual says its q-value corrects the MWU
#'   p-value, but that is WRONG for de-novo mode: verified by backing the
#'   denominator out of a Bonferroni run, `q / p_2dks` is constant at metilene's
#'   reported test count while `q / p_mwu` varies over ~1200 orders of
#'   magnitude. Defaulting to "2dks" keeps `q_bh` and `q_metilene` measuring the
#'   same underlying statistic, so they differ only by correction method and
#'   denominator. See `inst/scripts/diagnose_metilene_qvalues.R`.
#' @param min_effect Minimum `|mean_diff|` required of a DMR, applied in R
#'   AFTER q-value filtering. Default 0 (no filtering) preserves the historical
#'   behaviour, where the effect-size cut was done solely by metilene's `-d`.
#'   Set this together with `min_diff = 0` and `q_source = "BH"` to compute the
#'   FDR yourself over all tested segments.
#' @param mtc Multiple-testing correction metilene applies to its column-4
#'   q-value, passed through as `-c`: 1 = Bonferroni (metilene's default,
#'   retained here so existing results stay reproducible), 2 = Benjamini-Hochberg
#'   (FDR). Note that `mtc = 2` computes BH over ALL tested segments inside
#'   metilene, which is NOT the same as - and is statistically preferable to -
#'   the selection-biased `q_source = "BH"` recomputation.
#' @return data.table of significant DMRs with columns: chr, start, end,
#'   q_value (the chosen one), mean_diff, n_sites, p_mwu, q_metilene, q_bh.
call_dmrs_metilene <- function(split_data,
                               out_dir,
                               metilene_path = "metilene",
                               min_cpg = 5,
                               min_diff = 0.1,
                               metilene_max_dist = 1500,
                               q_cutoff = 0.05,
                               q_source = c("metilene", "BH"),
                               p_column = c("2dks", "mwu"),
                               min_effect = 0,
                               mtc = 1L) {
  q_source <- match.arg(q_source)
  p_column <- match.arg(p_column)

  mtc <- as.integer(mtc)
  if (length(mtc) != 1L || is.na(mtc) || !mtc %in% c(1L, 2L)) {
    stop("call_dmrs_metilene: `mtc` must be 1 (Bonferroni) or 2 ",
         "(Benjamini-Hochberg/FDR); got ", deparse(mtc), ".", call. = FALSE)
  }

  # The R-side BH is only honest over a complete set of tested segments. With
  # metilene's -d pre-filter in play, the .bed holds only high-difference
  # candidates, so p.adjust() sees a selected subset and under-corrects.
  if (q_source == "BH" && min_diff > 0) {
    warning("call_dmrs_metilene: q_source = 'BH' with min_diff = ", min_diff,
            " (> 0) is ANTI-CONSERVATIVE. metilene's -d pre-filters the ",
            "emitted segments, so the R-side BH corrects over a selected ",
            "subset. Use min_diff = 0 with min_effect = ",
            if (min_effect > 0) min_effect else "<your effect cut>",
            ", or q_source = 'metilene'.", call. = FALSE)
  }

  if (length(metilene_path) == 0L || !nzchar(metilene_path)) {
    stop("call_dmrs_metilene: metilene_path is empty/NULL. ",
         "Set config$metilene_path to the metilene binary path.", call. = FALSE)
  }

  # Prepare input files
  input_dir <- file.path(out_dir, "metilene_input")
  bedgraph_files <- prepare_metilene_input(split_data, input_dir)

  # Build metilene command
  # metilene expects a single input file with group columns interleaved
  # We need to merge group A and group B bedGraphs
  merged_input <- .merge_group_bedgraphs(
    bedgraph_files$groupA,
    bedgraph_files$groupB,
    split_data,
    file.path(input_dir, "merged_input.bed")
  )

  n_A <- length(split_data$PseudoA)
  n_B <- length(split_data$PseudoB)

  out_file <- file.path(out_dir, "metilene_dmrs.bed")

  # metilene reports its total test count on stderr ("Number of Tests: N").
  # That is the correct BH denominator, and it is much larger than the number of
  # segments emitted: on one comparison, 985949 tests vs 32162 emitted rows even
  # at -d 0. Capture the log rather than discarding it to /dev/null.
  log_file <- file.path(out_dir, "metilene_log.txt")

  # -c selects metilene's multiple-testing correction for its column-4 q-value:
  # 1 = Bonferroni (metilene's default), 2 = Benjamini-Hochberg (FDR). Passed
  # explicitly so the correction is recorded rather than left implicit.
  cmd <- sprintf(
    "%s -a PseudoA -b PseudoB -m %d -d %f -M %d -c %d %s > %s 2> %s",
    metilene_path, min_cpg, min_diff, metilene_max_dist, mtc,
    merged_input, out_file, log_file
  )

  message("Running metilene: ", cmd)
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("metilene exited with code ", exit_code,
            ". See log: ", log_file)
    if (file.exists(log_file)) {
      message(paste(utils::head(readLines(log_file, warn = FALSE), 20),
                    collapse = "\n"))
    }
  }

  n_tests <- .metilene_n_tests(log_file)

  # Parse metilene output
  empty <- data.table(
    chr = character(), start = integer(), end = integer(),
    q_value = numeric(), mean_diff = numeric(), n_sites = integer(),
    p_mwu = numeric(), q_metilene = numeric(), q_bh = numeric()
  )

  if (!file.exists(out_file) || file.size(out_file) == 0) {
    message("No DMRs called by metilene")
    return(empty)
  }

  dmrs <- fread(out_file, header = FALSE)
  # metilene DMR output columns:
  #   1 chr  2 start  3 end  4 q(metilene)  5 mean_diff  6 #CpGs
  #   7 p(MWU)  8 p(2D-KS)  9 mean_g1  10 mean_g2
  if (ncol(dmrs) >= 8) {
    setnames(dmrs, 1:8, c("chr", "start", "end", "q_metilene", "mean_diff",
                          "n_sites", "p_mwu", "p_2dks"))

    # BH over ALL tests metilene performed, not just the rows it emitted.
    # p.adjust(n =) treats the unreported tests as non-significant, which is the
    # correct conservative handling: metilene emits its best-scoring segments,
    # so the withheld tests are less significant than the reported ones.
    n_emitted <- nrow(dmrs)
    if (is.na(n_tests)) {
      warning("call_dmrs_metilene: could not read 'Number of Tests' from ",
              log_file, "; falling back to BH over the ", n_emitted,
              " emitted segments only. This UNDER-CORRECTS - metilene ",
              "typically tests far more segments than it emits.", call. = FALSE)
      n_bh <- n_emitted
    } else if (n_tests < n_emitted) {
      warning("call_dmrs_metilene: metilene reported ", n_tests, " tests but ",
              "emitted ", n_emitted, " segments; using the emitted count as ",
              "the BH denominator.", call. = FALSE)
      n_bh <- n_emitted
    } else {
      n_bh <- n_tests
    }
    p_vec <- if (p_column == "2dks") dmrs$p_2dks else dmrs$p_mwu
    dmrs[, q_bh := stats::p.adjust(p_vec, method = "BH", n = n_bh)]
  } else if (ncol(dmrs) >= 6) {
    warning("metilene output has <8 columns; only column-4 q available")
    setnames(dmrs, 1:6, c("chr", "start", "end", "q_metilene", "mean_diff",
                          "n_sites"))
    dmrs[, `:=`(p_mwu = NA_real_, q_bh = NA_real_)]
  } else {
    return(empty)
  }

  # Significance q-value. "BH" is the R-side value corrected over metilene's
  # full test count; "metilene" is its own column-4 q.
  dmrs[, q_value := if (q_source == "BH") q_bh else q_metilene]

  n_total <- nrow(dmrs)
  dmrs <- dmrs[is.finite(q_value) & q_value < q_cutoff]
  n_sig <- nrow(dmrs)

  # Effect-size filter as its own step, so significance and effect size stay
  # separable (see the min_diff = 0 + q_source = "BH" workflow above).
  if (min_effect > 0) {
    dmrs <- dmrs[is.finite(mean_diff) & abs(mean_diff) >= min_effect]
  }

  dmrs <- dmrs[, .(chr, start, end, q_value, mean_diff, n_sites,
                   p_mwu, q_metilene, q_bh)]
  message(sprintf(
    "metilene: %d candidate DMRs, %d significant (%s q < %.3g)%s",
    n_total, n_sig, q_source, q_cutoff,
    if (min_effect > 0)
      sprintf(", %d after |mean_diff| >= %.3g", nrow(dmrs), min_effect) else ""))
  return(dmrs)
}


#' Read metilene's total test count from its log
#'
#' metilene prints `Number of Tests: N` on stderr after segmenting. That N is
#' the correct denominator for any multiple-testing correction, and it is far
#' larger than the number of segments metilene emits: on one comparison, 985949
#' tests against 32162 emitted rows even with `-d 0`. Correcting over only the
#' emitted rows under-corrects by roughly that ratio.
#'
#' @param log_file Path to the captured metilene stderr log.
#' @return Integer test count, or `NA_integer_` if the log is missing or does
#'   not contain the line.
#' @keywords internal
.metilene_n_tests <- function(log_file) {
  if (is.null(log_file) || !file.exists(log_file)) return(NA_integer_)
  txt <- tryCatch(readLines(log_file, warn = FALSE),
                  error = function(e) character(0))
  if (!length(txt)) return(NA_integer_)

  hit <- regmatches(
    txt, regexpr("Number\\s+of\\s+Tests\\s*:\\s*[0-9]+", txt, ignore.case = TRUE))
  hit <- hit[nzchar(hit)]
  if (!length(hit)) return(NA_integer_)

  # Take the last occurrence in case metilene reports per-chromosome counts.
  n <- suppressWarnings(as.numeric(
    sub(".*:\\s*", "", hit[length(hit)])))
  if (!is.finite(n) || n < 1) return(NA_integer_)
  as.integer(n)
}


#' Merge two group bedGraphs into metilene's expected format
#' @keywords internal
.merge_group_bedgraphs <- function(path_A, path_B, split_data, out_path) {
  dt_A <- fread(path_A, header = FALSE)
  dt_B <- fread(path_B, header = FALSE)
  
  n_A <- length(split_data$PseudoA)
  n_B <- length(split_data$PseudoB)
  
  # Rename rate columns with group prefix to match metilene -a/-b args
  rate_cols_A <- paste0("V", 4:(3 + n_A))
  rate_cols_B <- paste0("V", 4:(3 + n_B))
  
  setnames(dt_A, rate_cols_A, paste0("PseudoA_", seq_len(n_A)))
  setnames(dt_B, rate_cols_B, paste0("PseudoB_", seq_len(n_B)))
  
  merged <- merge(
    dt_A, dt_B,
    by = c("V1", "V2", "V3"),
    sort = FALSE
  )
  
  # Rename coordinate columns
  setnames(merged, c("V1", "V2", "V3"), c("chr", "start", "end"))
  
  # Reorder: chr, start, end, all A reps, all B reps
  a_cols <- paste0("PseudoA_", seq_len(n_A))
  b_cols <- paste0("PseudoB_", seq_len(n_B))
  merged <- merged[, c("chr", "start", "end", a_cols, b_cols), with = FALSE]
  
  setorder(merged, chr, start)
  fwrite(merged, out_path, sep = "\t", col.names = TRUE)
  
  return(out_path)
}


# Master method runner ----------------------------------------------------

#' Run all normalization methods on pseudo-group data
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param methods Character vector of methods to run.
#'   Options: "raw", "downsampled", "SMFnorm", "ComBatMet"
#' @param SMFnorm_params Named list of SMFnorm parameters.
#' @return Named list of results, each with $method and $data.
run_all_methods <- function(pseudo_groups,
                            methods = c("raw", "downsampled", "SMFnorm", "ComBatMet"),
                            SMFnorm_params = list(
                              within_alpha = 0.3,
                              between_alpha = 0.9,
                              min_coverage = 10
                            )) {
  # NOTE: min_coverage is applied upstream at load time (prepare_wt_replicates /
  # load_real_data), before find_shared_sites, so every method here receives the
  # identical filtered, shared-site set. Do not filter again at this stage.
  if (length(methods) == 0L) {
    stop("run_all_methods: `methods` is empty/NULL. Set config$methods to a ",
         "non-empty subset of c('raw','downsampled','SMFnorm','ComBatMet').",
         call. = FALSE)
  }

  results <- list()

  for (m in methods) {
    message(sprintf("\n=== Running method: %s ===", m))
    t0 <- Sys.time()
    
    # Deep copy so each method gets fresh data
    pseudo_copy <- list(
      PseudoA = lapply(pseudo_groups$PseudoA, data.table::copy),
      PseudoB = lapply(pseudo_groups$PseudoB, data.table::copy),
      params = pseudo_groups$params
    )
    
    result <- tryCatch({
      switch(m,
             raw = run_raw(pseudo_copy),
             downsampled = run_downsampled(pseudo_copy),
             SMFnorm = do.call(run_SMFnorm, c(list(pseudo_groups = pseudo_copy),
                                                  SMFnorm_params)),
             ComBatMet = run_combatmet(pseudo_copy),
             stop("Unknown method: ", m)
      )
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", m, e$message))
      return(NULL)
    })
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("Method %s completed in %.1f seconds", m, elapsed))
    
    if (!is.null(result)) {
      result$elapsed_seconds <- elapsed
      results[[m]] <- result
    }
  }

  if (length(results) == 0L) {
    stop(sprintf("run_all_methods: every requested method failed (%s) - no ",
                 "results to return. See warnings() for the per-method error(s).",
                 paste(methods, collapse = ", ")), call. = FALSE)
  }

  return(results)
}


