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
                        min_coverage = 5,
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

#' Run ComBatMet normalization
#'
#' ComBatMet uses a ComBat-based approach on methylation beta values.
#' This wrapper prepares a matrix, runs ComBat, and converts back.
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @return Named list with method = "ComBatMet" and normalized split data.
run_combatmet <- function(pseudo_groups) {
  message("Method: ComBatMet (manual mean-only)")
  
  all_samples <- lapply(
    c(pseudo_groups$PseudoA, pseudo_groups$PseudoB),
    data.table::copy
  )
  names(all_samples) <- c(names(pseudo_groups$PseudoA), names(pseudo_groups$PseudoB))
  n_A <- length(pseudo_groups$PseudoA)
  nc <- length(all_samples)
  
  # Build rate matrix
  rate_matrix <- do.call(cbind, lapply(all_samples, function(dt) as.numeric(dt$rate)))
  colnames(rate_matrix) <- names(all_samples)
  nr <- nrow(rate_matrix)
  
  # Bound rates away from 0 and 1 (in place to preserve matrix)
  rate_matrix[rate_matrix < 1e-6] <- 1e-6
  rate_matrix[rate_matrix > (1 - 1e-6)] <- 1 - 1e-6
  
  # Logit transform (in place to preserve matrix)
  logit_rates <- rate_matrix
  logit_rates[] <- log(rate_matrix / (1 - rate_matrix))
  
  message("  logit_rates dim: ", nrow(logit_rates), " x ", ncol(logit_rates))
  
  # Mean-only batch correction
  grand_mean <- rowMeans(logit_rates)
  batch_A_mean <- rowMeans(logit_rates[, 1:n_A, drop = FALSE])
  batch_B_mean <- rowMeans(logit_rates[, (n_A+1):nc, drop = FALSE])
  
  adjusted <- logit_rates
  adjusted[, 1:n_A] <- adjusted[, 1:n_A] - batch_A_mean + grand_mean
  adjusted[, (n_A+1):nc] <- adjusted[, (n_A+1):nc] - batch_B_mean + grand_mean
  
  # Inverse logit
  adjusted_rates <- adjusted
  adjusted_rates[] <- 1 / (1 + exp(-adjusted))
  
  message("  adjusted_rates dim: ", nrow(adjusted_rates), " x ", ncol(adjusted_rates))
  
  result_samples <- lapply(seq_len(nc), function(i) {
    out <- copy(all_samples[[i]])
    out[, rate := adjusted_rates[, i]]
    out[, mc := as.integer(round(cov * rate))]
    return(out)
  })
  names(result_samples) <- names(all_samples)
  
  split_data <- list(
    PseudoA = result_samples[seq_len(n_A)],
    PseudoB = result_samples[seq(n_A + 1, nc)]
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
#'   \item{`q_metilene` (default)}{metilene's column-4 q-value, which is
#'     Benjamini-Hochberg corrected over ALL tested segments. Validated against
#'     the null (gives ~0 false positives where it should) and is the
#'     statistically valid choice for the standard `min_diff > 0` workflow.}
#'   \item{`q_bh`}{Benjamini-Hochberg recomputed in R from metilene's
#'     Mann-Whitney U p-value (column 7). NOTE: this is ANTI-CONSERVATIVE when
#'     `min_diff > 0`, because the reported DMRs are pre-filtered to
#'     high-difference candidates (selection bias). It is only valid when
#'     metilene is run with `min_diff = 0` so every tested segment is emitted.}
#' }
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
#' @return data.table of significant DMRs with columns: chr, start, end,
#'   q_value (the chosen one), mean_diff, n_sites, p_mwu, q_metilene, q_bh.
call_dmrs_metilene <- function(split_data,
                               out_dir,
                               metilene_path = "metilene",
                               min_cpg = 5,
                               min_diff = 0.1,
                               metilene_max_dist = 1500,
                               q_cutoff = 0.05,
                               q_source = c("metilene", "BH")) {
  q_source <- match.arg(q_source)

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

  cmd <- sprintf(
    "%s -a PseudoA -b PseudoB -m %d -d %f -M %d %s > %s 2>/dev/null",
    metilene_path, min_cpg, min_diff, metilene_max_dist, merged_input, out_file
  )

  message("Running metilene: ", cmd)
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("metilene exited with code ", exit_code)

    # Try without the 2>/dev/null to see errors
    cmd_debug <- sprintf(
      "%s -a PseudoA -b PseudoB -m %d -d %f -M %d %s",
      metilene_path, min_cpg, min_diff, metilene_max_dist, merged_input
    )
    system(cmd_debug)
  }

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
    dmrs[, q_bh := p.adjust(p_mwu, method = "BH")]
  } else if (ncol(dmrs) >= 6) {
    warning("metilene output has <8 columns; only column-4 q available")
    setnames(dmrs, 1:6, c("chr", "start", "end", "q_metilene", "mean_diff",
                          "n_sites"))
    dmrs[, `:=`(p_mwu = NA_real_, q_bh = NA_real_)]
  } else {
    return(empty)
  }

  # Significance q-value. Default = metilene's column-4 q (valid genome-wide
  # correction). "BH" uses the R-recomputed value, valid only with min_diff = 0.
  dmrs[, q_value := if (q_source == "BH") q_bh else q_metilene]

  n_total <- nrow(dmrs)
  dmrs <- dmrs[is.finite(q_value) & q_value < q_cutoff]
  dmrs <- dmrs[, .(chr, start, end, q_value, mean_diff, n_sites,
                   p_mwu, q_metilene, q_bh)]
  message(sprintf("metilene: %d candidate DMRs, %d significant (%s q < %.3g)",
                  n_total, nrow(dmrs), q_source, q_cutoff))
  return(dmrs)
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
                              min_coverage = 5
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


