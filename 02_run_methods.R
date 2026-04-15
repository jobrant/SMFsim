#!/usr/bin/env Rscript
# =============================================================================
# 02_run_methods.R
# Apply normalization methods to simulated data and call DMRs
# =============================================================================

library(data.table)

# ── Format helpers ───────────────────────────────────────────────────────────

#' Convert pseudo-group data to the format expected by MAPitNorm
#'
#' Builds a nested list with sample_metadata attribute, mimicking load_data() output.
#'
#' @param pseudo_groups List from create_pseudo_groups() with PseudoA and PseudoB.
#' @return List structured as load_data() output (flat sample list with metadata).
format_for_mapinorm <- function(pseudo_groups) {
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


# ── Method: Raw (no normalization) ───────────────────────────────────────────

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


# ── Method: Downsampling ─────────────────────────────────────────────────────

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


# ── Method: MAPitNorm ────────────────────────────────────────────────────────

#' Run MAPitNorm normalization pipeline
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param within_alpha Alpha for within-group rate normalization.
#' @param between_alpha Alpha for between-group rate normalization.
#' @param min_coverage Minimum coverage threshold.
#' @return Named list with method = "MAPitNorm" and normalized split data.
run_mapitnorm <- function(pseudo_groups,
                           within_alpha = 0.3,
                           between_alpha = 0.5,
                           min_coverage = 5) {
  message("Method: MAPitNorm")

  # Format as MAPitNorm expects
  formatted <- format_for_mapinorm(pseudo_groups)

  # Split by groups
  split_data <- MAPitNorm::split_by_groups(formatted)

  # Run full normalization
  normalized <- MAPitNorm::normalize_methylation_data(
    data_list = split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    coverage_between_groups = FALSE,
    rate_within_groups = TRUE,
    rate_between_groups = TRUE,
    within_alpha = within_alpha,
    between_alpha = between_alpha,
    min_coverage = min_coverage
  )

  return(list(method = "MAPitNorm", data = normalized))
}


# ── Method: ComBatMet (via external call) ─────────────────────────────────────

#' Run ComBatMet normalization
#'
#' ComBatMet uses a ComBat-based approach on methylation beta values.
#' This wrapper prepares a matrix, runs ComBat, and converts back.
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @return Named list with method = "ComBatMet" and normalized split data.
run_combatmet <- function(pseudo_groups) {
  message("Method: ComBatMet")

  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for ComBatMet. ",
         "Install with: BiocManager::install('sva')")
  }

  all_samples <- c(pseudo_groups$PseudoA, pseudo_groups$PseudoB)
  n_A <- length(pseudo_groups$PseudoA)

  # Build rate matrix (sites × samples)
  rate_matrix <- do.call(cbind, lapply(all_samples, function(dt) dt$rate))
  colnames(rate_matrix) <- names(all_samples)

  # Build batch variable — here we treat each sample as its own "batch"
  # (ComBatMet models per-sample efficiency as a batch effect)
  # Alternatively, treat the two pseudo-groups as batches
  batch <- ifelse(grepl("^PseudoA", names(all_samples)), 1, 2)

  # Biological group (same as batch in null simulation, but needed for spike-in)
  group <- factor(ifelse(grepl("^PseudoA", names(all_samples)), "A", "B"))
  mod <- model.matrix(~ group)

  # Apply ComBat
  # logit-transform rates (ComBat assumes ~normal data)
  # Add small offset to avoid log(0)
  eps <- 1e-6
  rate_matrix_bounded <- pmin(1 - eps, pmax(eps, rate_matrix))
  logit_rates <- log(rate_matrix_bounded / (1 - rate_matrix_bounded))

  adjusted <- tryCatch({
    sva::ComBat(
      dat = logit_rates,
      batch = batch,
      mod = mod
    )
  }, error = function(e) {
    message("ComBat failed: ", e$message)
    message("Falling back to ComBat with mean.only = TRUE")
    sva::ComBat(
      dat = logit_rates,
      batch = batch,
      mod = mod,
      mean.only = TRUE
    )
  })

  # Inverse logit to get back to [0, 1]
  adjusted_rates <- 1 / (1 + exp(-adjusted))

  # Rebuild sample data.tables
  result_samples <- lapply(seq_along(all_samples), function(i) {
    out <- copy(all_samples[[i]])
    out[, rate := adjusted_rates[, i]]
    out[, mc := as.integer(round(cov * rate))]
    return(out)
  })
  names(result_samples) <- names(all_samples)

  split_data <- list(
    PseudoA = result_samples[seq_len(n_A)],
    PseudoB = result_samples[seq(n_A + 1, length(result_samples))]
  )

  return(list(method = "ComBatMet", data = split_data))
}


# ── DMR calling with metilene ────────────────────────────────────────────────

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
#' @param split_data Named list with PseudoA and PseudoB.
#' @param out_dir Output directory.
#' @param metilene_path Path to metilene binary (default: assumes in PATH).
#' @param min_cpg Minimum number of CpGs/sites in a DMR (default: 5).
#' @param min_diff Minimum mean difference to report (default: 0.1).
#' @return data.table of called DMRs with columns: chr, start, end, q_value, mean_diff, n_sites
call_dmrs_metilene <- function(split_data,
                                out_dir,
                                metilene_path = "metilene",
                                min_cpg = 5,
                                min_diff = 0.1) {

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
    "%s -a PseudoA -b PseudoB -m %d -d %f %s > %s 2>/dev/null",
    metilene_path, min_cpg, min_diff, merged_input, out_file
  )

  message("Running metilene: ", cmd)
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("metilene exited with code ", exit_code)

    # Try without the 2>/dev/null to see errors
    cmd_debug <- sprintf(
      "%s -a PseudoA -b PseudoB -m %d -d %f %s",
      metilene_path, min_cpg, min_diff, merged_input
    )
    system(cmd_debug)
  }

  # Parse metilene output
  if (!file.exists(out_file) || file.size(out_file) == 0) {
    message("No DMRs called by metilene")
    return(data.table(
      chr = character(), start = integer(), end = integer(),
      q_value = numeric(), mean_diff = numeric(), n_sites = integer()
    ))
  }

  dmrs <- fread(out_file, header = FALSE)
  # metilene output: chr, start, end, q-value, mean_diff, #CpGs, ...
  if (ncol(dmrs) >= 6) {
    setnames(dmrs, 1:6, c("chr", "start", "end", "q_value", "mean_diff", "n_sites"))
    dmrs <- dmrs[, .(chr, start, end, q_value, mean_diff, n_sites)]
  }

  message(sprintf("metilene called %d DMRs", nrow(dmrs)))
  return(dmrs)
}


#' Merge two group bedGraphs into metilene's expected format
#' @keywords internal
.merge_group_bedgraphs <- function(path_A, path_B, split_data, out_path) {
  dt_A <- fread(path_A, header = FALSE)
  dt_B <- fread(path_B, header = FALSE)

  n_A <- length(split_data$PseudoA)
  n_B <- length(split_data$PseudoB)

  # Metilene format: chr start end groupA_rep1 groupA_rep2 ... groupB_rep1 groupB_rep2 ...
  # dt_A has: chr, start, end, rate_1, rate_2, ...
  # dt_B has: chr, start, end, rate_1, rate_2, ...
  # Merge on chr, start, end

  # Rename rate columns to avoid collision
  rate_cols_A <- paste0("V", 4:(3 + n_A))
  rate_cols_B <- paste0("V", 4:(3 + n_B))

  setnames(dt_A, rate_cols_A, paste0("A_", seq_len(n_A)))
  setnames(dt_B, rate_cols_B, paste0("B_", seq_len(n_B)))

  merged <- merge(
    dt_A, dt_B,
    by = c("V1", "V2", "V3"),
    sort = FALSE
  )

  # Reorder: chr, start, end, all A reps, all B reps
  a_cols <- paste0("A_", seq_len(n_A))
  b_cols <- paste0("B_", seq_len(n_B))
  merged <- merged[, c("V1", "V2", "V3", a_cols, b_cols), with = FALSE]

  setorder(merged, V1, V2)
  fwrite(merged, out_path, sep = "\t", col.names = FALSE)

  return(out_path)
}


# ── Master method runner ─────────────────────────────────────────────────────

#' Run all normalization methods on pseudo-group data
#'
#' @param pseudo_groups List from create_pseudo_groups().
#' @param methods Character vector of methods to run.
#'   Options: "raw", "downsampled", "MAPitNorm", "ComBatMet"
#' @param mapitnorm_params Named list of MAPitNorm parameters.
#' @return Named list of results, each with $method and $data.
run_all_methods <- function(pseudo_groups,
                             methods = c("raw", "downsampled", "MAPitNorm", "ComBatMet"),
                             mapitnorm_params = list(
                               within_alpha = 0.3,
                               between_alpha = 0.5,
                               min_coverage = 5
                             )) {
  results <- list()

  for (m in methods) {
    message(sprintf("\n=== Running method: %s ===", m))
    t0 <- Sys.time()

    result <- tryCatch({
      switch(m,
        raw = run_raw(pseudo_groups),
        downsampled = run_downsampled(pseudo_groups),
        MAPitNorm = do.call(run_mapitnorm, c(list(pseudo_groups = pseudo_groups),
                                              mapitnorm_params)),
        ComBatMet = run_combatmet(pseudo_groups),
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

  return(results)
}
