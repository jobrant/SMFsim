#!/usr/bin/env Rscript

# 03_evaluate.R
# Evaluate DMR calls: classify TP/FP/FN and compute performance metrics

library(data.table)

# DMR classification ------------------------------------------------------

#' Classify called DMRs against ground-truth spike-in regions
#'
#' A called DMR is a true positive (TP) if it overlaps a spike-in region by
#' at least `min_overlap_frac` of the smaller region's width.
#' Spike-in regions without any overlapping DMR call are false negatives (FN).
#' Called DMRs without any overlapping spike-in are false positives (FP).
#'
#' @param called_dmrs data.table of called DMRs (chr, start, end, ...).
#' @param truth_regions data.table of spike-in regions (region_id, chr, start, end, ...).
#'   If NULL, all called DMRs are treated as FP (null simulation).
#' @param min_overlap_bp Minimum overlap in base pairs to count as a match (default: 1).
#' @return List with:
#'   \item{classification}{data.table with each called DMR labeled TP or FP}
#'   \item{fn_regions}{data.table of spike-in regions not recovered (FN)}
#'   \item{summary}{Named vector: TP, FP, FN counts}
classify_dmrs <- function(called_dmrs, truth_regions = NULL, min_overlap_bp = 1) {

  # Null simulation: everything is FP

  if (is.null(truth_regions) || nrow(truth_regions) == 0) {
    n_called <- nrow(called_dmrs)
    classification <- copy(called_dmrs)
    if (n_called > 0) {
      classification[, label := "FP"]
      classification[, matched_region := NA_integer_]
    }
    return(list(
      classification = classification,
      fn_regions = data.table(),
      summary = c(TP = 0L, FP = n_called, FN = 0L)
    ))
  }

  # Spike-in simulation: check overlaps
  if (nrow(called_dmrs) == 0) {
    return(list(
      classification = data.table(),
      fn_regions = truth_regions,
      summary = c(TP = 0L, FP = 0L, FN = nrow(truth_regions))
    ))
  }

  # Use GenomicRanges for overlap detection if available, otherwise interval logic
  if (requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("IRanges", quietly = TRUE)) {
    result <- .classify_with_granges(called_dmrs, truth_regions, min_overlap_bp)
  } else {
    result <- .classify_with_intervals(called_dmrs, truth_regions, min_overlap_bp)
  }

  return(result)
}


#' Classify DMRs using GenomicRanges overlaps
#' @keywords internal
.classify_with_granges <- function(called_dmrs, truth_regions, min_overlap_bp) {
  called_gr <- GenomicRanges::GRanges(
    seqnames = called_dmrs$chr,
    ranges = IRanges::IRanges(start = called_dmrs$start, end = called_dmrs$end)
  )

  truth_gr <- GenomicRanges::GRanges(
    seqnames = truth_regions$chr,
    ranges = IRanges::IRanges(start = truth_regions$start, end = truth_regions$end)
  )

  # Find overlaps
  hits <- GenomicRanges::findOverlaps(called_gr, truth_gr,
                                       minoverlap = min_overlap_bp)

  # Classify called DMRs
  classification <- copy(called_dmrs)
  classification[, label := "FP"]
  classification[, matched_region := NA_integer_]

  tp_called_idx <- unique(S4Vectors::queryHits(hits))
  if (length(tp_called_idx) > 0) {
    classification[tp_called_idx, label := "TP"]

    # Record which truth region was matched
    for (i in tp_called_idx) {
      matched <- S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]
      classification[i, matched_region := truth_regions$region_id[matched[1]]]
    }
  }

  # Identify FN (truth regions not matched)
  matched_truth_idx <- unique(S4Vectors::subjectHits(hits))
  all_truth_idx <- seq_len(nrow(truth_regions))
  fn_idx <- setdiff(all_truth_idx, matched_truth_idx)
  fn_regions <- truth_regions[fn_idx]

  summary_vec <- c(
    TP = sum(classification$label == "TP"),
    FP = sum(classification$label == "FP"),
    FN = length(fn_idx)
  )

  return(list(
    classification = classification,
    fn_regions = fn_regions,
    summary = summary_vec
  ))
}


#' Classify DMRs using simple interval overlap (fallback)
#' @keywords internal
.classify_with_intervals <- function(called_dmrs, truth_regions, min_overlap_bp) {
  classification <- copy(called_dmrs)
  classification[, label := "FP"]
  classification[, matched_region := NA_integer_]

  truth_matched <- rep(FALSE, nrow(truth_regions))

  for (i in seq_len(nrow(called_dmrs))) {
    c_chr <- called_dmrs$chr[i]
    c_start <- called_dmrs$start[i]
    c_end <- called_dmrs$end[i]

    for (j in seq_len(nrow(truth_regions))) {
      if (truth_matched[j] && classification$label[i] == "TP") next

      t_chr <- truth_regions$chr[j]
      t_start <- truth_regions$start[j]
      t_end <- truth_regions$end[j]

      if (c_chr != t_chr) next

      overlap <- min(c_end, t_end) - max(c_start, t_start)
      if (overlap >= min_overlap_bp) {
        classification[i, label := "TP"]
        classification[i, matched_region := truth_regions$region_id[j]]
        truth_matched[j] <- TRUE
        break
      }
    }
  }

  fn_regions <- truth_regions[!truth_matched]

  summary_vec <- c(
    TP = sum(classification$label == "TP"),
    FP = sum(classification$label == "FP"),
    FN = sum(!truth_matched)
  )

  return(list(
    classification = classification,
    fn_regions = fn_regions,
    summary = summary_vec
  ))
}


# Performance metrics -----------------------------------------------------

#' Compute performance metrics from TP/FP/FN counts
#'
#' @param summary_vec Named vector with TP, FP, FN counts.
#' @return Named list of metrics: sensitivity, precision, specificity, F1, FDR.
compute_metrics <- function(summary_vec) {
  TP <- summary_vec["TP"]
  FP <- summary_vec["FP"]
  FN <- summary_vec["FN"]

  sensitivity <- if (TP + FN > 0) TP / (TP + FN) else NA_real_  # recall / TPR
  precision   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_  # PPV
  FDR         <- if (TP + FP > 0) FP / (TP + FP) else NA_real_  # false discovery rate
  F1          <- if (!is.na(sensitivity) && !is.na(precision) && (sensitivity + precision) > 0) {
    2 * sensitivity * precision / (sensitivity + precision)
  } else {
    NA_real_
  }

  return(list(
    TP = TP, FP = FP, FN = FN,
    sensitivity = round(sensitivity, 4),
    precision = round(precision, 4),
    FDR = round(FDR, 4),
    F1 = round(F1, 4),
    total_called = TP + FP
  ))
}


# Aggregation across scenarios --------------------------------------------

#' Compile results from multiple scenarios and methods into a summary table
#'
#' @param all_results List of lists. Outer: scenario names. Inner: method names.
#'   Each leaf should be the output of classify_dmrs().
#' @param scenario_labels Named character vector mapping scenario names to labels.
#' @return data.table with one row per scenario × method combination.
compile_results <- function(all_results, scenario_labels = NULL) {
  rows <- list()

  for (scenario in names(all_results)) {
    for (method in names(all_results[[scenario]])) {
      res <- all_results[[scenario]][[method]]

      if (is.null(res)) {
        rows[[length(rows) + 1]] <- data.table(
          scenario = scenario,
          method = method,
          TP = NA_integer_, FP = NA_integer_, FN = NA_integer_,
          sensitivity = NA_real_, precision = NA_real_,
          FDR = NA_real_, F1 = NA_real_, total_called = NA_integer_
        )
        next
      }

      metrics <- compute_metrics(res$summary)

      rows[[length(rows) + 1]] <- data.table(
        scenario = scenario,
        method = method,
        TP = metrics$TP,
        FP = metrics$FP,
        FN = metrics$FN,
        sensitivity = metrics$sensitivity,
        precision = metrics$precision,
        FDR = metrics$FDR,
        F1 = metrics$F1,
        total_called = metrics$total_called
      )
    }
  }

  summary_dt <- rbindlist(rows)

  if (!is.null(scenario_labels)) {
    summary_dt[, scenario_label := scenario_labels[scenario]]
  }

  return(summary_dt)
}


#' Compile results specifically for null simulations (FP-only)
#'
#' @param null_results Named list (by scenario) of named lists (by method),
#'   each containing classify_dmrs() output.
#' @return data.table with FP counts per scenario × method.
compile_null_results <- function(null_results) {
  rows <- list()

  for (scenario in names(null_results)) {
    for (method in names(null_results[[scenario]])) {
      res <- null_results[[scenario]][[method]]

      n_fp <- if (!is.null(res)) res$summary["FP"] else NA_integer_

      rows[[length(rows) + 1]] <- data.table(
        scenario = scenario,
        method = method,
        FP = n_fp
      )
    }
  }

  return(rbindlist(rows))
}
