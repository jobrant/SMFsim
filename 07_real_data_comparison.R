#!/usr/bin/env Rscript
# =============================================================================
# 07_real_data_comparison.R
# Compare normalization methods on real biological data (PC3 vs PrEC)
#
# Unlike the simulation, there is no ground truth here. Instead we evaluate:
#   1. DMR counts per method (overcorrection vs undercorrection)
#   2. DMR overlap between methods (Venn/UpSet)
#   3. Effect size distributions per method
#   4. Genomic feature enrichment of method-specific DMRs
#   5. (Optional) Concordance with orthogonal ATAC-seq data
# =============================================================================

library(data.table)
library(ggplot2)
library(MAPitNorm)

# Data loading ------------------------------------------------------------

#' Load and prepare real data for two-group comparison
#'
#' @param config List with data_dir, sample_sheet, and group names.
#' @param group_A Character: name of group A (e.g., "PrEC").
#' @param group_B Character: name of group B (e.g., "PC3").
#' @return List with:
#'   \item{group_A}{Named list of sample data.tables}
#'   \item{group_B}{Named list of sample data.tables}
#'   \item{group_names}{Character vector c(group_A, group_B)}
load_real_data <- function(config, group_A = "PrEC", group_B = "PC3") {
  message("\n", strrep("=", 60))
  message(sprintf("Loading real data: %s vs %s", group_A, group_B))
  message(strrep("=", 60))

  # Load all samples from both groups
  all_samples <- load_data(
    dir_path = config$data_dir,
    sample_sheet = config$sample_sheet,
    groups = c(group_A, group_B)
  )

  message(sprintf("Loaded %d total samples", length(all_samples)))

  # Find shared sites
  all_samples <- find_shared_sites(all_samples, filter = TRUE)
  message(sprintf("Shared sites: %d", nrow(all_samples[[1]])))

  # Split into groups based on sample metadata
  metadata <- attr(all_samples, "sample_metadata")

  idx_A <- which(metadata$group_id == group_A)
  idx_B <- which(metadata$group_id == group_B)

  samples_A <- all_samples[idx_A]
  samples_B <- all_samples[idx_B]

  names(samples_A) <- metadata$sample_id[idx_A]
  names(samples_B) <- metadata$sample_id[idx_B]

  message(sprintf("%s: %d samples", group_A, length(samples_A)))
  message(sprintf("%s: %d samples", group_B, length(samples_B)))

  return(list(
    group_A = samples_A,
    group_B = samples_B,
    group_names = c(group_A, group_B),
    metadata = metadata
  ))
}


# Method runners ----------------------------------------------------------

#' Format real data as pseudo_groups-style list for method runners
#'
#' Adapts real two-group data to the PseudoA/PseudoB structure expected by
#' the existing method runners, preserving real group names for metilene.
#'
#' @param real_data List from load_real_data().
#' @return List with GroupA, GroupB (using real names) compatible with methods.
format_real_data <- function(real_data) {
  group_A_name <- real_data$group_names[1]
  group_B_name <- real_data$group_names[2]

  result <- list()
  result[[group_A_name]] <- real_data$group_A
  result[[group_B_name]] <- real_data$group_B

  return(result)
}


#' Run all normalization methods on real two-group data
#'
#' @param real_data List from load_real_data().
#' @param methods Character vector of methods to run.
#' @param mapitnorm_params Named list of MAPitNorm parameters.
#' @return Named list of results per method, each with $data (split by group).
run_all_methods_real <- function(real_data,
                                 methods = c("raw", "downsampled",
                                             "MAPitNorm", "ComBatMet"),
                                 mapitnorm_params = list(
                                   within_alpha = 0.3,
                                   between_alpha = 0.5,
                                   min_coverage = 5
                                 )) {
  group_A_name <- real_data$group_names[1]
  group_B_name <- real_data$group_names[2]

  results <- list()

  for (m in methods) {
    message(sprintf("\n=== Running method: %s ===", m))
    t0 <- Sys.time()

    # Deep copy for each method
    data_A <- lapply(real_data$group_A, data.table::copy)
    data_B <- lapply(real_data$group_B, data.table::copy)

    result <- tryCatch({
      switch(m,
        raw = {
          message("Method: Raw (no normalization)")
          split_data <- list()
          split_data[[group_A_name]] <- data_A
          split_data[[group_B_name]] <- data_B
          list(method = "raw", data = split_data)
        },

        downsampled = {
          message("Method: Downsampling to minimum coverage")
          all_samples <- c(data_A, data_B)
          n_sites <- nrow(all_samples[[1]])
          min_cov <- rep(Inf, n_sites)
          for (s in all_samples) min_cov <- pmin(min_cov, s$cov)
          min_cov <- as.integer(pmax(0, min_cov))

          downsampled <- lapply(all_samples, function(dt) {
            out <- copy(dt)
            new_mc <- rhyper(
              nn = n_sites,
              m  = pmax(0L, as.integer(out$mc)),
              n  = pmax(0L, as.integer(out$cov - out$mc)),
              k  = min_cov
            )
            out[, `:=`(mc = new_mc, cov = min_cov,
                       rate = ifelse(min_cov > 0, new_mc / min_cov, 0))]
            return(out)
          })
          names(downsampled) <- names(all_samples)
          n_A <- length(data_A)

          split_data <- list()
          split_data[[group_A_name]] <- downsampled[seq_len(n_A)]
          split_data[[group_B_name]] <- downsampled[seq(n_A + 1, length(downsampled))]
          message(sprintf("Downsampled to median min coverage: %d", median(min_cov)))
          list(method = "downsampled", data = split_data)
        },

        MAPitNorm = {
          message("Method: MAPitNorm")
          # Build nested list with metadata for MAPitNorm
          all_samples <- c(data_A, data_B)
          sample_ids <- names(all_samples)
          group_ids <- c(rep(group_A_name, length(data_A)),
                         rep(group_B_name, length(data_B)))

          metadata <- data.table(
            group_id = group_ids,
            replicate = seq_along(sample_ids),
            sample_id = sample_ids,
            file_name = paste0(sample_ids, ".tsv.gz")
          )
          attr(all_samples, "sample_metadata") <- metadata

          split_data <- MAPitNorm::split_by_groups(all_samples)
          normalized <- MAPitNorm::normalize_methylation_data(
            data_list = split_data,
            do_coverage_norm = TRUE,
            normalize_rates = TRUE,
            coverage_between_groups = FALSE,
            rate_within_groups = TRUE,
            rate_between_groups = FALSE,
            within_alpha = mapitnorm_params$within_alpha,
            between_alpha = mapitnorm_params$between_alpha,
            min_coverage = mapitnorm_params$min_coverage
          )
          list(method = "MAPitNorm", data = normalized)
        },

        ComBatMet = {
          message("Method: ComBatMet (manual mean-only)")
          all_samples <- lapply(c(data_A, data_B), data.table::copy)
          names(all_samples) <- c(names(data_A), names(data_B))
          n_A <- length(data_A)
          nc <- length(all_samples)

          rate_matrix <- do.call(cbind, lapply(all_samples,
                                               function(dt) as.numeric(dt$rate)))
          colnames(rate_matrix) <- names(all_samples)

          rate_matrix[rate_matrix < 1e-6] <- 1e-6
          rate_matrix[rate_matrix > (1 - 1e-6)] <- 1 - 1e-6

          logit_rates <- rate_matrix
          logit_rates[] <- log(rate_matrix / (1 - rate_matrix))

          grand_mean <- rowMeans(logit_rates)
          batch_A_mean <- rowMeans(logit_rates[, 1:n_A, drop = FALSE])
          batch_B_mean <- rowMeans(logit_rates[, (n_A+1):nc, drop = FALSE])

          adjusted <- logit_rates
          adjusted[, 1:n_A] <- adjusted[, 1:n_A] - batch_A_mean + grand_mean
          adjusted[, (n_A+1):nc] <- adjusted[, (n_A+1):nc] - batch_B_mean + grand_mean

          adjusted_rates <- adjusted
          adjusted_rates[] <- 1 / (1 + exp(-adjusted))

          result_samples <- lapply(seq_len(nc), function(i) {
            out <- copy(all_samples[[i]])
            out[, rate := adjusted_rates[, i]]
            out[, mc := as.integer(round(cov * rate))]
            return(out)
          })
          names(result_samples) <- names(all_samples)

          split_data <- list()
          split_data[[group_A_name]] <- result_samples[seq_len(n_A)]
          split_data[[group_B_name]] <- result_samples[seq(n_A + 1, nc)]
          list(method = "ComBatMet", data = split_data)
        },

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


# DMR calling -------------------------------------------------------------

#' Call DMRs with metilene using real group names
#'
#' @param split_data Named list (by group) of named lists (by sample) of data.tables.
#' @param out_dir Output directory.
#' @param metilene_path Path to metilene binary.
#' @param min_cpg Minimum sites per DMR.
#' @param min_diff Minimum mean difference.
#' @return data.table of DMRs.
call_dmrs_real <- function(split_data,
                            out_dir,
                            metilene_path = "metilene",
                            min_cpg = 3,
                            min_diff = 0.1) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  group_names <- names(split_data)
  group_A_name <- group_names[1]
  group_B_name <- group_names[2]
  n_A <- length(split_data[[group_A_name]])
  n_B <- length(split_data[[group_B_name]])

  # Build merged bedGraph with proper header
  all_dts <- list()

  for (i in seq_len(n_A)) {
    dt <- split_data[[group_A_name]][[i]][, .(chr, pos, rate)]
    setnames(dt, "rate", paste0(group_A_name, "_", i))
    all_dts[[length(all_dts) + 1]] <- dt
  }
  for (i in seq_len(n_B)) {
    dt <- split_data[[group_B_name]][[i]][, .(chr, pos, rate)]
    setnames(dt, "rate", paste0(group_B_name, "_", i))
    all_dts[[length(all_dts) + 1]] <- dt
  }

  # Merge all on chr + pos
  merged <- Reduce(function(x, y) merge(x, y, by = c("chr", "pos")), all_dts)
  merged[, start := pos - 1L]
  merged[, end := pos]

  a_cols <- paste0(group_A_name, "_", seq_len(n_A))
  b_cols <- paste0(group_B_name, "_", seq_len(n_B))

  out_bed <- merged[, c("chr", "start", "end", a_cols, b_cols), with = FALSE]
  setorder(out_bed, chr, start)

  input_path <- file.path(out_dir, "merged_input.bed")
  fwrite(out_bed, input_path, sep = "\t", col.names = TRUE)

  # Run metilene
  out_file <- file.path(out_dir, "metilene_dmrs.bed")
  cmd <- sprintf(
    "%s -a %s -b %s -m %d -d %f %s > %s 2>/dev/null",
    metilene_path, group_A_name, group_B_name,
    min_cpg, min_diff, input_path, out_file
  )

  message("Running metilene: ", cmd)
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("metilene exited with code ", exit_code)
  }

  if (!file.exists(out_file) || file.size(out_file) == 0) {
    message("No DMRs called by metilene")
    return(data.table(
      chr = character(), start = integer(), end = integer(),
      q_value = numeric(), mean_diff = numeric(), n_sites = integer()
    ))
  }

  dmrs <- fread(out_file, header = FALSE)
  if (ncol(dmrs) >= 6) {
    setnames(dmrs, 1:6, c("chr", "start", "end", "q_value", "mean_diff", "n_sites"))
    dmrs <- dmrs[, .(chr, start, end, q_value, mean_diff, n_sites)]
  }

  message(sprintf("metilene called %d DMRs", nrow(dmrs)))
  return(dmrs)
}


# Comparison metrics ------------------------------------------------------

#' Compute pairwise DMR overlap between methods
#'
#' @param dmr_list Named list of DMR data.tables (one per method).
#' @param min_overlap_bp Minimum overlap to count as shared.
#' @return data.table with pairwise overlap counts.
compute_dmr_overlap <- function(dmr_list, min_overlap_bp = 1) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges required for overlap computation")
  }

  methods <- names(dmr_list)

  # Convert to GRanges
  gr_list <- lapply(dmr_list, function(dmrs) {
    if (nrow(dmrs) == 0) return(GenomicRanges::GRanges())
    GenomicRanges::GRanges(
      seqnames = dmrs$chr,
      ranges = IRanges::IRanges(start = dmrs$start, end = dmrs$end)
    )
  })

  # Pairwise overlaps
  overlap_dt <- data.table(
    method_A = character(), method_B = character(),
    n_A = integer(), n_B = integer(),
    shared = integer(), only_A = integer(), only_B = integer()
  )

  for (i in seq_along(methods)) {
    for (j in seq_along(methods)) {
      if (i >= j) next
      m_A <- methods[i]
      m_B <- methods[j]

      if (length(gr_list[[m_A]]) == 0 || length(gr_list[[m_B]]) == 0) {
        overlap_dt <- rbind(overlap_dt, data.table(
          method_A = m_A, method_B = m_B,
          n_A = length(gr_list[[m_A]]), n_B = length(gr_list[[m_B]]),
          shared = 0L,
          only_A = length(gr_list[[m_A]]),
          only_B = length(gr_list[[m_B]])
        ))
        next
      }

      hits <- GenomicRanges::findOverlaps(
        gr_list[[m_A]], gr_list[[m_B]],
        minoverlap = min_overlap_bp
      )

      n_shared_A <- length(unique(S4Vectors::queryHits(hits)))
      n_shared_B <- length(unique(S4Vectors::subjectHits(hits)))

      overlap_dt <- rbind(overlap_dt, data.table(
        method_A = m_A, method_B = m_B,
        n_A = length(gr_list[[m_A]]),
        n_B = length(gr_list[[m_B]]),
        shared = n_shared_A,
        only_A = length(gr_list[[m_A]]) - n_shared_A,
        only_B = length(gr_list[[m_B]]) - n_shared_B
      ))
    }
  }

  return(overlap_dt)
}


#' Classify DMRs as shared across methods or method-specific
#'
#' @param dmr_list Named list of DMR data.tables.
#' @param reference_method Method to use as reference (default: "raw").
#' @param min_overlap_bp Minimum overlap.
#' @return List with:
#'   \item{shared}{DMRs found by all methods}
#'   \item{method_specific}{Per-method DMRs not found by reference}
#'   \item{lost}{Reference DMRs lost by each method}
classify_method_dmrs <- function(dmr_list,
                                  reference_method = "raw",
                                  min_overlap_bp = 1) {
  if (!reference_method %in% names(dmr_list)) {
    stop("Reference method '", reference_method, "' not found in dmr_list")
  }

  ref_dmrs <- dmr_list[[reference_method]]
  if (nrow(ref_dmrs) == 0) {
    warning("Reference method has no DMRs")
    return(list(shared = data.table(), method_specific = list(), lost = list()))
  }

  ref_gr <- GenomicRanges::GRanges(
    seqnames = ref_dmrs$chr,
    ranges = IRanges::IRanges(start = ref_dmrs$start, end = ref_dmrs$end)
  )

  other_methods <- setdiff(names(dmr_list), reference_method)
  method_specific <- list()
  lost <- list()

  for (m in other_methods) {
    m_dmrs <- dmr_list[[m]]
    if (nrow(m_dmrs) == 0) {
      method_specific[[m]] <- data.table()
      lost[[m]] <- ref_dmrs
      next
    }

    m_gr <- GenomicRanges::GRanges(
      seqnames = m_dmrs$chr,
      ranges = IRanges::IRanges(start = m_dmrs$start, end = m_dmrs$end)
    )

    # DMRs in m not in reference → method-specific (potentially novel or FP)
    hits_m_to_ref <- GenomicRanges::findOverlaps(m_gr, ref_gr, minoverlap = min_overlap_bp)
    novel_idx <- setdiff(seq_len(nrow(m_dmrs)), unique(S4Vectors::queryHits(hits_m_to_ref)))
    method_specific[[m]] <- m_dmrs[novel_idx]

    # Reference DMRs not in m → lost by this method
    hits_ref_to_m <- GenomicRanges::findOverlaps(ref_gr, m_gr, minoverlap = min_overlap_bp)
    lost_idx <- setdiff(seq_len(nrow(ref_dmrs)), unique(S4Vectors::queryHits(hits_ref_to_m)))
    lost[[m]] <- ref_dmrs[lost_idx]
  }

  return(list(method_specific = method_specific, lost = lost))
}


# Visualization -----------------------------------------------------------

theme_manuscript <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size, color = "grey40"),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80")
    )
}

method_colors <- c(
  "raw"         = "#E41A1C",
  "downsampled" = "#377EB8",
  "MAPitNorm"   = "#4DAF4A",
  "ComBatMet"   = "#984EA3"
)


#' Bar plot of total DMR counts per method
plot_dmr_counts <- function(dmr_list, out_dir = NULL) {
  count_dt <- data.table(
    method = names(dmr_list),
    n_dmrs = sapply(dmr_list, nrow)
  )
  count_dt[, method := factor(method, levels = names(method_colors))]

  p <- ggplot(count_dt, aes(x = method, y = n_dmrs, fill = method)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n_dmrs), vjust = -0.5, size = 4) +
    scale_fill_manual(values = method_colors, guide = "none") +
    labs(
      title = "Total DMRs Called per Method (Real Data)",
      subtitle = "Cancer (PC3) vs Normal (PrEC)",
      x = NULL, y = "Number of DMRs"
    ) +
    theme_manuscript()

  if (!is.null(out_dir)) .save_figure(p, "fig_real_dmr_counts", out_dir)
  return(p)
}


#' Effect size distribution per method
plot_effect_size_distribution <- function(dmr_list, out_dir = NULL) {
  effect_dt <- rbindlist(lapply(names(dmr_list), function(m) {
    dt <- dmr_list[[m]]
    if (nrow(dt) == 0) return(data.table())
    data.table(method = m, mean_diff = dt$mean_diff)
  }))

  if (nrow(effect_dt) == 0) {
    message("No DMRs to plot effect sizes")
    return(NULL)
  }

  effect_dt[, method := factor(method, levels = names(method_colors))]

  p <- ggplot(effect_dt, aes(x = abs(mean_diff), fill = method)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ method, ncol = 2) +
    scale_fill_manual(values = method_colors, guide = "none") +
    labs(
      title = "Effect Size Distribution of Called DMRs",
      subtitle = "Absolute mean rate difference between groups",
      x = "|Mean Difference|", y = "Density"
    ) +
    theme_manuscript()

  if (!is.null(out_dir)) .save_figure(p, "fig_real_effect_sizes", out_dir)
  return(p)
}


#' Effect size distribution: lost vs retained DMRs
#'
#' Shows whether the DMRs lost by a method after normalization tend to be
#' small-effect (likely artifacts) or large-effect (likely biology).
plot_lost_dmr_effects <- function(dmr_list, classification, out_dir = NULL) {
  plots <- list()

  for (m in names(classification$lost)) {
    lost <- classification$lost[[m]]
    if (nrow(lost) == 0) next

    # Find retained (reference DMRs NOT lost)
    ref_dmrs <- dmr_list[["raw"]]
    ref_gr <- GenomicRanges::GRanges(
      seqnames = ref_dmrs$chr,
      ranges = IRanges::IRanges(start = ref_dmrs$start, end = ref_dmrs$end)
    )
    lost_gr <- GenomicRanges::GRanges(
      seqnames = lost$chr,
      ranges = IRanges::IRanges(start = lost$start, end = lost$end)
    )
    hits <- GenomicRanges::findOverlaps(ref_gr, lost_gr)
    retained_idx <- setdiff(seq_len(nrow(ref_dmrs)), unique(S4Vectors::queryHits(hits)))

    plot_dt <- rbind(
      data.table(status = "Retained", mean_diff = ref_dmrs$mean_diff[retained_idx]),
      data.table(status = "Lost", mean_diff = lost$mean_diff)
    )

    p <- ggplot(plot_dt, aes(x = abs(mean_diff), fill = status)) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(values = c("Retained" = "#4DAF4A", "Lost" = "#E41A1C")) +
      labs(
        title = sprintf("DMRs Lost by %s vs Raw", m),
        subtitle = "Small-effect lost DMRs suggest artifact removal",
        x = "|Mean Difference|", y = "Density", fill = NULL
      ) +
      theme_manuscript()

    if (!is.null(out_dir)) .save_figure(p, paste0("fig_real_lost_", tolower(m)), out_dir)
    plots[[m]] <- p
  }

  return(plots)
}


#' Overlap summary heatmap
plot_overlap_heatmap <- function(overlap_dt, out_dir = NULL) {
  # Build a symmetric matrix of shared counts
  methods <- unique(c(overlap_dt$method_A, overlap_dt$method_B))

  # Make full pairwise table including self-comparisons
  full_dt <- rbind(
    overlap_dt[, .(method_A, method_B, shared)],
    overlap_dt[, .(method_A = method_B, method_B = method_A, shared)]
  )

  full_dt[, method_A := factor(method_A, levels = names(method_colors))]
  full_dt[, method_B := factor(method_B, levels = names(method_colors))]

  p <- ggplot(full_dt, aes(x = method_A, y = method_B, fill = shared)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = shared), size = 4) +
    scale_fill_gradient(low = "#FFFFBF", high = "#4DAF4A", na.value = "grey90") +
    labs(
      title = "Pairwise DMR Overlap Between Methods",
      x = NULL, y = NULL, fill = "Shared\nDMRs"
    ) +
    theme_manuscript() +
    theme(legend.position = "right")

  if (!is.null(out_dir)) .save_figure(p, "fig_real_overlap_heatmap", out_dir)
  return(p)
}


.save_figure <- function(plot, name, out_dir, width = 8, height = 6) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, device = "pdf")
  ggsave(file.path(out_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, device = "png")
  message(sprintf("Saved: %s", name))
}


# Main pipeline -----------------------------------------------------------

#' Run the full real data comparison
#'
#' @param config Configuration list (data_dir, sample_sheet, etc.).
#' @param group_A Name of the normal/reference group.
#' @param group_B Name of the cancer/test group.
#' @param out_dir Output directory.
#' @return List of all results and figures.
run_real_data_comparison <- function(config,
                                     group_A = "PrEC",
                                     group_B = "PC3",
                                     out_dir = "results/real_data") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: Load data
  real_data <- load_real_data(config, group_A, group_B)

  # Step 2: Run all methods
  method_results <- run_all_methods_real(real_data)

  # Step 3: Call DMRs for each method
  message("\n", strrep("=", 60))
  message("Calling DMRs with metilene")
  message(strrep("=", 60))

  dmr_list <- list()
  for (m_name in names(method_results)) {
    message(sprintf("\nDMR calling for: %s", m_name))
    dmr_dir <- file.path(out_dir, "dmrs", m_name)

    dmrs <- tryCatch({
      call_dmrs_real(
        split_data = method_results[[m_name]]$data,
        out_dir = dmr_dir,
        metilene_path = config$metilene_path %||% "metilene",
        min_cpg = config$metilene_min_cpg %||% 3,
        min_diff = config$metilene_min_diff %||% 0.1
      )
    }, error = function(e) {
      warning(sprintf("DMR calling failed for %s: %s", m_name, e$message))
      data.table(chr = character(), start = integer(), end = integer(),
                 q_value = numeric(), mean_diff = numeric(), n_sites = integer())
    })

    dmr_list[[m_name]] <- dmrs
  }

  # Step 4: Summary statistics
  message("\n", strrep("=", 60))
  message("DMR Summary")
  message(strrep("=", 60))

  summary_dt <- data.table(
    method = names(dmr_list),
    n_dmrs = sapply(dmr_list, nrow),
    median_effect = sapply(dmr_list, function(dt) {
      if (nrow(dt) > 0) round(median(abs(dt$mean_diff)), 4) else NA_real_
    }),
    mean_effect = sapply(dmr_list, function(dt) {
      if (nrow(dt) > 0) round(mean(abs(dt$mean_diff)), 4) else NA_real_
    })
  )
  print(summary_dt)
  fwrite(summary_dt, file.path(out_dir, "dmr_summary.csv"))

  # Step 5: Overlap analysis
  message("\nComputing DMR overlaps...")
  overlap_dt <- compute_dmr_overlap(dmr_list)
  fwrite(overlap_dt, file.path(out_dir, "dmr_overlap.csv"))
  print(overlap_dt)

  # Step 6: Classify DMRs relative to raw
  classification <- classify_method_dmrs(dmr_list, reference_method = "raw")

  message("\nDMRs lost relative to raw:")
  for (m in names(classification$lost)) {
    message(sprintf("  %s: %d DMRs lost, %d method-specific",
                    m, nrow(classification$lost[[m]]),
                    nrow(classification$method_specific[[m]])))
  }

  # Step 7: Generate figures
  fig_dir <- file.path(out_dir, "figures")
  message("\nGenerating figures...")

  fig_counts <- plot_dmr_counts(dmr_list, fig_dir)
  fig_effects <- plot_effect_size_distribution(dmr_list, fig_dir)
  fig_overlap <- plot_overlap_heatmap(overlap_dt, fig_dir)
  fig_lost <- plot_lost_dmr_effects(dmr_list, classification, fig_dir)

  # Save everything
  saveRDS(list(
    dmr_list = dmr_list,
    summary = summary_dt,
    overlap = overlap_dt,
    classification = classification,
    config = config
  ), file.path(out_dir, "real_data_results.rds"))

  message("\n", strrep("=", 60))
  message("Real data comparison complete")
  message(sprintf("Results saved to: %s", out_dir))
  message(strrep("=", 60))

  return(invisible(list(
    dmr_list = dmr_list,
    summary = summary_dt,
    overlap = overlap_dt,
    classification = classification,
    figures = list(counts = fig_counts, effects = fig_effects,
                   overlap = fig_overlap, lost = fig_lost)
  )))
}

# Convenience null-coalescing operator ------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
