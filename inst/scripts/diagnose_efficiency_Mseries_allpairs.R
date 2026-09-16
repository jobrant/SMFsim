#!/usr/bin/env Rscript

# diagnose_efficiency_Mseries_allpairs.R
#
# Run the conversion-efficiency diagnostic across ALL pairs of the M-series
# (M1, M2, M3, M4) to find which pairs carry a systematic between-group
# efficiency offset (artifact candidate -- a uniform global shift, the
# "parallel-offset cloud") versus locus-specific biology (like M1 vs M2, which
# had no offset). Pairs with a genuine offset are the candidates for showing
# real-data correction efficacy; pairs without one are within-group/specificity
# cases.
#
# Replicate counts: M1 = 4, M2 = M3 = M4 = 2.
#
# For each pair this computes:
#   - per-sample global GpC rate (coverage-weighted = efficiency proxy)
#   - within-group SD vs between-group offset, and their ratio
#   - per-bin M-vs-M comparison: systematic offset vs locus-specific scatter
# and writes per-pair figures plus a MASTER summary table + ranking figure.
#
# Usage (from package root or one level up):
#   Rscript inst/scripts/diagnose_efficiency_Mseries_allpairs.R
#
# Interpretation of the master summary (one row per pair):
#   between_within_ratio  >> 1  => group-level offset dominates (ARTIFACT candidate)
#                          ~ 0.3 => within-group noise dominates (like M1/M2; no artifact)
#   offset_score = |median per-bin diff| / SD(per-bin diff)
#                          high  => systematic shift (artifact); ~0 => locus-specific (biology)

pkg_dir <- if (file.exists("DESCRIPTION")) "." else "SMFsim"
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  devtools::load_all(pkg_dir)
} else {
  library(SMFsim)
}
suppressWarnings(suppressMessages({
  library(data.table)
  library(ggplot2)
}))

config <- list(
  data_dir          = "data/allc",
  sample_sheet      = "data/sample_sheet.csv",
  output_dir        = "results/Mseries_efficiency",
  groups            = c("M1", "M2", "M3", "M4"),
  min_coverage      = 10,
  standard_chr_only = TRUE,
  chr_pattern       = "^(chr)?([0-9]{1,2}|X|Y)$",
  bin_size          = 1000L,
  min_sites_per_bin = 3L
)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load all M-series samples once ---------------------------------------
message("Loading M-series groups: ", paste(config$groups, collapse = ", "))
all_samples <- load_data(
  dir_path     = config$data_dir,
  sample_sheet = config$sample_sheet,
  groups       = config$groups,
  min_coverage = config$min_coverage
)
metadata <- attr(all_samples, "sample_metadata")
if (is.null(metadata)) stop("load_data did not attach sample_metadata.")

# Standard-chromosome filter (preserve metadata; lapply drops list attributes).
if (config$standard_chr_only) {
  n_before <- nrow(all_samples[[1]])
  all_samples <- lapply(all_samples, function(dt) dt[grepl(config$chr_pattern, chr)])
  attr(all_samples, "sample_metadata") <- metadata
  message(sprintf("Standard-chromosome filter: %d -> %d rows in sample 1",
                  n_before, nrow(all_samples[[1]])))
}

present <- intersect(config$groups, unique(metadata$group_id))
if (length(present) < 2) stop("Fewer than 2 of the requested groups are present.")
pairs <- utils::combn(present, 2, simplify = FALSE)
message(sprintf("Groups present: %s  ->  %d pairs",
                paste(present, collapse = ", "), length(pairs)))

# --- Per-pair diagnostic --------------------------------------------------
pair_diagnostic <- function(all_samples, metadata, gA, gB, config) {
  idx  <- which(metadata$group_id %in% c(gA, gB))
  sub  <- all_samples[idx]
  smeta <- metadata[idx]
  attr(sub, "sample_metadata") <- smeta

  # Shared sites within just this pair (maximises sites per comparison).
  sub <- find_shared_sites(sub, filter = TRUE)
  groups <- smeta$group_id            # order preserved by find_shared_sites
  n_shared <- nrow(sub[[1]])

  pair_label <- sprintf("%s_vs_%s", gA, gB)
  out_dir <- file.path(config$output_dir, pair_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Per-sample global (coverage-weighted) rate.
  global <- data.table(
    sample = names(sub),
    group  = groups,
    rate_weighted = vapply(sub, function(dt) sum(dt$mc) / sum(dt$cov), numeric(1)),
    mean_cov      = vapply(sub, function(dt) mean(dt$cov), numeric(1))
  )
  fwrite(global, file.path(out_dir, "per_sample_global_rate.csv"))

  gstat <- global[, .(within_sd = sd(rate_weighted),
                      group_mean = mean(rate_weighted)), by = group]
  mean_A <- gstat[group == gA]$group_mean
  mean_B <- gstat[group == gB]$group_mean
  between_offset <- abs(mean_A - mean_B)
  mean_within <- mean(gstat$within_sd, na.rm = TRUE)
  ratio <- between_offset / mean_within

  # Bin means and per-bin group comparison.
  bs <- config$bin_size
  bin_means <- rbindlist(lapply(seq_along(sub), function(i) {
    b <- sub[[i]][, .(bin_rate = sum(mc) / sum(cov), n = .N),
                  by = .(chr, bin = pos %/% bs)]
    b <- b[n >= config$min_sites_per_bin]
    b[, group := groups[i]]
    b
  }))
  group_bin <- bin_means[, .(grp_rate = mean(bin_rate)), by = .(chr, bin, group)]
  wide <- dcast(group_bin, chr + bin ~ group, value.var = "grp_rate")
  wide <- wide[stats::complete.cases(wide)]
  wide[, diff := get(gB) - get(gA)]
  diff_median <- median(wide$diff)
  diff_sd     <- sd(wide$diff)
  offset_score <- abs(diff_median) / diff_sd

  # Figures: global rate, hexbin, diff histogram.
  thm <- theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  p_rate <- ggplot(global, aes(x = group, y = rate_weighted, color = group)) +
    geom_point(size = 3.5) +
    geom_text(aes(label = sample), vjust = -1, size = 2.8, show.legend = FALSE) +
    labs(title = sprintf("Per-sample global GpC rate: %s vs %s", gA, gB),
         subtitle = sprintf("between/within ratio = %.2f", ratio),
         x = NULL, y = "Global rate") +
    thm + theme(legend.position = "none")

  p_hex <- ggplot(wide, aes(x = get(gA), y = get(gB))) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_fill_viridis_c(trans = "log10") +
    coord_equal() +
    labs(title = sprintf("Per-bin group means: %s vs %s", gA, gB),
         subtitle = "Off-diagonal parallel cloud = offset (artifact); on-diagonal scatter = biology",
         x = sprintf("%s bin mean", gA), y = sprintf("%s bin mean", gB), fill = "Bins") +
    thm

  p_hist <- ggplot(wide, aes(x = diff)) +
    geom_histogram(bins = 120, fill = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = diff_median, color = "red") +
    labs(title = sprintf("Per-bin difference (%s - %s)", gB, gA),
         subtitle = sprintf("median = %.4f (red); offset_score = %.2f", diff_median, offset_score),
         x = sprintf("%s - %s bin mean", gB, gA), y = "Bins") +
    thm

  ggsave(file.path(out_dir, "fig_global_rate.png"), p_rate, width = 6, height = 5, dpi = 150)
  ggsave(file.path(out_dir, "fig_bin_hexbin.png"),  p_hex,  width = 7, height = 6, dpi = 150)
  ggsave(file.path(out_dir, "fig_bin_diff_hist.png"), p_hist, width = 7, height = 5, dpi = 150)

  message(sprintf("  %s: n_shared=%d  ratio=%.2f  offset_score=%.2f  (rate %s=%.3f, %s=%.3f)",
                  pair_label, n_shared, ratio, offset_score, gA, mean_A, gB, mean_B))

  data.table(
    pair = pair_label, group_A = gA, group_B = gB,
    n_A = sum(groups == gA), n_B = sum(groups == gB),
    n_shared_sites = n_shared,
    rate_A = round(mean_A, 4), rate_B = round(mean_B, 4),
    within_sd_A = round(gstat[group == gA]$within_sd, 4),
    within_sd_B = round(gstat[group == gB]$within_sd, 4),
    between_offset = round(between_offset, 4),
    between_within_ratio = round(ratio, 3),
    bin_diff_median = round(diff_median, 4),
    bin_diff_sd = round(diff_sd, 4),
    offset_score = round(offset_score, 3)
  )
}

# --- Run all pairs --------------------------------------------------------
summary_rows <- lapply(pairs, function(p) {
  message(sprintf("\n=== Pair: %s vs %s ===", p[1], p[2]))
  tryCatch(pair_diagnostic(all_samples, metadata, p[1], p[2], config),
           error = function(e) {
             warning(sprintf("Pair %s vs %s failed: %s", p[1], p[2], e$message))
             NULL
           })
})
summary_dt <- rbindlist(Filter(Negate(is.null), summary_rows))
setorder(summary_dt, -between_within_ratio)
fwrite(summary_dt, file.path(config$output_dir, "efficiency_summary_all_pairs.csv"))

message("\n--- Efficiency summary (ranked by between/within ratio) ---")
print(summary_dt[, .(pair, between_within_ratio, offset_score,
                     rate_A, rate_B, between_offset)])

# --- Master ranking figure ------------------------------------------------
summary_dt[, pair := factor(pair, levels = pair)]
p_master <- ggplot(summary_dt, aes(x = pair, y = between_within_ratio)) +
  geom_col(fill = "#377EB8", width = 0.65) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.2f", between_within_ratio)), vjust = -0.4, size = 3.5) +
  labs(title = "Between-group efficiency offset across M-series pairs",
       subtitle = "Ratio >> 1 (above red line) = systematic offset = artifact candidate; ~0.3 = no offset (M1/M2-like)",
       x = NULL, y = "between / within global-rate ratio") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())
ggsave(file.path(config$output_dir, "fig_offset_ranking_all_pairs.png"),
       p_master, width = 9, height = 6, dpi = 150)
ggsave(file.path(config$output_dir, "fig_offset_ranking_all_pairs.pdf"),
       p_master, width = 9, height = 6)

message("\nDone. Per-pair folders + master summary in: ",
        normalizePath(config$output_dir, winslash = "/"))
message("\nNext step: any pair with between_within_ratio >> 1 (and high offset_score)")
message("is an artifact candidate -- run inst/scripts/run_real_data_parameter_test.R")
message("on that pair (set group_A/group_B) to show between-group correction efficacy.")
