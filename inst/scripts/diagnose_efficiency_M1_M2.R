#!/usr/bin/env Rscript

# diagnose_efficiency_M1_M2.R
#
# Diagnostic for the M1 vs M2 (MCF10) real-data comparison. Assesses:
#   1. Per-sample global GpC methylation rate (a conversion-efficiency proxy)
#      and the within-group vs between-group spread.
#   2. The distribution of bin-mean rates per replicate (the within-replicate
#      variation that motivated trying between-group normalization).
#   3. Whether the M1-vs-M2 difference is a UNIFORM GLOBAL OFFSET (efficiency
#      artifact, correctable) or LOCUS-SPECIFIC (transformation biology, not).
#
# A true efficiency artifact scales every site -> shifts the whole distribution.
# Real biology moves specific loci -> scatters off the M1=M2 diagonal.
#
# Usage (from package root or one level up):
#   Rscript inst/scripts/diagnose_efficiency_M1_M2.R
# Reuses an existing `real_data` object if one is already in the session.

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
  output_dir        = "results/real_data_M1_M2/efficiency_diagnostic",
  min_coverage      = 10,
  bin_size          = 1000L,   # genomic bin width (bp)
  min_sites_per_bin = 3L       # drop sparse, noisy bins
)
group_A <- "M1"
group_B <- "M2"

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse real_data if already loaded interactively; otherwise load it.
if (!exists("real_data")) {
  message("No 'real_data' in environment; loading from ", config$data_dir)
  real_data <- load_real_data(config, group_A, group_B)
} else {
  message("Reusing existing 'real_data' from the session.")
}

all_samples <- c(real_data$group_A, real_data$group_B)
sample_ids  <- names(all_samples)
groups      <- c(rep(group_A, length(real_data$group_A)),
                 rep(group_B, length(real_data$group_B)))

group_colors <- c("#1B9E77", "#D95F02")
names(group_colors) <- c(group_A, group_B)

# --- 1. Per-sample global rate (efficiency proxy) -------------------------
global <- data.table(
  sample = sample_ids,
  group  = groups,
  rate_unweighted = vapply(all_samples, function(dt) mean(dt$rate), numeric(1)),
  rate_weighted   = vapply(all_samples, function(dt) sum(dt$mc) / sum(dt$cov), numeric(1)),
  mean_cov        = vapply(all_samples, function(dt) mean(dt$cov), numeric(1)),
  n_sites         = vapply(all_samples, nrow, integer(1))
)
message("\n--- Per-sample global GpC rate (coverage-weighted = efficiency proxy) ---")
print(global)
fwrite(global, file.path(config$output_dir, "per_sample_global_rate.csv"))

# Within-group spread vs between-group offset, on the weighted global rate.
grp_stats <- global[, .(within_sd = sd(rate_weighted),
                        group_mean = mean(rate_weighted)), by = group]
between_diff <- abs(diff(grp_stats$group_mean))
mean_within  <- mean(grp_stats$within_sd, na.rm = TRUE)
message(sprintf("\nWithin-group SD of global rate: %s",
                paste(sprintf("%s=%.4f", grp_stats$group, grp_stats$within_sd),
                      collapse = "  ")))
message(sprintf("Between-group mean offset:      %.4f", between_diff))
message(sprintf("Between/within ratio:          %.2f  (>>1 => group-level offset; ~1 => within-group noise dominates)",
                between_diff / mean_within))

# --- 2. Bin-mean rates per sample -----------------------------------------
bs <- config$bin_size
message(sprintf("\nComputing %d-bp bin means (min %d sites/bin)...",
                bs, config$min_sites_per_bin))
bin_means <- rbindlist(lapply(seq_along(all_samples), function(i) {
  dt <- all_samples[[i]]
  b  <- dt[, .(bin_rate = sum(mc) / sum(cov), n = .N),
           by = .(chr, bin = pos %/% bs)]
  b  <- b[n >= config$min_sites_per_bin]
  b[, `:=`(sample = sample_ids[i], group = groups[i])]
  b
}))

# --- 3. Per-bin group means: global offset vs locus-specific --------------
group_bin <- bin_means[, .(grp_rate = mean(bin_rate)), by = .(chr, bin, group)]
wide <- dcast(group_bin, chr + bin ~ group, value.var = "grp_rate")
wide <- wide[stats::complete.cases(wide)]   # bins covered in both groups
wide[, diff := get(group_B) - get(group_A)]
message(sprintf("\nBins shared by both groups: %d", nrow(wide)))
message(sprintf("Median per-bin (%s - %s): %.4f   Mean: %.4f   SD: %.4f",
                group_B, group_A, median(wide$diff), mean(wide$diff), sd(wide$diff)))
message("  (median far from 0 => systematic global offset; centered ~0 with wide SD => locus-specific)")

# --- Figures --------------------------------------------------------------
theme_diag <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p1 <- ggplot(global, aes(x = group, y = rate_weighted, color = group)) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_color_manual(values = group_colors, guide = "none") +
  labs(title = "Per-sample global GpC rate (efficiency proxy)",
       subtitle = "Coverage-weighted genome-wide mean; vertical spread within a group = within-group efficiency variation",
       x = NULL, y = "Global methylation rate") +
  theme_diag

p2 <- ggplot(bin_means, aes(x = bin_rate, group = sample, color = group)) +
  geom_density(linewidth = 0.6, alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ group, ncol = 1) +
  labs(title = sprintf("Distribution of %d-bp bin-mean rates per replicate", bs),
       subtitle = "Each curve = one replicate. Spread within a panel = within-group variation; shift between panels = between-group",
       x = "Bin mean rate", y = "Density", color = NULL) +
  theme_diag

p3 <- ggplot(wide, aes(x = get(group_A), y = get(group_B))) +
  geom_bin2d(bins = 80) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10") +
  coord_equal() +
  labs(title = "Per-bin group means: M1 vs M2",
       subtitle = "Cloud parallel to but off the diagonal = uniform global offset (artifact-like); broad scatter = locus-specific (biology)",
       x = sprintf("%s bin mean", group_A), y = sprintf("%s bin mean", group_B),
       fill = "Bins") +
  theme_diag

p4 <- ggplot(wide, aes(x = diff)) +
  geom_histogram(bins = 120, fill = "grey50") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_vline(xintercept = median(wide$diff), color = "red") +
  labs(title = sprintf("Per-bin difference (%s - %s)", group_B, group_A),
       subtitle = "Red = median. Far from 0 = systematic offset (efficiency); centered at 0 with wide spread = locus-specific biology",
       x = sprintf("%s - %s bin mean", group_B, group_A), y = "Bins") +
  theme_diag

save_fig <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(config$output_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 150)
  ggsave(file.path(config$output_dir, paste0(name, ".pdf")), p,
         width = w, height = h)
  message("Saved: ", name)
}

# ggrepel is optional; fall back to plain labels if unavailable.
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  p1 <- ggplot(global, aes(x = group, y = rate_weighted, color = group)) +
    geom_point(size = 3.5) +
    geom_text(aes(label = sample), vjust = -1, size = 3, show.legend = FALSE) +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(title = "Per-sample global GpC rate (efficiency proxy)",
         subtitle = "Coverage-weighted genome-wide mean",
         x = NULL, y = "Global methylation rate") +
    theme_diag
}

save_fig(p1, "fig_per_sample_global_rate")
save_fig(p2, "fig_bin_mean_distributions", h = 7)
save_fig(p3, "fig_bin_means_M1_vs_M2")
save_fig(p4, "fig_bin_difference_histogram")

message("\nDiagnostic complete. Figures and CSV in: ",
        normalizePath(config$output_dir, winslash = "/"))
