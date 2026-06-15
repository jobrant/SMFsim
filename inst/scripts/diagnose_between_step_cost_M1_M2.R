#!/usr/bin/env Rscript

# diagnose_between_step_cost_M1_M2.R
#
# Decomposes what SMFnorm removes on the M1 vs M2 comparison into:
#   (1) within-group cleanup  = raw DMRs dropped by within-group correction
#       (rate_between_groups = FALSE). These are calls within-correction judged
#       unstable -- the legitimate noise cleanup.
#   (2) between-step cost      = DMRs that SURVIVED within-group correction but
#       were dropped only once between-group correction was added
#       (rate_between_groups = TRUE). These are NOT within-group noise (they
#       passed that filter); their removal is the between-step flattening
#       group-mean differences.
#
# If the between-step cost set has LARGER effects than the within cleanup set,
# the between-group step is removing real (locus-specific) biology, not artifact
# -- consistent with the efficiency diagnostic showing no between-group offset
# on this pair. If it is uniformly tiny, the cost is small.
#
# Reads the saved real_data_results.rds from both runs; identifies which is
# within-only vs both from the rate_between_groups flag stored in each, so it
# does not depend on folder names.
#
# Usage (from package root or one level up):
#   Rscript inst/scripts/diagnose_between_step_cost_M1_M2.R

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
  results_dir = "results/real_data_M1_M2",
  output_dir  = "results/real_data_M1_M2/between_step_cost",
  min_overlap_bp = 1L,
  # For the within-group variance test (needs per-site data):
  data_dir          = "data/allc",
  sample_sheet      = "data/sample_sheet.csv",
  min_coverage      = 10,
  group_A           = "M1",
  group_B           = "M2",
  near_threshold_max = 0.15   # focus the variance test on the fragile tail
)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Locate and classify the two runs -------------------------------------
rds_files <- list.files(config$results_dir, pattern = "real_data_results\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if (length(rds_files) < 2) {
  stop("Need both a within-only and a both run under ", config$results_dir,
       " (found ", length(rds_files), " real_data_results.rds files).")
}

runs <- lapply(rds_files, readRDS)
# Each saved object stores config = list(config=, mode=, SMFnorm_params=).
rbg <- vapply(runs, function(r) {
  isTRUE(r$config$SMFnorm_params$rate_between_groups)
}, logical(1))
modes <- vapply(runs, function(r) r$config$mode %||% NA_character_, character(1))

if (!any(rbg) || !any(!rbg)) {
  stop("Could not find one within-only (rate_between_groups=FALSE) and one ",
       "both (TRUE) run. Modes found: ", paste(modes, collapse = ", "))
}

within_run <- runs[[which(!rbg)[1]]]
both_run   <- runs[[which(rbg)[1]]]
message(sprintf("within-only run: mode='%s'  |  both run: mode='%s' (between_alpha=%s)",
                within_run$config$mode, both_run$config$mode,
                both_run$config$SMFnorm_params$between_alpha %||% "?"))

R <- both_run$dmr_list$raw            # raw is identical across modes
W <- within_run$dmr_list$SMFnorm      # within-only SMFnorm
B <- both_run$dmr_list$SMFnorm        # both (within + between) SMFnorm
message(sprintf("raw=%d  within-only SMFnorm=%d  both SMFnorm=%d",
                nrow(R), nrow(W), nrow(B)))

# --- Overlap helper -------------------------------------------------------
overlaps_any <- function(q, r) {
  if (nrow(q) == 0) return(logical(0))
  if (nrow(r) == 0) return(rep(FALSE, nrow(q)))
  qg <- GenomicRanges::GRanges(q$chr, IRanges::IRanges(q$start, q$end))
  rg <- GenomicRanges::GRanges(r$chr, IRanges::IRanges(r$start, r$end))
  GenomicRanges::countOverlaps(qg, rg, minoverlap = config$min_overlap_bp) > 0
}

# Per-DMR within-group replicate variance: for each DMR, compute each replicate's
# region methylation (sum(mc)/sum(cov) over sites in the DMR), then the SD across
# replicates within each group, combined across groups. High value = replicates
# disagree = noise-driven call.
compute_within_var <- function(real_data, dmrs) {
  ref <- dmrs[, .(dmr_id, chr, start, end)]
  data.table::setkey(ref, chr, start, end)

  all_samples <- c(real_data$group_A, real_data$group_B)
  groups <- c(rep(real_data$group_names[1], length(real_data$group_A)),
              rep(real_data$group_names[2], length(real_data$group_B)))

  per_sample <- rbindlist(lapply(seq_along(all_samples), function(i) {
    s <- all_samples[[i]][, .(chr, start = pos, end = pos, mc, cov)]
    data.table::setkey(s, chr, start, end)
    ov <- data.table::foverlaps(s, ref, type = "within", nomatch = 0L)
    if (nrow(ov) == 0) return(NULL)
    agg <- ov[, .(region_rate = sum(mc) / sum(cov)), by = dmr_id]
    agg[, group := groups[i]]
    agg
  }))

  # SD across replicates within each group, then average the per-group SDs.
  grp_sd <- per_sample[, .(grp_sd = sd(region_rate)), by = .(dmr_id, group)]
  grp_sd[, .(within_sd = mean(grp_sd, na.rm = TRUE)), by = dmr_id]
}

# --- Decompose raw DMRs ---------------------------------------------------
in_W <- overlaps_any(R, W)
within_cleanup <- R[!in_W]            # raw DMRs removed by within-group step
survived_within <- R[in_W]            # raw DMRs kept by within-group step

in_B_sw <- overlaps_any(survived_within, B)
between_cost <- survived_within[!in_B_sw]     # removed only by the between step
between_retained <- survived_within[in_B_sw]  # survived both steps

message(sprintf("\nWithin-group cleanup (raw dropped by within step):   %d  (median |eff| %.4f)",
                nrow(within_cleanup), median(abs(within_cleanup$mean_diff))))
message(sprintf("Between-step cost (survived within, dropped by both): %d  (median |eff| %.4f)",
                nrow(between_cost), median(abs(between_cost$mean_diff))))

frac_large <- function(x, thr = 0.2) mean(abs(x) >= thr)
message(sprintf("Fraction with |eff| >= 0.2  --  within cleanup: %.1f%%   between cost: %.1f%%",
                100 * frac_large(within_cleanup$mean_diff),
                100 * frac_large(between_cost$mean_diff)))
message("  (between cost >= within cleanup on effect size => between step is removing real biology, not noise)")

# --- Save tables ----------------------------------------------------------
fwrite(within_cleanup, file.path(config$output_dir, "within_group_cleanup_dmrs.csv"))
fwrite(between_cost,   file.path(config$output_dir, "between_step_cost_dmrs.csv"))

# --- Figure ---------------------------------------------------------------
plot_dt <- rbindlist(list(
  data.table(set = "Within-group cleanup", eff = abs(within_cleanup$mean_diff)),
  data.table(set = "Between-step cost",    eff = abs(between_cost$mean_diff))
))

p <- ggplot(plot_dt, aes(x = eff, fill = set)) +
  geom_density(alpha = 0.55) +
  scale_fill_manual(values = c("Within-group cleanup" = "#377EB8",
                               "Between-step cost"     = "#E41A1C")) +
  labs(
    title = "What each normalization step removes (M1 vs M2)",
    subtitle = "Between-step cost shifted toward larger effects = real biology removed, not artifact",
    x = "|Mean Difference| of removed DMRs", y = "Density", fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(config$output_dir, "fig_removal_decomposition.png"), p,
       width = 8, height = 6, dpi = 150)
ggsave(file.path(config$output_dir, "fig_removal_decomposition.pdf"), p,
       width = 8, height = 6)

# --- Within-group variance enrichment (small-FP vs small-real) ------------
# Among near-threshold raw DMRs, do the ones the between step REMOVES have higher
# within-group replicate variance (noise-driven, FP-like) than the ones it
# RETAINS? Effect size can't separate small-real from small-false; variance can.
near <- config$near_threshold_max
fate <- rbindlist(list(
  data.table(within_cleanup[abs(mean_diff) <= near, .(chr, start, end)],
             fate = "Removed (within step)"),
  data.table(between_cost[abs(mean_diff) <= near, .(chr, start, end)],
             fate = "Removed (between step)"),
  data.table(between_retained[abs(mean_diff) <= near, .(chr, start, end)],
             fate = "Retained (both)")
))
fate[, dmr_id := .I]

# Variance needs per-site rates; reuse real_data if present, else load it.
if (!exists("real_data")) {
  message("\nLoading real_data for the variance test...")
  real_data <- load_real_data(
    list(data_dir = config$data_dir, sample_sheet = config$sample_sheet,
         min_coverage = config$min_coverage),
    group_A = config$group_A, group_B = config$group_B)
}

within_var <- compute_within_var(real_data, fate)
fate <- merge(fate, within_var, by = "dmr_id", all.x = TRUE)
fate <- fate[is.finite(within_sd)]

med <- fate[, .(n = .N, median_within_sd = round(median(within_sd), 4)), by = fate]
message("\n--- Within-group variance by fate (near-threshold DMRs) ---")
print(med)

wt <- stats::wilcox.test(
  fate[fate == "Removed (between step)"]$within_sd,
  fate[fate == "Retained (both)"]$within_sd
)
message(sprintf("\nWilcoxon (removed-between vs retained): p = %.3g", wt$p.value))
message("  removed > retained in within-group SD => between step removes the noisier (FP-like) calls;")
message("  no difference => removed near-threshold DMRs look like the retained ones (small real biology).")

fate[, fate := factor(fate, levels = c("Removed (within step)",
                                       "Removed (between step)", "Retained (both)"))]
pv <- ggplot(fate, aes(x = fate, y = within_sd, fill = fate)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Removed (within step)"  = "#377EB8",
                               "Removed (between step)"  = "#E41A1C",
                               "Retained (both)"         = "#4DAF4A"),
                    guide = "none") +
  labs(title = "Within-group replicate variance of near-threshold DMRs by fate",
       subtitle = sprintf("|effect| <= %.2f. Higher variance among REMOVED = noise-driven (FP) cleanup",
                          near),
       x = NULL, y = "Within-group SD of region rate") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(config$output_dir, "fig_within_variance_by_fate.png"), pv,
       width = 8, height = 6, dpi = 150)
ggsave(file.path(config$output_dir, "fig_within_variance_by_fate.pdf"), pv,
       width = 8, height = 6)
fwrite(fate, file.path(config$output_dir, "near_threshold_dmr_fate_variance.csv"))

message("\nDone. Tables + figures in: ",
        normalizePath(config$output_dir, winslash = "/"))
