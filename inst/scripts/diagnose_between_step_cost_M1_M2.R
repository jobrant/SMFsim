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
  min_overlap_bp = 1L
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

# --- Decompose raw DMRs ---------------------------------------------------
in_W <- overlaps_any(R, W)
within_cleanup <- R[!in_W]            # raw DMRs removed by within-group step
survived_within <- R[in_W]            # raw DMRs kept by within-group step

in_B_sw <- overlaps_any(survived_within, B)
between_cost <- survived_within[!in_B_sw]   # removed only by the between step

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

message("\nDone. Tables + figure in: ",
        normalizePath(config$output_dir, winslash = "/"))
