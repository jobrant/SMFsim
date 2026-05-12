#!/usr/bin/env Rscript

# 05_plot_results.R
# Generate publication-quality figures from simulation results

library(data.table)
library(ggplot2)

# Plotting theme ----------------------------------------------------------

theme_manuscript <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80")
    )
}

method_colors <- c(
  "raw"         = "#E41A1C",
  "downsampled" = "#377EB8",
  "SMFnorm"   = "#4DAF4A",
  "ComBatMet"   = "#984EA3"
)

method_labels <- c(
  "raw"         = "Raw (no normalization)",
  "downsampled" = "Downsampled",
  "SMFnorm"   = "SMFnorm",
  "ComBatMet"   = "ComBatMet"
)


# Figure 1: Null simulation — false positive counts -----------------------

#' Bar plot of false positive DMR counts across scenarios and methods
#'
#' @param null_dt data.table from compile_null_results()
#' @param out_dir Output directory for figure files.
#' @return ggplot object
plot_null_fp <- function(null_dt, out_dir = NULL) {
  null_dt[, method := factor(method, levels = names(method_colors))]

  # Order scenarios by severity
  scenario_order <- c("mild", "moderate", "severe")
  null_dt[, scenario := factor(scenario, levels = scenario_order)]

  p <- ggplot(null_dt, aes(x = scenario, y = FP, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = method_colors, labels = method_labels) +
    labs(
      title = "False Positive DMRs Under Null (No True Differences)",
      x = "Enzyme Efficiency Scenario",
      y = "Number of False Positive DMRs",
      fill = "Method"
    ) +
    theme_manuscript() +
    # Add count labels on bars
    geom_text(aes(label = FP),
              position = position_dodge(width = 0.8),
              vjust = -0.5, size = 3.5)

  if (!is.null(out_dir)) {
    .save_figure(p, "fig1_null_fp", out_dir, width = 8, height = 5)
  }

  return(p)
}


# Figure 2: Sensitivity curves --------------------------------------------

#' Line plot of sensitivity vs effect size, faceted by efficiency scenario
#'
#' @param spikein_dt data.table from compile_results() for spike-in sim.
#' @param out_dir Output directory.
#' @return ggplot object
plot_sensitivity_curves <- function(spikein_dt, out_dir = NULL) {
  spikein_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  spikein_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  p <- ggplot(spikein_dt, aes(x = effect_size, y = sensitivity,
                               color = method, shape = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_wrap(~ efficiency_scenario,
               labeller = labeller(efficiency_scenario = c(
                 mild = "Mild Efficiency Bias",
                 moderate = "Moderate Efficiency Bias",
                 severe = "Severe Efficiency Bias"
               ))) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_shape_manual(values = c(16, 17, 15, 18), labels = method_labels) +
    scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.30)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Sensitivity to Detect Spike-In DMRs",
      x = "Effect Size (Rate Difference)",
      y = "Sensitivity (TP / (TP + FN))",
      color = "Method", shape = "Method"
    ) +
    theme_manuscript()

  if (!is.null(out_dir)) {
    .save_figure(p, "fig2_sensitivity", out_dir, width = 10, height = 5)
  }

  return(p)
}


# Figure 3: F1 score heatmap ----------------------------------------------

#' Heatmap of F1 scores across effect sizes and methods
#'
#' @param spikein_dt data.table from compile_results().
#' @param out_dir Output directory.
#' @return ggplot object
plot_f1_heatmap <- function(spikein_dt, out_dir = NULL) {
  spikein_dt[, method := factor(method, levels = rev(names(method_colors)))]
  spikein_dt[, effect_label := factor(paste0(effect_size * 100, "%"))]

  scenario_order <- c("mild", "moderate", "severe")
  spikein_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  p <- ggplot(spikein_dt, aes(x = effect_label, y = method, fill = F1)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", F1)), size = 3.2, color = "black") +
    facet_wrap(~ efficiency_scenario, ncol = 1,
               labeller = labeller(efficiency_scenario = c(
                 mild = "Mild", moderate = "Moderate", severe = "Severe"
               ))) +
    scale_fill_gradient2(low = "#E41A1C", mid = "#FFFFBF", high = "#4DAF4A",
                         midpoint = 0.5, limits = c(0, 1),
                         na.value = "grey90") +
    labs(
      title = "F1 Score Across Methods and Conditions",
      x = "Effect Size",
      y = "Method",
      fill = "F1 Score"
    ) +
    theme_manuscript() +
    theme(legend.position = "right")

  if (!is.null(out_dir)) {
    .save_figure(p, "fig3_f1_heatmap", out_dir, width = 8, height = 8)
  }

  return(p)
}


# Figure 4: Precision-Recall tradeoff -------------------------------------

#' Precision vs sensitivity scatter, one point per scenario × effect × method
#'
#' @param spikein_dt data.table from compile_results().
#' @param out_dir Output directory.
#' @return ggplot object
plot_precision_recall <- function(spikein_dt, out_dir = NULL) {
  spikein_dt[, method := factor(method, levels = names(method_colors))]

  p <- ggplot(spikein_dt[!is.na(sensitivity) & !is.na(precision)],
              aes(x = sensitivity, y = precision, color = method)) +
    geom_point(aes(size = effect_size), alpha = 0.7) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_size_continuous(name = "Effect Size",
                          breaks = c(0.05, 0.10, 0.20, 0.30),
                          range = c(1.5, 5)) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Precision vs Sensitivity Trade-off",
      x = "Sensitivity (Recall)",
      y = "Precision",
      color = "Method"
    ) +
    theme_manuscript()

  if (!is.null(out_dir)) {
    .save_figure(p, "fig4_precision_recall", out_dir, width = 8, height = 6)
  }

  return(p)
}


# Figure 5: FDR comparison ------------------------------------------------

#' FDR across methods and effect sizes
#'
#' @param spikein_dt data.table from compile_results().
#' @param out_dir Output directory.
#' @return ggplot object
plot_fdr <- function(spikein_dt, out_dir = NULL) {
  spikein_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  spikein_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  p <- ggplot(spikein_dt, aes(x = effect_size, y = FDR,
                               color = method, shape = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
    facet_wrap(~ efficiency_scenario,
               labeller = labeller(efficiency_scenario = c(
                 mild = "Mild", moderate = "Moderate", severe = "Severe"
               ))) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_shape_manual(values = c(16, 17, 15, 18), labels = method_labels) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "False Discovery Rate Across Conditions",
      x = "Effect Size",
      y = "FDR (FP / (TP + FP))",
      color = "Method", shape = "Method"
    ) +
    annotate("text", x = 0.28, y = 0.08, label = "FDR = 0.05",
             color = "grey50", size = 3) +
    theme_manuscript()

  if (!is.null(out_dir)) {
    .save_figure(p, "fig5_fdr", out_dir, width = 10, height = 5)
  }

  return(p)
}


# Helper ------------------------------------------------------------------

#' Save figure as PDF and PNG
#' @keywords internal
.save_figure <- function(plot, name, out_dir, width = 8, height = 6) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))

  ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
  ggsave(png_path, plot, width = width, height = height, dpi = 300, device = "png")

  message(sprintf("Saved: %s (.pdf + .png)", name))
}


# Generate all figures ----------------------------------------------------

#' Generate all manuscript figures from saved simulation results
#'
#' @param results_file Path to all_results.rds from run_simulation.
#' @param out_dir Output directory for figures.
generate_all_figures <- function(results_file, out_dir = "results/figures") {
  results <- readRDS(results_file)

  message("Generating figures...")

  fig1 <- plot_null_fp(results$null, out_dir)
  fig2 <- plot_sensitivity_curves(results$spikein, out_dir)
  fig3 <- plot_f1_heatmap(results$spikein, out_dir)
  fig4 <- plot_precision_recall(results$spikein, out_dir)
  fig5 <- plot_fdr(results$spikein, out_dir)

  message("All figures generated.")

  return(invisible(list(fig1 = fig1, fig2 = fig2, fig3 = fig3,
                        fig4 = fig4, fig5 = fig5)))
}
