#!/usr/bin/env Rscript

# 06_additional_figures.R
# Additional figures to better illustrate the precision-sensitivity tradeoff


# Plotting theme and shared aesthetics ------------------------------------
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
  "SMFnorm"   = "#4DAF4A",
  "ComBatMet"   = "#984EA3"
)

method_labels <- c(
  "raw"         = "Raw",
  "downsampled" = "Downsampled",
  "SMFnorm"   = "SMFnorm",
  "ComBatMet"   = "ComBatMet"
)

scenario_labels <- c(
  mild = "Mild Efficiency Bias",
  moderate = "Moderate Efficiency Bias",
  severe = "Severe Efficiency Bias"
)


# Figure 6: Side-by-side sensitivity + precision panels -------------------

#' Combined sensitivity and precision panels
#'
#' @param spikein_dt data.table of spike-in results
#' @param out_dir Output directory
#' @return ggplot object
plot_sensitivity_precision_panels <- function(spikein_dt, out_dir = NULL) {

  # Exclude ComBatMet if it's just returning raw data
  plot_dt <- copy(spikein_dt)
  plot_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  plot_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  # Reshape to long format for faceting by metric
  sens_dt <- plot_dt[, .(method, effect_size, efficiency_scenario, value = sensitivity)]
  sens_dt[, metric := "Sensitivity (Recall)"]

  prec_dt <- plot_dt[, .(method, effect_size, efficiency_scenario, value = precision)]
  prec_dt[, metric := "Precision (1 - FDR)"]

  long_dt <- rbind(sens_dt, prec_dt)
  long_dt[, metric := factor(metric, levels = c("Sensitivity (Recall)", "Precision (1 - FDR)"))]

  p <- ggplot(long_dt[!is.na(value)],
              aes(x = effect_size, y = value, color = method, shape = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_grid(metric ~ efficiency_scenario,
               labeller = labeller(efficiency_scenario = scenario_labels)) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_shape_manual(values = c(16, 17, 15, 18), labels = method_labels) +
    scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.30)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Sensitivity vs Precision Across Normalization Methods",
      subtitle = "SMFnorm trades modest sensitivity for near-perfect precision under efficiency artifacts",
      x = "Effect Size (Rate Difference)",
      y = NULL,
      color = "Method", shape = "Method"
    ) +
    theme_manuscript() +
    theme(strip.text.y = element_text(angle = 0))

  if (!is.null(out_dir)) {
    .save_figure(p, "fig6_sensitivity_precision_panels", out_dir, width = 12, height = 7)
  }

  return(p)
}


# Figure 7: Stacked bar — TP vs FP breakdown ------------------------------

#' Stacked bar chart showing composition of DMR calls (TP vs FP)
#'
#' @param spikein_dt data.table of spike-in results
#' @param out_dir Output directory
#' @return ggplot object
plot_call_breakdown <- function(spikein_dt, out_dir = NULL) {

  plot_dt <- copy(spikein_dt)
  plot_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  plot_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  # Reshape TP and FP to long
  tp_dt <- plot_dt[, .(method, effect_size, efficiency_scenario, count = TP, type = "True Positive")]
  fp_dt <- plot_dt[, .(method, effect_size, efficiency_scenario, count = FP, type = "False Positive")]
  long_dt <- rbind(tp_dt, fp_dt)
  long_dt[, type := factor(type, levels = c("False Positive", "True Positive"))]
  long_dt[, effect_label := paste0(effect_size * 100, "%")]
  long_dt[, effect_label := factor(effect_label,
                                    levels = paste0(sort(unique(effect_size)) * 100, "%"))]

  p <- ggplot(long_dt, aes(x = effect_label, y = count, fill = type)) +
    geom_col(position = "stack", width = 0.7) +
    facet_grid(efficiency_scenario ~ method,
               labeller = labeller(
                 efficiency_scenario = scenario_labels,
                 method = method_labels
               ),
               scales = "free_y") +
    scale_fill_manual(values = c("True Positive" = "#4DAF4A", "False Positive" = "#E41A1C")) +
    labs(
      title = "Composition of DMR Calls by Method",
      subtitle = "SMFnorm calls are almost exclusively true positives; other methods are dominated by false positives",
      x = "Effect Size",
      y = "Number of DMRs Called",
      fill = NULL
    ) +
    theme_manuscript() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(size = 10)
    )

  if (!is.null(out_dir)) {
    .save_figure(p, "fig7_call_breakdown", out_dir, width = 14, height = 9)
  }

  return(p)
}


# Figure 8: "What a researcher actually sees" — calls table ---------------

#' Focused comparison table at a single effect size
#'
#' @param spikein_dt data.table of spike-in results
#' @param target_effect Effect size to highlight (default 0.3)
#' @param out_dir Output directory
#' @return ggplot object
plot_researcher_view <- function(spikein_dt, target_effect = 0.3, out_dir = NULL) {

  plot_dt <- spikein_dt[effect_size == target_effect]
  plot_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  plot_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  # Build a "researcher's view" — for each call, what % is real?
  plot_dt[, pct_real := ifelse(total_called > 0, TP / total_called * 100, NA)]

  # Lollipop chart: total calls with color indicating purity
  p <- ggplot(plot_dt, aes(x = method, y = total_called)) +
    geom_segment(aes(xend = method, y = 0, yend = total_called), color = "grey60", linewidth = 0.8) +
    geom_point(aes(color = pct_real), size = 5) +
    geom_text(aes(label = paste0("TP=", TP, "\nFP=", FP)),
              vjust = -0.8, size = 3, lineheight = 0.85) +
    facet_wrap(~ efficiency_scenario, ncol = 3,
               labeller = labeller(efficiency_scenario = scenario_labels)) +
    scale_color_gradient2(
      low = "#E41A1C", mid = "#FFFFBF", high = "#4DAF4A",
      midpoint = 50, limits = c(0, 100),
      na.value = "grey80",
      name = "% True Positives"
    ) +
    scale_x_discrete(labels = method_labels) +
    labs(
      title = paste0("What a Researcher Sees: DMR Calls at Effect Size = ", target_effect),
      subtitle = "Point color shows the fraction of calls that are genuine; SMFnorm's calls are 100% real",
      x = NULL,
      y = "Total DMRs Called"
    ) +
    theme_manuscript() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )

  if (!is.null(out_dir)) {
    .save_figure(p, "fig8_researcher_view", out_dir, width = 12, height = 6)
  }

  return(p)
}


# Figure 9: F1 line plot (clearer than heatmap for this data) -------------

#' F1 score as line plot, easier to read than the heatmap
#'
#' @param spikein_dt data.table of spike-in results
#' @param out_dir Output directory
#' @return ggplot object
plot_f1_lines <- function(spikein_dt, out_dir = NULL) {

  plot_dt <- copy(spikein_dt)
  plot_dt[, method := factor(method, levels = names(method_colors))]

  scenario_order <- c("mild", "moderate", "severe")
  plot_dt[, efficiency_scenario := factor(efficiency_scenario, levels = scenario_order)]

  p <- ggplot(plot_dt[!is.na(F1)],
              aes(x = effect_size, y = F1, color = method, shape = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_wrap(~ efficiency_scenario,
               labeller = labeller(efficiency_scenario = scenario_labels)) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_shape_manual(values = c(16, 17, 15, 18), labels = method_labels) +
    scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.30)) +
    scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.15)) +
    labs(
      title = "F1 Score: Balancing Sensitivity and Precision",
      subtitle = "Under moderate-to-severe distortion, SMFnorm achieves the best balance",
      x = "Effect Size (Rate Difference)",
      y = "F1 Score",
      color = "Method", shape = "Method"
    ) +
    theme_manuscript()

  if (!is.null(out_dir)) {
    .save_figure(p, "fig9_f1_lines", out_dir, width = 10, height = 5)
  }

  return(p)
}


# Helper ------------------------------------------------------------------

.save_figure <- function(plot, name, out_dir, width = 8, height = 6) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))

  ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
  ggsave(png_path, plot, width = width, height = height, dpi = 300, device = "png")

  message(sprintf("Saved: %s (.pdf + .png)", name))
}


# Generate all additional figures -----------------------------------------

#' Generate all additional figures
#'
#' @param results_file Path to all_results.rds or the spikein CSV.
#' @param out_dir Output directory for figures.
#' @export
generate_additional_figures <- function(results_file, out_dir = "results/figures") {

  if (grepl("\\.rds$", results_file)) {
    results <- readRDS(results_file)
    spikein_dt <- results$spikein
  } else {
    spikein_dt <- fread(results_file)
  }

  # Ensure column names are correct
  if (!"efficiency_scenario" %in% names(spikein_dt)) {
    # Parse from the scenario column if needed
    spikein_dt[, c("efficiency_scenario", "effect_info") :=
                  tstrsplit(scenario, "_effect", fixed = TRUE)]
    spikein_dt[, effect_size := as.numeric(effect_info)]
    spikein_dt[, effect_info := NULL]
  }

  message("Generating additional figures...")

  fig6 <- plot_sensitivity_precision_panels(spikein_dt, out_dir)
  fig7 <- plot_call_breakdown(spikein_dt, out_dir)
  fig8 <- plot_researcher_view(spikein_dt, target_effect = 0.3, out_dir)
  fig9 <- plot_f1_lines(spikein_dt, out_dir)

  message("All additional figures generated.")

  return(invisible(list(fig6 = fig6, fig7 = fig7, fig8 = fig8, fig9 = fig9)))
}
