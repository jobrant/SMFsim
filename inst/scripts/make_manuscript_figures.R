#!/usr/bin/env Rscript

# make_manuscript_figures.R
#
# Manuscript figures and tables from the across-seed benchmark results.
#
# Reads the per-seed tables that summarize_seeds.R and
# diagnose_spikein_calls.R write under
#   <base>/<model>/all_seeds/<block>/
# for the blocks control_within, bias_parametric, broad_between and
# broad_within, and writes to <out_dir>, per efficiency model:
#
#   fig_null_false_positives_<model>     false-positive calls in the null runs
#   fig_spikein_sensitivity_<model>      sensitivity vs effect size, bias scenarios
#   fig_spikein_sensitivity_within_<model>  the same, within-group scenarios
#   fig_broad_tradeoff_<model>           raw vs within-only vs within + between
#
# each as .png (300 dpi) and .pdf, and one table per figure with both models
# (mean, sd, min and max across seeds):
#
#   table_null_false_positives.csv
#   table_spikein_sensitivity.csv
#   table_broad_tradeoff.csv
#
# Every figure shows each seed as a small point and the mean across seeds as
# a large one (or a line); with five seeds that is more honest than an error
# bar.
#
# Usage, from the manuscript root:
#   Rscript SMFsim/inst/scripts/make_manuscript_figures.R [base_dir] [out_dir]
#
# base_dir defaults to results/bias_experiment, out_dir to
# results/figures/manuscript.

suppressMessages({
    library(data.table)
    library(ggplot2)
})
if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("This script needs the patchwork package.", call. = FALSE)
}
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
base    <- if (length(args) >= 1) args[1] else "results/bias_experiment"
out_dir <- if (length(args) >= 2) args[2] else "results/figures/manuscript"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

models <- c("multiplicative", "logit")
models <- models[dir.exists(file.path(base, models, "all_seeds"))]
if (!length(models)) stop("No <model>/all_seeds under ", base, call. = FALSE)


# Labels, colours, theme ----------------------------------------------------

method_labels <- c(raw = "Raw", downsampled = "Downsampled", SMFnorm = "SMFnorm",
                   ComBatMet_balanced = "ComBat-met (balanced)",
                   ComBatMet_median = "ComBat-met (median)")
# Fixed colour and shape per method, the same in every figure.
method_colors <- setNames(c("#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac"),
                          method_labels)
method_shapes <- setNames(c(16, 17, 15, 18, 4), method_labels)

arm_labels <- c(raw = "Raw", within = "SMFnorm, within-group only",
                both = "SMFnorm, within + between")
arm_colors <- setNames(c("#257de0", "#00b894", "#6a4ce2"), arm_labels)
arm_shapes <- setNames(c(16, 15, 17), arm_labels)

scenario_labels <- c(
    mild = "Mild", moderate = "Moderate", severe = "Severe",
    aligned_moderate = "Aligned,\nmoderate", aligned_strong = "Aligned,\nstrong",
    imbalanced_moderate = "Imbalanced,\nmoderate",
    imbalanced_strong = "Imbalanced,\nstrong",
    broad_10 = "10%", broad_30 = "30%", broad_50 = "50%",
    broad_30_aligned = "30%, with\nefficiency bias")

.lab <- function(x, labels) factor(labels[as.character(x)], levels = labels)

ink <- "#191919"; ink2 <- "#666666"; grid_col <- "#e6e6e6"
theme_ms <- function(base_size = 9) {
    theme_minimal(base_size = base_size) +
        theme(text = element_text(colour = ink),
              axis.text = element_text(colour = ink2),
              axis.title = element_text(colour = ink),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
              axis.line.x = element_line(colour = ink2, linewidth = 0.3),
              axis.ticks.x = element_line(colour = ink2, linewidth = 0.3),
              strip.text = element_text(colour = ink, face = "bold", hjust = 0),
              plot.title = element_text(face = "bold", size = rel(1.05)),
              plot.title.position = "plot",
              plot.subtitle = element_text(colour = ink2),
              legend.position = "bottom", legend.title = element_blank(),
              legend.key.height = unit(10, "pt"),
              plot.background = element_rect(fill = "white", colour = NA))
}

.save <- function(p, name, width, height) {
    for (ext in c("png", "pdf")) {
        ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = width,
               height = height, dpi = 300, bg = "white")
    }
    message("  wrote ", name, " (.png, .pdf)")
}

.read <- function(model, block, file) {
    f <- file.path(base, model, "all_seeds", block, file)
    if (!file.exists(f)) return(NULL)
    fread(f)
}

.across <- function(dt, value, by) {
    dt[, .(n_seeds = .N, mean = mean(.v), sd = stats::sd(.v),
           min = min(.v), max = max(.v)), by = by, env = list(.v = value)]
}


# 1. Null false positives ---------------------------------------------------

null_block_labels <- c(
    control_within  = "Efficiency differs within groups",
    bias_parametric = "Efficiency differs between groups")

.null_data <- function(model) {
    rbindlist(lapply(names(null_block_labels), function(b) {
        d <- .read(model, b, "null_by_seed.csv")
        if (is.null(d)) return(NULL)
        d[, .(model = model, block = b, seed, scenario, method, FP)]
    }))
}

plot_null <- function(model) {
    d <- .null_data(model)
    if (!nrow(d)) return(NULL)
    d[, `:=`(m = .lab(method, method_labels), s = .lab(scenario, scenario_labels),
             panel = .lab(block, null_block_labels))]
    mean_d <- d[, .(FP = mean(FP)), by = .(panel, s, m)]
    dodge <- position_dodge(width = 0.7)
    p <- ggplot(d, aes(s, FP, colour = m, shape = m, group = m)) +
        geom_point(position = dodge, size = 0.9, alpha = 0.45) +
        geom_point(data = mean_d, position = dodge, size = 2.6) +
        facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
        scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                           breaks = c(0, 10, 100, 1000, 10000),
                           labels = scales::comma) +
        scale_colour_manual(values = method_colors, drop = FALSE) +
        scale_shape_manual(values = method_shapes, drop = FALSE) +
        labs(x = NULL, y = "False-positive calls (log scale)",
             title = "False-positive calls when the groups do not differ",
             subtitle = paste0("Null simulation, ", model, " efficiency model. ",
                               "Large symbols: mean of ", uniqueN(d$seed),
                               " seeds; small: each seed.")) +
        theme_ms()
    .save(p, paste0("fig_null_false_positives_", model), 7.2, 3.8)
    .across(d, "FP", c("model", "block", "scenario", "method"))
}


# 2. Spike-in sensitivity ---------------------------------------------------

measure_labels <- c(
    called      = "Called: q < 0.05 and |difference| >= 0.1",
    significant = "Significant: q < 0.05")

# Per seed: the share of planted regions called, and the share significant
# (called, or significant with an effect under the 0.1 cut).
.sens_data <- function(model, block) {
    st <- .read(model, block, "spikein_region_status_by_seed.csv")
    if (!is.null(st)) {
        d <- st[, .(called = mean(status == "called"),
                    significant = mean(status %in% c("called", "sig_below_cut"))),
                by = .(seed, scenario, effect_size, method)]
    } else {
        d <- .read(model, block, "spikein_by_seed.csv")
        if (is.null(d)) return(NULL)
        if (!"sensitivity_significant" %in% names(d)) d[, sensitivity_significant := NA_real_]
        d <- d[, .(seed, scenario = efficiency_scenario, effect_size, method,
                   called = sensitivity, significant = sensitivity_significant)]
    }
    d <- melt(d, measure.vars = c("called", "significant"),
              variable.name = "measure", value.name = "sensitivity")
    d[, `:=`(model = model, block = block)]
    d[]
}

plot_sensitivity <- function(model, block, name, title) {
    d <- .sens_data(model, block)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[!is.na(sensitivity)]
    d[, `:=`(m = .lab(method, method_labels),
             s = .lab(scenario, gsub("\n", " ", scenario_labels)),
             meas = .lab(measure, measure_labels))]
    mean_d <- d[, .(sensitivity = mean(sensitivity)), by = .(meas, s, m, effect_size)]
    p <- ggplot(d, aes(effect_size, sensitivity, colour = m, shape = m, group = m)) +
        geom_point(size = 0.7, alpha = 0.35) +
        geom_line(data = mean_d, linewidth = 0.6) +
        geom_point(data = mean_d, size = 1.9) +
        facet_grid(meas ~ s, labeller = labeller(meas = label_wrap_gen(24))) +
        scale_x_continuous(breaks = sort(unique(d$effect_size)),
                           labels = function(x) format(x, drop0trailing = TRUE),
                           expand = expansion(mult = 0.08)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        scale_colour_manual(values = method_colors, drop = FALSE) +
        scale_shape_manual(values = method_shapes, drop = FALSE) +
        labs(x = "Planted difference in methylation rate",
             y = "Share of planted regions detected", title = title,
             subtitle = paste0("Spike-in simulation, ", model, " efficiency model. ",
                               "Lines: mean of ", uniqueN(d$seed),
                               " seeds; small symbols: each seed.")) +
        theme_ms() +
        theme(strip.text.y = element_text(face = "plain", angle = 270, hjust = 0.5),
              panel.spacing = unit(12, "pt"))
    .save(p, name, 7.2, 4.6)
    .across(d, "sensitivity", c("model", "block", "scenario", "effect_size",
                                "method", "measure"))
}


# 3. Broad-difference trade-off ---------------------------------------------

.broad_arms <- function(model, file) {
    b <- .read(model, "broad_between", file)
    w <- .read(model, "broad_within", file)
    if (is.null(b) || is.null(w)) return(NULL)
    rbind(b[method == "raw"][, arm := "raw"],
          w[method == "SMFnorm"][, arm := "within"],
          b[method == "SMFnorm"][, arm := "both"], fill = TRUE)
}

plot_broad <- function(model, effect = 0.15) {
    n <- .broad_arms(model, "null_by_seed.csv")
    s <- .broad_arms(model, "spikein_by_seed.csv")
    if (is.null(n) || is.null(s)) return(NULL)
    s <- s[abs(effect_size - effect) < 1e-8]
    d <- rbind(
        n[, .(seed, scenario, arm, measure = "retained", value = domain_retained)],
        n[, .(seed, scenario, arm, measure = "off_domain", value = as.numeric(FP))],
        s[, .(seed, scenario = efficiency_scenario, arm, measure = "sensitivity",
              value = sensitivity)])
    d[, `:=`(a = .lab(arm, arm_labels), s = .lab(scenario, scenario_labels))]
    mean_d <- d[, .(value = mean(value)), by = .(measure, s, a)]
    dodge <- position_dodge(width = 0.55)

    panel <- function(meas, ylab, title, log = FALSE) {
        p <- ggplot(d[measure == meas], aes(s, value, colour = a, shape = a, group = a)) +
            geom_point(position = dodge, size = 0.9, alpha = 0.45) +
            geom_point(data = mean_d[measure == meas], position = dodge, size = 2.6) +
            scale_colour_manual(values = arm_colors, drop = FALSE) +
            scale_shape_manual(values = arm_shapes, drop = FALSE) +
            labs(x = "Share of genome that truly differs", y = ylab, title = title) +
            theme_ms() + theme(plot.title = element_text(size = rel(0.95)))
        if (log) {
            p + scale_y_log10(labels = scales::comma)
        } else {
            p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
        }
    }
    p <- panel("retained", "Share of domain windows with a call",
               "Real broad difference kept") +
        panel("off_domain", "Calls outside the domains (log scale)",
              "False calls elsewhere", log = TRUE) +
        panel("sensitivity", "Share of planted regions called",
              sprintf("Focal regions found (difference %.2f)", effect)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = "The between-group step removes real broad differences",
            subtitle = paste0("Group B raised by 0.15 over 1 Mb domains, ", model,
                              " efficiency model. Large symbols: mean of ",
                              uniqueN(d$seed), " seeds; small: each seed."),
            theme = theme_ms()) &
        theme(legend.position = "bottom")
    .save(p, paste0("fig_broad_tradeoff_", model), 9.5, 4)
    d[, model := model]
    .across(d, "value", c("model", "scenario", "arm", "measure"))
}


# Run -----------------------------------------------------------------------

tabs <- list(null = list(), sens = list(), broad = list())
for (model in models) {
    message("\n== ", model)
    tabs$null[[model]] <- plot_null(model)
    tabs$sens[[paste(model, "bias")]] <- plot_sensitivity(
        model, "bias_parametric", paste0("fig_spikein_sensitivity_", model),
        "Detecting planted differences under a between-group efficiency bias")
    tabs$sens[[paste(model, "within")]] <- plot_sensitivity(
        model, "control_within", paste0("fig_spikein_sensitivity_within_", model),
        "Detecting planted differences when efficiency differs within groups")
    tabs$broad[[model]] <- plot_broad(model)
}

.write <- function(x, name) {
    x <- rbindlist(x, fill = TRUE)
    if (!nrow(x)) return(invisible())
    num <- c("mean", "sd", "min", "max")
    x[, (num) := lapply(.SD, signif, 4), .SDcols = num]
    fwrite(x, file.path(out_dir, name))
    message("  wrote ", name)
}
message("\n== tables")
.write(tabs$null,  "table_null_false_positives.csv")
.write(tabs$sens,  "table_spikein_sensitivity.csv")
.write(tabs$broad, "table_broad_tradeoff.csv")
