#!/usr/bin/env Rscript

# retabulate_metilene.R
#
# Recompute DMR significance directly from existing metilene_dmrs.bed files and
# re-tabulate FP / TP / FN metrics, WITHOUT re-running metilene. Use this to
# correct a run made before R-side q-value filtering was added (metilene's own
# column-4 q-value is unreliable).
#
# For each DMR it recomputes a Benjamini-Hochberg q-value from metilene's
# Mann-Whitney U p-value (column 7), correcting across the candidate DMRs within
# each comparison, then keeps DMRs with q < Q_CUTOFF.
#
# Walks <ROOT>/<block>/ and writes, per block (originals are left untouched):
#   null_results_qBH.csv      scenario, method, FP
#   spikein_results_qBH.csv   scenario, effect_size, method, TP/FP/FN + metrics
#
# Usage (from the package root):
#   Rscript inst/scripts/retabulate_metilene.R
#   # or override the root / cutoff at the top of this file

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".")
} else {
    library(SMFsim)
}
suppressWarnings(suppressMessages(library(data.table)))

# === EDIT IF NEEDED ======================================================
ROOT     <- "results/bias_experiment"
Q_CUTOFF <- 0.05
# =========================================================================

# Internal package functions (available under load_all; fall back to namespace).
.classify <- if (exists("classify_dmrs")) classify_dmrs else
    getFromNamespace("classify_dmrs", "SMFsim")
.metrics  <- if (exists("compute_metrics")) compute_metrics else
    getFromNamespace("compute_metrics", "SMFsim")

# Read a metilene_dmrs.bed, recompute BH q from the MWU p-value (col 7), and
# return only the significant DMRs.
read_significant_dmrs <- function(bed, q_cutoff) {
    empty <- data.table(chr = character(), start = integer(), end = integer(),
                        q_value = numeric(), mean_diff = numeric(),
                        n_sites = integer())
    if (!file.exists(bed) || file.size(bed) == 0) return(empty)
    d <- tryCatch(fread(bed, header = FALSE), error = function(e) NULL)
    if (is.null(d) || ncol(d) < 8) return(empty)
    setnames(d, 1:8, c("chr", "start", "end", "q_metilene", "mean_diff",
                       "n_sites", "p_mwu", "p_2dks"))
    d[, q_value := p.adjust(p_mwu, method = "BH")]
    d[is.finite(q_value) & q_value < q_cutoff,
      .(chr, start, end, q_value, mean_diff, n_sites)]
}

dirs_of <- function(path) {
    if (!dir.exists(path)) return(character(0))
    list.dirs(path, recursive = FALSE)
}

# --- Per-block re-tabulation ---------------------------------------------

retab_null <- function(block_dir) {
    null_dir <- file.path(block_dir, "null")
    rows <- list()
    for (sc in dirs_of(null_dir)) {
        for (m in dirs_of(sc)) {
            dmrs <- read_significant_dmrs(
                file.path(m, "metilene_dmrs.bed"), Q_CUTOFF)
            rows[[length(rows) + 1]] <- data.table(
                scenario = basename(sc), method = basename(m),
                FP = nrow(dmrs))
        }
    }
    if (!length(rows)) return(NULL)
    rbindlist(rows)
}

retab_spikein <- function(block_dir) {
    spk_dir <- file.path(block_dir, "spikein")
    regions_file <- file.path(block_dir, "spikein_regions.csv")
    if (!dir.exists(spk_dir)) return(NULL)
    if (!file.exists(regions_file)) {
        message("  [skip spikein] no spikein_regions.csv in ", block_dir)
        return(NULL)
    }
    regions <- fread(regions_file)
    rows <- list()
    for (sc in dirs_of(spk_dir)) {
        for (eff in dirs_of(sc)) {
            effect_size <- suppressWarnings(
                as.numeric(sub("^effect_", "", basename(eff))))
            for (m in dirs_of(eff)) {
                dmrs <- read_significant_dmrs(
                    file.path(m, "metilene_dmrs.bed"), Q_CUTOFF)
                ev  <- .classify(dmrs, regions)
                met <- .metrics(ev$summary)
                rows[[length(rows) + 1]] <- data.table(
                    scenario = basename(sc), effect_size = effect_size,
                    method = basename(m),
                    TP = as.integer(met$TP), FP = as.integer(met$FP),
                    FN = as.integer(met$FN),
                    sensitivity = as.numeric(met$sensitivity),
                    precision = as.numeric(met$precision),
                    FDR = as.numeric(met$FDR), F1 = as.numeric(met$F1),
                    total_called = as.integer(met$total_called))
            }
        }
    }
    if (!length(rows)) return(NULL)
    rbindlist(rows)
}

# --- Run over all blocks -------------------------------------------------

if (!dir.exists(ROOT)) stop("ROOT not found: ", ROOT)
blocks <- dirs_of(ROOT)
if (!length(blocks)) stop("No block directories under: ", ROOT)

message(sprintf("Re-tabulating with BH q < %.3g over: %s", Q_CUTOFF, ROOT))

for (block in blocks) {
    bn <- basename(block)
    message("\n=== block: ", bn, " ===")

    null_dt <- retab_null(block)
    if (!is.null(null_dt)) {
        out <- file.path(block, "null_results_qBH.csv")
        fwrite(null_dt, out)
        message("  null -> ", out)
        print(null_dt[order(scenario, method)])
    }

    spk_dt <- retab_spikein(block)
    if (!is.null(spk_dt)) {
        out <- file.path(block, "spikein_results_qBH.csv")
        fwrite(spk_dt, out)
        message("  spikein -> ", out)
        print(spk_dt[order(scenario, effect_size, method),
                     .(scenario, effect_size, method, TP, FP, FN,
                       sensitivity, precision, F1)])
    }
}

message("\nDone. Corrected summaries written as *_qBH.csv beside the originals.")
