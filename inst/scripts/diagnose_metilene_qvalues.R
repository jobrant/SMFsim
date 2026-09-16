#!/usr/bin/env Rscript

# diagnose_metilene_qvalues.R
#
# Settles empirically -- against your own metilene output, no source reading --
# three things the documentation and the source disagree about:
#
#   1. WHICH p-value does metilene's column-4 q correct?
#      The manual says Mann-Whitney U (column 7). The source reportedly
#      corrects the 2D-KS value (column 8) in de-novo mode.
#   2. WHAT denominator does it use?
#      Not the emitted row count. Back it out exactly.
#   3. IS `-d` applied after q-values are computed?
#      If so, q-values for a given region are identical at -d 0 and -d 0.1,
#      and the emitted row count never entered the correction at all.
#
# HOW IT WORKS
#   Under Bonferroni (-c 1), q = min(1, p * numberTests) for every uncapped row.
#   So q/p is a CONSTANT equal to numberTests -- but only for the p-value
#   metilene actually corrected. The other column gives noise. Whichever column
#   yields a single clean integer answers questions 1 and 2 simultaneously.
#
# USAGE
#   Produce a Bonferroni run at -d 0 (and optionally a -d 0.1 run for check 3):
#     metilene -a A -b B -m 10 -M 300 -d 0   -c 1 merged.bed > d0_bonf.bed
#     metilene -a A -b B -m 10 -M 300 -d 0.1 -c 1 merged.bed > d01_bonf.bed
#   then:
#     Rscript inst/scripts/diagnose_metilene_qvalues.R d0_bonf.bed [d01_bonf.bed]

suppressWarnings(suppressMessages(library(data.table)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: diagnose_metilene_qvalues.R <bonferroni.bed> [filtered.bed]")
}
f_all <- args[1]
f_filt <- if (length(args) >= 2) args[2] else NA_character_

COLS <- c("chr", "start", "stop", "q_metilene", "meandiff",
          "nCpG", "p_mwu", "p_ks", "mean_g1", "mean_g2")

read_metilene <- function(path) {
    d <- fread(path, header = FALSE)
    if (ncol(d) < 8) stop("Expected >=8 columns in ", path, ", got ", ncol(d))
    setnames(d, seq_len(min(ncol(d), length(COLS))), COLS[seq_len(min(ncol(d), length(COLS)))])
    d[]
}

d <- read_metilene(f_all)
cat("Rows in", f_all, ":", nrow(d), "\n\n")

# --- 1 + 2: which p-value, and what denominator ---------------------------
cat("=== Which p-value does metilene correct, and with what denominator? ===\n")
cat("Under -c 1 (Bonferroni), q/p is constant == numberTests for the p-value\n")
cat("metilene actually used. Rows with q == 1 are capped and excluded.\n\n")

uncapped <- d[is.finite(q_metilene) & q_metilene < 1]
cat("Uncapped rows:", nrow(uncapped), "of", nrow(d), "\n\n")

if (nrow(uncapped) == 0) {
    cat("!! Every q is capped at 1 -- cannot back out the denominator.\n")
    cat("   Re-run on a dataset with stronger signal, or with -c 1.\n")
} else {
    # Judge constancy on RELATIVE spread, not on the count of distinct rounded
    # values. metilene prints p-values at limited precision, so a genuinely
    # constant ratio near 1e6 still rounds to dozens of distinct integers --
    # what matters is that (max-min)/median is tiny.
    for (pcol in c("p_mwu", "p_ks")) {
        r <- uncapped[[pcol]]
        ok <- is.finite(r) & r > 0
        if (!any(ok)) { cat(pcol, ": no usable p-values\n"); next }
        ratio <- uncapped$q_metilene[ok] / r[ok]
        med <- stats::median(ratio)
        rel <- (max(ratio) - min(ratio)) / med

        cat(sprintf("%-6s ratio q/p: min=%.6g  max=%.6g  median=%.6g\n",
                    pcol, min(ratio), max(ratio), med))
        cat(sprintf("        relative spread = %.3g\n", rel))
        if (is.finite(rel) && rel < 1e-3) {
            cat(sprintf("        --> CONSTANT.\n"))
            cat(sprintf("        ==> metilene corrected %s, numberTests = %d\n",
                        pcol, as.integer(round(med))))
        } else {
            cat("        (varies by orders of magnitude -- NOT this column)\n")
        }
        cat("\n")
    }
    cat("Compare the constant above with the 'Number of Tests: N' line metilene\n")
    cat("prints to stderr. They should agree.\n\n")
}

# --- 3: is -d applied after q-values are computed? ------------------------
if (!is.na(f_filt) && file.exists(f_filt)) {
    cat("=== Is -d applied AFTER q-values are computed? ===\n")
    f <- read_metilene(f_filt)
    cat("Rows:", nrow(d), "(unfiltered) vs", nrow(f), "(filtered)\n")

    key <- c("chr", "start", "stop")
    m <- merge(d[, c(key, "q_metilene"), with = FALSE],
               f[, c(key, "q_metilene"), with = FALSE],
               by = key, suffixes = c("_all", "_filt"))
    cat("Regions present in both:", nrow(m), "of", nrow(f), "filtered rows\n")

    if (nrow(m) > 0) {
        same <- isTRUE(all.equal(m$q_metilene_all, m$q_metilene_filt,
                                 tolerance = 1e-8))
        cat("q-values identical for shared regions:", same, "\n")
        if (same) {
            cat("  ==> -d filters OUTPUT ONLY. The emitted row count never\n")
            cat("      entered the correction; the denominator is unaffected by -d.\n")
        } else {
            cat("  ==> q-values DIFFER, so -d changes the correction itself.\n")
            cat("      Max |difference|: ",
                max(abs(m$q_metilene_all - m$q_metilene_filt)), "\n")
        }
    }
    cat("\n")
}

# --- Practical consequence ------------------------------------------------
cat("=== Denominator options, for reference ===\n")
cat("emitted rows in this file      :", nrow(d), "\n")
cat("(use the backed-out numberTests above as the BH denominator:\n")
cat("   p.adjust(p, 'BH', n = numberTests)  )\n\n")
cat("NOTE: correcting the denominator does NOT fix p-value calibration.\n")
cat("metilene flattens CpG x replicate into the rank test, so a 40-CpG 5v5\n")
cat("region is tested as n=200 vs 200 rather than 5 vs 5, and boundaries were\n")
cat("chosen by segmentation to maximise the same signal being tested. Validate\n")
cat("against a TRUE NULL (expect ~0 calls) before trusting any absolute FDR.\n")
