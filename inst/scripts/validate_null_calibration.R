#!/usr/bin/env Rscript

# validate_null_calibration.R
#
# RUN THIS BEFORE ANY SWEEP OR HEADLINE EXPERIMENT.
#
# It answers three questions, in the order they can invalidate downstream work:
#
#   1. Does the new q-value plumbing actually work on THIS metilene build?
#      call_dmrs_metilene() reads metilene's "Number of Tests: N" line off
#      stderr and uses it as the BH denominator. If that line is worded
#      differently in your build, the code falls back to the emitted row count
#      and UNDER-CORRECTS by ~30x -- with only a warning. This script checks the
#      captured logs directly and fails loudly if the count was not parsed.
#
#   2. Is the FDR calibrated? The control_within scenarios are a TRUE NULL:
#      within-group efficiency variation only, matched group means, no spike-in.
#      Every method should call ~0 DMRs. Anything else means the significance
#      path is anti-conservative and every sweep built on it is worthless.
#
#   3. Does the rewritten ComBatMet run end-to-end on real data? These scenarios
#      have overlapping efficiency ranges, so the efficiency-stratum batch
#      crosses the group boundary and ComBat is applicable (unlike the aligned_*
#      scenarios, where it correctly refuses).
#
# Usage (from the package root):
#   Rscript inst/scripts/validate_null_calibration.R
#
# Start with FAST_MODE <- TRUE (chr1 only, ~20-40 min) to catch plumbing errors
# cheaply, then re-run with FAST_MODE <- FALSE for the real calibration verdict.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".")
} else {
    library(SMFsim)
}
suppressWarnings(suppressMessages(library(data.table)))

# Internal package function (visible under load_all; fall back to the namespace
# when SMFsim is installed rather than loaded from source).
.n_tests_of <- if (exists(".metilene_n_tests")) .metilene_n_tests else
    getFromNamespace(".metilene_n_tests", "SMFsim")

# === EDIT THESE FOR YOUR ENVIRONMENT =====================================
config <- parse_args()
config$data_dir      <- "../m-series-data/data/allc"
config$sample_sheet  <- "../m-series-data/data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "results/null_calibration"

FAST_MODE <- TRUE   # TRUE = chr1 only. Do this first, then set FALSE.
# =========================================================================

# True null: within-group efficiency variation, matched group means, no
# between-group artifact, so nothing for any method to legitimately call.
config$sim_mode            <- "parametric"
config$dispersion_s        <- 26
config$rate_between_groups <- FALSE
config$scenarios           <- c("mild", "moderate", "severe")
config$methods             <- c("raw", "downsampled", "SMFnorm", "ComBatMet")
config$seed                <- 42

if (FAST_MODE) {
    config$chr_pattern <- "^(chr)?1$"
    base_output        <- "results/null_calibration_fast"
    message("FAST_MODE: chr1 only -> ", base_output)
}
config$output_dir <- base_output
dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# --- Record the significance path actually in force ----------------------
message("\n", strrep("=", 64))
message("Significance settings (inherited from parse_args)")
message(strrep("=", 64))
for (k in c("metilene_min_diff", "metilene_min_effect", "metilene_qval",
            "metilene_q_source", "metilene_p_column", "metilene_mtc",
            "metilene_min_cpg", "min_coverage", "within_alpha", "dispersion_s")) {
    message(sprintf("  %-20s = %s", k, format(config[[k]])))
}

has_combat <- requireNamespace("ComBatMet", quietly = TRUE)
message("\n  ComBatMet installed   = ", has_combat)
if (!has_combat) {
    message("  !! ComBatMet missing - it will be skipped with a warning.")
    message("     Install: remotes::install_github('JmWangBio/ComBatMet')")
}

# --- Run the null --------------------------------------------------------
wt_reps <- prepare_wt_replicates(config)
null_dt <- run_null_simulation(wt_reps, config)

# --- CHECK 1: was metilene's test count actually parsed? -----------------
message("\n", strrep("=", 64))
message("CHECK 1: did we read metilene's 'Number of Tests' line?")
message(strrep("=", 64))

logs <- list.files(base_output, pattern = "^metilene_log\\.txt$",
                   recursive = TRUE, full.names = TRUE)
if (length(logs) == 0) {
    message("  !! No metilene_log.txt found under ", base_output)
    message("     Either no DMR calling ran, or stderr capture is not working.")
    check1 <- FALSE
} else {
    counts <- vapply(logs, .n_tests_of, integer(1))
    n_ok <- sum(!is.na(counts))
    message(sprintf("  logs found: %d,  test count parsed in: %d",
                    length(logs), n_ok))
    if (n_ok > 0) {
        message(sprintf("  parsed counts: min %d, median %d, max %d",
                        min(counts, na.rm = TRUE),
                        as.integer(stats::median(counts, na.rm = TRUE)),
                        max(counts, na.rm = TRUE)))
    }
    check1 <- n_ok == length(logs)
    if (!check1) {
        message("  !! FAILED for ", length(logs) - n_ok, " log(s).")
        message("     The BH denominator silently fell back to the emitted row")
        message("     count, which under-corrects by roughly 30x.")
        message("     Inspect one log and adjust .metilene_n_tests():")
        message("       ", logs[which(is.na(counts))[1]])
    } else {
        message("  OK - every metilene run reported a usable test count.")
    }
}

# --- CHECK 2: is the null clean? -----------------------------------------
message("\n", strrep("=", 64))
message("CHECK 2: false positives under a TRUE NULL (expect ~0)")
message(strrep("=", 64))

fp_file <- file.path(base_output, "null_simulation_results.csv")
if (file.exists(fp_file)) null_dt <- fread(fp_file)

if (!is.data.table(null_dt) || !"FP" %in% names(null_dt)) {
    message("  !! No usable null results to summarise.")
    check2 <- NA
} else {
    setorder(null_dt, scenario, method)
    print(null_dt[, .(scenario, method, FP)])

    by_method <- null_dt[, .(total_FP = sum(FP, na.rm = TRUE),
                             max_FP = max(FP, na.rm = TRUE)), by = method]
    message("\n  Totals by method:")
    print(by_method)

    # A handful of calls genome-wide is tolerable; hundreds is not.
    TOL <- 10L
    worst <- by_method[which.max(total_FP)]
    check2 <- worst$total_FP <= TOL
    if (check2) {
        message(sprintf(
            "\n  OK - worst method (%s) has %d FP across all null scenarios.",
            worst$method, worst$total_FP))
    } else {
        message(sprintf(
            "\n  !! FAILED - %s produced %d false positives in a TRUE NULL.",
            worst$method, worst$total_FP))
        message("     The significance path is anti-conservative. Do NOT run the")
        message("     alpha grid until this is resolved. First suspects:")
        message("       - metilene_p_column: should be '2dks' for de-novo mode")
        message("       - the BH denominator (see CHECK 1)")
        message("       - fall back to metilene_q_source = 'metilene' to compare")
    }
}

# --- CHECK 3: did ComBatMet run? -----------------------------------------
message("\n", strrep("=", 64))
message("CHECK 3: did the rewritten ComBatMet run end-to-end?")
message(strrep("=", 64))

if (!is.data.table(null_dt) || !"method" %in% names(null_dt)) {
    check3 <- NA
    message("  (no results to check)")
} else {
    check3 <- "ComBatMet" %in% unique(null_dt$method)
    if (check3) {
        message("  OK - ComBatMet produced results (it is present in the output).")
        message("  NOTE: these scenarios have OVERLAPPING efficiency ranges, so")
        message("  batch crosses the group boundary and ComBat is applicable.")
        message("  In the aligned_* scenarios it will correctly refuse to run.")
    } else {
        message("  !! ComBatMet is ABSENT from the results - it failed or was")
        message("     skipped. Check warnings() above for the reason.")
    }
}

# --- Verdict -------------------------------------------------------------
message("\n", strrep("=", 64))
message("VERDICT")
message(strrep("=", 64))
fmt <- function(x) if (isTRUE(x)) "PASS" else if (isFALSE(x)) "FAIL" else "N/A"
message("  1. metilene test count parsed : ", fmt(check1))
message("  2. null is clean              : ", fmt(check2))
message("  3. ComBatMet ran              : ", fmt(check3))

if (isTRUE(check1) && isTRUE(check2) && isTRUE(check3)) {
    if (FAST_MODE) {
        message("\n  All checks pass on chr1. Now set FAST_MODE <- FALSE and")
        message("  re-run for the genome-wide calibration verdict.")
    } else {
        message("\n  All checks pass genome-wide. Cleared to run the alpha grid.")
    }
} else {
    message("\n  Resolve the failures above BEFORE running the alpha grid.")
}
