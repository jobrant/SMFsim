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
#      The check gates on the UN-NORMALIZED arms (raw, downsampled), which have
#      nothing to legitimately call and so must sit at ~0; anything else means
#      the significance path is anti-conservative and every sweep built on it
#      is worthless. The normalization arms (SMFnorm, ComBatMet) are REPORTED
#      but not gated -- their null FP counts depend on within_alpha and on how
#      far the efficiency batch split is confounded with group, which makes
#      them results rather than plumbing failures. See the CHECK 2 block.
#
#   3. Did every requested method actually run? run_all_methods() drops a
#      method that errored, so a failed arm is indistinguishable from one that
#      was never requested, and its row is simply missing from the FP table --
#      which reads as if it scored zero. This compares the output against
#      config$methods and fails if any arm is absent.
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
config$data_dir      <- "data/allc"
config$sample_sheet  <- "data/sample_sheet.csv"
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
base_output          <- "/orange/jobrant/yuki/smfnorm-sim-manuscript/results/null_calibration"

FAST_MODE <- FALSE   # TRUE = chr1 only. Do this first, then set FALSE.
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
    base_output        <- "/orange/jobrant/yuki/smfnorm-sim-manuscript/results/null_calibration_fast"
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

# --- CHECK 2: is the SIGNIFICANCE PATH calibrated? -----------------------
#
# This check gates on the UN-NORMALIZED arms only (raw, downsampled). They
# apply no normalization, so under a matched-mean null they have nothing to
# call: any DMR they produce is a defect in the significance path itself --
# the BH denominator, the p-value column, or the q source. That is the thing
# this script exists to validate, and it is the thing that must be correct
# before a grid is worth running.
#
# The normalization arms are deliberately NOT gated. Their null FP counts are
# properties of the methods, i.e. results (Figure 1), not plumbing failures:
#
#   - SMFnorm's count is a function of within_alpha. Shrinking within-group
#     variance lowers the noise floor, so residual artifact clears
#     significance. At within_alpha = 0.3 severe gave 104 FP; at 0.5 it gave 4.
#   - ComBatMet's count is a function of how far the efficiency batch split is
#     confounded with group. In `severe` the median split lands 1 high / 3 low
#     in A and 3 high / 1 low in B, i.e. batch nearly encodes group, and
#     removing it removes group signal.
#
# Gating on those blocked the alpha grid in the Aug 2026 run, where the real
# causes were a scenario-construction bug (rep_len() recycling skewed the
# group means, injecting a 0.075 between-group artifact into a "true null")
# plus an unlocked within_alpha inherited from the parse_args default. Neither
# was the significance path, which raw and downsampled showed was clean by
# sitting at exactly 0 throughout.
message("\n", strrep("=", 64))
message("CHECK 2: significance path under a TRUE NULL")
message(strrep("=", 64))

GATE_METHODS <- c("raw", "downsampled")
TOL <- 10L

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

    gated    <- by_method[method %in% GATE_METHODS]
    reported <- by_method[!method %in% GATE_METHODS]

    message(sprintf("\n  GATED (un-normalized; expect ~0, tolerance %d):", TOL))
    if (nrow(gated) == 0L) {
        message("    !! none of ", paste(GATE_METHODS, collapse = ", "),
                " are present in the results.")
    } else {
        print(gated)
    }

    if (nrow(reported)) {
        message("\n  REPORTED (normalization arms; not gated - these are results):")
        print(reported)
    }

    # Every gated arm must be present AND within tolerance. A missing gated arm
    # is a failure, not a pass: run_all_methods() drops a method that errored,
    # so absence means the check never ran.
    missing_gated <- setdiff(intersect(GATE_METHODS, config$methods),
                             by_method$method)
    if (length(missing_gated)) {
        message("\n  !! FAILED - gated method(s) absent from the results: ",
                paste(missing_gated, collapse = ", "))
        message("     run_all_methods() drops a method that errored, so this")
        message("     check did not actually run. See warnings() above.")
        check2 <- FALSE
    } else {
        worst   <- gated[which.max(total_FP)]
        check2  <- worst$total_FP <= TOL
        if (check2) {
            message(sprintf(
                "\n  OK - worst un-normalized arm (%s) has %d FP across all null",
                worst$method, worst$total_FP))
            message("       scenarios. The significance path is calibrated.")
        } else {
            message(sprintf(
                "\n  !! FAILED - %s produced %d false positives in a TRUE NULL",
                worst$method, worst$total_FP))
            message("     WITHOUT any normalization applied, so the significance")
            message("     path itself is anti-conservative. First suspects:")
            message("       - metilene_p_column: should be '2dks' for de-novo mode")
            message("       - the BH denominator (see CHECK 1)")
            message("       - fall back to metilene_q_source = 'metilene' to compare")
        }
    }

    if (nrow(reported)) {
        message("\n  Interpret the reported arms against the run's parameters,")
        message(sprintf("  not against zero (this run: within_alpha = %s).",
                        format(config$within_alpha)))
    }
}

# --- CHECK 3: did EVERY requested method run? ----------------------------
#
# run_all_methods() wraps each method in tryCatch, emits warning() on failure
# and returns NULL, which is then dropped from the results list. Under Rscript
# warning() is deferred to the end of the log, so a failed arm looks exactly
# like one that was never requested -- and if the deferred warnings are never
# flushed ("There were 41 warnings"), the reason is lost entirely.
#
# That is how the 2026-09-02 headline run completed "successfully" with SMFnorm
# absent from both blocks and ComBatMet absent from bias_parametric, having
# logged "Method SMFnorm completed in 9.5 seconds" (that message prints BEFORE
# the is.null() guard). Checking only ComBatMet, as this check used to, would
# not have caught it. Compare against config$methods, not against a hardcoded
# name.
message("\n", strrep("=", 64))
message("CHECK 3: did every requested method produce results?")
message(strrep("=", 64))

if (!is.data.table(null_dt) || !"method" %in% names(null_dt)) {
    check3 <- NA
    message("  (no results to check)")
} else {
    requested <- config$methods
    present   <- unique(null_dt$method)
    absent    <- setdiff(requested, present)

    message("  requested : ", paste(requested, collapse = ", "))
    message("  present   : ", paste(present, collapse = ", "))

    check3 <- length(absent) == 0L
    if (check3) {
        message("  OK - all ", length(requested), " requested methods are present.")
        if ("ComBatMet" %in% present) {
            message("  NOTE: these scenarios have OVERLAPPING efficiency ranges, so")
            message("  batch crosses the group boundary and ComBat is applicable.")
            message("  In the aligned_* scenarios it will correctly refuse to run.")
        }
    } else {
        message("  !! FAILED - ABSENT from the results: ",
                paste(absent, collapse = ", "))
        message("     Each was requested but produced nothing, so it failed and")
        message("     was silently dropped by run_all_methods(). Do not read the")
        message("     FP table as if that arm had scored zero - it did not run.")
        message("     The reason is in the deferred warnings; if the log ends in")
        message("     'There were N warnings', re-run with")
        message("       options(warn = 1)")
        message("     to interleave them with the block that produced them.")
    }
}

# --- Verdict -------------------------------------------------------------
message("\n", strrep("=", 64))
message("VERDICT")
message(strrep("=", 64))
fmt <- function(x) if (isTRUE(x)) "PASS" else if (isFALSE(x)) "FAIL" else "N/A"
message("  1. metilene test count parsed : ", fmt(check1))
message("  2. significance path clean    : ", fmt(check2))
message("  3. all methods ran            : ", fmt(check3))

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
