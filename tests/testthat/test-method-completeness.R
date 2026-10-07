# Regression tests for the run_all_methods() completeness backstop and for the
# matched-mean guard on the within-group scenarios.
#
# Background, both from the 2026-09-02 headline run:
#
#   1. run_all_methods() wrapped each method in tryCatch, warned on failure and
#      returned NULL, which was then dropped. The "Method X completed in N
#      seconds" message printed BEFORE the is.null() guard, so a failed method
#      still logged as completed. Under Rscript the warnings were deferred and
#      never flushed ("There were 41 warnings"), so the run finished with
#      SMFnorm -- the method the package exists to evaluate -- absent from every
#      block, and its rows simply missing from the FP table as if it had scored
#      zero.
#
#   2. The within-group scenarios (mild/moderate/severe) defined 3 efficiencies
#      per group but the M-series has 4 replicates. rep_len() duplicated the
#      FIRST element, shifting each group mean by a different amount. `severe`
#      acquired a 0.075 between-group efficiency difference with the sign
#      REVERSED relative to the authored design -- so the "true null" used to
#      validate FDR calibration contained a real between-group artifact.

.mk_simple_groups <- function(n_rep = 3, n_sites = 40) {
    d <- data.table::data.table(
        chr = "chr1", pos = seq_len(n_sites), strand = "+",
        site = seq_len(n_sites), mc = 5L, cov = 10L, rate = 0.5)
    mk <- function(prefix) stats::setNames(
        lapply(seq_len(n_rep), function(i) data.table::copy(d)),
        paste0(prefix, "_", seq_len(n_rep)))
    list(PseudoA = mk("PseudoA"), PseudoB = mk("PseudoB"),
         params = list(efficiency_A = rep(0.9, n_rep),
                       efficiency_B = rep(0.8, n_rep)))
}


test_that("run_all_methods returns every requested method on success", {
    res <- run_all_methods(.mk_simple_groups(),
                           methods = c("raw", "downsampled"))

    expect_setequal(names(res), c("raw", "downsampled"))
    expect_length(attr(res, "method_failures"), 0L)
})


test_that("run_all_methods errors rather than silently dropping a method", {
    expect_error(
        run_all_methods(.mk_simple_groups(), methods = c("raw", "bogus")),
        "produced no results")
})


test_that("the completeness error names the method and its reason", {
    err <- tryCatch(
        run_all_methods(.mk_simple_groups(), methods = c("raw", "bogus")),
        error = function(e) conditionMessage(e))

    # The reason must survive into the error. Relying on warning() lost it
    # entirely under Rscript.
    expect_match(err, "bogus")
    expect_match(err, "Unknown method: bogus")
    # And it must not read as a zero score.
    expect_match(err, "ABSENT from the output, not zero")
})


test_that("SMFSIM_ALLOW_METHOD_FAILURE=1 downgrades the error to a warning", {
    # Set directly rather than via withr, which is not a declared dependency.
    Sys.setenv(SMFSIM_ALLOW_METHOD_FAILURE = "1")
    on.exit(Sys.unsetenv("SMFSIM_ALLOW_METHOD_FAILURE"), add = TRUE)

    expect_warning(
        res <- run_all_methods(.mk_simple_groups(),
                               methods = c("raw", "bogus")),
        "produced no results")

    expect_named(res, "raw")
    expect_identical(attr(res, "method_failures")[["bogus"]],
                     "Unknown method: bogus")
})


test_that("within-group scenarios have exactly matched group means", {
    sc <- get_efficiency_scenarios()

    for (nm in c("mild", "moderate", "severe")) {
        s <- sc[[nm]]
        expect_true(isTRUE(s$matched_means),
                    info = paste(nm, "must declare matched_means"))
        expect_equal(mean(s$efficiency_A), mean(s$efficiency_B),
                     tolerance = 1e-12,
                     info = paste(nm, "group means must match exactly"))
    }
})


test_that("recycling an explicit vector preserves the authored mean", {
    # A length-4 scenario vector applied to a 3-replicate dataset (PrEC) must
    # still realize the authored mean. Plain rep_len() does not.
    sc <- get_efficiency_scenarios()

    for (nm in c("mild", "moderate", "severe")) {
        s <- sc[[nm]]
        for (n_reps in c(3L, 5L, 6L)) {
            a <- .resolve_efficiencies(s$efficiency_A, n_reps, "A")
            b <- .resolve_efficiencies(s$efficiency_B, n_reps, "B")

            expect_length(a, n_reps)
            expect_equal(mean(a), mean(s$efficiency_A), tolerance = 1e-12,
                         info = paste(nm, "group A at n_reps =", n_reps))
            expect_equal(mean(a), mean(b), tolerance = 1e-12,
                         info = paste(nm, "still matched at n_reps =", n_reps))
        }
    }
})


test_that("a matched-mean scenario that drifts is an error, not a warning", {
    reps <- stats::setNames(
        lapply(seq_len(4), function(i) data.table::data.table(
            chr = "chr1", pos = 1:40, strand = "+", site = 1:40,
            mc = 5L, cov = 10L, rate = 0.5)),
        paste0("R", seq_len(4)))

    # The pre-fix severe spec: 3 values, unequal authored means.
    expect_error(
        create_pseudo_groups(reps,
                             efficiency_A = c(0.55, 0.95, 0.75),
                             efficiency_B = c(0.90, 0.50, 0.80),
                             mode = "parametric", seed = 42,
                             dispersion_s = 26, matched_means = TRUE),
        "matched_means = TRUE")

    # Editing one group only must also be caught -- this is the case a guard
    # that inferred intent from the values would have silently skipped.
    sc <- get_efficiency_scenarios()
    expect_error(
        create_pseudo_groups(reps,
                             efficiency_A = c(0.55, 0.95, 0.75, 0.80),
                             efficiency_B = sc$severe$efficiency_B,
                             mode = "parametric", seed = 42,
                             dispersion_s = 26, matched_means = TRUE),
        "matched_means = TRUE")

    # The corrected vectors must pass.
    expect_no_error(
        create_pseudo_groups(reps,
                             efficiency_A = sc$severe$efficiency_A,
                             efficiency_B = sc$severe$efficiency_B,
                             mode = "parametric", seed = 42,
                             dispersion_s = 26, matched_means = TRUE))
})


test_that("range-based bias scenarios are not subject to the matched guard", {
    reps <- stats::setNames(
        lapply(seq_len(4), function(i) data.table::data.table(
            chr = "chr1", pos = 1:40, strand = "+", site = 1:40,
            mc = 5L, cov = 10L, rate = 0.5)),
        paste0("R", seq_len(4)))

    sc <- get_efficiency_scenarios()
    expect_false(isTRUE(sc$aligned_strong$matched_means))
    expect_no_error(
        create_pseudo_groups(reps,
                             efficiency_A = sc$aligned_strong$efficiency_A,
                             efficiency_B = sc$aligned_strong$efficiency_B,
                             mode = "parametric", seed = 42,
                             dispersion_s = 26))
})


# --- Inapplicable vs failed ----------------------------------------------
#
# ComBatMet correctly refuses on the aligned_* scenarios: their efficiency
# ranges do not overlap, so the median batch split reproduces the group labels
# and batch is perfectly confounded with group. The completeness backstop
# originally treated that refusal as a failure and killed the headline run, so
# the backstop was disabled wholesale with SMFSIM_ALLOW_METHOD_FAILURE=1 --
# which would also have waved through a genuine SMFnorm crash. These tests pin
# the distinction.

.mk_confounded <- function(aligned = TRUE, n = 300) {
    set.seed(1)
    base <- stats::runif(n, 0.15, 0.55)
    mk <- function(effs, prefix) stats::setNames(lapply(effs, function(e) {
        r <- pmin(pmax(base * e + stats::rnorm(n, 0, .02), 0.01), 0.99)
        data.table::data.table(chr = "chr1", pos = seq_len(n), strand = "+",
                               site = seq_len(n), cov = 20L, rate = r,
                               mc = as.integer(round(20 * r)))
    }), paste0(prefix, "_", seq_along(effs)))

    if (aligned) {   # non-overlapping, like aligned_strong
        eA <- c(0.95, 0.90, 0.88, 0.92); eB <- c(0.60, 0.57, 0.63, 0.58)
    } else {         # overlapping, like imbalanced_strong
        eA <- c(0.93, 0.70, 0.88, 0.75); eB <- c(0.80, 0.62, 0.85, 0.66)
    }
    list(PseudoA = mk(eA, "PseudoA"), PseudoB = mk(eB, "PseudoB"),
         params = list(efficiency_A = eA, efficiency_B = eB))
}


test_that(".method_inapplicable is an error subclass, so unhandled it still stops", {
    cond <- .method_inapplicable("nope", method = "X")
    expect_s3_class(cond, "smfsim_method_inapplicable")
    expect_s3_class(cond, "error")
    expect_error(stop(cond), "X: nope")
})


test_that("run_combatmet signals INAPPLICABLE, not a plain error, under confounding", {
    skip_if_not_installed("ComBatMet")
    expect_error(run_combatmet(.mk_confounded(aligned = TRUE)),
                 class = "smfsim_method_inapplicable")
})


test_that("an inapplicable method is skipped and recorded, and the run continues", {
    skip_if_not_installed("ComBatMet")
    res <- run_all_methods(.mk_confounded(aligned = TRUE),
                           methods = c("raw", "SMFnorm", "ComBatMet"))

    expect_setequal(names(res), c("raw", "SMFnorm"))
    expect_named(attr(res, "method_inapplicable"), "ComBatMet")
    expect_match(attr(res, "method_inapplicable")[["ComBatMet"]],
                 "confounded with group")
    expect_length(attr(res, "method_failures"), 0L)
})


test_that("the same method runs normally when batch crosses the group boundary", {
    skip_if_not_installed("ComBatMet")
    res <- run_all_methods(.mk_confounded(aligned = FALSE),
                           methods = c("raw", "ComBatMet"))

    expect_setequal(names(res), c("raw", "ComBatMet"))
    expect_length(attr(res, "method_inapplicable"), 0L)
})


test_that("an inapplicable arm does not mask a genuine failure elsewhere", {
    skip_if_not_installed("ComBatMet")
    # ComBatMet inapplicable AND a real failure in the same scenario: the run
    # must still abort, and must name only the real failure as missing.
    err <- tryCatch(
        run_all_methods(.mk_confounded(aligned = TRUE),
                        methods = c("raw", "ComBatMet", "bogus")),
        error = function(e) conditionMessage(e))

    expect_match(err, "1 of 3 requested method\\(s\\) produced no results")
    expect_match(err, "bogus")
    expect_no_match(err, "- ComBatMet")
})


test_that(".applicability_rows distinguishes applied / inapplicable / failed", {
    res <- structure(
        list(raw = list(), SMFnorm = list()),
        method_inapplicable = c(ComBatMet = "confounded"),
        method_failures     = c(bogus = "Unknown method: bogus"))

    dt <- .applicability_rows("aligned_strong",
                              c("raw", "SMFnorm", "ComBatMet", "bogus"), res)

    expect_identical(dt$status,
                     c("applied", "applied", "inapplicable", "failed"))
    expect_identical(dt$reason[3], "confounded")
    expect_identical(unique(dt$scenario), "aligned_strong")
})
