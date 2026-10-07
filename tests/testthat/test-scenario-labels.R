# Regression tests for scenario ordering and labelling.
#
# Background: every scenario-faceted figure once hardcoded
#   scenario_order <- c("mild", "moderate", "severe")
#   factor(efficiency_scenario, levels = scenario_order)
# When the between-group bias scenarios (aligned_*/imbalanced_*) were added to
# get_efficiency_scenarios(), factor() silently mapped every row to NA, so all
# facets collapsed into a single empty "NA" panel. factor() does not warn on
# unmatched levels, so this was invisible until the figures were inspected.
# These tests pin the behaviour that replaced it.

test_that("every canonical scenario survives .scenario_factor() without NA", {
    scenarios <- names(get_efficiency_scenarios())
    expect_gt(length(scenarios), 0)

    f <- .scenario_factor(scenarios)

    expect_false(any(is.na(f)))
    expect_setequal(levels(f), scenarios)
})

test_that("the between-group bias scenarios round-trip (the original bug)", {
    bias <- c("aligned_moderate", "aligned_strong",
              "imbalanced_moderate", "imbalanced_strong")

    # These must be real scenarios, not just strings this test invented.
    expect_true(all(bias %in% names(get_efficiency_scenarios())))

    f <- .scenario_factor(bias)

    expect_equal(sum(is.na(f)), 0)
    expect_equal(as.character(f), bias)
})

test_that("levels follow get_efficiency_scenarios() order, not input order", {
    canonical <- names(get_efficiency_scenarios())
    shuffled <- rev(canonical)

    f <- .scenario_factor(shuffled)

    expect_equal(levels(f), canonical)
})

test_that("levels are restricted to the values actually present", {
    canonical <- names(get_efficiency_scenarios())
    subset <- canonical[c(1, length(canonical))]

    f <- .scenario_factor(subset)

    expect_equal(levels(f), canonical[canonical %in% subset])
    expect_false(any(is.na(f)))
})

test_that("unknown scenarios are kept and warned about, never silently NA", {
    x <- c("aligned_strong", "some_future_scenario")

    expect_warning(f <- .scenario_factor(x), "Unrecognised")

    # The whole point: an unrecognised value must still plot.
    expect_false(any(is.na(f)))
    expect_true("some_future_scenario" %in% levels(f))
    # Known scenarios sort ahead of unknown ones.
    expect_equal(levels(f), c("aligned_strong", "some_future_scenario"))
})

test_that("NA input stays NA without error or warning", {
    expect_silent(f <- .scenario_factor(c("mild", NA)))
    expect_true(is.na(f[2]))
    expect_false(is.na(f[1]))
})

test_that("scenario_labels covers every scenario in get_efficiency_scenarios()", {
    # Guards against adding a scenario without a display label, which would
    # otherwise surface as a raw snake_case name in a manuscript figure.
    expect_true(all(names(get_efficiency_scenarios()) %in% names(scenario_labels)))
})

test_that(".scenario_label() maps known scenarios and falls back to the raw name", {
    expect_equal(.scenario_label("aligned_strong"), "Aligned Bias (strong)")

    # Fallback keeps an unlabelled scenario readable instead of returning NA.
    expect_equal(.scenario_label("some_future_scenario"), "some_future_scenario")

    mixed <- .scenario_label(c("mild", "some_future_scenario"))
    expect_equal(mixed, c("Mild Within-Group", "some_future_scenario"))
    expect_false(any(is.na(mixed)))
})

test_that(".scenario_label() accepts a factor, as ggplot labellers pass one", {
    f <- .scenario_factor(c("aligned_strong", "aligned_moderate"))

    expect_equal(.scenario_label(f),
                 c("Aligned Bias (strong)", "Aligned Bias (moderate)"))
})


# --- imbalanced_* batch structure ----------------------------------------
#
# The imbalanced scenarios exist to give ComBatMet a batch that CROSSES the
# group boundary. When they were [min, max] ranges sampled with runif(), the
# ranges overlapped so little that all of A usually landed above all of B
# (P = 0.99 for strong, 0.65 for moderate), the median split reproduced the
# group labels, and ComBatMet was refused as inapplicable in every scenario of
# the bias block. These tests pin the explicit vectors that replaced them.

.median_split <- function(sc, n_reps = 4L) {
    a <- .resolve_efficiencies(sc$efficiency_A, n_reps, "A")
    b <- .resolve_efficiencies(sc$efficiency_B, n_reps, "B")
    .efficiency_batch(list(efficiency_A = a, efficiency_B = b), n_reps, n_reps)
}

test_that("imbalanced_* are explicit vectors, not sampled ranges", {
    sc <- get_efficiency_scenarios()
    for (nm in c("imbalanced_moderate", "imbalanced_strong")) {
        expect_length(sc[[nm]]$efficiency_A, 4L)
        expect_length(sc[[nm]]$efficiency_B, 4L)
    }
})

test_that("imbalanced_* hit their stated mean gap exactly", {
    sc <- get_efficiency_scenarios()
    gap <- function(s) mean(s$efficiency_A) - mean(s$efficiency_B)
    expect_equal(gap(sc$imbalanced_moderate), 0.13, tolerance = 1e-12)
    expect_equal(gap(sc$imbalanced_strong),   0.22, tolerance = 1e-12)
})

test_that("imbalanced_* median split crosses groups: A = HHHL, B = HLLL", {
    sc <- get_efficiency_scenarios()
    for (nm in c("imbalanced_moderate", "imbalanced_strong")) {
        batch <- as.character(.median_split(sc[[nm]]))
        expect_identical(sum(batch[1:4] == "high_eff"), 3L, info = nm)
        expect_identical(sum(batch[5:8] == "high_eff"), 1L, info = nm)
    }
})

test_that("aligned_* remain perfectly confounded for any draw", {
    # Non-overlapping ranges: confounded by construction, so ComBatMet is
    # correctly inapplicable there whatever the seed.
    sc <- get_efficiency_scenarios()
    for (nm in c("aligned_moderate", "aligned_strong")) {
        for (seed in 1:20) {
            set.seed(seed)
            batch <- as.character(.median_split(sc[[nm]]))
            expect_true(all(batch[1:4] == "high_eff") &&
                        all(batch[5:8] == "low_eff"),
                        info = paste(nm, "seed", seed))
        }
    }
})
