# Regression tests for the ComBatMet comparator.
#
# Background: run_combatmet once hand-rolled a "mean-only" adjustment that
# subtracted each group's own mean and added the grand mean. That forces every
# group mean to the grand mean at every site, so the A-vs-B difference is
# identically zero and no DMR can ever be called -- the method scored 0
# sensitivity in every scenario, including the controls. It was also not
# ComBat-met at all. These tests pin the replacement.

.mk_pseudo_groups <- function(eff_A, eff_B, n = 200, effect_size = 0.25,
                              n_spiked = 50, seed = 1) {
    set.seed(seed)
    base <- stats::runif(n, 0.15, 0.55)
    effect <- c(rep(effect_size, n_spiked), rep(0, n - n_spiked))

    mk <- function(rate) {
        cov <- rep(20L, n)
        data.table::data.table(
            chr = "chr1", pos = seq_len(n), cov = cov,
            rate = rate, mc = as.integer(round(cov * rate)))
    }
    clamp <- function(x) pmin(pmax(x, 0.01), 0.99)

    A <- lapply(eff_A, function(e) mk(clamp(base * e + stats::rnorm(n, 0, .02))))
    B <- lapply(eff_B, function(e) mk(clamp((base + effect) * e + stats::rnorm(n, 0, .02))))
    names(A) <- paste0("PseudoA_", seq_along(eff_A))
    names(B) <- paste0("PseudoB_", seq_along(eff_B))

    list(PseudoA = A, PseudoB = B,
         params = list(efficiency_A = eff_A, efficiency_B = eff_B))
}

.group_diff <- function(split_data) {
    gm <- function(lst) rowMeans(do.call(cbind, lapply(lst, function(d) d$rate)))
    gm(split_data$PseudoB) - gm(split_data$PseudoA)
}


test_that(".efficiency_batch splits replicates at the median efficiency", {
    b <- .efficiency_batch(
        list(efficiency_A = c(0.72, 0.83, 0.95),
             efficiency_B = c(0.50, 0.61, 0.73)), 3, 3)

    expect_s3_class(b, "factor")
    expect_length(b, 6)
    expect_setequal(levels(b), c("low_eff", "high_eff"))
    # Overlapping ranges: the split must cross the group boundary.
    expect_true(length(unique(b[1:3])) > 1 || length(unique(b[4:6])) > 1)
})

test_that(".efficiency_batch errors when efficiencies do not match replicates", {
    expect_error(
        .efficiency_batch(list(efficiency_A = c(0.8, 0.9),
                               efficiency_B = c(0.6, 0.7)), 3, 3),
        "replicates")
})

test_that("run_combatmet aborts when batch is confounded with group", {
    # Non-overlapping efficiency ranges (the aligned_* scenarios): the median
    # split reproduces the group labels exactly, so ComBat cannot separate
    # artifact from biology. Must fail loudly, not return centred data.
    pg <- .mk_pseudo_groups(c(0.85, 0.90, 0.95), c(0.55, 0.60, 0.65))

    expect_error(run_combatmet(pg), "confounded with group")
})

test_that("run_combatmet preserves biological signal when batch is separable", {
    skip_if_not_installed("ComBatMet")

    # Overlapping efficiency ranges (the imbalanced_* scenarios).
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))
    res <- suppressMessages(run_combatmet(pg))

    expect_equal(res$method, "ComBatMet")
    expect_length(res$data$PseudoA, 3)
    expect_length(res$data$PseudoB, 3)

    after <- .group_diff(res$data)

    # THE REGRESSION: the old implementation drove every group difference to
    # ~0 (max |diff| ~2.5e-3, far under metilene's 0.1 min_diff), so no site
    # was ever callable. A working comparator must leave real signal standing.
    expect_gt(max(abs(after)), 0.1)
    expect_gt(sum(abs(after) > 0.1), 0)

    # Spiked sites must retain more signal than null sites.
    expect_gt(mean(abs(after[1:50])), mean(abs(after[51:200])))
})

test_that("run_combatmet accepts an explicit batch and keeps sample order", {
    skip_if_not_installed("ComBatMet")

    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))
    res <- suppressMessages(
        run_combatmet(pg, batch = factor(c("b1", "b2", "b1", "b2", "b1", "b2"))))

    expect_named(res$data$PseudoA, paste0("PseudoA_", 1:3))
    expect_named(res$data$PseudoB, paste0("PseudoB_", 1:3))
    # Coverage is untouched; only rate (and the derived mc) are adjusted.
    expect_equal(res$data$PseudoA[[1]]$cov, pg$PseudoA[[1]]$cov)
    expect_equal(res$data$PseudoA[[1]]$pos, pg$PseudoA[[1]]$pos)
})

test_that("run_combatmet rejects a batch with only one level", {
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))

    expect_error(
        run_combatmet(pg, batch = factor(rep("only_one", 6))),
        "confounded with group")
})
