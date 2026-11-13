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
             efficiency_B = c(0.50, 0.61, 0.73)), 3, 3, rule = "median")

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
    skip_if_not_installed("ComBatMet")
    pg <- .mk_pseudo_groups(c(0.85, 0.90, 0.95), c(0.55, 0.60, 0.65))

    expect_error(run_combatmet(pg, batch_rule = "median"),
                 "confounded with group")
})

test_that("run_combatmet preserves biological signal when batch is separable", {
    skip_if_not_installed("ComBatMet")

    # Overlapping efficiency ranges (the imbalanced_* scenarios).
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))
    res <- suppressMessages(run_combatmet(pg, batch_rule = "median"))

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
    skip_if_not_installed("ComBatMet")
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))

    expect_error(
        run_combatmet(pg, batch = factor(rep("only_one", 6))),
        "confounded with group")
})


# --- Batch labelling rules -----------------------------------------------
#
# ComBat needs a discrete batch label; SMFsim derives one from the simulated
# efficiencies. The original median split of all replicates is confounded with
# group whenever the groups differ in efficiency, and in the matched-mean
# `severe` null it produced 7,791 false-positive DMRs against 11 with balanced
# labels (diagnose_combatmet_batch_labels.R, 2026-11-13). Both rules are kept
# as separate method arms; balanced is the default.

.batch_counts <- function(batch, n_A, n_B) {
    b <- as.character(batch)
    c(A_high = sum(b[seq_len(n_A)] == "high_eff"),
      B_high = sum(b[n_A + seq_len(n_B)] == "high_eff"))
}

test_that("the default batch rule is balanced", {
    expect_identical(eval(formals(.efficiency_batch)$rule)[1], "balanced")
    expect_identical(eval(formals(run_combatmet)$batch_rule)[1], "balanced")
})

test_that("balanced labels put the more efficient half of EACH group in high_eff", {
    eA <- c(0.55, 0.95, 0.75, 0.71); eB <- c(0.90, 0.50, 0.80, 0.76)   # severe
    b <- .efficiency_batch(list(efficiency_A = eA, efficiency_B = eB), 4, 4,
                           rule = "balanced")

    expect_identical(as.character(b),
                     c("low_eff", "high_eff", "high_eff", "low_eff",
                       "high_eff", "low_eff", "high_eff", "low_eff"))
    expect_identical(levels(b), c("low_eff", "high_eff"))
})

test_that("balanced labels hold equal numbers of each group in every scenario", {
    sc <- get_efficiency_scenarios()
    for (nm in names(sc)) {
        set.seed(1)
        a <- suppressMessages(.resolve_efficiencies(sc[[nm]]$efficiency_A, 4L, "A"))
        b <- suppressMessages(.resolve_efficiencies(sc[[nm]]$efficiency_B, 4L, "B"))
        batch <- .efficiency_batch(list(efficiency_A = a, efficiency_B = b), 4, 4,
                                   rule = "balanced")
        expect_identical(unname(.batch_counts(batch, 4, 4)), c(2L, 2L), info = nm)
    }
})

test_that("balanced labels are never confounded, including where median labels are", {
    eA <- c(0.95, 0.90, 0.88, 0.92); eB <- c(0.60, 0.57, 0.63, 0.58)   # aligned
    p <- list(efficiency_A = eA, efficiency_B = eB)

    med <- .efficiency_batch(p, 4, 4, rule = "median")
    bal <- .efficiency_batch(p, 4, 4, rule = "balanced")

    expect_identical(unname(.batch_counts(med, 4, 4)), c(4L, 0L))
    expect_identical(unname(.batch_counts(bal, 4, 4)), c(2L, 2L))
})

test_that("balanced labels work with odd and tied replicate efficiencies", {
    b3 <- .efficiency_batch(list(efficiency_A = c(0.9, 0.8, 0.7),
                                 efficiency_B = c(0.6, 0.5, 0.4)), 3, 3,
                            rule = "balanced")
    expect_identical(unname(.batch_counts(b3, 3, 3)), c(1L, 1L))

    tied <- .efficiency_batch(list(efficiency_A = rep(0.8, 4),
                                   efficiency_B = rep(0.8, 4)), 4, 4,
                              rule = "balanced")
    expect_identical(unname(.batch_counts(tied, 4, 4)), c(2L, 2L))
})

test_that("an unknown batch rule is an error", {
    expect_error(.efficiency_batch(list(efficiency_A = c(0.8, 0.9),
                                        efficiency_B = c(0.6, 0.7)), 2, 2,
                                   rule = "random"))
})

# Matched group means with the `severe` efficiencies: no true group difference.
.mk_severe_null <- function(n = 400, seed = 9) {
    set.seed(seed)
    base <- stats::runif(n, 0.15, 0.55)
    eA <- c(0.55, 0.95, 0.75, 0.71); eB <- c(0.90, 0.50, 0.80, 0.76)
    mk <- function(effs, prefix) stats::setNames(lapply(effs, function(e) {
        r <- pmin(pmax(base * e + stats::rnorm(n, 0, .02), 0.01), 0.99)
        data.table::data.table(chr = "chr1", pos = seq_len(n), strand = "+",
                               site = seq_len(n), cov = 20L, rate = r,
                               mc = as.integer(round(20 * r)))
    }), paste0(prefix, "_", seq_along(effs)))
    list(PseudoA = mk(eA, "PseudoA"), PseudoB = mk(eB, "PseudoB"),
         params = list(efficiency_A = eA, efficiency_B = eB))
}

test_that("the two ComBatMet arms use different labels and give different data", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_severe_null()
    res <- suppressMessages(run_all_methods(
        pg, methods = c("ComBatMet_balanced", "ComBatMet_median")))

    expect_named(res, c("ComBatMet_balanced", "ComBatMet_median"))
    expect_identical(res$ComBatMet_balanced$method, "ComBatMet_balanced")
    expect_identical(res$ComBatMet_median$method, "ComBatMet_median")
    expect_false(isTRUE(all.equal(res$ComBatMet_balanced$data$PseudoA[[1]]$rate,
                                  res$ComBatMet_median$data$PseudoA[[1]]$rate)))
})

test_that("plain ComBatMet is the balanced arm", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_severe_null()
    res <- suppressMessages(run_all_methods(
        pg, methods = c("ComBatMet", "ComBatMet_balanced")))

    expect_identical(res$ComBatMet$method, "ComBatMet")
    for (g in c("PseudoA", "PseudoB")) {
        for (i in seq_along(res$ComBatMet$data[[g]])) {
            expect_equal(res$ComBatMet$data[[g]][[i]]$rate,
                         res$ComBatMet_balanced$data[[g]][[i]]$rate)
        }
    }
})

test_that("median labels create a group gap in a matched-mean null; balanced do not", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_severe_null()
    gap <- function(d) abs(mean(.group_diff(d)))
    raw_gap <- gap(pg)

    bal <- suppressMessages(run_combatmet(pg, batch_rule = "balanced"))
    med <- suppressMessages(run_combatmet(pg, batch_rule = "median"))

    expect_lt(raw_gap, 0.005)                # groups start with matched means
    expect_lt(gap(bal$data), 0.01)
    expect_gt(gap(med$data), 5 * gap(bal$data))
    expect_gt(gap(med$data), 0.02)
})

test_that("balanced runs where median is inapplicable (aligned efficiencies)", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_pseudo_groups(c(0.85, 0.90, 0.95, 0.88), c(0.55, 0.60, 0.65, 0.58))
    res <- suppressMessages(run_all_methods(
        pg, methods = c("raw", "ComBatMet_balanced", "ComBatMet_median")))

    expect_setequal(names(res), c("raw", "ComBatMet_balanced"))
    expect_named(attr(res, "method_inapplicable"), "ComBatMet_median")
    expect_length(attr(res, "method_failures"), 0L)
})

test_that("a caller-supplied batch with an unused level is treated as one level", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))
    one <- factor(rep("low_eff", 6), levels = c("low_eff", "high_eff"))

    expect_error(run_combatmet(pg, batch = one),
                 class = "smfsim_method_inapplicable")
})

test_that("a batch of fewer than three samples is inapplicable, not a crash", {
    skip_if_not_installed("ComBatMet")
    # Three replicates per group: balanced labels give a batch of two samples,
    # which ComBat_met cannot fit ("subscript out of bounds").
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95), c(0.50, 0.61, 0.73))

    expect_error(run_combatmet(pg, batch_rule = "balanced"),
                 class = "smfsim_method_inapplicable")
    expect_error(run_combatmet(pg, batch_rule = "balanced"),
                 "at least 3 samples in every batch")

    res <- suppressMessages(run_all_methods(
        pg, methods = c("raw", "ComBatMet_balanced", "ComBatMet_median")))
    expect_setequal(names(res), c("raw", "ComBatMet_median"))
    expect_named(attr(res, "method_inapplicable"), "ComBatMet_balanced")
})

test_that("balanced labels run with four replicates per group", {
    skip_if_not_installed("ComBatMet")
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95, 0.80), c(0.50, 0.61, 0.73, 0.60))

    res <- suppressMessages(run_combatmet(pg, batch_rule = "balanced"))

    expect_length(res$data$PseudoA, 4)
    expect_length(res$data$PseudoB, 4)
    expect_false(isTRUE(all.equal(res$data$PseudoA[[1]]$rate,
                                  pg$PseudoA[[1]]$rate)))
})

test_that("balanced labels leave a group-aligned efficiency gap in place", {
    skip_if_not_installed("ComBatMet")
    # Group A is more efficient than group B (mean 0.825 vs 0.61). Balanced
    # batches are independent of group, so ComBat has no batch term that can
    # absorb the gap: unspiked sites keep most of their artifactual difference.
    pg <- .mk_pseudo_groups(c(0.72, 0.83, 0.95, 0.80), c(0.50, 0.61, 0.73, 0.60))
    null_sites <- 51:200
    raw_gap <- mean(.group_diff(pg)[null_sites])

    res <- suppressMessages(run_combatmet(pg, batch_rule = "balanced"))
    adj_gap <- mean(.group_diff(res$data)[null_sites])

    expect_lt(raw_gap, -0.05)                 # B lower: the artifact is there
    expect_lt(adj_gap, 0.7 * raw_gap)         # and at least 70% of it remains
})

test_that("every method arm has a colour, a label and a point shape", {
    arms <- c("raw", "downsampled", "SMFnorm", "ComBatMet",
              "ComBatMet_balanced", "ComBatMet_median")
    expect_true(all(arms %in% names(method_colors)))
    expect_true(all(arms %in% names(method_labels)))
    expect_true(all(arms %in% names(method_shapes)))
    expect_false(method_colors[["ComBatMet_balanced"]] ==
                     method_colors[["ComBatMet_median"]])
    expect_false(method_labels[["ComBatMet_balanced"]] ==
                     method_labels[["ComBatMet_median"]])
})
