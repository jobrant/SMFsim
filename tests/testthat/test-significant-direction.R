# Significance-only sensitivity and the planted direction of spike-in regions.
#
# A planted effect of 0.10 is observed as 0.10 x efficiency, so at low
# efficiency a region can have a significant segment that the 0.1 effect cut
# removes. Reported sensitivity then falls with efficiency although detection
# does not (aligned_strong, effect 0.10: 0.13 called against 0.37 significant).
# The direction of each planted effect was also not saved, so a detected region
# could not be told from a call an efficiency bias produced.

skip_if_not_installed("GenomicRanges")

.sd_truth <- function() data.table::data.table(
    region_id = 1:4, chr = c("chr1", "chr1", "chr1", "chr2"),
    start = c(1000L, 5000L, 9000L, 1000L),
    end   = c(2000L, 6000L, 10000L, 2000L),
    density_stratum = c("low", "low", "high", "high"))


test_that("call_dmrs_metilene keeps the significant segments it cuts on effect size", {
    # Two segments with tiny p-values: one of 0.15 and one of 0.06; a third is
    # not significant.
    bed <- c("chr1\t100\t400\t1\t0.15\t10\t1e-12\t1e-12\t0.4\t0.25",
             "chr1\t600\t900\t1\t-0.06\t10\t1e-12\t1e-12\t0.2\t0.26",
             "chr1\t1100\t1400\t1\t0.30\t10\t0.5\t0.5\t0.5\t0.2")
    fake <- tempfile("fake_metilene_")
    writeLines(c("#!/bin/sh", "echo 'Number of Tests: 100' >&2",
                 sprintf("printf '%s\\n'", gsub("\t", "\\\\t", bed))), fake)
    Sys.chmod(fake, "0755")
    d <- function(rate) data.table::data.table(
        chr = "chr1", pos = seq(100, by = 50, length.out = 30), strand = "+",
        site = 1:30, cov = 20L, rate = rate, mc = as.integer(round(20 * rate)))
    dat <- list(PseudoA = list(PseudoA_1 = d(0.2), PseudoA_2 = d(0.25)),
                PseudoB = list(PseudoB_1 = d(0.3), PseudoB_2 = d(0.35)))

    res <- suppressMessages(call_dmrs_metilene(
        dat, tempfile("dmr_"), fake, min_diff = 0, q_source = "BH",
        min_effect = 0.1, mtc = 2))

    expect_identical(res$start, 100L)
    sig <- attr(res, "significant")
    expect_identical(sig$start, c(100L, 600L))
    expect_identical(names(sig), names(res))
})


test_that("compute_metrics reports significance-only sensitivity when supplied", {
    m <- compute_metrics(c(TP = 1, FP = 0, FN = 3, n_regions = 4,
                           regions_detected = 1, regions_significant = 3))
    expect_equal(m$sensitivity, 0.25)
    expect_equal(m$regions_significant, 3)
    expect_equal(m$sensitivity_significant, 0.75)

    # Older summaries have no significant count.
    m0 <- compute_metrics(c(TP = 1, FP = 0, FN = 3, n_regions = 4,
                            regions_detected = 1))
    expect_true(is.na(m0$sensitivity_significant))
    expect_equal(m0$sensitivity, 0.25)
})


test_that("compiled tables carry the significance-only columns", {
    calls <- data.table::data.table(chr = "chr1", start = 1500L, end = 1800L,
                                    mean_diff = 0.2)
    r <- classify_dmrs(calls, .sd_truth())
    r$significant_regions <- c(1L, 3L, 4L)
    r$summary["regions_significant"] <- 3
    plain <- classify_dmrs(calls, .sd_truth())

    dt <- compile_results(list(k = list(SMFnorm = r, raw = plain)))
    expect_equal(dt[method == "SMFnorm", sensitivity_significant], 0.75)
    expect_identical(dt[method == "SMFnorm", regions_significant], 3L)
    expect_true(is.na(dt[method == "raw", sensitivity_significant]))

    by <- compile_density_results(list(k = list(SMFnorm = r, raw = plain)),
                                  .sd_truth())
    expect_equal(by[method == "SMFnorm" & density_stratum == "low",
                    c(regions_detected, regions_significant)], c(1, 1))
    expect_equal(by[method == "SMFnorm" & density_stratum == "high",
                    sensitivity_significant], 1)
    expect_true(all(is.na(by[method == "raw", sensitivity_significant])))
})


test_that("spikein_direction matches what inject_spikein_parametric plants", {
    set.seed(3)
    ref <- data.table::data.table(chr = "chr1", pos = seq(10L, by = 10L, length.out = 400))
    # Regions at low, high and two middling baseline rates, plus one without sites.
    p <- rep(c(0.1, 0.9, 0.45, 0.55), each = 100) + runif(400, -0.02, 0.02)
    regions <- data.table::data.table(
        region_id = 1:5, chr = c(rep("chr1", 4), "chr2"),
        start = c(10L, 1010L, 2010L, 3010L, 10L),
        end   = c(1000L, 2000L, 3000L, 4000L, 1000L),
        target_group = c("A", "B", "B", "A", "A"))

    d <- spikein_direction(ref, p, regions)
    expect_identical(d$direction, c("increase", "decrease", "increase", "decrease", NA))
    # Sign of the planted A - B difference.
    expect_identical(d$planted_sign, c(1L, 1L, -1L, -1L, NA))
    expect_equal(d$baseline_rate[1], mean(p[1:100]))

    inj <- suppressMessages(inject_spikein_parametric(
        ref, p, regions, effect_size = 0.1, target_group = regions$target_group))
    planted <- vapply(1:4, function(i) {
        idx <- ref$pos >= regions$start[i] & ref$pos <= regions$end[i]
        sign(mean(inj$p_A[idx] - inj$p_B[idx]))
    }, numeric(1))
    expect_equal(planted, d$planted_sign[1:4])
    truth_dir <- unique(inj$truth[, .(region_id, direction)])
    expect_identical(truth_dir$direction, d$direction[1:4])
})
