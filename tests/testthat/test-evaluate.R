# Tests for DMR classification and metrics (03_evaluate.R).
#
# Every reported sensitivity, precision and false-positive count comes from
# this code, and it had no tests: labelling true positives as false, or
# computing sensitivity from call counts instead of regions, passed the suite.
# There were also two classifiers that scored the same calls differently,
# chosen by which packages happened to be installed.

skip_if_not_installed("GenomicRanges")

.truth <- function() data.table::data.table(
    region_id = 1:4, chr = c("chr1", "chr1", "chr1", "chr2"),
    start = c(1000L, 5000L, 9000L, 1000L),
    end   = c(2000L, 6000L, 10000L, 2000L),
    density_stratum = c("low", "low", "high", "high"))

.calls <- function(chr, start, end) data.table::data.table(
    chr = chr, start = as.integer(start), end = as.integer(end),
    mean_diff = 0.2)


test_that("calls overlapping a truth region are TP, others FP, missed regions FN", {
    calls <- .calls(c("chr1", "chr1", "chr1"),
                    c(1500, 3000, 9500), c(1800, 3500, 9900))
    r <- classify_dmrs(calls, .truth())

    expect_identical(r$classification$label, c("TP", "FP", "TP"))
    expect_identical(r$classification$matched_region, c(1L, NA, 3L))
    expect_equal(unname(r$summary[c("TP", "FP", "FN")]), c(2, 1, 2))
    expect_identical(sort(r$fn_regions$region_id), c(2L, 4L))
    expect_identical(r$detected_regions, c(1L, 3L))
})


test_that("the same coordinates on another chromosome do not match", {
    r <- classify_dmrs(.calls("chr3", 1500, 1800), .truth())

    expect_identical(r$classification$label, "FP")
    expect_length(r$detected_regions, 0L)
})


test_that("a one-base overlap counts and an adjacent call does not", {
    touching <- classify_dmrs(.calls("chr1", 2000, 2500), .truth())   # shares base 2000
    adjacent <- classify_dmrs(.calls("chr1", 2001, 2500), .truth())

    expect_identical(touching$classification$label, "TP")
    expect_identical(adjacent$classification$label, "FP")
})


test_that("min_overlap_bp is the number of shared bases required", {
    call <- .calls("chr1", 1991, 2500)        # shares 1991..2000 = 10 bases
    expect_identical(classify_dmrs(call, .truth(), min_overlap_bp = 10)$classification$label, "TP")
    expect_identical(classify_dmrs(call, .truth(), min_overlap_bp = 11)$classification$label, "FP")
})


test_that("several calls in one region are several TP calls but one detected region", {
    calls <- .calls(rep("chr1", 3), c(1000, 1300, 1700), c(1200, 1500, 1900))
    r <- classify_dmrs(calls, .truth())

    expect_equal(unname(r$summary["TP"]), 3)
    expect_equal(unname(r$summary["regions_detected"]), 1)
    expect_equal(unname(r$summary["FN"]), 3)

    m <- compute_metrics(r$summary)
    expect_equal(m$sensitivity, 1 / 4)        # region-level, not 3 / 4
    expect_equal(m$precision, 1)              # call-level
})


test_that("one call spanning two regions detects both", {
    r <- classify_dmrs(.calls("chr1", 1500, 5500), .truth())

    expect_equal(unname(r$summary["TP"]), 1)
    expect_identical(r$detected_regions, c(1L, 2L))
    expect_equal(unname(r$summary["regions_detected"]), 2)
})


test_that("with no truth regions every call is a false positive", {
    calls <- .calls(c("chr1", "chr2"), c(10, 10), c(50, 50))
    r <- classify_dmrs(calls, truth_regions = NULL)

    expect_identical(r$classification$label, c("FP", "FP"))
    expect_equal(unname(r$summary[c("TP", "FP", "FN")]), c(0, 2, 0))
    expect_length(r$detected_regions, 0L)
})


test_that("with no calls every truth region is a false negative", {
    empty <- .calls(character(), integer(), integer())
    r <- classify_dmrs(empty, .truth())

    expect_equal(unname(r$summary[c("TP", "FP", "FN")]), c(0, 0, 4))
    expect_equal(unname(r$summary["regions_detected"]), 0)
    expect_length(r$detected_regions, 0L)
    expect_equal(compute_metrics(r$summary)$sensitivity, 0)
})


test_that("compute_metrics gives the stated formulas", {
    m <- compute_metrics(c(TP = 6, FP = 2, FN = 5, n_regions = 10,
                           regions_detected = 5))

    expect_equal(m$sensitivity, 0.5)                       # 5 of 10 regions
    expect_equal(m$precision, 0.75)                        # 6 of 8 calls
    expect_equal(m$FDR, 0.25)
    expect_equal(m$F1, round(2 * 0.5 * 0.75 / (0.5 + 0.75), 4))
    expect_equal(m$total_called, 8)
})


test_that("compute_metrics returns NA precision when nothing is called", {
    m <- compute_metrics(c(TP = 0, FP = 0, FN = 4, n_regions = 4,
                           regions_detected = 0))

    expect_equal(m$sensitivity, 0)
    expect_true(is.na(m$precision))
    expect_true(is.na(m$FDR))
    expect_true(is.na(m$F1))
})


test_that("compile_results reports one row per scenario and method", {
    r1 <- classify_dmrs(.calls(c("chr1", "chr1"), c(1500, 3000), c(1800, 3500)), .truth())
    r2 <- classify_dmrs(.calls("chr1", 3000, 3500), .truth())
    dt <- compile_results(list(sc_effect0.2 = list(raw = r1, SMFnorm = r2)))

    expect_equal(nrow(dt), 2L)
    expect_equal(dt[method == "raw", c(TP, FP, regions_detected)], c(1, 1, 1))
    expect_equal(dt[method == "raw", sensitivity], 0.25)
    expect_equal(dt[method == "SMFnorm", c(TP, FP, regions_detected)], c(0, 1, 0))
    expect_equal(dt[method == "SMFnorm", sensitivity], 0)
})


test_that("compile_null_results counts false positives", {
    r <- classify_dmrs(.calls(c("chr1", "chr2"), c(10, 10), c(50, 50)), NULL)
    dt <- compile_null_results(list(mild = list(raw = r)))

    expect_equal(dt$FP, 2)
    expect_identical(dt$scenario, "mild")
    expect_identical(dt$method, "raw")
})


test_that("compile_density_results gives region-level sensitivity per stratum", {
    # Detects region 1 (low) and regions 3 and 4 (high).
    calls <- .calls(c("chr1", "chr1", "chr2"), c(1500, 9500, 1500), c(1800, 9900, 1800))
    r <- classify_dmrs(calls, .truth())
    dt <- compile_density_results(list(mild_effect0.2 = list(raw = r)), .truth())

    expect_equal(nrow(dt), 2L)
    expect_equal(dt[density_stratum == "low", c(n_regions, regions_detected)], c(2, 1))
    expect_equal(dt[density_stratum == "low", sensitivity], 0.5)
    expect_equal(dt[density_stratum == "high", sensitivity], 1)
})


test_that("compile_density_results is empty when regions carry no stratum", {
    truth <- .truth()[, density_stratum := NULL]
    r <- classify_dmrs(.calls("chr1", 1500, 1800), truth)

    expect_equal(nrow(compile_density_results(list(k = list(raw = r)), truth)), 0L)
})
