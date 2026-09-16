# Tests for the metilene multiple-testing-correction option.
#
# Background: the metilene command never passed -c, so it silently used
# metilene's DEFAULT correction (Bonferroni) while the roxygen docs claimed
# Benjamini-Hochberg. The correction is now explicit and validated, so the
# choice is recorded rather than implied.
#
# `mtc` is validated before any file I/O, so these run without a metilene
# binary present.

test_that("mtc rejects values other than 1 or 2", {
    expect_error(
        call_dmrs_metilene(split_data = NULL, out_dir = tempdir(), mtc = 3),
        "must be 1 \\(Bonferroni\\) or 2")

    expect_error(
        call_dmrs_metilene(split_data = NULL, out_dir = tempdir(), mtc = 0),
        "must be 1 \\(Bonferroni\\) or 2")
})

test_that("mtc rejects non-scalar and empty input", {
    expect_error(
        call_dmrs_metilene(split_data = NULL, out_dir = tempdir(), mtc = c(1, 2)),
        "must be 1 \\(Bonferroni\\) or 2")

    expect_error(
        call_dmrs_metilene(split_data = NULL, out_dir = tempdir(), mtc = NULL),
        "must be 1 \\(Bonferroni\\) or 2")
})

test_that("mtc validation happens before metilene_path is even checked", {
    # Ordering matters: a bad mtc must not be masked by an unrelated error, and
    # must not create input files on disk first.
    expect_error(
        call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                           metilene_path = "", mtc = 99),
        "must be 1 \\(Bonferroni\\) or 2")
})

test_that("valid mtc values pass validation and reach the path check", {
    # With a valid mtc, the next guard (empty metilene_path) is what fires.
    for (v in c(1, 2)) {
        expect_error(
            call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                               metilene_path = "", mtc = v),
            "metilene_path is empty")
    }
})

test_that("parse_args exposes metilene_mtc as a valid correction code", {
    cfg <- parse_args()
    expect_true(cfg$metilene_mtc %in% c(1, 2))
})
