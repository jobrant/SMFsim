# Tests for reading metilene's total test count.
#
# Background: metilene emits far fewer segments than it tests -- on one real
# comparison, 985949 tests against 32162 emitted rows at -d 0 (and 568 at
# -d 0.1). Correcting over the emitted rows instead of the reported test count
# under-corrects by ~30x and ~1700x respectively. The count is printed on
# stderr, which the pipeline used to discard via 2>/dev/null.

.write_log <- function(lines) {
    f <- tempfile(fileext = ".txt")
    writeLines(lines, f)
    f
}

test_that("the reported test count is read from a metilene log", {
    f <- .write_log(c("metilene starting",
                      "Number of Tests: 985949",
                      "done"))
    expect_equal(.metilene_n_tests(f), 985949L)
})

test_that("parsing tolerates spacing and case variation", {
    expect_equal(.metilene_n_tests(.write_log("number of tests:12345")), 12345L)
    expect_equal(.metilene_n_tests(.write_log("Number  of  Tests  :  678")), 678L)
})

test_that("the last count wins when several are reported", {
    # metilene may report per-chromosome counts before a total.
    f <- .write_log(c("Number of Tests: 100",
                      "Number of Tests: 250"))
    expect_equal(.metilene_n_tests(f), 250L)
})

test_that("a missing, empty, or uninformative log yields NA", {
    expect_true(is.na(.metilene_n_tests(file.path(tempdir(), "nope.txt"))))
    expect_true(is.na(.metilene_n_tests(.write_log(character(0)))))
    expect_true(is.na(.metilene_n_tests(.write_log("no count here"))))
    expect_true(is.na(.metilene_n_tests(NULL)))
})

test_that("a zero or malformed count is rejected rather than trusted", {
    expect_true(is.na(.metilene_n_tests(.write_log("Number of Tests: 0"))))
    expect_true(is.na(.metilene_n_tests(.write_log("Number of Tests: abc"))))
})

test_that("p.adjust with the full test count is far more conservative", {
    # The behaviour the fix depends on: same p-values, different denominator.
    set.seed(1)
    p <- stats::runif(32162, 0, 0.02)

    q_emitted <- stats::p.adjust(p, "BH")
    q_full    <- stats::p.adjust(p, "BH", n = 985949)

    expect_true(all(q_full >= q_emitted))
    expect_gt(sum(q_emitted < 0.05), sum(q_full < 0.05))
})
