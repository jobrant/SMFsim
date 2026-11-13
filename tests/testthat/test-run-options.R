# Options a run is configured with from outside the scripts: where the input
# lives, and SMFnorm's bin count.

test_that("parse_args takes the input location from the environment", {
    # Set directly rather than via withr, which is not a declared dependency.
    vars <- c("SMFSIM_DATA_DIR", "SMFSIM_SAMPLE_SHEET")
    old <- Sys.getenv(vars, unset = NA)
    on.exit({
        Sys.unsetenv(vars[is.na(old)])
        if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    }, add = TRUE)

    Sys.unsetenv(vars)
    cfg <- parse_args(character(0))
    expect_identical(cfg$data_dir, "data/allc_methylseq4.2_clip5")
    expect_identical(cfg$sample_sheet, "data/sample_sheet.csv")

    Sys.setenv(SMFSIM_DATA_DIR = "d/allc", SMFSIM_SAMPLE_SHEET = "d/sheet.csv")
    cfg <- parse_args(character(0))
    expect_identical(cfg$data_dir, "d/allc")
    expect_identical(cfg$sample_sheet, "d/sheet.csv")
    # The command line still wins.
    expect_identical(parse_args(c("--data_dir", "x"))$data_dir, "x")
})


test_that("parse_args reads and validates max_quantiles", {
    expect_identical(parse_args(character(0))$max_quantiles, 50L)
    expect_identical(parse_args(c("--max_quantiles", "200"))$max_quantiles, 200L)
    expect_error(parse_args(c("--max_quantiles", "0")), ">= 1")
    expect_error(parse_args(c("--max_quantiles", "many")), "whole number")
})


# A between-group efficiency gap on rates spread over (0, 1): the number of
# bins changes the normalized values, so the option must reach SMFnorm.
test_that("max_quantiles reaches SMFnorm through run_all_methods", {
    set.seed(21)
    n <- 20000
    p <- stats::rbeta(n, 1, 3)
    mk <- function(eff) {
        cov <- rep(40L, n)
        mc <- stats::rbinom(n, cov, pmin(p * eff, 1))
        data.table::data.table(chr = "chr1", pos = seq_len(n) * 10L, strand = "+",
                               site = "GCH", mc = mc, cov = cov, rate = mc / cov)
    }
    pg <- list(PseudoA = list(PseudoA_1 = mk(1.0), PseudoA_2 = mk(0.95)),
               PseudoB = list(PseudoB_1 = mk(0.7), PseudoB_2 = mk(0.65)),
               params = list())
    run <- function(mq) suppressMessages(run_all_methods(
        pg, methods = "SMFnorm",
        SMFnorm_params = list(min_coverage = 10, max_quantiles = mq)))$SMFnorm$data

    a <- run(5); b <- run(5); c <- run(10)
    rate1 <- function(x) x[[1]][[1]]$rate
    expect_identical(rate1(a), rate1(b))
    expect_false(isTRUE(all.equal(rate1(a), rate1(c))))
})


test_that("parse_args accepts a chromosome pattern and defaults to none", {
    cfg <- parse_args(character(0))
    expect_true("chr_pattern" %in% names(cfg))
    expect_null(cfg$chr_pattern)
    expect_identical(parse_args(c("--chr_pattern", "^(chr)?1$"))$chr_pattern,
                     "^(chr)?1$")
})
