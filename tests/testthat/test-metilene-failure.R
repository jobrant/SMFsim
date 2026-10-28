# A metilene run that does not complete must stop the analysis.
#
# call_dmrs_metilene() used to warn on a non-zero exit and return an empty
# table, and both simulation loops caught errors and substituted an empty
# table too. A crash therefore scored as "0 DMRs": a perfect null and zero
# sensitivity in a spike-in. Under Rscript the warning was deferred to the end
# of the log. These tests use a stand-in metilene script to pin the behaviour.

.fake_metilene <- function(body) {
    f <- tempfile("fake_metilene_")
    writeLines(c("#!/bin/sh", body), f)
    Sys.chmod(f, "0755")
    f
}

.tiny_split <- function(n = 30) {
    d <- function(rate) data.table::data.table(
        chr = "chr1", pos = seq(100, by = 50, length.out = n), strand = "+",
        site = seq_len(n), cov = 20L, rate = rate,
        mc = as.integer(round(20 * rate)))
    list(PseudoA = list(PseudoA_1 = d(0.2), PseudoA_2 = d(0.25)),
         PseudoB = list(PseudoB_1 = d(0.3), PseudoB_2 = d(0.35)))
}

.call <- function(metilene) {
    call_dmrs_metilene(.tiny_split(), tempfile("dmr_"), metilene_path = metilene,
                       min_cpg = 5, min_diff = 0, q_source = "BH",
                       min_effect = 0.1, mtc = 2)
}


test_that("a non-zero metilene exit is an error, not 0 DMRs", {
    crash <- .fake_metilene(c("echo 'segmenting chr1' >&2", "exit 1"))
    expect_error(.call(crash), "metilene did not complete \\(exit status")
})


test_that("a metilene log without 'Number of Tests' is an error", {
    # Exits 0 but never reaches the end of its run.
    truncated <- .fake_metilene(c("echo 'segmenting chr1' >&2", "exit 0"))
    expect_error(.call(truncated), "no 'Number of Tests' line")
})


test_that("the error reports the command, the log and its last lines", {
    crash <- .fake_metilene(c("echo 'out of memory somewhere' >&2", "exit 137"))
    err <- tryCatch(.call(crash), error = conditionMessage)
    expect_match(err, "exit status 137")
    expect_match(err, "out of memory somewhere")
    expect_match(err, "metilene_log.txt")
})


test_that("a completed run with no DMRs still returns an empty table", {
    ok <- .fake_metilene(c("echo 'Number of Tests: 12' >&2", "exit 0"))
    res <- .call(ok)
    expect_s3_class(res, "data.table")
    expect_equal(nrow(res), 0L)
})
