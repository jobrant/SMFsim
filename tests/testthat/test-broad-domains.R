# Broad, genuine group differences (2026-10-28).
#
# The broad_* scenarios raise group B's true rate over 1 Mb domains covering a
# fraction of the genome. They measure what SMFnorm's between-group step costs:
# its per-bin offset absorbs roughly fraction x effect and pushes the opposite
# shift onto everything else. Calls inside a domain are genuine, so they must
# not count as false positives; they are scored as domain retention.

.broad_sites <- function() {
    data.table::data.table(
        chr = rep(c("chr1", "chr2"), each = 4000),
        pos = c(seq(1, by = 2500, length.out = 4000),
                seq(1, by = 2500, length.out = 4000)))
}


test_that("domains cover the requested fraction and avoid excluded regions", {
    s <- .broad_sites()                       # 10 Mb per chromosome -> 20 windows
    ex <- data.table::data.table(chr = "chr1", start = 2.5e6, end = 3.5e6)
    d <- .select_broad_domains(s, 0.3, 1e6, exclude = ex, seed = 4)
    expect_equal(nrow(d), round(0.3 * 18))    # 20 windows minus 2 excluded
    expect_false(any(d$chr == "chr1" & d$win %in% c(2L, 3L)))
    expect_true(all(d$end - d$start + 1 == 1e6))
})


test_that("the broad shift raises the true rate only inside domains", {
    s <- .broad_sites()
    d <- .select_broad_domains(s, 0.2, 1e6, seed = 2)
    p <- rep(0.2, nrow(s))
    p2 <- .apply_broad_shift(s, p, d, 0.15)
    inside <- .in_domains(s, d)
    expect_true(any(inside) && !all(inside))
    expect_equal(p2[inside], rep(0.35, sum(inside)))
    expect_equal(p2[!inside], p[!inside])
    expect_equal(.apply_broad_shift(s, rep(0.95, nrow(s)), d, 0.15)[inside],
                 rep(1, sum(inside)))         # clamped
})


test_that("no domains means no site is inside one", {
    s <- .broad_sites()
    expect_false(any(.in_domains(s, NULL)))
    expect_false(any(.in_domains(s, s[0, .(chr, start = pos, end = pos, win = pos)])))
})


test_that("calls inside or near domains are not false positives", {
    dom <- data.table::data.table(chr = "chr1", win = 1L, start = 1e6, end = 2e6 - 1)
    win <- data.table::data.table(chr = "chr1", start = c(1e6, 1.5e6),
                                  end = c(1e6 + 1999, 1.5e6 + 1999))
    data.table::setkey(win, chr, start, end)
    calls <- data.table::data.table(
        chr = "chr1",
        start = c(1.5e6, 2e6 + 500, 5e6),     # inside / within margin / outside
        end   = c(1.5e6 + 300, 2e6 + 900, 5e6 + 300),
        mean_diff = c(-0.2, -0.15, 0.12))
    r <- .split_domain_calls(calls, dom, win, margin_bp = 2000)
    expect_equal(r$n_domain_calls, 2L)
    expect_equal(nrow(r$outside), 1L)
    expect_equal(r$outside$start, 5e6)
    expect_equal(r$retained, 0.5)             # 1 of 2 domain windows hit
})


test_that("no calls gives zero retention and no false positives", {
    dom <- data.table::data.table(chr = "chr1", win = 1L, start = 1e6, end = 2e6 - 1)
    win <- data.table::data.table(chr = "chr1", start = 1e6, end = 1e6 + 1999)
    data.table::setkey(win, chr, start, end)
    empty <- data.table::data.table(chr = character(), start = numeric(),
                                    end = numeric(), mean_diff = numeric())
    r <- .split_domain_calls(empty, dom, win)
    expect_equal(r$n_domain_calls, 0L)
    expect_equal(r$retained, 0)
    expect_equal(nrow(r$outside), 0L)
})


test_that("broad scenarios are defined with matched or aligned efficiencies", {
    sc <- get_efficiency_scenarios()
    for (nm in c("broad_10", "broad_30", "broad_50", "broad_30_aligned")) {
        expect_false(is.null(sc[[nm]]$broad), info = nm)
        expect_equal(sc[[nm]]$broad$effect, 0.15, info = nm)
    }
    for (nm in c("broad_10", "broad_30", "broad_50")) {
        expect_true(isTRUE(sc[[nm]]$matched_means), info = nm)
    }
    expect_equal(c(sc$broad_10$broad$fraction, sc$broad_30$broad$fraction,
                   sc$broad_50$broad$fraction), c(0.1, 0.3, 0.5))
})


test_that("create_pseudo_groups applies a broad difference to group B only", {
    s <- .broad_sites()[1:3000]
    reps <- stats::setNames(lapply(1:4, function(i) {
        data.table::data.table(chr = s$chr, pos = s$pos, strand = "+",
                               site = seq_len(nrow(s)), mc = 6L, cov = 30L,
                               rate = 0.2)
    }), paste0("R", 1:4))
    d <- .select_broad_domains(reps[[1]], 0.3, 1e6, seed = 3)
    inside <- .in_domains(reps[[1]], d)
    pg <- suppressMessages(create_pseudo_groups(
        reps, rep(0.9, 4), rep(0.9, 4), mode = "parametric", seed = 1,
        dispersion_s = 128, broad = list(domains = d, effect = 0.15)))
    gm <- function(g, idx) mean(sapply(pg[[g]], function(dt) mean(dt$rate[idx])))
    expect_equal(gm("PseudoB", inside) - gm("PseudoA", inside), 0.15 * 0.9,
                 tolerance = 0.02)
    expect_equal(gm("PseudoB", !inside) - gm("PseudoA", !inside), 0,
                 tolerance = 0.01)
    expect_error(suppressMessages(create_pseudo_groups(
        reps, rep(0.9, 4), rep(0.9, 4), mode = "clone",
        broad = list(domains = d, effect = 0.15))), "parametric")
})
