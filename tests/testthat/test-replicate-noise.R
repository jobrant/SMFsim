# Structured replicate noise and two-group spike-ins (2026-10-28).
#
# SMFsim used to draw every site's replicate noise independently and to put
# every spike-in effect in group B. Real M1 replicates share noise between
# sites on the same sequenced fragments (residual correlation ~0.34 within
# 50 bp), tissue replicates also vary by region, and in a real experiment an
# efficiency artifact does not always oppose the biology. These tests pin the
# new behaviour and check that the original path is unchanged when it is off.

.noise_sites <- function(n = 20000, seed = 3) {
    set.seed(seed)
    pos <- cumsum(sample(5:40, n, TRUE))
    data.table::data.table(chr = "chr1", pos = pos, strand = "+",
                           site = seq_len(n), cov = 40L,
                           mc = 0L, rate = 0)
}

.neighbour_cor <- function(reps, lo = 1, hi = 50) {
    # Correlation of replicate residuals between site pairs lo..hi bp apart,
    # pairing each site with sites up to 40 positions downstream.
    R <- do.call(cbind, lapply(reps, `[[`, "rate"))
    Z <- R - rowMeans(R)
    pos <- reps[[1]]$pos
    n <- length(pos)
    pairs <- do.call(rbind, lapply(1:40, function(lag) {
        i <- seq_len(n - lag)
        d <- pos[i + lag] - pos[i]
        keep <- d >= lo & d <= hi
        cbind(i[keep], i[keep] + lag)
    }))
    stopifnot(nrow(pairs) > 100)
    mean(vapply(seq_len(ncol(Z)), function(j) {
        cor(Z[pairs[, 1], j], Z[pairs[, 2], j])
    }, numeric(1)))
}

.noise <- function(rho = 0, bio_sd = 0) {
    list(fragment_bp = 50L, fragment_rho = rho, bio_region_bp = 1000L,
         bio_region_sd = bio_sd)
}


test_that(".noise_from_config is NULL when both components are off", {
    cfg <- parse_args(character(0))
    cfg$dispersion_s <- 26
    cfg$fragment_rho <- 0
    cfg$bio_region_sd <- 0
    expect_null(.noise_from_config(cfg))
})


test_that("with noise = NULL the draw is exactly the original one", {
    tmpl <- .noise_sites()
    p <- stats::runif(nrow(tmpl), 0.05, 0.6)
    a <- .draw_parametric_rep(tmpl, p, 0.8, 11, 26, noise = NULL)
    set.seed(11)
    expected <- .betabinom_counts(p * 0.8, tmpl$cov, 26)
    expect_identical(a$mc, expected)
})


test_that("shared-fragment noise correlates neighbours and keeps the mean", {
    tmpl <- .noise_sites()
    set.seed(5)
    p <- stats::runif(nrow(tmpl), 0.1, 0.5)
    draw <- function(rho) lapply(1:4, function(i) {
        .draw_parametric_rep(tmpl, p, 0.9, 100 + i, 26,
                             noise = if (rho > 0) .noise(rho) else NULL)
    })
    indep <- draw(0)
    shared <- draw(0.7)

    expect_lt(abs(.neighbour_cor(indep)), 0.03)
    expect_gt(.neighbour_cor(shared), .neighbour_cor(indep) + 0.05)
    # Same marginal: mean rate unchanged
    expect_equal(mean(sapply(shared, function(d) mean(d$rate))),
                 mean(sapply(indep, function(d) mean(d$rate))), tolerance = 0.01)
    # Nothing shared between sites farther apart than a fragment block
    expect_lt(abs(.neighbour_cor(shared, 200, 600)), 0.05)
})


test_that("region-level biological variation correlates sites ~1 kb apart", {
    tmpl <- .noise_sites(n = 30000)
    set.seed(6)
    p <- stats::runif(nrow(tmpl), 0.1, 0.5)
    reps <- lapply(1:4, function(i) {
        .draw_parametric_rep(tmpl, p, 0.9, 200 + i, 26,
                             noise = .noise(0, bio_sd = 0.5))
    })
    # Sites 100-400 bp apart usually share a 1 kb block.
    expect_gt(.neighbour_cor(reps, 100, 400), 0.1)
})


test_that("spike-in effects go to the group they are assigned", {
    tmpl <- .noise_sites(n = 2000)
    p <- rep(0.2, nrow(tmpl))
    regions <- data.table::data.table(
        region_id = 1:2, chr = "chr1",
        start = c(tmpl$pos[100], tmpl$pos[1000]),
        end   = c(tmpl$pos[150], tmpl$pos[1050]))
    inj <- inject_spikein_parametric(tmpl, p, regions, effect_size = 0.1,
                                     target_group = c("A", "B"))
    in1 <- 100:150
    in2 <- 1000:1050
    expect_equal(inj$p_A[in1], rep(0.3, length(in1)))
    expect_equal(inj$p_B[in1], rep(0.2, length(in1)))
    expect_equal(inj$p_A[in2], rep(0.2, length(in2)))
    expect_equal(inj$p_B[in2], rep(0.3, length(in2)))
    expect_identical(inj$p, inj$p_B)
    expect_setequal(unique(inj$truth$group), c("A", "B"))
})


test_that("default target_group keeps every effect in group B", {
    tmpl <- .noise_sites(n = 500)
    p <- rep(0.2, nrow(tmpl))
    regions <- data.table::data.table(region_id = 1, chr = "chr1",
                                      start = tmpl$pos[10], end = tmpl$pos[60])
    inj <- inject_spikein_parametric(tmpl, p, regions, effect_size = 0.1)
    expect_identical(inj$p_A, p)
    expect_equal(inj$p_B[10:60], rep(0.3, 51))
})


test_that("parse_args validates the new options", {
    expect_error(parse_args(c("--spikein_groups", "A")), "'both' or 'B'")
    expect_error(parse_args(c("--fragment_rho", "1")), "\\[0, 1\\)")
    expect_error(parse_args(c("--bio_region_sd", "-0.1")), ">= 0")
    cfg <- parse_args(c("--fragment_rho", "0.3", "--bio_region_sd", "0.25",
                        "--fragment_bp", "60", "--spikein_groups", "B"))
    expect_identical(cfg$fragment_rho, 0.3)
    expect_identical(cfg$bio_region_sd, 0.25)
    expect_identical(cfg$fragment_bp, 60L)
    expect_identical(cfg$spikein_groups, "B")
})


# ── Efficiency model ─────────────────────────────────────────────────────────
#
# "multiplicative" (default; fits most real replicate pairs) scales every site
# by the efficiency; "logit" applies a shift on the logit scale with the same
# genome-wide loss (fits the most extreme M1 pair, N3 vs N2, better).

test_that("multiplicative model is rate * efficiency", {
    p <- c(0, 0.1, 0.5, 0.9, 1)
    expect_equal(.apply_efficiency(p, 0.6, "multiplicative"), p * 0.6)
})

test_that("logit model matches the genome-wide mean loss of the efficiency", {
    set.seed(8)
    p <- c(stats::rbeta(50000, 1.2, 4), rep(0, 500), rep(1, 50))
    for (e in c(0.9, 0.7, 0.55)) {
        expect_equal(mean(.apply_efficiency(p, e, "logit")), e * mean(p),
                     tolerance = 1e-3, info = paste("efficiency", e))
    }
})

test_that("logit model is a constant shift on the logit scale", {
    p <- seq(0.02, 0.98, by = 0.02)
    q <- .apply_efficiency(p, 0.6, "logit")
    d <- stats::qlogis(q) - stats::qlogis(p)
    expect_lt(diff(range(d)), 1e-8)
    expect_lt(d[1], 0)
    # ...so, unlike the multiplicative model, the ratio is not constant
    expect_gt(diff(range(q / p)), 0.1)
})

test_that("logit model leaves rates of 0 and 1, and efficiency 1, unchanged", {
    p <- c(0, 0.3, 1)
    expect_equal(.apply_efficiency(p, 1, "logit"), p)
    out <- .apply_efficiency(c(p, stats::runif(100, 0.1, 0.5)), 0.7, "logit")
    expect_equal(out[c(1, 3)], c(0, 1))
})

test_that("default draw path is unchanged by the efficiency-model argument", {
    tmpl <- .noise_sites(n = 2000)
    p <- stats::runif(nrow(tmpl), 0.05, 0.6)
    a <- .draw_parametric_rep(tmpl, p, 0.8, 21, 26)
    set.seed(21)
    expect_identical(a$mc, .betabinom_counts(p * 0.8, tmpl$cov, 26))
})

test_that("parse_args defaults to the multiplicative model and validates it", {
    expect_identical(parse_args(character(0))$efficiency_model, "multiplicative")
    expect_error(parse_args(c("--efficiency_model", "linear")),
                 "'logit' or 'multiplicative'")
})
