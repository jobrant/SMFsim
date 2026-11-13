# Spike-in placement across site densities, and repeated seeds.
#
# Spike-ins used to go in the densest regions genome-wide (median 156 sites
# against the 10 metilene needs), where every method detected nearly all of
# them, and every result came from seed 42. Regions are now sampled evenly
# from density strata, sensitivity is also reported per stratum, and the
# simulation can be repeated over seeds.

# Sites on two chromosomes in 2 kb windows holding 12 to 120 sites each.
.density_sites <- function(n_windows = 180, seed = 5) {
    set.seed(seed)
    n_per <- sample(12:120, n_windows, replace = TRUE)
    data.table::rbindlist(lapply(seq_len(n_windows), function(w) {
        chr <- if (w <= n_windows / 2) "chr1" else "chr2"
        base <- (w - 1L) * 2000L
        data.table::data.table(chr = chr,
                               pos = base + sort(sample(1:1999, n_per[w])))
    }))
}

.sel <- function(sites, ...) suppressMessages(
    select_stratified_spikein_regions(sites, ...))


test_that("regions are split evenly across the density strata", {
    r <- .sel(.density_sites(), n_regions = 60, seed = 1)

    expect_equal(nrow(r), 60L)
    expect_identical(as.vector(table(r$density_stratum)[c("low", "medium", "high")]),
                     c(20L, 20L, 20L))
    expect_identical(r$region_id, 1:60)
})


test_that("strata are ordered by site count and do not overlap in density", {
    r <- .sel(.density_sites(), n_regions = 60, seed = 1)
    rng <- r[, .(lo = min(n_sites), hi = max(n_sites)), by = density_stratum]
    data.table::setkey(rng, density_stratum)

    expect_lte(rng["low", hi], rng["medium", lo])
    expect_lte(rng["medium", hi], rng["high", lo])
    # The point of the design: low-density regions are far sparser than the
    # densest ones.
    expect_lt(median(r[density_stratum == "low", n_sites]),
              0.5 * median(r[density_stratum == "high", n_sites]))
})


test_that("region coordinates and site counts match the data", {
    sites <- .density_sites()
    r <- .sel(sites, n_regions = 30, seed = 2)

    for (i in seq_len(nrow(r))) {
        inside <- sites[chr == r$chr[i] & pos >= r$start[i] & pos <= r$end[i]]
        expect_equal(nrow(inside), r$n_sites[i])
        expect_lt(r$end[i] - r$start[i], 2000)
    }
})


test_that("regions never overlap", {
    r <- .sel(.density_sites(), n_regions = 90, seed = 3)
    data.table::setorder(r, chr, start)
    gap <- r[, .(d = start[-1] - end[-.N]), by = chr]

    expect_true(all(gap$d > 0))
})


test_that("eligibility limits are respected", {
    r <- .sel(.density_sites(), n_regions = 30, min_sites = 40,
              max_median_spacing = 30, seed = 4)

    expect_true(all(r$n_sites >= 40))
    expect_true(all(r$median_spacing <= 30))
})


test_that("the same seed gives the same regions and another seed different ones", {
    sites <- .density_sites()
    a <- .sel(sites, n_regions = 30, seed = 11)
    b <- .sel(sites, n_regions = 30, seed = 11)
    c <- .sel(sites, n_regions = 30, seed = 12)

    expect_identical(a, b)
    expect_false(identical(a[, .(chr, start)], c[, .(chr, start)]))
})


test_that("asking for more regions than there are eligible windows is an error", {
    expect_error(.sel(.density_sites(n_windows = 20), n_regions = 50),
                 "eligible windows")
})


test_that("other numbers of strata are labelled D1..Dk", {
    r <- .sel(.density_sites(), n_regions = 20, n_strata = 4, seed = 1)

    expect_setequal(unique(r$density_stratum), paste0("D", 1:4))
    expect_true(all(table(r$density_stratum) == 5L))
})


test_that("parse_args defaults to stratified placement and validates the options", {
    cfg <- parse_args(character(0))
    expect_identical(cfg$spikein_density, "stratified")
    expect_identical(cfg$spikein_density_strata, 3L)

    expect_identical(parse_args(c("--spikein_density", "densest"))$spikein_density,
                     "densest")
    expect_identical(parse_args(c("--spikein_density_strata", "5"))$spikein_density_strata,
                     5L)
    expect_error(parse_args(c("--spikein_density", "random")),
                 "'stratified' or 'densest'")
    expect_error(parse_args(c("--spikein_density_strata", "0")), ">= 1")
})


# --- Across seeds ------------------------------------------------------------

test_that("summarize_across_seeds gives mean, sd, min, max and the seed count", {
    dt <- data.table::data.table(
        seed = c(1, 2, 3, 1, 2),
        scenario = c("a", "a", "a", "b", "b"),
        method = "raw",
        FP = c(10, 20, 30, 4, NA))
    s <- summarize_across_seeds(dt, c("scenario", "method"), "FP")
    data.table::setkey(s, scenario)

    expect_equal(s["a", n_seeds], 3L)
    expect_equal(s["a", FP], 20)
    expect_equal(s["a", FP_sd], 10)
    expect_equal(s["a", c(FP_min, FP_max)], c(10, 30))
    # One usable value: a mean but no spread.
    expect_equal(s["b", n_seeds], 2L)
    expect_equal(s["b", FP], 4)
    expect_true(is.na(s["b", FP_sd]))
})


test_that("combine_seed_dirs stacks seed directories and reports missing ones", {
    base <- tempfile("seeds_")
    for (s in c(7, 8)) {
        d <- file.path(base, paste0("seed_", s), "blockA")
        dir.create(d, recursive = TRUE)
        data.table::fwrite(data.table::data.table(scenario = "mild", method = "raw",
                                                  FP = s), file.path(d, "null_simulation_results.csv"))
    }
    data.table::fwrite(
        data.table::data.table(efficiency_scenario = "mild", effect_size = 0.2,
                               method = "raw", sensitivity = 0.8),
        file.path(base, "seed_7", "blockA", "spikein_simulation_results.csv"))

    expect_warning(res <- suppressMessages(combine_seed_dirs(base, block = "blockA")),
                   "missing for seed\\(s\\) 8")

    expect_identical(res$by_seed$null$seed, c(7L, 8L))
    expect_equal(res$across_seeds$null$FP, 7.5)
    expect_equal(res$across_seeds$null$n_seeds, 2L)
    expect_equal(res$across_seeds$spikein$n_seeds, 1L)
    out <- file.path(base, "all_seeds", "blockA")
    expect_true(file.exists(file.path(out, "null_by_seed.csv")))
    expect_true(file.exists(file.path(out, "null_across_seeds.csv")))
})


test_that("combine_seed_dirs errors when there are no seed directories", {
    base <- tempfile("seeds_")
    dir.create(base)
    expect_error(combine_seed_dirs(base), "no seed_<n> directories")
})


test_that("run_simulation_seeds rejects repeated or missing seeds", {
    cfg <- parse_args(character(0))
    expect_error(run_simulation_seeds(list(), cfg, c(1, 1)), "distinct whole numbers")
    expect_error(run_simulation_seeds(list(), cfg, integer(0)), "distinct whole numbers")
})


# End to end with a stand-in metilene that reports one significant segment
# covering the whole chromosome: every spike-in region is detected, and the
# null has exactly one false positive.
test_that("the simulation runs over seeds and reports sensitivity by density", {
    skip_if_not_installed("GenomicRanges")
    skip_on_os("windows")

    sites <- .density_sites(n_windows = 120, seed = 9)[chr == "chr1"]
    set.seed(10)
    p <- stats::rbeta(nrow(sites), 1.5, 4)
    reps <- stats::setNames(lapply(1:4, function(i) {
        cov <- rep(30L, nrow(sites))
        mc <- stats::rbinom(nrow(sites), cov, p)
        data.table::data.table(chr = sites$chr, pos = sites$pos, strand = "+",
                               site = "GCH", mc = mc, cov = cov, rate = mc / cov)
    }), paste0("R", 1:4))

    fake <- tempfile("fake_metilene_")
    writeLines(c("#!/bin/sh",
                 "printf 'chr1\\t0\\t999999999\\t1e-20\\t0.25\\t500\\t1e-20\\t1e-20\\t0.4\\t0.15\\n'",
                 "echo 'Number of Tests: 50' >&2"), fake)
    Sys.chmod(fake, "0755")

    cfg <- parse_args(character(0))
    cfg$metilene_path <- fake
    cfg$output_dir <- tempfile("seedrun_")
    cfg$scenarios <- "mild"
    cfg$methods <- "raw"
    cfg$effect_sizes <- 0.2
    cfg$n_spikein_regions <- 30

    res <- suppressMessages(run_simulation_seeds(reps, cfg, seeds = c(1, 2)))

    # Each seed has its own directory, config and regions.
    for (s in 1:2) {
        d <- file.path(cfg$output_dir, paste0("seed_", s))
        expect_true(file.exists(file.path(d, "config.rds")))
        expect_identical(readRDS(file.path(d, "config.rds"))$seed, s)
        expect_true(file.exists(file.path(d, "spikein_by_density.csv")))
    }
    reg1 <- data.table::fread(file.path(cfg$output_dir, "seed_1", "spikein_regions.csv"))
    reg2 <- data.table::fread(file.path(cfg$output_dir, "seed_2", "spikein_regions.csv"))
    expect_equal(as.vector(table(reg1$density_stratum)), c(10L, 10L, 10L))
    expect_false(identical(reg1$start, reg2$start))

    # Stacked and summarised across seeds.
    expect_identical(sort(unique(res$by_seed$spikein$seed)), 1:2)
    expect_equal(res$across_seeds$null$n_seeds, 2L)
    expect_equal(res$across_seeds$null$FP, 1)
    expect_equal(res$across_seeds$spikein$sensitivity, 1)
    dens <- res$across_seeds$spikein_by_density
    expect_setequal(dens$density_stratum, c("low", "medium", "high"))
    expect_true(all(dens$sensitivity == 1))
    expect_true(all(dens$n_seeds == 2L))
    expect_true(file.exists(file.path(cfg$output_dir, "spikein_by_density_across_seeds.csv")))
})


test_that("the densest placement is still available and labelled as such", {
    skip_if_not_installed("GenomicRanges")
    skip_on_os("windows")

    sites <- .density_sites(n_windows = 120, seed = 9)[chr == "chr1"]
    set.seed(10)
    reps <- stats::setNames(lapply(1:4, function(i) {
        mc <- stats::rbinom(nrow(sites), 30L, 0.3)
        data.table::data.table(chr = sites$chr, pos = sites$pos, strand = "+",
                               site = "GCH", mc = mc, cov = 30L, rate = mc / 30)
    }), paste0("R", 1:4))
    fake <- tempfile("fake_metilene_")
    writeLines(c("#!/bin/sh", "echo 'Number of Tests: 50' >&2"), fake)
    Sys.chmod(fake, "0755")

    cfg <- parse_args(c("--spikein_density", "densest"))
    cfg$metilene_path <- fake
    cfg$output_dir <- tempfile("dense_")
    cfg$scenarios <- "mild"
    cfg$methods <- "raw"
    cfg$effect_sizes <- 0.2
    cfg$n_spikein_regions <- 10

    out <- suppressMessages(run_spikein_simulation(reps, cfg))
    reg <- data.table::fread(file.path(cfg$output_dir, "spikein_regions.csv"))

    expect_identical(unique(reg$density_stratum), "densest")
    expect_equal(out$sensitivity, 0)             # the stand-in called nothing
    expect_identical(unique(attr(out, "by_density")$density_stratum), "densest")
})
