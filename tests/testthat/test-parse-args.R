# Tests for command-line config overrides.
#
# Background: every key not explicitly typed used to fall through to a raw
# character assignment, so `--within_alpha 0.5` produced the STRING "0.5". That
# fails as "non-numeric argument to binary operator" deep inside SMFnorm --
# potentially hours into a multi-day job. Unknown keys were silently ignored, so
# a typo ran the whole job at defaults with no indication anything was wrong.

test_that("numeric options are parsed as numbers and are usable in arithmetic", {
    cfg <- parse_args(c("--within_alpha", "0.5", "--between_alpha", "0.7"))

    expect_type(cfg$within_alpha, "double")
    expect_equal(cfg$within_alpha, 0.5)
    expect_equal(cfg$between_alpha, 0.7)

    # The exact operation that failed when these came through as characters.
    expect_equal(cfg$within_alpha * 2, 1.0)
})

test_that("integer-typed options are parsed as whole numbers", {
    cfg <- parse_args(c("--min_coverage", "10", "--metilene_min_cpg", "8",
                        "--seed", "99", "--metilene_mtc", "2"))
    expect_equal(cfg$min_coverage, 10L)
    expect_equal(cfg$metilene_min_cpg, 8L)
    expect_equal(cfg$seed, 99L)
    expect_equal(cfg$metilene_mtc, 2L)
})

test_that("flags accept the usual spellings and stay logical", {
    for (v in c("true", "TRUE", "T", "yes", "1")) {
        expect_true(parse_args(c("--rate_between_groups", v))$rate_between_groups,
                    info = v)
    }
    for (v in c("false", "FALSE", "F", "no", "0")) {
        expect_false(parse_args(c("--rate_between_groups", v))$rate_between_groups,
                     info = v)
    }
})

test_that("comma-separated options split correctly", {
    cfg <- parse_args(c("--methods", "raw,SMFnorm",
                        "--effect_sizes", "0.1,0.2,0.3"))
    expect_equal(cfg$methods, c("raw", "SMFnorm"))
    expect_equal(cfg$effect_sizes, c(0.1, 0.2, 0.3))
    expect_type(cfg$effect_sizes, "double")
})

test_that("a non-numeric value for a numeric option is rejected", {
    expect_error(parse_args(c("--within_alpha", "abc")), "expects a number")
    expect_error(parse_args(c("--min_coverage", "ten")), "expects a whole number")
})

test_that("an unknown option is rejected rather than silently ignored", {
    # A hyphen/underscore typo previously ran the entire job at defaults.
    expect_error(parse_args(c("--within-alpha", "0.5")), "unknown option")
    expect_error(parse_args(c("--nonsense", "1")), "unknown option")
})

test_that("a bad flag value and a dangling option are rejected", {
    expect_error(parse_args(c("--rate_between_groups", "maybe")),
                 "expects true/false")
    expect_error(parse_args(c("--within_alpha")), "no value")
})

test_that("defaults are untouched when no arguments are given", {
    cfg <- parse_args(character(0))
    expect_type(cfg$within_alpha, "double")
    expect_type(cfg$between_alpha, "double")
    expect_equal(cfg$metilene_q_source, "BH")
    expect_equal(cfg$metilene_p_column, "2dks")
})

test_that("config keys the runner relies on are all present", {
    # A missing key silently becomes NULL and hits a %||% fallback - which is
    # how the PrEC config-rebuild produced a run with no methods set.
    cfg <- parse_args(character(0))
    required <- c("methods", "scenarios", "effect_sizes", "within_alpha",
                  "between_alpha", "min_coverage", "dispersion_s", "sim_mode",
                  "metilene_min_diff", "metilene_min_effect", "metilene_qval",
                  "metilene_q_source", "metilene_p_column", "metilene_mtc",
                  "rate_between_groups")
    for (k in required) {
        expect_false(is.null(cfg[[k]]), info = paste("missing config key:", k))
    }
})
