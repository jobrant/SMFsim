# Tests for q-value sourcing and the separated effect-size filter.
#
# Background: `q_source = "BH"` recomputes Benjamini-Hochberg in R from
# metilene's MWU p-value, but only over the rows metilene emitted. With
# metilene's -d pre-filter active those rows are a high-difference SUBSET, so
# p.adjust() under-corrects -- it produced 819 false positives in a true null
# where the column-4 q gave 0. The honest version emits every segment
# (min_diff = 0), corrects across all of them, and applies the effect-size cut
# separately via min_effect.

test_that("BH with a non-zero min_diff warns about selection bias", {
    expect_warning(
        try(call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                               metilene_path = "", min_diff = 0.1,
                               q_source = "BH"),
            silent = TRUE),
        "ANTI-CONSERVATIVE")
})

test_that("BH with min_diff = 0 does not warn", {
    # The valid configuration: every segment emitted, so the BH denominator is
    # the full set of tested segments.
    expect_no_warning(
        try(call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                               metilene_path = "", min_diff = 0,
                               q_source = "BH", min_effect = 0.1),
            silent = TRUE))
})

test_that("the default q_source does not warn at any min_diff", {
    expect_no_warning(
        try(call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                               metilene_path = "", min_diff = 0.1),
            silent = TRUE))
})

test_that("the default config is a self-consistent own-FDR workflow", {
    cfg <- parse_args()

    # metilene emits every tested segment...
    expect_equal(cfg$metilene_min_diff, 0)
    # ...BH is computed in R across all of them...
    expect_equal(cfg$metilene_q_source, "BH")
    # ...and effect size is applied separately.
    expect_gt(cfg$metilene_min_effect, 0)
})

test_that("the default p-value column is the one metilene actually corrects", {
    # Verified empirically (inst/scripts/diagnose_metilene_qvalues.R): backing
    # the denominator out of a Bonferroni run gives q/p_2dks constant at
    # metilene's reported test count (985949, relative spread ~1e-4), while
    # q/p_mwu ranges over ~1200 orders of magnitude. metilene's manual claims
    # MWU and is wrong for de-novo mode. Using "mwu" would make q_bh and
    # q_metilene incomparable quantities.
    cfg <- parse_args()
    expect_equal(cfg$metilene_p_column, "2dks")

    expect_equal(formals(call_dmrs_metilene)$p_column[[2]], "2dks")
})

test_that("the default config does not trip the selection-bias warning", {
    # Guards the combination as a whole: if min_diff is ever raised above 0
    # while q_source stays "BH", the BH denominator silently becomes a
    # pre-selected subset again. This test fails if that happens.
    cfg <- parse_args()

    expect_no_warning(
        try(call_dmrs_metilene(split_data = NULL, out_dir = tempdir(),
                               metilene_path = "",
                               min_diff = cfg$metilene_min_diff,
                               q_source = cfg$metilene_q_source,
                               min_effect = cfg$metilene_min_effect),
            silent = TRUE))
})
