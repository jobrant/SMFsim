#!/usr/bin/env Rscript

# diagnose_spikein_calls.R
#
# What happens to each spike-in region, and what the false calls look like,
# from the saved metilene output of a multi-seed bias experiment. No metilene
# run and no simulation: it reads <seed>/<block>/spikein/<scenario>/
# effect_<e>/<method>/metilene_dmrs.bed (every tested segment that metilene
# emitted) and its log.
#
# Two questions it answers:
#
# 1. Why does sensitivity fall at the smallest effect in low-efficiency
#    scenarios? Each region is put in one class:
#      called         a significant segment with |difference| >= the effect cut
#      sig_below_cut  a significant segment, all below the effect cut
#      not_sig        segments overlap it, none significant
#      no_segment     no emitted segment overlaps it
#    A planted effect is observed as effect x efficiency, so at effect 0.10
#    many regions are significant and still under a 0.1 cut. Regions are also
#    split by the sign of the planted group A - group B difference, which
#    shows whether a method detects effects only where an efficiency bias
#    adds to them.
#
# 2. Where do the false calls come from? Their number, size (sites, width,
#    difference) and how they respond to stricter rules (lower q, more sites,
#    larger difference), next to the sensitivity each rule costs.
#
# Expects  <base>/seed_<n>/<block>/  as written by run_bias_experiment.R, and
# writes to <base>/all_seeds/<block>/:
#   spikein_region_status_by_seed.csv   one row per seed x scenario x effect x
#                                       method x region
#   spikein_region_status.csv           class fractions, pooled over seeds
#   spikein_region_status_by_sign.csv   the same, split by planted sign
#   spikein_false_calls_by_seed.csv     false-call counts and sizes per run
#   spikein_false_calls.csv             mean over seeds
#   spikein_filter_rules.csv            false calls and sensitivity under
#                                       stricter rules, mean over seeds
#
# Runs made before spikein_regions.csv recorded the planted direction: the
# direction is recomputed from the source data (it depends only on the pooled
# rate of each region) and written beside the original as
# spikein_regions_direction.csv. This loads the source replicates once and
# needs the data paths in the run's saved config to resolve, so run from the
# directory the experiment was run from. SMFSIM_SKIP_DIRECTION=1 skips it.
#
# Usage, from the manuscript root:
#   Rscript SMFsim/inst/scripts/diagnose_spikein_calls.R results/bias_experiment/multiplicative
#   SMFSIM_METHODS=raw,SMFnorm Rscript SMFsim/inst/scripts/diagnose_spikein_calls.R <base>

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
suppressWarnings(suppressMessages({
    library(data.table)
    library(GenomicRanges)
}))

.internal <- function(name) {
    if (exists(name, mode = "function")) get(name, mode = "function") else
        utils::getFromNamespace(name, "SMFsim")
}
.n_tests <- .internal(".metilene_n_tests")
.pooled  <- .internal(".pooled_rate")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: diagnose_spikein_calls.R <base_dir>", call. = FALSE)
base <- args[1]
if (!dir.exists(base)) stop("No such directory: ", base, call. = FALSE)
only_methods <- Sys.getenv("SMFSIM_METHODS", "")
only_methods <- if (nzchar(only_methods)) strsplit(only_methods, ",")[[1]] else NULL
skip_direction <- identical(Sys.getenv("SMFSIM_SKIP_DIRECTION"), "1")

seed_dirs <- list.dirs(base, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("^seed_-?[0-9]+$", basename(seed_dirs))]
if (!length(seed_dirs)) stop("No seed_<n> directories under ", base, call. = FALSE)
seeds <- as.integer(sub("^seed_", "", basename(seed_dirs)))
seed_dirs <- seed_dirs[order(seeds)]
seeds <- sort(seeds)

subdirs <- function(path) list.dirs(path, recursive = FALSE, full.names = FALSE)
blocks <- unique(unlist(lapply(seed_dirs, function(d) {
    b <- subdirs(d)
    b[file.exists(file.path(d, b, "spikein_regions.csv"))]
})))
if (!length(blocks)) stop("No spike-in results under ", base, call. = FALSE)

as_gr <- function(x) GRanges(x$chr, IRanges(x$start, x$end))

# Every emitted segment, with the BH q-value over all of metilene's tests.
read_segments <- function(dir, p_column) {
    d <- fread(file.path(dir, "metilene_dmrs.bed"), header = FALSE,
               colClasses = list(character = 1))
    setnames(d, 1:10, c("chr", "start", "end", "q_metilene", "diff", "n_sites",
                        "p_mwu", "p_2dks", "mean_A", "mean_B"))
    n <- .n_tests(file.path(dir, "metilene_log.txt"))
    p <- if (identical(p_column, "mwu")) d$p_mwu else d$p_2dks
    d[, q := stats::p.adjust(p, method = "BH", n = max(n, .N, na.rm = TRUE))]
    d[, width := end - start]
    attr(d, "n_tests") <- n
    d
}

# Source replicates, loaded once: the planted direction does not depend on the
# seed, only on the regions.
wt_cache <- NULL
planted_direction <- function(block_dir, regions) {
    if ("planted_sign" %in% names(regions)) return(regions)
    side <- file.path(block_dir, "spikein_regions_direction.csv")
    if (file.exists(side)) {
        return(fread(side, colClasses = list(character = "chr")))
    }
    if (skip_direction) return(regions)
    if (is.null(wt_cache)) {
        rds <- list.files(block_dir, pattern = "_results\\.rds$", full.names = TRUE)
        if (!length(rds)) {
            message("  no saved config in ", block_dir, "; planted direction skipped")
            return(regions)
        }
        cfg <- readRDS(rds[1])$config
        message("  loading source replicates to recompute the planted direction")
        wt <- tryCatch(suppressMessages(prepare_wt_replicates(cfg)),
                       error = function(e) {
                           message("  could not load them (", conditionMessage(e),
                                   "); planted direction skipped")
                           NULL
                       })
        if (is.null(wt)) {
            skip_direction <<- TRUE
            return(regions)
        }
        wt_cache <<- list(ref = wt[[1]][, .(chr = as.character(chr), pos)],
                          p = .pooled(wt))
    }
    planted <- spikein_direction(wt_cache$ref, wt_cache$p, regions)
    regions <- cbind(regions, planted)
    fwrite(regions, side)
    regions
}

rules <- function(d, q_cut, eff_cut) list(
    current      = d$q < q_cut & abs(d$diff) >= eff_cut,
    q_0.01       = d$q < 0.01 & abs(d$diff) >= eff_cut,
    sites_15     = d$q < q_cut & abs(d$diff) >= eff_cut & d$n_sites >= 15,
    sites_20     = d$q < q_cut & abs(d$diff) >= eff_cut & d$n_sites >= 20,
    diff_plus_02 = d$q < q_cut & abs(d$diff) >= eff_cut + 0.02,
    significant  = d$q < q_cut)

for (b in blocks) {
    message("\n== Block: ", b)
    status <- list(); false_calls <- list(); rule_rows <- list()

    for (i in seq_along(seed_dirs)) {
        block_dir <- file.path(seed_dirs[i], b)
        reg_file <- file.path(block_dir, "spikein_regions.csv")
        if (!file.exists(reg_file)) next
        regions <- planted_direction(
            block_dir, fread(reg_file, colClasses = list(character = "chr")))
        if (!"planted_sign" %in% names(regions)) regions[, planted_sign := NA_integer_]
        if (!"density_stratum" %in% names(regions)) regions[, density_stratum := NA_character_]
        reg_gr <- as_gr(regions)

        rds <- list.files(block_dir, pattern = "_results\\.rds$", full.names = TRUE)
        cfg <- if (length(rds)) readRDS(rds[1])$config else list()
        q_cut   <- cfg$metilene_qval %||% 0.05
        eff_cut <- cfg$metilene_min_effect %||% 0.1
        p_col   <- cfg$metilene_p_column %||% "2dks"
        message(sprintf("  seed %d (q < %g, |difference| >= %g)", seeds[i], q_cut, eff_cut))

        for (sc in subdirs(file.path(block_dir, "spikein"))) {
            for (ef_dir in subdirs(file.path(block_dir, "spikein", sc))) {
                effect <- as.numeric(sub("^effect_", "", ef_dir))
                for (m in subdirs(file.path(block_dir, "spikein", sc, ef_dir))) {
                    if (!is.null(only_methods) && !m %in% only_methods) next
                    dir <- file.path(block_dir, "spikein", sc, ef_dir, m)
                    if (!file.exists(file.path(dir, "metilene_dmrs.bed")) ||
                        file.size(file.path(dir, "metilene_dmrs.bed")) == 0) next
                    d <- read_segments(dir, p_col)
                    key <- data.table(seed = seeds[i], scenario = sc,
                                      effect_size = effect, method = m)

                    h <- findOverlaps(as_gr(d), reg_gr)
                    o <- data.table(r = subjectHits(h), q = d$q[queryHits(h)],
                                    diff = d$diff[queryHits(h)])
                    per <- o[, {
                        sig <- q < q_cut
                        # The segment that decides the class: the largest
                        # significant one, else the most significant one.
                        best <- if (any(sig)) diff[sig][which.max(abs(diff[sig]))] else
                            diff[which.min(q)]
                        .(sig = any(sig), called = any(sig & abs(diff) >= eff_cut),
                          diff = best)
                    }, by = r]
                    st <- rep("no_segment", nrow(regions))
                    st[per$r] <- ifelse(per$called, "called",
                                        ifelse(per$sig, "sig_below_cut", "not_sig"))
                    obs <- rep(NA_real_, nrow(regions))
                    obs[per$r] <- per$diff
                    status[[length(status) + 1L]] <- cbind(key, data.table(
                        region_id = regions$region_id,
                        density_stratum = regions$density_stratum,
                        planted_sign = regions$planted_sign,
                        status = st, observed_diff = obs))

                    on <- seq_len(nrow(d)) %in% queryHits(h)
                    region_of <- rep(NA_integer_, nrow(d))
                    region_of[queryHits(h)] <- subjectHits(h)
                    rl <- rules(d, q_cut, eff_cut)
                    f <- d[rl$current & !on]
                    false_calls[[length(false_calls) + 1L]] <- cbind(key, data.table(
                        n_tests = attr(d, "n_tests"),
                        significant_off_region = sum(rl$significant & !on),
                        false_calls = nrow(f),
                        true_calls = sum(rl$current & on),
                        median_sites = as.numeric(stats::median(f$n_sites)),
                        frac_min_sites = if (nrow(f)) mean(f$n_sites == min(d$n_sites)) else NA_real_,
                        median_width_bp = as.numeric(stats::median(f$width)),
                        median_abs_diff = stats::median(abs(f$diff)),
                        frac_within_0.02_of_cut = if (nrow(f)) mean(abs(f$diff) < eff_cut + 0.02) else NA_real_,
                        frac_A_higher = if (nrow(f)) mean(f$diff > 0) else NA_real_))
                    rule_rows[[length(rule_rows) + 1L]] <- rbindlist(lapply(names(rl), function(nm) {
                        cbind(key, data.table(
                            rule = nm, false_calls = sum(rl[[nm]] & !on),
                            sensitivity = uniqueN(region_of[rl[[nm]] & on]) / nrow(regions)))
                    }))
                }
            }
        }
    }
    if (!length(status)) next
    status <- rbindlist(status); false_calls <- rbindlist(false_calls)
    rule_rows <- rbindlist(rule_rows)
    out_dir <- file.path(base, "all_seeds", b)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    classes <- c("called", "sig_below_cut", "not_sig", "no_segment")

    frac <- function(by) {
        x <- status[, .N, by = c(by, "status")]
        x[, fraction := round(N / sum(N), 4), by = by]
        w <- dcast(x, stats::as.formula(paste(paste(by, collapse = " + "), "~ status")),
                   value.var = "fraction", fill = 0)
        for (cl in setdiff(classes, names(w))) w[, (cl) := 0]
        n <- status[, .(n_regions = .N, n_seeds = uniqueN(seed)), by = by]
        w <- merge(n, w, by = by)
        w[, significant := round(called + sig_below_cut, 4)]
        setcolorder(w, c(by, "n_seeds", "n_regions", "called", "sig_below_cut",
                         "significant", "not_sig", "no_segment"))
        w[]
    }
    by0 <- c("scenario", "effect_size", "method")
    region_status <- frac(by0)
    med <- status[status %in% c("called", "sig_below_cut"),
                  .(median_abs_observed_diff = round(stats::median(abs(observed_diff)), 4)),
                  by = by0]
    region_status <- merge(region_status, med, by = by0, all.x = TRUE)

    fwrite(status, file.path(out_dir, "spikein_region_status_by_seed.csv"))
    fwrite(region_status, file.path(out_dir, "spikein_region_status.csv"))
    if (!all(is.na(status$planted_sign))) {
        fwrite(frac(c(by0, "planted_sign")),
               file.path(out_dir, "spikein_region_status_by_sign.csv"))
    }
    num <- setdiff(names(false_calls), c("seed", by0))
    fwrite(false_calls, file.path(out_dir, "spikein_false_calls_by_seed.csv"))
    fwrite(false_calls[, c(list(n_seeds = uniqueN(seed)),
                           lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 3))),
                       by = by0, .SDcols = num],
           file.path(out_dir, "spikein_false_calls.csv"))
    fwrite(rule_rows[, .(n_seeds = uniqueN(seed),
                         false_calls = round(mean(false_calls), 1),
                         sensitivity = round(mean(sensitivity), 4)),
                     by = c(by0, "rule")],
           file.path(out_dir, "spikein_filter_rules.csv"))
    message("  written to ", out_dir)
}
