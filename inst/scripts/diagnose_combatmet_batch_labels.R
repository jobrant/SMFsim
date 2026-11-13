#!/usr/bin/env Rscript

# diagnose_combatmet_batch_labels.R
#
# Are ComBatMet's null false positives a ComBat_met problem, or an artifact of
# how SMFsim labels batches?
#
# In the severe null (matched group means, so every DMR is false) ComBatMet
# calls 7,791 DMRs against 17-31 for the other methods, and its corrected data
# show a group gap the input does not have (mean rate A 0.149 vs B 0.108;
# raw 0.1205 / 0.1206). SMFsim makes ComBat's batch label by splitting the
# eight per-replicate efficiencies at their median (.efficiency_batch()). In
# severe that puts 3 of A's 4 replicates in low_eff and 1 of B's, and within
# each batch A's replicates are the more efficient ones (low: A 0.67 vs B
# 0.50; high: A 0.95 vs B 0.82). With group in the model and protected, ComBat
# reads that within-batch difference as biology and keeps it.
#
# This script rebuilds each matched-mean null exactly as
# validate_null_calibration.R does (same settings, same seed) and runs
# ComBatMet under three batch labellings:
#   median    the current rule (should reproduce the null calibration count)
#   balanced  two levels, the top half of EACH group by efficiency = high_eff
#   paired    four levels, batch k = the k-th most efficient replicate of each
#             group: one A and one B per batch, as when an experiment processes
#             one sample per group in each batch. ComBat_met errored on this
#             design on chr21 (2026-11-13); an error is recorded, not fatal.
# For each it records the within-batch group efficiency gap (the suspected
# cause), the group mean rates after correction, and the false-positive DMRs.
#
# If balanced / paired remove the gap and the false positives, the median
# split is what fails, not ComBat_met.
#
# Usage, from the manuscript root (genome-wide ~15 min per ComBatMet run, so
# ~2.5 h for the default 3 scenarios x 3 labellings; chr1 ~10x faster):
#   Rscript SMFsim/inst/scripts/diagnose_combatmet_batch_labels.R
#
# Environment variables (all optional):
#   CBM_SCENARIOS  comma-separated scenarios      (default mild,moderate,severe)
#   CBM_CHR        chromosome regex, e.g. '^(chr)?1$' (default: all standard)
#   CBM_OUT        output directory  (default results/combatmet_batch_labels[_chr])
#   SMFSIM_OVERWRITE  1 to write into a non-empty output directory
#
# Output: combatmet_batch_labels.csv, one row per scenario x labelling.
# metilene_input/ is deleted after each run.

if (requireNamespace("devtools", quietly = TRUE) && file.exists("SMFsim/DESCRIPTION")) {
    suppressMessages(devtools::load_all("SMFsim/", quiet = TRUE))
} else {
    library(SMFsim)
}
library(data.table)
.internal <- function(f) if (exists(f)) get(f) else getFromNamespace(f, "SMFsim")
.n_tests_of     <- .internal(".metilene_n_tests")
.noise_of       <- .internal(".noise_from_config")
.median_batch   <- .internal(".efficiency_batch")

# --- Settings: match validate_null_calibration.R ------------------------------
scen_names  <- strsplit(Sys.getenv("CBM_SCENARIOS", "mild,moderate,severe"), ",")[[1]]
chr_pattern <- Sys.getenv("CBM_CHR", "")
out <- Sys.getenv("CBM_OUT", file.path(
    "results", paste0("combatmet_batch_labels",
                      if (nzchar(chr_pattern)) paste0("_", gsub("[^A-Za-z0-9]", "", chr_pattern)))))

config <- parse_args(character(0))
config$metilene_path <- "/apps/metilene/0.2.8/metilene"
config$wt_group_id   <- "M1"
config$sim_mode      <- "parametric"
config$dispersion_s  <- 128
config$seed          <- 42
if (nzchar(chr_pattern)) config$chr_pattern <- chr_pattern

scenarios <- get_efficiency_scenarios()
bad <- setdiff(scen_names, names(scenarios))
if (length(bad)) stop("Unknown scenario(s): ", paste(bad, collapse = ", "), call. = FALSE)
not_null <- scen_names[!vapply(scenarios[scen_names], function(s) isTRUE(s$matched_means),
                               logical(1))]
if (length(not_null)) {
    stop("Not matched-mean nulls (their groups differ in efficiency by design, so ",
         "relabelling would change the question): ", paste(not_null, collapse = ", "),
         call. = FALSE)
}

message(sprintf("ComBatMet batch-label diagnostic: %s, %s -> %s",
                paste(scen_names, collapse = "/"),
                if (nzchar(chr_pattern)) chr_pattern else "all standard chromosomes", out))

if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 &&
    !identical(Sys.getenv("SMFSIM_OVERWRITE"), "1")) {
    stop("Output directory is not empty: ", normalizePath(out, mustWork = FALSE),
         "\nArchive it, set CBM_OUT, or re-run with SMFSIM_OVERWRITE=1.", call. = FALSE)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# --- Labellings ----------------------------------------------------------------
# Each returns a factor, PseudoA replicates first then PseudoB.
labellings <- list(
    median   = function(eA, eB) .median_batch(list(efficiency_A = eA, efficiency_B = eB),
                                              length(eA), length(eB), rule = "median"),
    balanced = function(eA, eB) {
        top <- function(e) rank(-e, ties.method = "first") <= length(e) %/% 2
        factor(ifelse(c(top(eA), top(eB)), "high_eff", "low_eff"),
               levels = c("low_eff", "high_eff"))
    },
    paired   = function(eA, eB) {
        stopifnot(length(eA) == length(eB))
        factor(paste0("b", c(rank(-eA, ties.method = "first"),
                             rank(-eB, ties.method = "first"))))
    }
)

# Mean over batches of (A's mean efficiency - B's) within the batch, weighted
# by batch size; batches holding only one group contribute nothing. This is
# the group difference ComBat sees once batch is accounted for.
within_batch_gap <- function(eff, batch, group) {
    d <- data.table(eff, batch, group)
    per <- d[, .(gap = if (uniqueN(group) == 2)
                     mean(eff[group == "A"]) - mean(eff[group == "B"]) else NA_real_,
                 n = .N), by = batch]
    per[!is.na(gap), sum(gap * n) / sum(n)]
}

# --- Run -----------------------------------------------------------------------
reps <- prepare_wt_replicates(config)
res <- list()
for (sc_name in scen_names) {
    sc <- scenarios[[sc_name]]
    message("\n", strrep("=", 64), "\n", sc$label, "\n", strrep("=", 64))
    pseudo <- create_pseudo_groups(
        replicates = reps, efficiency_A = sc$efficiency_A, efficiency_B = sc$efficiency_B,
        mode = config$sim_mode, seed = config$seed, dispersion_s = config$dispersion_s,
        matched_means = isTRUE(sc$matched_means), noise = .noise_of(config),
        efficiency_model = config$efficiency_model %||% "multiplicative")
    eA <- pseudo$params$efficiency_A
    eB <- pseudo$params$efficiency_B
    group <- rep(c("A", "B"), c(length(eA), length(eB)))
    raw_means <- vapply(c(pseudo$PseudoA, pseudo$PseudoB), function(dt) mean(dt$rate), numeric(1))

    for (lab in names(labellings)) {
        batch <- labellings[[lab]](eA, eB)
        gap_eff <- within_batch_gap(c(eA, eB), batch, group)
        message(sprintf("\n--- %s / %s: within-batch efficiency gap A-B %+.3f ---",
                        sc_name, lab, gap_eff))
        t0 <- Sys.time()
        # ComBat_met itself can fail on some designs (paired, four batches of
        # two, gave "subscript out of bounds" in mu_hat_mat on chr21). Record it
        # and go on; the other labellings still answer the question.
        adj <- tryCatch(run_combatmet(pseudo, batch = batch)$data, error = function(e) e)
        if (inherits(adj, "error")) {
            message("  ComBat_met failed: ", conditionMessage(adj))
            res[[length(res) + 1]] <- data.table(
                scenario = sc_name, labelling = lab,
                batch = paste(as.character(batch), collapse = ","),
                efficiency = paste(round(c(eA, eB), 2), collapse = ","),
                within_batch_eff_gap = round(gap_eff, 4),
                raw_mean_A = mean(raw_means[group == "A"]),
                raw_mean_B = mean(raw_means[group == "B"]),
                status = paste("ComBat_met error:", conditionMessage(adj)))
            fwrite(rbindlist(res, fill = TRUE), file.path(out, "combatmet_batch_labels.csv"))
            next
        }
        m <- vapply(c(adj$PseudoA, adj$PseudoB), function(dt) mean(dt$rate), numeric(1))

        dir <- file.path(out, sc_name, lab)
        dmrs <- call_dmrs_metilene(
            adj, dir, config$metilene_path,
            min_cpg = config$metilene_min_cpg, min_diff = config$metilene_min_diff,
            metilene_max_dist = config$metilene_max_dist,
            q_cutoff = config$metilene_qval, q_source = config$metilene_q_source,
            p_column = config$metilene_p_column, min_effect = config$metilene_min_effect,
            mtc = config$metilene_mtc)
        # An empty table from a crashed metilene would read as a perfect null.
        n_tests <- .n_tests_of(file.path(dir, "metilene_log.txt"))
        if (is.na(n_tests)) {
            stop("metilene did not complete for ", sc_name, " / ", lab, ". See ", dir,
                 call. = FALSE)
        }
        unlink(file.path(dir, "metilene_input"), recursive = TRUE)

        a_higher <- if (nrow(dmrs) > 0 && "mean_diff" %in% names(dmrs))
            mean(dmrs$mean_diff > 0) else NA_real_
        res[[length(res) + 1]] <- data.table(
            scenario = sc_name, labelling = lab,
            batch = paste(as.character(batch), collapse = ","),
            efficiency = paste(round(c(eA, eB), 2), collapse = ","),
            within_batch_eff_gap = round(gap_eff, 4),
            raw_mean_A = mean(raw_means[group == "A"]), raw_mean_B = mean(raw_means[group == "B"]),
            adj_mean_A = mean(m[group == "A"]), adj_mean_B = mean(m[group == "B"]),
            FP = nrow(dmrs), frac_A_higher = a_higher, n_tests = n_tests,
            minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
            status = "ok")
        message(sprintf("  group means after: A %.4f  B %.4f  |  %d FP DMRs",
                        mean(m[group == "A"]), mean(m[group == "B"]), nrow(dmrs)))
        fwrite(rbindlist(res, fill = TRUE), file.path(out, "combatmet_batch_labels.csv"))
        rm(adj); invisible(gc())
    }
}

res <- rbindlist(res, fill = TRUE)
res[, adj_gap := adj_mean_A - adj_mean_B]
message("\n", strrep("=", 64), "\nComBatMet NULL FALSE POSITIVES BY BATCH LABELLING\n",
        strrep("=", 64))
print(res[, .(scenario, labelling, within_batch_eff_gap,
              adj_gap = round(adj_gap, 4), FP, frac_A_higher = round(frac_A_higher, 2),
              status = substr(status, 1, 40))])
message("\nResults: ", normalizePath(file.path(out, "combatmet_batch_labels.csv")))
