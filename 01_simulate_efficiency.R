#!/usr/bin/env Rscript
# =============================================================================
# 01_simulate_efficiency.R
# Core simulation functions for enzyme efficiency distortion and DMR spike-ins
# =============================================================================

library(data.table)

# Enzyme efficiency simulation --------------------------------------------

#' Simulate enzyme efficiency distortion on a single sample
#'
#' For each site, resample methylated counts from a Binomial distribution:
#'   new_mc ~ Binomial(original_mc, efficiency)
#' Coverage stays constant (efficiency affects methylation labeling, not
#' sequencing depth).
#'
#' @param dt data.table with columns: chr, pos, strand, site, uniqueID, mc, cov, rate
#' @param efficiency Numeric scalar in (0, 1]. The simulated enzyme efficiency.
#'   Values < 1 reduce observed methylation proportional to accessibility.
#' @param seed Optional integer seed for reproducibility.
#' @return data.table with distorted mc and recalculated rate (cov unchanged).
simulate_efficiency <- function(dt, efficiency, seed = NULL) {
  stopifnot(
    is.data.table(dt),
    all(c("mc", "cov", "rate") %in% names(dt)),
    is.numeric(efficiency), length(efficiency) == 1,
    efficiency > 0, efficiency <= 1
  )

  if (!is.null(seed)) set.seed(seed)

  out <- copy(dt)

  # Binomial resampling: each methylated read is retained with prob = efficiency

  # mc must be integer for rbinom; round in case of prior normalization
  mc_int <- pmax(0L, as.integer(round(out$mc)))
  out[, mc := rbinom(.N, size = mc_int, prob = efficiency)]
  out[, rate := ifelse(cov > 0, mc / cov, 0)]

  return(out)
}


#' Create pseudo-groups with differential enzyme efficiency
#'
#' Takes a list of replicate data.tables (e.g., WT or normal replicates) and
#' creates two pseudo-groups by applying different efficiency ranges.
#'
#' Two modes are supported:
#' \describe{
#'   \item{mode = "clone"}{(Default) All replicates are used for BOTH groups.
#'     Each replicate is copied and distorted independently for group A and B,
#'     yielding a balanced N vs N design. This is the recommended approach when
#'     you have few replicates (e.g., 3), as it maximizes statistical power for
#'     downstream DMR calling while maintaining the null (identical biology in
#'     both groups, differing only by simulated efficiency).}
#'   \item{mode = "split"}{Replicates are split in half: first half goes to
#'     group A, second half to group B. Each replicate appears in only one group.
#'     Requires ≥4 replicates for a balanced design (≥2 per group).}
#' }
#'
#' @param replicates Named list of data.tables (≥2 replicates).
#' @param efficiency_A Numeric vector of length 2: [min, max] efficiency for group A.
#' @param efficiency_B Numeric vector of length 2: [min, max] efficiency for group B.
#' @param mode Character: "clone" (default) or "split". See Details.
#' @param seed Base seed for reproducibility.
#' @return List with elements:
#'   \item{PseudoA}{Named list of distorted data.tables}
#'   \item{PseudoB}{Named list of distorted data.tables}
#'   \item{params}{List recording simulation parameters}
create_pseudo_groups <- function(replicates,
                                 efficiency_A = c(0.95, 1.0),
                                 efficiency_B = c(0.70, 0.80),
                                 mode = c("clone", "split"),
                                 seed = 42) {
  mode <- match.arg(mode)
  stopifnot(is.list(replicates), length(replicates) >= 2)

  n_reps <- length(replicates)
  rep_names <- names(replicates)
  if (is.null(rep_names)) rep_names <- paste0("Rep", seq_len(n_reps))

  set.seed(seed)

  if (mode == "clone") {
    # Clone mode: all reps used for both groups (N vs N)
    # Each replicate is independently distorted for A and B using different
    # efficiency values and different random seeds, so the two copies diverge
    # only due to the efficiency artifact.
    eff_vals_A <- runif(n_reps, min = efficiency_A[1], max = efficiency_A[2])
    eff_vals_B <- runif(n_reps, min = efficiency_B[1], max = efficiency_B[2])

    message(sprintf("Clone mode: %d reps → %d vs %d design", n_reps, n_reps, n_reps))
    message(sprintf("PseudoA efficiencies: %s",
                    paste(round(eff_vals_A, 3), collapse = ", ")))
    message(sprintf("PseudoB efficiencies: %s",
                    paste(round(eff_vals_B, 3), collapse = ", ")))

    pseudo_A <- mapply(function(dt, eff, i) {
      simulate_efficiency(dt, eff, seed = seed + i)
    }, replicates, eff_vals_A, seq_along(replicates), SIMPLIFY = FALSE)
    names(pseudo_A) <- paste0("PseudoA_", rep_names)

    pseudo_B <- mapply(function(dt, eff, i) {
      simulate_efficiency(dt, eff, seed = seed + 1000 + i)
    }, replicates, eff_vals_B, seq_along(replicates), SIMPLIFY = FALSE)
    names(pseudo_B) <- paste0("PseudoB_", rep_names)

    params <- list(
      mode = "clone",
      n_source_reps = n_reps,
      n_reps_A = n_reps, n_reps_B = n_reps,
      efficiency_A = eff_vals_A,
      efficiency_B = eff_vals_B,
      seed = seed
    )

  } else {
    # Split mode: reps divided between groups
    n_A <- ceiling(n_reps / 2)
    n_B <- n_reps - n_A

    if (n_B < 2) {
      warning(sprintf(
        "Split mode with %d reps gives %d vs %d — consider mode='clone' for better power.",
        n_reps, n_A, n_B
      ))
    }

    idx_A <- seq_len(n_A)
    idx_B <- seq(n_A + 1, n_reps)

    eff_vals_A <- runif(n_A, min = efficiency_A[1], max = efficiency_A[2])
    eff_vals_B <- runif(n_B, min = efficiency_B[1], max = efficiency_B[2])

    message(sprintf("Split mode: %d reps → %d vs %d design", n_reps, n_A, n_B))
    message(sprintf("PseudoA efficiencies: %s",
                    paste(round(eff_vals_A, 3), collapse = ", ")))
    message(sprintf("PseudoB efficiencies: %s",
                    paste(round(eff_vals_B, 3), collapse = ", ")))

    pseudo_A <- mapply(function(dt, eff, i) {
      simulate_efficiency(dt, eff, seed = seed + i)
    }, replicates[idx_A], eff_vals_A, seq_along(idx_A), SIMPLIFY = FALSE)
    names(pseudo_A) <- paste0("PseudoA_", rep_names[idx_A])

    pseudo_B <- mapply(function(dt, eff, i) {
      simulate_efficiency(dt, eff, seed = seed + 1000 + i)
    }, replicates[idx_B], eff_vals_B, seq_along(idx_B), SIMPLIFY = FALSE)
    names(pseudo_B) <- paste0("PseudoB_", rep_names[idx_B])

    params <- list(
      mode = "split",
      n_source_reps = n_reps,
      n_reps_A = n_A, n_reps_B = n_B,
      efficiency_A = eff_vals_A,
      efficiency_B = eff_vals_B,
      seed = seed
    )
  }

  return(list(PseudoA = pseudo_A, PseudoB = pseudo_B, params = params))
}


# Spike-in DMR injection --------------------------------------------------

#' Select random genomic regions for spike-in DMRs
#'
#' Samples contiguous regions from the data to serve as ground-truth DMRs.
#'
#' @param dt Reference data.table (any sample; needs chr, pos columns).
#' @param n_regions Number of DMR regions to select.
#' @param region_width_bp Width of each region in base pairs.
#' @param min_sites Minimum number of sites a region must contain to be valid.
#' @param seed Random seed.
#' @return data.table with columns: region_id, chr, start, end, n_sites
select_spikein_regions <- function(dt,
                                    n_regions = 300,
                                    region_width_bp = 500,
                                    min_sites = 5,
                                    seed = 42) {
  set.seed(seed)

  # Get unique chromosomes and their site ranges
  chr_ranges <- dt[, .(min_pos = min(pos), max_pos = max(pos), n = .N), by = chr]
  # Only use chromosomes with enough span

  chr_ranges <- chr_ranges[max_pos - min_pos > region_width_bp * 2]

  regions <- data.table(
    region_id = integer(),
    chr = character(),
    start = integer(),
    end = integer(),
    n_sites = integer()
  )

  attempts <- 0
  max_attempts <- n_regions * 20

  while (nrow(regions) < n_regions && attempts < max_attempts) {
    attempts <- attempts + 1

    # Pick a random chromosome weighted by number of sites
    chr_pick <- sample(chr_ranges$chr, 1, prob = chr_ranges$n)
    chr_info <- chr_ranges[chr == chr_pick]

    # Pick a random start position
    start_pos <- sample(chr_info$min_pos:(chr_info$max_pos - region_width_bp), 1)
    end_pos <- start_pos + region_width_bp

    # Count sites in this region
    n_sites_in_region <- dt[chr == chr_pick & pos >= start_pos & pos <= end_pos, .N]

    if (n_sites_in_region < min_sites) next

    # Check for overlap with existing regions
    if (nrow(regions) > 0) {
      overlaps <- regions[chr == chr_pick & start <= end_pos & end >= start_pos]
      if (nrow(overlaps) > 0) next
    }

    regions <- rbind(regions, data.table(
      region_id = nrow(regions) + 1L,
      chr = chr_pick,
      start = start_pos,
      end = end_pos,
      n_sites = n_sites_in_region
    ))
  }

  if (nrow(regions) < n_regions) {
    warning(sprintf("Only found %d non-overlapping regions (requested %d)",
                    nrow(regions), n_regions))
  }

  message(sprintf("Selected %d spike-in regions (median %d sites/region)",
                  nrow(regions), median(regions$n_sites)))

  return(regions)
}


#' Inject spike-in DMRs into a sample
#'
#' For sites within spike-in regions, shift the methylation rate by a specified
#' effect size, then resample mc from Binomial(cov, new_rate).
#'
#' @param dt data.table of methylation data for one sample.
#' @param regions data.table of spike-in regions (from select_spikein_regions).
#' @param effect_size Numeric scalar: amount to shift rates (positive = increase).
#' @param direction Character: "increase", "decrease", or "both" (random per region).
#' @param seed Random seed.
#' @return List with:
#'   \item{data}{Modified data.table}
#'   \item{truth}{data.table mapping region_id to injected effect per site}
inject_spikein_dmrs <- function(dt,
                                 regions,
                                 effect_size = 0.15,
                                 direction = "both",
                                 seed = NULL) {
  stopifnot(
    is.data.table(dt),
    is.data.table(regions),
    direction %in% c("increase", "decrease", "both")
  )

  if (!is.null(seed)) set.seed(seed)

  out <- copy(dt)

  # Track ground truth
  truth_list <- list()

  for (i in seq_len(nrow(regions))) {
    r <- regions[i]

    # Find sites in this region
    site_idx <- which(out$chr == r$chr & out$pos >= r$start & out$pos <= r$end)

    if (length(site_idx) == 0) next

    # Determine direction for this region
    if (direction == "both") {
      dir_sign <- sample(c(1, -1), 1)
    } else if (direction == "increase") {
      dir_sign <- 1
    } else {
      dir_sign <- -1
    }

    # Shift rates
    old_rates <- out$rate[site_idx]
    new_rates <- pmin(1, pmax(0, old_rates + dir_sign * effect_size))

    # Resample mc from Binomial(cov, new_rate)
    covs <- out$cov[site_idx]
    new_mc <- rbinom(length(site_idx), size = as.integer(covs), prob = new_rates)

    out[site_idx, `:=`(mc = new_mc, rate = ifelse(cov > 0, new_mc / cov, 0))]

    # Record truth
    truth_list[[i]] <- data.table(
      region_id = r$region_id,
      chr = r$chr,
      pos = out$pos[site_idx],
      original_rate = old_rates,
      injected_rate = new_rates,
      actual_effect = dir_sign * effect_size,
      direction = ifelse(dir_sign == 1, "increase", "decrease")
    )
  }

  truth <- rbindlist(truth_list)
  message(sprintf("Injected DMRs in %d regions (%d sites affected, effect = %.2f)",
                  nrow(regions), nrow(truth), effect_size))

  return(list(data = out, truth = truth))
}


# Scenario definitions ----------------------------------------------------

#' Get predefined efficiency scenarios
#'
#' @return Named list of scenarios, each with efficiency_A and efficiency_B ranges.
get_efficiency_scenarios <- function() {
  list(
    mild = list(
      label = "Mild (0.90 vs 1.0)",
      efficiency_A = c(0.95, 1.00),
      efficiency_B = c(0.85, 0.90)
    ),
    moderate = list(
      label = "Moderate (0.75 vs 1.0)",
      efficiency_A = c(0.95, 1.00),
      efficiency_B = c(0.70, 0.80)
    ),
    severe = list(
      label = "Severe (0.55 vs 1.0)",
      efficiency_A = c(0.95, 1.00),
      efficiency_B = c(0.50, 0.60)
    )
  )
}

#' Get predefined effect sizes for spike-in simulation
#'
#' @return Numeric vector of effect sizes to test.
get_effect_sizes <- function() {
  c(0.05, 0.10, 0.15, 0.20, 0.30)
}
