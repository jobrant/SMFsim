#!/usr/bin/env Rscript

# 01_simulate_efficiency.R
# Core simulation functions for enzyme efficiency distortion and DMR spike-ins
#
# Efficiency scenarios model WITHIN-GROUP replicate-to-replicate variation,
# reflecting real experimental conditions where replicates are processed on
# different days with varying enzyme efficiencies. Both groups span similar
# efficiency ranges with no systematic between-group bias.

# Enzyme efficiency simulation --------------------------------------------

#' Simulate enzyme efficiency distortion on a single sample
#'
#' For each site, resample methylated counts from a Binomial distribution:
#'   new_mc ~ Binomial(original_mc, efficiency)
#' Coverage stays constant (efficiency affects methylation labeling, not
#' sequencing depth).
#'
#' @param dt data.table with columns: chr, pos, strand, site, mc, cov, rate
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
  mc_new <- rbinom(nrow(out), 
                   size = pmax(0L, as.integer(round(out$mc))), 
                   prob = efficiency)
  set(out, j = "mc", value = mc_new)
  set(out, j = "rate", value = ifelse(out$cov > 0, out$mc / out$cov, 0))
  
  return(out)
}


# Pseudo-group construction -----------------------------------------------

#' Resolve efficiency input to a per-replicate vector
#'
#' Handles both explicit vectors (per-replicate values) and min/max ranges
#' (legacy behavior, uniform sampling).
#'
#' @param eff Numeric vector: length n_reps (explicit) or length 2 (range).
#' @param n_reps Number of replicates.
#' @param label Group label for messages.
#' @return Numeric vector of length n_reps.
#' @keywords internal
.resolve_efficiencies <- function(eff, n_reps, label) {
  stopifnot(is.numeric(eff), all(eff > 0), all(eff <= 1))
  
  if (length(eff) == n_reps) {
    # Explicit per-replicate values
    return(eff)
  } else if (length(eff) == 2) {
    # Legacy range: sample uniformly
    message(sprintf("Group %s: sampling %d efficiencies from [%.2f, %.2f]",
                    label, n_reps, eff[1], eff[2]))
    return(runif(n_reps, min = eff[1], max = eff[2]))
  } else {
    stop(sprintf(
      "efficiency_%s must be length %d (per-replicate) or length 2 (min/max range), got %d",
      label, n_reps, length(eff)
    ))
  }
}


#' Create pseudo-groups with per-replicate enzyme efficiency variation
#'
#' Takes a list of replicate data.tables and creates two pseudo-groups by
#' applying replicate-specific efficiency values. This models real experimental
#' conditions where individual replicates are processed on different days or
#' batches and have different enzyme efficiencies.
#'
#' Two modes are supported:
#' \describe{
#'   \item{mode = "clone"}{(Default) All replicates are used for BOTH groups.
#'     Each replicate is copied and distorted independently for group A and B,
#'     yielding a balanced N vs N design.}
#'   \item{mode = "split"}{Replicates are split in half: first half goes to
#'     group A, second half to group B. Requires ≥4 replicates for balanced design.}
#' }
#'
#' @param replicates Named list of data.tables (≥2 replicates).
#' @param efficiency_A Numeric vector: per-replicate efficiencies for group A.
#'   Can be length n_reps (explicit per-replicate) or length 2 (min/max range).
#' @param efficiency_B Same as efficiency_A but for group B.
#' @param mode Character: "clone" (default) or "split".
#' @param seed Base seed for reproducibility.
#' @return List with elements:
#'   \item{PseudoA}{Named list of distorted data.tables}
#'   \item{PseudoB}{Named list of distorted data.tables}
#'   \item{params}{List recording simulation parameters}
create_pseudo_groups <- function(replicates,
                                 efficiency_A = c(0.95, 0.70, 0.85),
                                 efficiency_B = c(0.90, 0.65, 0.80),
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
    eff_vals_A <- .resolve_efficiencies(efficiency_A, n_reps, "A")
    eff_vals_B <- .resolve_efficiencies(efficiency_B, n_reps, "B")
    
    message(sprintf("Clone mode: %d reps → %d vs %d design", n_reps, n_reps, n_reps))
    message(sprintf("PseudoA efficiencies: %s",
                    paste(round(eff_vals_A, 3), collapse = ", ")))
    message(sprintf("PseudoB efficiencies: %s",
                    paste(round(eff_vals_B, 3), collapse = ", ")))
    message(sprintf("PseudoA efficiency range: %.3f (sd=%.3f)",
                    diff(range(eff_vals_A)), sd(eff_vals_A)))
    message(sprintf("PseudoB efficiency range: %.3f (sd=%.3f)",
                    diff(range(eff_vals_B)), sd(eff_vals_B)))
    message(sprintf("Between-group mean difference: %.3f",
                    abs(mean(eff_vals_A) - mean(eff_vals_B))))
    
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
      mean_eff_A = mean(eff_vals_A),
      mean_eff_B = mean(eff_vals_B),
      sd_eff_A = sd(eff_vals_A),
      sd_eff_B = sd(eff_vals_B),
      seed = seed
    )
    
  } else {
    # Split mode: reps divided between groups
    n_A <- ceiling(n_reps / 2)
    n_B <- n_reps - n_A
    
    if (n_B < 2) {
      warning(sprintf(
        "Split mode with %d reps gives %d vs %d — consider mode = 'clone' for better power.",
        n_reps, n_A, n_B
      ))
    }
    
    idx_A <- seq_len(n_A)
    idx_B <- seq(n_A + 1, n_reps)
    
    eff_vals_A <- .resolve_efficiencies(efficiency_A, n_A, "A")
    eff_vals_B <- .resolve_efficiencies(efficiency_B, n_B, "B")
    
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


# Spike-in region selection -----------------------------------------------

#' Select random genomic regions for spike-in DMRs (legacy)
#'
#' Samples contiguous regions from the data to serve as ground-truth DMRs.
#' For sparse data, prefer select_dense_spikein_regions() instead.
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
                                   seed = 123) {
  set.seed(seed)
  
  # Get unique chromosomes and their site ranges
  chr_ranges <- dt[, .(min_pos = min(pos), max_pos = max(pos), n = .N), by = chr]
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


#' Select spike-in regions at naturally dense site clusters
#'
#' Scans the data for regions where sites are closely spaced, suitable for
#' DMR calling by tools like metilene. Uses an efficient two-pointer algorithm.
#'
#' @param dt Reference data.table (any sample; needs chr, pos columns).
#' @param n_regions Number of spike-in regions to select.
#' @param region_width_bp Maximum width of each region in base pairs.
#' @param min_sites Minimum sites within a region to be eligible.
#' @param max_median_spacing Maximum median inter-site spacing in bp.
#' @param seed Random seed.
#' @return data.table with columns: region_id, chr, start, end, n_sites, median_spacing
#' @noRd
select_dense_spikein_regions <- function(dt,
                                         n_regions = 300,
                                         region_width_bp = 2000,
                                         min_sites = 10,
                                         max_median_spacing = 200,
                                         seed = 123) {
  set.seed(seed)
  
  chr_list <- unique(dt$chr)
  candidates <- list()
  
  message("Scanning for dense regions...")
  t0 <- Sys.time()
  
  for (ch in chr_list) {
    chr_dt <- dt[chr == ch]
    if (nrow(chr_dt) < min_sites) next
    
    setorder(chr_dt, pos)
    positions <- chr_dt$pos
    n <- length(positions)
    
    # Two-pointer scan
    right <- 1L
    left <- 1L
    
    while (left <= n - min_sites + 1L) {
      # Advance right pointer to furthest site within the window
      while (right < n && (positions[right + 1L] - positions[left]) <= region_width_bp) {
        right <- right + 1L
      }
      
      n_in_window <- right - left + 1L
      
      if (n_in_window >= min_sites) {
        spacings <- diff(positions[left:right])
        med_spacing <- median(spacings)
        
        if (med_spacing <= max_median_spacing) {
          candidates[[length(candidates) + 1L]] <- data.table(
            chr = ch,
            start = positions[left],
            end = positions[right],
            n_sites = n_in_window,
            median_spacing = med_spacing
          )
          
          # Skip ahead past this region to avoid heavy overlap
          left <- left + max(1L, as.integer(n_in_window / 2))
          next
        }
      }
      
      left <- left + 1L
    }
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  candidates_dt <- rbindlist(candidates)
  message(sprintf("Found %d candidate dense regions in %.1f seconds",
                  nrow(candidates_dt), elapsed))
  
  if (nrow(candidates_dt) == 0) {
    stop("No dense regions found. Try increasing region_width_bp or decreasing min_sites.")
  }
  
  # Greedy non-overlapping selection: pick densest first
  setorder(candidates_dt, -n_sites)
  
  selected <- data.table()
  
  for (i in seq_len(nrow(candidates_dt))) {
    if (nrow(selected) >= n_regions) break
    
    r <- candidates_dt[i]
    
    if (nrow(selected) > 0) {
      overlaps <- selected[chr == r$chr & start <= r$end & end >= r$start]
      if (nrow(overlaps) > 0) next
    }
    
    selected <- rbind(selected, r)
  }
  
  selected[, region_id := seq_len(.N)]
  
  # Fix data.table internal reference after rbind loop
  selected <- selected[]
  
  message(sprintf("Selected %d dense regions (median %d sites, median spacing %.0f bp)",
                  nrow(selected), median(selected$n_sites),
                  median(selected$median_spacing)))
  
  return(selected)
}


# Spike-in DMR injection --------------------------------------------------

#' Inject spike-in DMRs into a sample
#'
#' For sites within spike-in regions, shift the methylation rate by a specified
#' effect size, then resample mc from Binomial(cov, new_rate).
#'
#' Direction is chosen based on the mean baseline rate of each region to
#' maximize the achievable effect. Regions with high baseline rates receive
#' decreasing effects; regions with low baseline rates receive increasing
#' effects. This avoids wasting injections on clamped sites that cannot
#' move in the chosen direction.
#'
#' @param dt data.table of methylation data for one sample.
#' @param regions data.table of spike-in regions (from select_spikein_regions
#'   or select_dense_spikein_regions).
#' @param effect_size Numeric scalar: amount to shift rates.
#' @param direction Character: "auto" (default, rate-aware), "increase",
#'   "decrease", or "both" (random per region, legacy behavior).
#' @param low_rate_threshold Numeric: regions with mean rate below this get
#'   "increase" direction. Default 0.4.
#' @param high_rate_threshold Numeric: regions with mean rate above this get
#'   "decrease" direction. Default 0.6.
#' @param seed Random seed.
#' @return List with:
#'   \item{data}{Modified data.table}
#'   \item{truth}{data.table mapping region_id to injected effect per site}
#'   \item{stats}{Summary of achievable effects across regions}
#' @noRd
inject_spikein_dmrs <- function(dt,
                                regions,
                                effect_size = 0.15,
                                direction = "auto",
                                low_rate_threshold = 0.4,
                                high_rate_threshold = 0.6,
                                seed = NULL) {
  stopifnot(
    is.data.table(dt),
    is.data.table(regions),
    direction %in% c("auto", "increase", "decrease", "both")
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- copy(dt)
  
  # Track ground truth
  truth_list <- list()
  n_injected <- 0L
  n_skipped <- 0L
  achievable_effects <- numeric(nrow(regions))
  
  for (i in seq_len(nrow(regions))) {
    r <- regions[i]
    
    # Find sites in this region
    site_idx <- which(out$chr == r$chr & out$pos >= r$start & out$pos <= r$end)
    
    if (length(site_idx) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }
    
    old_rates <- out$rate[site_idx]
    mean_rate <- mean(old_rates)
    
    # Determine direction
    dir_sign <- if (direction == "auto") {
      if (mean_rate < low_rate_threshold) {
        1L   # room to increase
      } else if (mean_rate > high_rate_threshold) {
        -1L  # room to decrease
      } else {
        # Middle range: pick direction with most room to move
        room_up   <- 1 - mean_rate
        room_down <- mean_rate
        if (room_up >= room_down) 1L else -1L
      }
    } else if (direction == "increase") {
      1L
    } else if (direction == "decrease") {
      -1L
    } else {
      # "both" — legacy random behavior
      sample(c(1L, -1L), 1L)
    }
    
    # Compute new rates
    new_rates <- pmin(1, pmax(0, old_rates + dir_sign * effect_size))
    
    # Track achievable effect (actual mean shift after clamping)
    actual_mean_effect <- mean(new_rates - old_rates)
    achievable_effects[i] <- actual_mean_effect
    
    # Resample mc from Binomial(cov, new_rate)
    covs <- out$cov[site_idx]
    new_mc <- rbinom(length(site_idx), size = as.integer(covs), prob = new_rates)
    set(out, i = site_idx, j = "mc", value = new_mc)
    set(out, i = site_idx, j = "rate", 
        value = ifelse(out$cov[site_idx] > 0, new_mc / out$cov[site_idx], 0))
    
    truth_list[[i]] <- data.table(
      region_id     = r$region_id,
      chr           = r$chr,
      pos           = out$pos[site_idx],
      original_rate = old_rates,
      injected_rate = new_rates,
      actual_effect = new_rates - old_rates,
      direction     = ifelse(dir_sign == 1L, "increase", "decrease"),
      mean_baseline = mean_rate
    )
    
    n_injected <- n_injected + 1L
  }
  
  truth <- rbindlist(truth_list)
  
  # Summary statistics
  achievable_effects <- achievable_effects[achievable_effects != 0]
  stats <- list(
    n_regions_injected = n_injected,
    n_regions_skipped  = n_skipped,
    n_sites_affected   = nrow(truth),
    mean_achievable_effect = mean(abs(achievable_effects)),
    median_achievable_effect = median(abs(achievable_effects)),
    pct_regions_gt_half_effect = mean(abs(achievable_effects) >= effect_size / 2) * 100
  )
  
  message(sprintf(
    "Injected DMRs in %d regions (%d sites, effect=%.2f, median achievable=%.3f, %.0f%% >%.2f)",
    n_injected, nrow(truth), effect_size,
    stats$median_achievable_effect,
    stats$pct_regions_gt_half_effect,
    effect_size / 2
  ))
  
  return(list(data = out, truth = truth, stats = stats))
}



# Scenario definitions ----------------------------------------------------

#' Get predefined efficiency scenarios
#'
#' Scenarios model WITHIN-GROUP replicate-to-replicate efficiency variation,
#' reflecting real experimental conditions where replicates are processed on
#' different days with varying enzyme activity. Both groups span similar
#' efficiency ranges; there is no systematic between-group bias. The artifact
#' is per-replicate, not per-group.
#'
#' Each scenario specifies explicit per-replicate efficiencies for 3 replicates
#' per group. Group means are approximately matched to avoid confounding
#' within-group variation with between-group differences.
#'
#' @return Named list of scenarios, each with:
#'   \item{label}{Human-readable description}
#'   \item{efficiency_A}{Numeric vector of 3 per-replicate efficiencies for group A}
#'   \item{efficiency_B}{Numeric vector of 3 per-replicate efficiencies for group B}
get_efficiency_scenarios <- function() {
  list(
    mild = list(
      label = "Mild within-group variation (sd ~ 0.04)",
      # Group A: mean ~ 0.90, range 0.85-0.95
      efficiency_A = c(0.95, 0.88, 0.87),
      # Group B: mean ~ 0.90, range 0.86-0.93
      efficiency_B = c(0.86, 0.93, 0.90)
    ),
    moderate = list(
      label = "Moderate within-group variation (sd ~ 0.10)",
      # Group A: mean ~ 0.83, range 0.70-0.95
      efficiency_A = c(0.95, 0.70, 0.85),
      # Group B: mean ~ 0.83, range 0.65-0.95
      efficiency_B = c(0.90, 0.65, 0.95)
    ),
    severe = list(
      label = "Severe within-group variation (sd ~ 0.17)",
      # Group A: mean ~ 0.75, range 0.55-0.95
      efficiency_A = c(0.55, 0.95, 0.75),
      # Group B: mean ~ 0.73, range 0.50-0.90
      efficiency_B = c(0.90, 0.50, 0.80)
    )
  )
}


#' Get predefined effect sizes for spike-in simulation
#'
#' @return Numeric vector of effect sizes to test.
get_effect_sizes <- function() {
  c(0.10, 0.15, 0.20, 0.30)
}


