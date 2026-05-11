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
          # Jump left pointer to halfway through the current window
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
  
  # Fix data.table internal reference after rbind loop (this nearly drove me mad)
  selected <- selected[]
  
  message(sprintf("Selected %d dense regions (median %d sites, median spacing %.0f bp)",
                  nrow(selected), median(selected$n_sites),
                  median(selected$median_spacing)))
  
  return(selected)
}
