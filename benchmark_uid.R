#!/usr/bin/env Rscript
# =============================================================================
# benchmark_uid.R
# Quick benchmark to compare string uniqueID vs numeric uid performance
#
# Run this after loading your data to see the expected speedup.
# Usage:
#   source("benchmark_uid.R")
#   benchmark_uid(wt_reps[[1]])
# =============================================================================

library(data.table)

# Source the uid helper functions
# source("path/to/MAPitNorm/R/utils.R")  # or wherever uid functions live

benchmark_uid <- function(dt, n_trials = 3) {
  dt <- copy(dt)  # don't modify original

  cat("Benchmark: String uniqueID vs Numeric uid\n")
  cat(sprintf("Dataset: %s rows\n", format(nrow(dt), big.mark = ",")))
  cat(strrep("=", 50), "\n\n")

  # ── 1. Creation speed ──
  cat("1. ID CREATION\n")

  # String (old method)
  string_times <- replicate(n_trials, {
    dt_copy <- copy(dt)
    if ("uniqueID" %in% names(dt_copy)) dt_copy[, uniqueID := NULL]
    system.time(dt_copy[, uniqueID := paste(chr, pos, site, sep = "_")])["elapsed"]
  })

  # Numeric (new method)
  numeric_times <- replicate(n_trials, {
    dt_copy <- copy(dt)
    if ("uid" %in% names(dt_copy)) dt_copy[, uid := NULL]
    system.time(.create_uid(dt_copy))["elapsed"]
  })

  cat(sprintf("  String paste():  %.2f sec (mean of %d trials)\n",
              mean(string_times), n_trials))
  cat(sprintf("  Numeric uid:     %.2f sec (mean of %d trials)\n",
              mean(numeric_times), n_trials))
  cat(sprintf("  Speedup:         %.1fx\n\n", mean(string_times) / mean(numeric_times)))

  # ── 2. Memory usage ──
  cat("2. MEMORY USAGE\n")

  dt_string <- copy(dt)
  dt_string[, uniqueID := paste(chr, pos, site, sep = "_")]
  string_size <- object.size(dt_string$uniqueID)

  dt_numeric <- copy(dt)
  .create_uid(dt_numeric)
  numeric_size <- object.size(dt_numeric$uid)

  cat(sprintf("  String column:   %s\n", format(string_size, units = "MB")))
  cat(sprintf("  Numeric column:  %s\n", format(numeric_size, units = "MB")))
  cat(sprintf("  Memory savings:  %.1fx\n\n", as.numeric(string_size) / as.numeric(numeric_size)))

  # ── 3. setkey speed ──
  cat("3. SETKEY SPEED\n")

  key_string_times <- replicate(n_trials, {
    dt_copy <- copy(dt_string)
    system.time(setkey(dt_copy, uniqueID))["elapsed"]
  })

  key_numeric_times <- replicate(n_trials, {
    dt_copy <- copy(dt_numeric)
    system.time(setkey(dt_copy, uid))["elapsed"]
  })

  cat(sprintf("  setkey(string):  %.2f sec\n", mean(key_string_times)))
  cat(sprintf("  setkey(numeric): %.2f sec\n", mean(key_numeric_times)))
  cat(sprintf("  Speedup:         %.1fx\n\n", mean(key_string_times) / mean(key_numeric_times)))

  # ── 4. Intersection speed (simulates find_shared_sites) ──
  cat("4. INTERSECTION SPEED (simulates find_shared_sites)\n")

  # Create a second sample with 90% overlap
  n <- nrow(dt)
  keep_idx <- sort(sample(n, floor(n * 0.9)))

  # String intersection
  ids_1 <- dt_string$uniqueID
  ids_2 <- dt_string$uniqueID[keep_idx]

  intersect_string_times <- replicate(n_trials, {
    system.time(intersect(ids_1, ids_2))["elapsed"]
  })

  # Numeric intersection
  uids_1 <- dt_numeric$uid
  uids_2 <- dt_numeric$uid[keep_idx]

  intersect_numeric_times <- replicate(n_trials, {
    system.time(intersect(uids_1, uids_2))["elapsed"]
  })

  cat(sprintf("  intersect(string):  %.2f sec\n", mean(intersect_string_times)))
  cat(sprintf("  intersect(numeric): %.2f sec\n", mean(intersect_numeric_times)))
  cat(sprintf("  Speedup:            %.1fx\n\n", mean(intersect_string_times) / mean(intersect_numeric_times)))

  # ── 5. Key-based subsetting speed ──
  cat("5. KEY-BASED SUBSETTING\n")

  shared_string <- intersect(ids_1, ids_2)
  shared_numeric <- intersect(uids_1, uids_2)

  setkey(dt_string, uniqueID)
  setkey(dt_numeric, uid)

  subset_string_times <- replicate(n_trials, {
    system.time(dt_string[shared_string])["elapsed"]
  })

  shared_dt <- data.table(uid = shared_numeric, key = "uid")
  subset_numeric_times <- replicate(n_trials, {
    system.time(dt_numeric[shared_dt, nomatch = NULL])["elapsed"]
  })

  cat(sprintf("  subset(string key):  %.2f sec\n", mean(subset_string_times)))
  cat(sprintf("  subset(numeric key): %.2f sec\n", mean(subset_numeric_times)))
  cat(sprintf("  Speedup:             %.1fx\n\n", mean(subset_string_times) / mean(subset_numeric_times)))

  # ── Summary ──
  cat(strrep("=", 50), "\n")
  cat("SUMMARY: Estimated total speedup for load + find_shared_sites\n")
  total_string <- mean(string_times) + mean(key_string_times) + mean(intersect_string_times) + mean(subset_string_times)
  total_numeric <- mean(numeric_times) + mean(key_numeric_times) + mean(intersect_numeric_times) + mean(subset_numeric_times)
  cat(sprintf("  String pipeline:  %.2f sec\n", total_string))
  cat(sprintf("  Numeric pipeline: %.2f sec\n", total_numeric))
  cat(sprintf("  Overall speedup:  %.1fx\n", total_string / total_numeric))

  invisible(list(
    creation = c(string = mean(string_times), numeric = mean(numeric_times)),
    memory = c(string = as.numeric(string_size), numeric = as.numeric(numeric_size)),
    setkey = c(string = mean(key_string_times), numeric = mean(key_numeric_times)),
    intersect = c(string = mean(intersect_string_times), numeric = mean(intersect_numeric_times)),
    subset = c(string = mean(subset_string_times), numeric = mean(subset_numeric_times))
  ))
}

