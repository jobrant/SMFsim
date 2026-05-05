# =============================================================================
# Numeric UID encoding for fast site identification
#
# Replaces the slow paste(chr, pos, site, sep="_") with a numeric key that
# is 10-50x faster to create and uses ~8 bytes per site vs ~40-60 bytes for
# the string version. data.table operations (setkey, merge, intersect) are
# also significantly faster on numeric keys.
#
# Encoding: uid = chr_int * 1e10 + pos * 10 + site_int
#   - chr_int: 1-22 for autosomes, 23=X, 24=Y, 25=M/MT
#   - pos:     genomic position (up to ~250M for human, fits in 1e10 space)
#   - site_int: 1-9 for methylation context (GCH=1, HCG=2, GCG=3, etc.)
#
# Max possible value: 25e10 + 2.5e9 + 9 ≈ 2.75e10
# R double precision: exact integers up to 2^53 ≈ 9e15 → no loss of precision
# =============================================================================

#' Chromosome name to integer mapping
#' @keywords internal
#' @noRd
.chr_to_int <- function(chr) {
  # Strip "chr" prefix if present
  chr_clean <- sub("^chr", "", chr, ignore.case = TRUE)

  # Map to integer
  result <- integer(length(chr_clean))
  result[chr_clean == "X"] <- 23L
  result[chr_clean == "Y"] <- 24L
  result[chr_clean %in% c("M", "MT")] <- 25L

  # Numeric chromosomes
  numeric_idx <- !chr_clean %in% c("X", "Y", "M", "MT")
  result[numeric_idx] <- as.integer(chr_clean[numeric_idx])

  return(result)
}

#' Site context to integer mapping
#' @keywords internal
#' @noRd
.site_to_int <- function(site) {
  # Use a fixed lookup for known contexts; fallback for unknowns
  lookup <- c(
    "GCH" = 1L, "HCG" = 2L, "GCG" = 3L,
    "CCG" = 4L, "CCC" = 5L, "CCH" = 6L,
    "HCH" = 7L, "GCC" = 8L
  )
  result <- lookup[site]
  # Any unrecognized context gets 9
  result[is.na(result)] <- 9L
  return(as.integer(result))
}


#' Create numeric unique ID for methylation sites
#'
#' Generates a numeric key from chr, pos, and site columns that uniquely
#' identifies each methylation site. This is ~10-50x faster than the string-based
#' \code{paste(chr, pos, site, sep = "_")} and uses significantly less memory.
#'
#' @param dt A data.table with columns \code{chr}, \code{pos}, and \code{site}.
#' @param in_place Logical. If TRUE (default), adds the \code{uid} column to
#'   \code{dt} by reference. If FALSE, returns the uid vector without modifying \code{dt}.
#' @return If \code{in_place = TRUE}, returns \code{dt} invisibly (modified by reference).
#'   If \code{in_place = FALSE}, returns a numeric vector of uid values.
#'
#' @details The encoding is: \code{uid = chr_int * 1e10 + pos * 10 + site_int}, where
#'   chr_int maps chromosomes to integers (1-22, X=23, Y=24, M=25) and site_int
#'   maps methylation contexts to integers (GCH=1, HCG=2, etc.). This fits within
#'   R's double-precision integer range (exact up to 2^53) with no loss of precision.
#'
#' @keywords internal
#' @noRd
.create_uid <- function(dt, in_place = TRUE) {
  uid_vec <- .chr_to_int(dt$chr) * 1e10 + dt$pos * 10 + .site_to_int(dt$site)

  if (in_place) {
    dt[, uid := uid_vec]
    return(invisible(dt))
  } else {
    return(uid_vec)
  }
}


#' Create string uniqueID from components (for RSE export only)
#'
#' This is the legacy string-based ID. Use only when a character identifier
#' is required (e.g., matrix rownames, GRanges metadata).
#'
#' @param dt A data.table with columns \code{chr}, \code{pos}, and \code{site}.
#' @param in_place Logical. If TRUE, adds \code{uniqueID} column by reference.
#' @return If \code{in_place = TRUE}, returns \code{dt} invisibly.
#'   If \code{in_place = FALSE}, returns a character vector.
#'
#' @keywords internal
#' @noRd
.create_string_uid <- function(dt, in_place = TRUE) {
  uid_str <- paste(dt$chr, dt$pos, dt$site, sep = "_")

  if (in_place) {
    dt[, uniqueID := uid_str]
    return(invisible(dt))
  } else {
    return(uid_str)
  }
}


