#!/usr/bin/env Rscript
# =============================================================================
# select_promoter_spikein_regions.R
# Select spike-in DMR regions at gene promoters where GCH site density is high
# =============================================================================

library(data.table)

#' Get promoter coordinates from a TxDb annotation
#'
#' @param txdb A TxDb object, or NULL to use the default hg38 knownGene.
#' @param upstream Bases upstream of TSS (default: 1500).
#' @param downstream Bases downstream of TSS (default: 500).
#' @param chr_style Character: "UCSC" (chr1) or "bare" (1). Set to match your data.
#' @return data.table with columns: gene_id, chr, start, end, strand
get_promoter_coords <- function(txdb = NULL,
                                upstream = 1500,
                                downstream = 500,
                                chr_style = "bare") {
  
  if (is.null(txdb)) {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop("Install TxDb.Hsapiens.UCSC.hg38.knownGene:\n",
           "  BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  
  # Get gene-level coordinates (one entry per gene, not per transcript)
  genes_gr <- GenomicFeatures::genes(txdb)
  
  # Get promoter regions
  promoters_gr <- GenomicRanges::promoters(genes_gr,
                                           upstream = upstream,
                                           downstream = downstream)
  
  # Convert to data.table
  promoter_dt <- data.table(
    gene_id = names(promoters_gr),
    chr = as.character(GenomicRanges::seqnames(promoters_gr)),
    start = GenomicRanges::start(promoters_gr),
    end = GenomicRanges::end(promoters_gr),
    strand = as.character(GenomicRanges::strand(promoters_gr))
  )
  
  # Strip chr prefix if needed
  if (chr_style == "bare") {
    promoter_dt[, chr := sub("^chr", "", chr)]
  } else if (chr_style == "UCSC" && !grepl("^chr", promoter_dt$chr[1])) {
    promoter_dt[, chr := paste0("chr", chr)]
  }
  
  # Keep only standard chromosomes
  standard_chr <- c(as.character(1:22), "X", "Y")
  if (chr_style == "UCSC") standard_chr <- paste0("chr", standard_chr)
  promoter_dt <- promoter_dt[chr %in% standard_chr]
  
  message(sprintf("Retrieved %d gene promoters (-%d/+%d around TSS)",
                  nrow(promoter_dt), upstream, downstream))
  
  return(promoter_dt)
}


#' Select spike-in regions at promoters with high GCH density
#'
#' Intersects promoter coordinates with actual GCH sites in the data,
#' then selects promoters that have enough sites for metilene to call DMRs.
#'
#' @param dt Reference data.table (any sample; needs chr, pos columns).
#' @param promoters data.table from get_promoter_coords(), or NULL to fetch automatically.
#' @param n_regions Number of spike-in regions to select.
#' @param min_sites Minimum GCH sites within a promoter to be eligible.
#' @param min_density Minimum sites per kb within the promoter (optional filter).
#' @param max_spacing Maximum median inter-site spacing in bp (optional filter).
#' @param seed Random seed.
#' @param chr_style Character: "bare" or "UCSC", must match your data.
#' @return data.table with columns: region_id, gene_id, chr, start, end,
#'   n_sites, median_spacing, density_per_kb
select_promoter_spikein_regions <- function(dt,
                                            promoters = NULL,
                                            n_regions = 300,
                                            min_sites = 10,
                                            min_density = NULL,
                                            max_spacing = 300,
                                            seed = 123,
                                            chr_style = "bare") {
  set.seed(seed)
  
  # Get promoter coordinates if not provided
  if (is.null(promoters)) {
    promoters <- get_promoter_coords(chr_style = chr_style)
  }
  
  message("Evaluating GCH site density at promoters...")
  
  # For each promoter, count sites and compute spacing
  promoter_stats <- lapply(seq_len(nrow(promoters)), function(i) {
    p <- promoters[i]
    sites <- dt[chr == p$chr & pos >= p$start & pos <= p$end]
    
    if (nrow(sites) < 2) {
      return(data.table(
        gene_id = p$gene_id, chr = p$chr,
        start = p$start, end = p$end,
        n_sites = nrow(sites),
        median_spacing = NA_real_,
        max_spacing = NA_real_,
        density_per_kb = nrow(sites) / ((p$end - p$start) / 1000)
      ))
    }
    
    spacings <- diff(sort(sites$pos))
    
    data.table(
      gene_id = p$gene_id, chr = p$chr,
      start = p$start, end = p$end,
      n_sites = nrow(sites),
      median_spacing = median(spacings),
      max_spacing = max(spacings),
      density_per_kb = nrow(sites) / ((p$end - p$start) / 1000)
    )
  })
  
  stats_dt <- rbindlist(promoter_stats)
  
  # Filter for eligible promoters
  eligible <- stats_dt[n_sites >= min_sites]
  
  if (!is.null(max_spacing)) {
    eligible <- eligible[!is.na(median_spacing) & median_spacing <= max_spacing]
  }
  
  if (!is.null(min_density)) {
    eligible <- eligible[density_per_kb >= min_density]
  }
  
  message(sprintf("Eligible promoters: %d of %d (min_sites=%d, max_spacing=%s)",
                  nrow(eligible), nrow(stats_dt), min_sites,
                  ifelse(is.null(max_spacing), "none", as.character(max_spacing))))
  
  if (nrow(eligible) == 0) {
    stop("No promoters meet the selection criteria. ",
         "Try lowering min_sites or increasing max_spacing.")
  }
  
  # Report density stats for eligible promoters
  message(sprintf("Eligible promoter stats: median %d sites, median spacing %.0f bp, median %.1f sites/kb",
                  median(eligible$n_sites),
                  median(eligible$median_spacing, na.rm = TRUE),
                  median(eligible$density_per_kb)))
  
  # Sample n_regions from eligible promoters
  if (nrow(eligible) < n_regions) {
    warning(sprintf("Only %d eligible promoters available (requested %d). Using all.",
                    nrow(eligible), n_regions))
    selected <- eligible
  } else {
    selected <- eligible[sample(.N, n_regions)]
  }
  
  # Add region_id
  selected[, region_id := seq_len(.N)]
  setcolorder(selected, c("region_id", "gene_id", "chr", "start", "end",
                          "n_sites", "median_spacing", "density_per_kb"))
  
  message(sprintf("Selected %d promoter spike-in regions (median %d sites/region, median spacing %.0f bp)",
                  nrow(selected),
                  median(selected$n_sites),
                  median(selected$median_spacing, na.rm = TRUE)))
  
  return(selected)
}


#' Quick diagnostic: compare promoter vs random region density
#'
#' @param dt Reference data.table.
#' @param promoters data.table from get_promoter_coords().
#' @param n_random Number of random regions to sample for comparison.
#' @param region_width Width of random regions (should match promoter width).
#' @return Prints comparison and returns a list with both distributions.
compare_promoter_vs_random_density <- function(dt,
                                               promoters = NULL,
                                               n_random = 1000,
                                               region_width = 2000) {
  if (is.null(promoters)) {
    promoters <- get_promoter_coords(chr_style = ifelse(
      grepl("^chr", dt$chr[1]), "UCSC", "bare"
    ))
  }
  
  # Promoter density
  promo_sites <- sapply(seq_len(min(nrow(promoters), n_random)), function(i) {
    p <- promoters[i]
    dt[chr == p$chr & pos >= p$start & pos <= p$end, .N]
  })
  
  # Random region density
  chr_ranges <- dt[, .(min_pos = min(pos), max_pos = max(pos)), by = chr]
  chr_ranges <- chr_ranges[max_pos - min_pos > region_width * 2]
  
  random_sites <- sapply(seq_len(n_random), function(i) {
    chr_pick <- sample(chr_ranges$chr, 1)
    chr_info <- chr_ranges[chr == chr_pick]
    start_pos <- sample(chr_info$min_pos:(chr_info$max_pos - region_width), 1)
    dt[chr == chr_pick & pos >= start_pos & pos <= start_pos + region_width, .N]
  })
  
  message(sprintf("Promoter regions: median %d sites (IQR: %d-%d)",
                  median(promo_sites),
                  quantile(promo_sites, 0.25),
                  quantile(promo_sites, 0.75)))
  message(sprintf("Random regions:   median %d sites (IQR: %d-%d)",
                  median(random_sites),
                  quantile(random_sites, 0.25),
                  quantile(random_sites, 0.75)))
  message(sprintf("Promoters are %.1fx denser on average",
                  mean(promo_sites) / mean(random_sites)))
  
  return(invisible(list(promoter = promo_sites, random = random_sites)))
}