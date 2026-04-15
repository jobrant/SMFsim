# MAPitNorm Simulation Framework

Benchmarking framework for evaluating methylation normalization methods using
simulated enzyme efficiency artifacts and spike-in DMRs.

## Test Dataset

NOMe-seq data from Taberlay/Statham et al. (PMID 24916973, GEO: GSE57498).
Four human cell lines with the following replicate structure:

| Cell line | Type             | SRA accessions                    | Reps |
|-----------|------------------|-----------------------------------|------|
| HMEC      | Normal breast    | SRR1282202, SRR1282203            | 2    |
| MCF7      | Breast cancer    | SRR1282204, SRR1282205            | 2    |
| **PrEC**  | **Normal prostate** | **SRR1282206–SRR1282208**      | **3** |
| **PC3**   | **Prostate cancer** | **SRR1282209–SRR1282211**      | **3** |

We use the 3-replicate groups for simulation:
- **PrEC** (normal) as the source for null and spike-in simulations
- **PC3 vs PrEC** as a real biological comparison for validation

### Upstream processing (FASTQ → allc)

Raw FASTQs must be processed before MAPitNorm can use them:

```bash
# 1. Trim adapters
trim_galore --paired SRR1282206_1.fastq.gz SRR1282206_2.fastq.gz

# 2. Align with bisulfite-aware aligner (bwameth or bismark)
bwameth.py --reference hg38.fa \
  SRR1282206_1_val_1.fq.gz SRR1282206_2_val_2.fq.gz \
  | samtools sort -o SRR1282206.sorted.bam

# 3. Mark duplicates
picard MarkDuplicates I=SRR1282206.sorted.bam O=SRR1282206.dedup.bam ...

# 4. Extract methylation calls (allc format)
# Using MethylDackel or allcools
MethylDackel extract --CHG --CHH hg38.fa SRR1282206.dedup.bam
# Then convert to allc format: chr, pos, strand, site, mc, cov
```

## Overview

### Layer 1 — Null Simulation (False Positive Assessment)
- Takes PrEC replicates and creates two pseudo-groups with **identical biology**
- Uses **clone mode**: all 3 reps are copied to both groups (3 vs 3 design),
  each copy distorted with a different simulated enzyme efficiency
- Every DMR called is a **false positive** — compare FP counts across methods

### Layer 2 — Spike-in Simulation (Sensitivity Assessment)
- Same clone-mode design, but injects known DMRs into PseudoB copies
  **before** efficiency scaling
- Classifies called DMRs as TP (overlaps spike-in), FP, or FN
- Reports sensitivity, precision, FDR, F1 across methods and effect sizes

## File Structure

```
simulation/
├── R/
│   ├── 01_simulate_efficiency.R   # Core: efficiency distortion + spike-in injection
│   ├── 02_run_methods.R           # Runners: raw, downsampled, MAPitNorm, ComBatMet
│   ├── 03_evaluate.R              # DMR classification (TP/FP/FN) + metrics
│   ├── 04_run_simulation.R        # Main orchestration script
│   └── 05_plot_results.R          # Manuscript-quality figures
├── data/                          # Input data (not tracked)
│   ├── allc/                      # allc methylation files
│   └── sample_sheet.csv           # Sample metadata
└── results/                       # Output (not tracked)
    └── simulation/
        ├── null/                  # Null simulation outputs
        ├── spikein/               # Spike-in simulation outputs
        ├── null_simulation_results.csv
        ├── spikein_simulation_results.csv
        └── all_results.rds
```

## Prerequisites

### R packages
```r
# Core
install.packages(c("data.table", "ggplot2"))

# Bioconductor
BiocManager::install(c("sva", "GenomicRanges", "IRanges", "S4Vectors"))

# MAPitNorm (your package)
devtools::install("path/to/MAPitNorm")
```

### External tools
- **metilene** (DMR caller): https://www.bioinf.uni-leipzig.de/Software/metilene/
- **bwameth** or **bismark** (bisulfite alignment)
- **MethylDackel** or **allcools** (methylation calling)

## Quick Start

### Interactive (in R)
```r
setwd("path/to/simulation")

source("R/01_simulate_efficiency.R")
source("R/02_run_methods.R")
source("R/03_evaluate.R")
source("R/04_run_simulation.R")

# Load PrEC data (3 replicates)
config <- parse_args()
config$data_dir <- "data/allc"
config$sample_sheet <- "data/sample_sheet.csv"
config$wt_group_id <- "PrEC"

wt_reps <- prepare_wt_replicates(config)

# Quick test with one scenario (clone mode: 3 vs 3)
pseudo <- create_pseudo_groups(
  replicates = wt_reps,
  efficiency_A = c(0.95, 1.0),
  efficiency_B = c(0.70, 0.80),
  mode = "clone"              # all 3 reps → both groups
)

# Run MAPitNorm only
result <- run_mapitnorm(pseudo)
```

### Command line (full pipeline)
```bash
Rscript R/04_run_simulation.R \
  --data_dir data/allc \
  --sample_sheet data/sample_sheet.csv \
  --output_dir results/simulation \
  --wt_group_id PrEC \
  --scenarios mild,moderate,severe \
  --seed 42
```

### Generate figures
```r
source("R/05_plot_results.R")
generate_all_figures("results/simulation/all_results.rds", "results/figures")
```

## Simulation Parameters

### Pseudo-group construction modes
| Mode    | Design     | When to use                                   |
|---------|------------|-----------------------------------------------|
| `clone` | N vs N     | Default. All reps in both groups. Best power.  |
| `split` | ⌈N/2⌉ vs ⌊N/2⌋ | When replicate independence matters.    |

### Efficiency scenarios
| Scenario | Group A (reference) | Group B (distorted) |
|----------|--------------------|--------------------|
| Mild     | 0.95–1.00          | 0.85–0.90          |
| Moderate | 0.95–1.00          | 0.70–0.80          |
| Severe   | 0.95–1.00          | 0.50–0.60          |

### Spike-in effect sizes
0.05, 0.10, 0.15, 0.20, 0.30 (absolute rate difference)

### Methods compared
1. **Raw** — no normalization
2. **Downsampled** — hypergeometric downsampling to minimum coverage
3. **MAPitNorm** — quantile-bin shrinkage (this package)
4. **ComBatMet** — ComBat on logit-transformed rates

## Key Design Decisions

- **Clone mode for 3-replicate data**: With only 3 reps, splitting gives 2 vs 1
  which is underpowered for metilene. Clone mode gives 3 vs 3 while preserving
  the null (both groups have the same underlying biology).
- **Coverage stays constant during efficiency simulation**: enzyme efficiency affects
  methylation labeling, not sequencing depth. Only `mc` is resampled.
- **Spike-ins injected BEFORE efficiency scaling**: this models real biology where
  true differences exist prior to the assay.
- **Same spike-in regions across all scenarios**: ensures comparability.
- **metilene for DMR calling**: well-established, fast, method-agnostic.
