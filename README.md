# NOMe-Simulator-9000

Benchmarking framework for evaluating normalization methods for exogenous
methyltransferase data (MAPit, NOMe-seq). Simulates enzyme efficiency artifacts
and spike-in DMRs to assess false positive rates, sensitivity, and precision
across normalization approaches.

Developed as a companion analysis project for [MAPitNorm](https://github.com/jobrant/MAPitNorm).

## Motivation

Exogenous methyltransferase assays (MAPit, NOMe-seq) probe chromatin
accessibility by treating genomic DNA with GpC methyltransferase (M.CviPI).
Variation in enzyme efficiency between samples introduces systematic artifacts
that can be mistaken for biological differences. Existing normalization
approaches, including downsampling and ComBat-based batch correction, are
not designed for this problem and can either amplify false positives or
eliminate genuine biological signal.

This framework provides a rigorous, simulation-based evaluation of how well
different normalization methods separate technical artifacts from real biology.

## Simulation Design

### Layer 1 — Null Simulation (False Positive Assessment)

Takes biological replicates from a single condition and creates two
pseudo-groups with **identical biology** using clone mode: all replicates are
copied to both groups (N vs N design), with each copy distorted by a different
simulated enzyme efficiency via binomial resampling of methylated counts.
Coverage is unchanged (enzyme efficiency affects methylation labeling, not
sequencing depth). Every DMR called under this design is a false positive.

### Layer 2 — Spike-in Simulation (Sensitivity Assessment)

Same clone-mode design, but known DMRs are injected into one pseudo-group
**before** efficiency scaling — modeling real biological differences that exist
prior to the assay. Called DMRs are classified as true positives (overlap a
spike-in region), false positives (no overlap), or false negatives (spike-in
regions missed). Reports sensitivity, precision, FDR, and F1 across methods
and effect sizes.

### Methods Compared

| Method | Approach | Implementation |
|--------|----------|----------------|
| Raw | No normalization (baseline) | Direct pass-through |
| Downsampled | Hypergeometric downsampling to minimum per-site coverage | `02_run_methods.R` |
| MAPitNorm | Quantile-bin shrinkage normalization | [MAPitNorm package](https://github.com/jobrant/MAPitNorm) |
| ComBatMet | Mean-only batch correction on logit-transformed rates | Manual implementation (sva-equivalent) |

### Efficiency Scenarios

| Scenario | Group A (reference) | Group B (distorted) |
|----------|---------------------|---------------------|
| Mild     | 0.95–1.00           | 0.85–0.90           |
| Moderate | 0.95–1.00           | 0.70–0.80           |
| Severe   | 0.95–1.00           | 0.50–0.60           |

### Spike-in Effect Sizes

0.05, 0.10, 0.15, 0.20, 0.30 (absolute rate difference)

## Test Dataset

NOMe-seq data from Taberlay/Statham et al. (PMID 24916973, GEO: GSE57498).

| Cell line | Type              | SRA accessions          | Reps |
|-----------|-------------------|-------------------------|------|
| HMEC      | Normal breast     | SRR1282202–SRR1282203   | 2    |
| MCF7      | Breast cancer     | SRR1282204–SRR1282205   | 2    |
| **PrEC**  | Normal prostate   | SRR1282206–SRR1282208   | **3** |
| **PC3**   | Prostate cancer   | SRR1282209–SRR1282211   | **3** |

PrEC (3 replicates) serves as the source for null and spike-in simulations.
PC3 vs PrEC provides a real biological comparison for validation.

### Data Preprocessing

Raw FASTQs are processed with the [nf-core/methylseq](https://nf-co.re/methylseq)
pipeline (Bismark aligner) followed by methylation calling to generate allc-format
files (e.g. methylpy or allcools). Pre-filtering for standard chromosomes and minimum coverage (≥5x) is
recommended before running the simulation:

```bash
for f in data/allc/PrEC*_GCH.tsv.gz; do
    out="${f%.tsv.gz}_cov5_std.tsv.gz"
    zcat "$f" | awk '$6 >= 5 && $1 ~ /^[0-9]+$|^[XYM]$/' | gzip > "$out"
done
```

## File Structure

```
nome-simulator-9000/
├── R/
│   ├── 01_simulate_efficiency.R   # Efficiency distortion + spike-in DMR injection
│   ├── 02_run_methods.R           # Method runners + metilene DMR calling
│   ├── 03_evaluate.R              # DMR classification (TP/FP/FN) + metrics
│   ├── 04_run_simulation.R        # Main orchestration script (CLI + interactive)
│   ├── 05_plot_results.R          # Primary manuscript figures
│   ├── 06_additional_figures.R    # Precision-sensitivity panels, call breakdown
│   └── uid_functions.R            # Fast numeric site identifiers
├── data/                          # Input data (not tracked)
│   ├── allc/                      # allc methylation files
│   └── sample_sheet.csv           # Sample metadata
├── results/                       # Output (not tracked)
│   └── simulation/
│       ├── null/
│       ├── spikein/
│       ├── null_simulation_results.csv
│       ├── spikein_simulation_results.csv
│       └── all_results.rds
└── README.md
```

## Prerequisites

### R packages

```r
install.packages(c("data.table", "ggplot2"))
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"))
devtools::install_github("jobrant/MAPitNorm")
```

### External tools

- [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/) — DMR caller

## Usage

### Interactive (R)

```r
source("R/01_simulate_efficiency.R")
source("R/02_run_methods.R")
source("R/03_evaluate.R")
source("R/04_run_simulation.R")

# Configure
config <- parse_args()
config$data_dir <- "data/allc"
config$sample_sheet <- "data/sample_sheet.csv"
config$wt_group_id <- "PrEC"
config$output_dir <- "results/simulation"
config$region_width_bp <- 10000
config$n_spikein_regions <- 300
config$metilene_min_cpg <- 3
config$seed <- 42

# Load and prepare data
wt_reps <- prepare_wt_replicates(config)

# Run simulations
null_results <- run_null_simulation(wt_reps, config)
spikein_results <- run_spikein_simulation(wt_reps, config)

# Generate figures
source("R/05_plot_results.R")
source("R/06_additional_figures.R")
generate_all_figures("results/simulation/all_results.rds", "results/figures")
generate_additional_figures(
  "results/simulation/spikein_simulation_results.csv",
  "results/figures"
)
```

### Command line

```bash
Rscript R/04_run_simulation.R \
  --data_dir data/allc \
  --sample_sheet data/sample_sheet.csv \
  --output_dir results/simulation \
  --wt_group_id PrEC \
  --scenarios mild,moderate,severe \
  --region_width_bp 10000 \
  --n_spikein_regions 300 \
  --metilene_min_cpg 3 \
  --seed 42
```

## Key Design Decisions

- **Clone mode for pseudo-group construction**: With ≤3 replicates, splitting
  gives an underpowered design. Clone mode copies all replicates to both groups
  (N vs N), maximizing statistical power while maintaining the null (identical
  biology, different simulated efficiency).

- **Coverage unchanged during efficiency simulation**: Enzyme efficiency affects
  whether accessible GpC sites are methylated, not whether the DNA is sequenced.
  Only methylated counts (mc) are resampled; coverage stays constant.

- **Spike-ins injected before efficiency scaling**: Models real biology where
  true differences exist prior to the assay.

- **Numeric site identifiers**: Sites are encoded as `chr_int * 1e10 + pos * 10
  + site_int` instead of string concatenation, providing ~20x faster loading
  and significantly reduced memory usage for genome-wide data.

- **Manual ComBat implementation**: The sva::ComBat function has compatibility
  issues with current R versions. The manual mean-only implementation is
  equivalent and more transparent. Note that ComBat is given an idealized setup 
  in this benchmark — the enzyme efficiency difference aligns perfectly with the 
  batch structure. In real experiments, efficiency variation is per-sample and 
  can occur within batches, making it invisible to batch-correction approaches 
  that rely on known group labels.

## Sample Sheet Format

```csv
group_id,replicate,sample_id,file_name
PrEC,1,PrECA,PrECA_cov5_std_GCH.tsv.gz
PrEC,2,PrECB,PrECB_cov5_std_GCH.tsv.gz
PrEC,3,PrECC,PrECC_cov5_std_GCH.tsv.gz
```

## Figures Generated

| Figure | Description |
|--------|-------------|
| `fig1_null_fp` | False positive DMR counts under null (bar chart) |
| `fig2_sensitivity` | Sensitivity curves by effect size, faceted by scenario |
| `fig3_f1_heatmap` | F1 score heatmap across methods and conditions |
| `fig4_precision_recall` | Precision vs sensitivity scatter |
| `fig5_fdr` | FDR curves by effect size |
| `fig6_sensitivity_precision_panels` | Side-by-side sensitivity and precision panels |
| `fig7_call_breakdown` | Stacked TP/FP composition of DMR calls |
| `fig8_researcher_view` | Lollipop chart showing call purity at a single effect size |
| `fig9_f1_lines` | F1 score line plot (alternative to heatmap) |

## License

MIT