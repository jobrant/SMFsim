# SMFsim: Single Molecule Footprinting Simulation Framework

Benchmarking framework for evaluating normalization methods for exogenous
methyltransferase data (MAPit, NOMe-seq, Fiber-seq). Simulates enzyme efficiency
artifacts and spike-in DMRs to assess false positive rates, sensitivity, and
precision across normalization approaches.

Developed as a companion analysis package for
[SMFnorm](https://github.com/jobrant/SMFnorm).

## Motivation

Exogenous methyltransferase assays (MAPit, NOMe-seq, Fiber-seq) probe chromatin
accessibility by treating genomic DNA with GpC methyltransferase (M.CviPI).
Variation in enzyme efficiency between replicates introduces systematic artifacts
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
simulated per-replicate enzyme efficiency via binomial resampling of methylated
counts. Coverage is unchanged (enzyme efficiency affects methylation labeling,
not sequencing depth). Every DMR called under this design is a false positive.

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
| SMFnorm | Quantile-bin shrinkage normalization | [SMFnorm package](https://github.com/jobrant/SMFnorm) |
| ComBatMet | Mean-only batch correction on logit-transformed rates | Manual implementation (sva-equivalent) |

### Efficiency Scenarios

Efficiency variation is modeled as **within-group, per-replicate** variation —
reflecting real experiments where replicates processed on different days may
have different enzyme efficiencies. Both pseudo-groups span overlapping
efficiency ranges with no systematic between-group bias.

| Scenario | Group A (per replicate) | Group B (per replicate) | Within-group SD |
|----------|-------------------------|-------------------------|-----------------|
| Mild     | 0.95, 0.88, 0.87        | 0.86, 0.93, 0.90        | ≈ 0.04          |
| Moderate | 0.95, 0.70, 0.85        | 0.90, 0.65, 0.95        | ≈ 0.10          |
| Severe   | 0.55, 0.95, 0.75        | 0.90, 0.50, 0.80        | ≈ 0.17          |

### Spike-in Effect Sizes

0.10, 0.15, 0.20, 0.30 (absolute rate difference)

## Test Dataset

NOMe-seq data from Taberlay/Statham et al. (PMID 24916973, GEO: GSE57498).

| Cell line | Type              | SRA accessions          | Reps used |
|-----------|-------------------|-------------------------|-----------|
| **PrEC**  | Normal prostate   | SRR1282206–SRR1282208   | 3         |
| **PC3**   | Prostate cancer   | SRR1282209–SRR1282211   | 2*        |

\* PC3 replicate PCC was excluded due to low coverage (~2.0M sites at ≥5x vs
7.9–8.9M for other replicates). Exclusion of PCC reduced shared sites from
78K to 24K across all six samples.

PrEC (3 replicates) serves as the source for null and spike-in simulations.
PC3 vs PrEC (2 vs 3 replicates, ~78K shared sites) provides a real biological
comparison for validation.

### Data Preprocessing

Raw FASTQs are processed with the [nf-core/methylseq](https://nf-co.re/methylseq)
pipeline (Bismark aligner) followed by methylation calling to generate
allc-format files (e.g., methylpy or allcools). Pre-filtering for standard
chromosomes and minimum coverage (≥5x) is recommended before running the
simulation:

```bash
for f in data/allc/PrEC*_GCH.tsv.gz; do
    out="${f%.tsv.gz}_cov5_std.tsv.gz"
    zcat "$f" | awk '$6 >= 5 && $1 ~ /^[0-9]+$|^[XYM]$/' | gzip > "$out"
done
```

## Installation

```r
# Install dependencies
install.packages(c("data.table", "ggplot2"))
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"))
devtools::install_github("jobrant/SMFnorm")

To compare:

devtools::install_github("JmWangBio/ComBatMet")

# Install SMFsim
devtools::install_github("jobrant/SMFsim")
```

### External tools

- [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/) — DMR caller (v0.2.8)

## File Structure

```
SMFsim/
├── R/
│   ├── SMFsim-package.R           # Package-level imports and documentation
│   ├── 01_simulate_efficiency.R   # Efficiency distortion + spike-in DMR injection
│   ├── 02_run_methods.R           # Method runners + metilene DMR calling
│   ├── 03_evaluate.R              # DMR classification (TP/FP/FN) + metrics
│   ├── 04_run_simulation.R        # Main orchestration functions
│   ├── 05_plot_results.R          # Primary manuscript figures
│   ├── 06_additional_figures.R    # Precision-sensitivity panels, call breakdown
│   └── 07_real_data_comparison.R  # PC3 vs PrEC real data validation
├── inst/
│   └── scripts/
│       └── simulation_testing.R   # Example script for running the simulation
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
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

## Usage

After installing the package, all functions are available via `library(SMFsim)`.
A complete worked example is provided in `inst/scripts/simulation_testing.R`.
The core workflow is:

```r
library(SMFsim)

# Configure
config <- parse_args()
config$data_dir        <- "data/allc"
config$sample_sheet    <- "data/sample_sheet.csv"
config$wt_group_id     <- "PrEC"
config$output_dir      <- "results/simulation"
config$n_spikein_regions <- 300
config$region_width_bp <- 2000
config$metilene_min_cpg <- 10
config$seed            <- 42

# Load and prepare data
wt_reps <- prepare_wt_replicates(config)

# Run simulations
null_results    <- run_null_simulation(wt_reps, config)
spikein_results <- run_spikein_simulation(wt_reps, config)

# Generate figures
generate_all_figures(
  "results/simulation/all_results.rds",
  "results/figures"
)
generate_additional_figures(
  "results/simulation/spikein_simulation_results.csv",
  "results/figures"
)
```

### Real data comparison (PC3 vs PrEC)

```r
library(SMFsim)

config <- parse_args()
config$data_dir     <- "data/allc"
config$sample_sheet <- "data/sample_sheet.csv"
config$output_dir   <- "results/real_data"
config$metilene_min_cpg <- 5   # use lower threshold for sparser real data

run_real_data_comparison(config, group_A = "PrEC", group_B = "PC3")
```

## Key Design Decisions

- **Within-group variance model**: Efficiency variation is modeled per-replicate
  within each group, not as a systematic between-group difference. This reflects
  real experimental conditions and tests the normalization mode users will
  actually apply (`rate_between_groups = FALSE` in SMFnorm).

- **Clone mode for pseudo-group construction**: With ≤3 replicates, splitting
  gives an underpowered design. Clone mode copies all replicates to both groups
  (N vs N), maximizing statistical power while maintaining the null (identical
  biology, different simulated efficiency per replicate).

- **Coverage unchanged during efficiency simulation**: Enzyme efficiency affects
  whether accessible GpC sites are methylated, not whether the DNA is sequenced.
  Only methylated counts (mc) are resampled; coverage stays constant.

- **Spike-ins injected before efficiency scaling**: Models real biology where
  true differences exist prior to the assay.

- **Dense region selection for spike-in placement**: Random genomic windows
  have too few GCH sites for reliable DMR calling (median ~875 bp spacing
  genome-wide). Spike-in regions are selected by a two-pointer sliding window
  algorithm that finds naturally dense clusters (median ~15 sites, 21–35 bp
  spacing).

- **Numeric site identifiers**: Sites are encoded as `chr_int * 1e10 + pos * 10
  + site_int` instead of string concatenation, providing ~20x faster loading
  and significantly reduced memory usage for genome-wide data.

- **Manual ComBatMet implementation**: The `sva::ComBat` function has
  compatibility issues with current R versions (matrix dimension dropping in
  `log()`/`pmin()`/`pmax()`, batch/group confounding). The manual mean-only
  implementation is equivalent and more transparent. Note that ComBatMet is
  given an idealized setup in this benchmark — the efficiency structure aligns
  perfectly with the declared batch labels. In real experiments, efficiency
  variation is per-sample and can occur within batches, making it invisible
  to batch-correction approaches that rely on known group labels.

## Sample Sheet Format

```csv
group_id,replicate,sample_id,file_name
PrEC,1,PrECA,PrECA_cov5_std_GCH.tsv.gz
PrEC,2,PrECB,PrECB_cov5_std_GCH.tsv.gz
PrEC,3,PrECC,PrECC_cov5_std_GCH.tsv.gz
```

## Figures Generated

### Simulation figures

| Figure | Function | Description |
|--------|----------|-------------|
| `fig1_null_fp` | `plot_null_fp()` | False positive DMR counts under null (bar chart) |
| `fig2_sensitivity` | `plot_sensitivity_curves()` | Sensitivity curves by effect size, faceted by scenario |
| `fig3_f1_heatmap` | `plot_f1_heatmap()` | F1 score heatmap across methods and conditions |
| `fig4_precision_recall` | `plot_precision_recall()` | Precision vs sensitivity scatter |
| `fig5_fdr` | `plot_fdr()` | FDR curves by effect size |
| `fig6_sensitivity_precision_panels` | `plot_sensitivity_precision_panels()` | Side-by-side sensitivity and precision panels |
| `fig7_call_breakdown` | `plot_call_breakdown()` | Stacked TP/FP composition of DMR calls |
| `fig8_researcher_view` | `plot_researcher_view()` | Lollipop chart showing call purity at a single effect size |
| `fig9_f1_lines` | `plot_f1_lines()` | F1 score line plot (alternative to heatmap) |

### Real data figures

| Figure | Function | Description |
|--------|----------|-------------|
| `fig_real_dmr_counts` | `plot_dmr_counts()` | DMR counts per method |
| `fig_real_effect_sizes` | `plot_effect_size_distribution()` | Effect size distributions per method |
| `fig_real_overlap_heatmap` | `plot_overlap_heatmap()` | Pairwise DMR overlap between methods |
| `fig_real_lost_smfnorm` | `plot_lost_dmr_effects()` | Density: effect sizes of DMRs lost vs retained by SMFnorm |
| `fig_real_lost_boxplot_smfnorm` | `plot_lost_dmr_effects()` | Boxplot: effect sizes of DMRs lost vs retained by SMFnorm |

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.

