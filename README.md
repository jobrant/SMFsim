# SMFsim: Single Molecule Footprinting Simulation Framework

SMFsim is a benchmarking simulation package for exogenous methyltransferase
assays (MAPit, NOMe-seq, Fiber-seq). It generates controlled null and spike-in
experiments to evaluate how normalization methods separate technical enzyme
efficiency artifacts from real biological signal.

It is designed as a companion package for [SMFnorm](https://github.com/jobrant/SMFnorm)
and the associated manuscript, but it can also be used as a general simulation
framework for benchmarking SMF-style normalization methods.

## Core Design

SMFsim produces two complementary simulation layers:

- **Null simulation** — identical biology in both pseudo-groups, with only
  simulated per-replicate enzyme efficiency variation. Every called DMR in this
  design is a false positive.
- **Spike-in simulation** — a known set of DMRs is injected into one pseudo-group
  before efficiency scaling, creating a true biological signal that later
  undergoes the same efficiency noise model.

The package supports multiple normalization methods and explicit SMFnorm
parameter tuning for the two meaningful correction regimes:

- `rate_between_groups = FALSE` for within-group variance scenarios
- `rate_between_groups = TRUE` for systematic between-group bias scenarios

## Noise and Sampling Model

### Parametric pseudo-group construction

The main simulation model is `sim_mode = "parametric"`. In this mode, pseudo-
replicates are generated from a pooled per-site true methylation rate rather
than by cloning raw counts.

For each site:

1. A pooled true rate `p_i = sum(mc) / sum(cov)` is computed across source
   replicates.
2. Each pseudo-replicate reuses a real per-site coverage `cov_ij` from the
   source data.
3. The effective mean rate is `mu_ij = p_i * efficiency_j`.
4. If `dispersion_s` is finite, a latent rate is drawn from
   `Beta(mu_ij * s, (1 - mu_ij) * s)` and then counts are drawn from
   `Binomial(cov_ij, latent_rate)`.
5. If `dispersion_s = Inf`, the model reduces to a pure binomial draw.

This parametric design ensures the null is properly calibrated because the two
groups share no realized noise from the source replicates.

### Overdispersion

- `dispersion_s = Inf` uses a binomial model.
- Finite `dispersion_s` injects beta-binomial overdispersion.

In the current code, the `run_alpha_sweep.R` script uses `dispersion_s = 26` for
between-group bias tuning, while the within-group tuning script
(`run_within_alpha_sweep.R`) and the bias comparison script
(`run_method_comparison_bias.R`) currently use `dispersion_s = Inf`.

## Simulation Scenarios

### Within-group scenarios

These scenarios model per-replicate efficiency variation within each pseudo-
group, with no systematic mean bias between groups.

- `mild`
- `moderate`
- `severe`

These correspond to increasingly wide efficiency ranges and replicate-to-
replicate variability.

### Between-group bias scenarios

These scenarios introduce systematic differences in the pseudo-group means.
Current code uses:

- `aligned_moderate`
- `aligned_strong`
- `imbalanced_moderate`
- `imbalanced_strong`

These scenario labels are used by `run_alpha_sweep.R` and
`run_method_comparison_bias.R` to evaluate how SMFnorm performs when a
true between-group efficiency artifact is present.

## Methods Compared

SMFsim currently compares these methods:

- `raw` — no normalization
- `downsampled` — coverage normalization by downsampling to the minimum
  coverage across samples at each site
- `SMFnorm` — the SMFnorm normalization pipeline with `within_alpha` and
  `between_alpha` control
- `ComBatMet` — manual mean-only batch correction on logit-transformed rates

## Installation

```r
install.packages(c("data.table", "ggplot2"))
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"))
devtools::install_github("jobrant/SMFnorm")
devtools::install_github("jobrant/SMFsim")
```

### Required external tools

- `metilene` — DMR caller (usually v0.2.8 in current scripts)

## Quick Start

```r
library(SMFsim)

config <- parse_args()
config$data_dir        <- "data/allc"
config$sample_sheet    <- "data/sample_sheet.csv"
config$wt_group_id     <- "PrEC"
config$output_dir      <- "results/simulation"
config$n_spikein_regions <- 300
config$region_width_bp <- 2000
config$metilene_min_cpg <- 10
config$seed            <- 42

wt_reps <- prepare_wt_replicates(config)
null_results    <- run_null_simulation(wt_reps, config)
spikein_results <- run_spikein_simulation(wt_reps, config)
```

## Recommended default config

For the current manuscript-style benchmarking, the following defaults are a good starting point:

```r
config <- parse_args()
config$data_dir          <- "data/allc"
config$sample_sheet      <- "data/sample_sheet.csv"
config$wt_group_id       <- "PrEC"
config$output_dir        <- "results/simulation"
config$n_spikein_regions <- 300
config$region_width_bp   <- 2000
config$metilene_min_cpg  <- 10
config$metilene_qval     <- 0.05
config$seed              <- 42
config$sim_mode          <- "parametric"
config$dispersion_s      <- Inf
config$rate_between_groups <- FALSE
config$within_alpha      <- 0.3
config$between_alpha     <- 0.8
config$methods           <- c("raw", "downsampled", "SMFnorm", "ComBatMet")
```

Adjust `dispersion_s`, `within_alpha`, `between_alpha`, and `rate_between_groups`
for the specific tuning run you are performing.

## Example Scripts

### Within-alpha tuning

```bash
Rscript inst/scripts/run_within_alpha_sweep.R
```

This script sweeps `within_alpha` across the within-group scenarios
(`mild`, `moderate`, `severe`) with `rate_between_groups = FALSE`.
It writes results under `results/within_alpha_sweep/`.

### Between-alpha tuning

```bash
Rscript inst/scripts/run_alpha_sweep.R
```

This script sweeps `between_alpha` across the biased scenarios with
`rate_between_groups = TRUE` and is used to tune SMFnorm when a systematic
between-group artifact is present.

### Bias comparison run

```bash
Rscript inst/scripts/run_method_comparison_bias.R
```

This script runs the full method comparison on the current between-group bias
scenarios and saves both results and figures.

### Real data comparison

```r
library(SMFsim)
config <- parse_args()
config$data_dir     <- "data/allc"
config$sample_sheet <- "data/sample_sheet.csv"
config$output_dir   <- "results/real_data"
config$metilene_min_cpg <- 5
run_real_data_comparison(config, group_A = "PrEC", group_B = "PC3")
```

## Data and Input Format

The package assumes a sample sheet with columns:

```csv
group_id,replicate,sample_id,file_name
PrEC,1,PrECA,PrECA_cov5_std_GCH.tsv.gz
PrEC,2,PrECB,PrECB_cov5_std_GCH.tsv.gz
PrEC,3,PrECC,PrECC_cov5_std_GCH.tsv.gz
```

Input files are expected in allc-like TSV format with columns:
`chr`, `pos`, `strand`, `site`, `mc`, `cov`.

## File Structure

```
SMFsim/
├── R/
│   ├── SMFsim-package.R
│   ├── 01_simulate_efficiency.R
│   ├── 02_run_methods.R
│   ├── 03_evaluate.R
│   ├── 04_run_simulation.R
│   ├── 05_plot_results.R
│   ├── 06_additional_figures.R
│   └── 07_real_data_comparison.R
├── inst/
│   └── scripts/
│       ├── simulation_testing.R
│       ├── run_alpha_sweep.R
│       ├── run_within_alpha_sweep.R
│       └── run_method_comparison_bias.R
├── data/
├── results/
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

## Key Concepts

- `sim_mode = "parametric"`: the default recommended mode for null and
  spike-in simulations.
- `dispersion_s`: finite values inject beta-binomial overdispersion; `Inf`
  reduces to simple binomial sampling.
- `rate_between_groups`: `TRUE` corrects systematic between-group artifacts;
  `FALSE` corrects only within-group rate variation.
- `within_alpha` and `between_alpha`: SMFnorm tuning parameters that control
  shrinkage strength for within-group and between-group rate normalization.

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.
