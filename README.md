# Anonymous implementation for submission 21076


This repository contains R code for computing distance-profile centrality
scores and reproducing the simulation studies, methodological illustrations,
and real-data applications reported in the accompanying manuscript.

## Repository layout

```text
.
├── README.md
├── reproduce_all.R                 # main simulation/figure workflow
├── R/                              # method, generators, simulation, plotting
├── config/                         # simulation profiles
├── scripts/                        # staged simulation and plotting entry points
├── tests/                          # reduced end-to-end test
├── figures/manuscript/             # reference manuscript figures
├── results/manuscript/             # clean summaries used by the plot stage
├── score_comparison/               # one-dimensional five-method comparison
│   └── README.md
└── dataIllustrations/              # four real-data applications
    ├── README.md
    ├── run_all.R
    ├── data/
    └── results/                    # reference numerical and figure outputs
```

The root simulation workflow and `score_comparison/` share the authoritative
implementation in `R/method.R`. The real-data applications use the separate
object-type implementation in `dataIllustrations/tolReg.R`.

## Software requirements

The workflows were tested with R 4.5.1 and these direct package versions:

| Package | Tested version | Used by |
|---|---:|---|
| `ggplot2` | 4.0.1 | simulations and real data |
| `cowplot` | 1.2.0 | simulation figure layouts |
| `MASS` | 7.3-65 | multivariate normal generation |
| `KernSmooth` | 2.23-26 | distance-profile smoothing |
| `scales` | 1.4.0 | simulation plot labels |
| `patchwork` | 1.3.2 | real-data figure layouts |
| `pracma` | 2.4.4 | colon transport-map interpolation |
| `dbscan` | 1.2.5 | taxi-network clustering |

Install the direct packages from an R session with:

```r
install.packages(c(
  "ggplot2", "cowplot", "MASS", "KernSmooth", "scales",
  "patchwork", "pracma", "dbscan"
))
```

Base or recommended packages used by the code include `parallel`, `stats`,
`utils`, `graphics`, `grDevices`, and `grid`. No package lockfile is supplied;
the table records the tested environment rather than promising byte-identical
package resolution.

## Quick start

Run commands in this README from the repository root. The reduced test checks
the method identities and exercises all three main simulation families,
summary generation, and plotting:

```sh
Rscript --vanilla tests/smoke_test.R
```

It writes ignored diagnostic output to `results/smoke/` and
`figures/smoke/`. A successful run ends with
`All smoke tests passed.` On a recent laptop it completes in under one minute
and remains below 1 GB of memory.

## Reproducibility

The root configuration is `config/simulation_settings.R`.

- Figures 1 and 2 use seed 2024.
- Each Monte Carlo replication resets its seed to the replication number.
- The RNG uses Mersenne-Twister, Inversion normal generation, and Rejection
  sampling.
- Results therefore do not depend on worker scheduling.
- The `manuscript` profile contains the full settings; `smoke` keeps all
  scenario types with smaller samples, grids, and replication counts.

Cached replications are reused only when their scientific settings and RNG
configuration match. Add `--force` to a simulation command to recompute them.

## Simulation studies

### Main simulation and figure workflow

The root workflow covers Euclidean Gaussian mixtures, tangent-noise objects on
the sphere, and normal distributions in univariate 2-Wasserstein space.

| Entry point | Purpose |
|---|---|
| `scripts/01_reproduce_figures_1_2.R` | methodological illustrations |
| `scripts/02_run_figures_3_4_simulations.R` | Monte Carlo simulation and clean summaries |
| `scripts/03_plot_figures_3_4.R` | plots from saved clean summaries |
| `reproduce_all.R` | all three stages in order |

A reduced standalone pipeline is:

```sh
Rscript --vanilla reproduce_all.R --profile=smoke --workers=2 --force
```

The saved full-profile summaries bundled in `results/manuscript/` can be
replotted without rerunning the Monte Carlo simulation:

```sh
Rscript --vanilla scripts/03_plot_figures_3_4.R --profile=manuscript
```

The full workflow is:

```sh
Rscript --vanilla reproduce_all.R --profile=manuscript --workers=4
```

The full profile uses 200 replications for each simulation family, separate
training and calibration samples of sizes 200, 400, and 800, and 2,000 test
objects per setting. Runtime is hardware- and worker-dependent and is expected
to be measured in hours rather than minutes; the largest runs require several
GB of available memory per concurrent worker. Replication caches below
`results/manuscript/cache/` allow interrupted runs to resume and are excluded
from version control.

Numerical summaries are written to `results/<profile>/`; PDFs are written to
`figures/<profile>/`. The full-profile reference summaries and figures are
deliberately included, while replicate caches and smoke outputs are not.

### One-dimensional score comparison

`score_comparison/` compares the proposed weighted C-score, transport rank,
profile transport cost, metric spatial depth, and metric halfspace depth. Its
pilot configuration runs the full population/simulation/aggregation/table
pipeline at reduced replication counts:

```sh
Rscript --vanilla score_comparison/run_all.R \
  --config=config/pilot_config.R \
  --output-dir=results/pilot \
  --workers=4
```

The pilot normally completes in under one minute. Full baseline and
independent-seed configurations use 200 replications in each of six scenarios
and produce replication-level results, Monte Carlo summaries, geometry
diagnostics, validation checks, and Markdown tables. See
[`score_comparison/README.md`](score_comparison/README.md) for definitions,
configuration details, stage commands, and output files.

## Real-data illustrations

`dataIllustrations/` contains applications to colon gene expression, Tampa
minimum temperatures, MNIST handwriting, and Manhattan taxi transport
networks. The temperature source CSV and the processed taxi-network matrices
are included. The colon and MNIST source files must be downloaded separately
because their source pages do not state explicit redistribution terms.

Run the two applications whose inputs are bundled with the repository:

```sh
Rscript --vanilla dataIllustrations/run_all.R \
  --analyses=temperature,taxi
```

With all documented inputs in place, omit `--analyses` to run all four
applications. Each script writes numerical CSV files and a PDF or PNG below
`dataIllustrations/results/<analysis>/`. Reference outputs are included for
comparison. See
[`dataIllustrations/README.md`](dataIllustrations/README.md) for exact input
paths, data access, citations, preprocessing, individual commands, expected
outputs, and limitations.

## Manuscript-output mapping

| Manuscript item | Output |
|---|---|
| Figure `fig:dp-1` | `figures/manuscript/Fig2_new_DP.pdf` |
| Figure `fig:dp-3` | `figures/manuscript/illu_2d_2plots.pdf` |
| Figure `fig:sim-2d` | `figures/manuscript/mix_2d_1.pdf`, `mix_2d_heatmap_2x3.pdf` |
| Figure `fig:sim-S2W2` | `figures/manuscript/S2_W2_lines_2x2.pdf`, `S2W2_heat_2_by_1.pdf` |
| Figure `fig:colon` | `dataIllustrations/results/colon/colon.pdf` |
| Figure `fig:temp` | `dataIllustrations/results/temperature/temperature.pdf` |
| Figure `fig:mnist` | `dataIllustrations/results/mnist/mnist.png` |
| Figure `fig:taxi` | `dataIllustrations/results/taxi/taxi_profile_mds.pdf` |
| Table `tab:taxi` | `dataIllustrations/results/taxi/taxi_dbscan_clusters.csv` |

## Data and licensing notes

Third-party data remain subject to their source terms and citation
requirements. NOAA/NCEI data are U.S. federal data; the bundled temperature
subset is attributed in the real-data README. The processed taxi object is
derived from publicly downloadable NYC TLC records and contains aggregate
graph Laplacians rather than trip-level records; the source terms are linked in
the real-data README. Colon and MNIST raw files are not distributed.

No repository-wide software license has been selected. The absence of a
license file should not be interpreted as permission beyond applicable law or
the separate terms attached to third-party data.

## Troubleshooting

- If a package check fails, install the missing package listed in the error and
  restart the command in a fresh R session.
- On Windows, `parallel::mclapply()` falls back to sequential execution, so
  `--workers` does not provide multicore speedup.
- If a real-data command reports a missing input, reproduce the exact directory
  tree in `dataIllustrations/README.md`; application paths are resolved from
  the script location, not the current working directory.
- If a saved simulation summary is rejected, its configuration does not match
  the requested profile. Recreate it with the matching profile rather than
  bypassing the validation.
