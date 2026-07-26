# One-dimensional five-method score comparison

This self-contained R module compares five centrality-score rows under the
canonical Euclidean metric on the real line:

1. the proposed weighted C-score;
2. transport rank, which is the $\tau=0$ C-score ordering;
3. unconditional profile transport cost;
4. Virta's main metric spatial depth; and
5. metric halfspace (Tukey) depth.

The module generates deterministic population diagnostics, paired Monte Carlo
results, exact or interpolated prediction-set geometry, independent recalculations,
validation checks, and Markdown tables. Numerical tables are generated only
from saved CSV results.

Only base and recommended R packages are needed: `stats`, `utils`, and
`parallel`.

## Score definitions and conventions

For a candidate $\omega$ with distance-profile CDF $F_\omega$, the proposed
score is

$$
R_\tau(\omega)=\operatorname{expit}\left[
E\int_0^1\{F_X^{-1}(u)-F_\omega^{-1}(u)\}e^{-\tau u}\,du
\right].
$$

Transport rank is the same construction at $\tau=0$. The published
leave-one-out training-anchor convention and the common full-reference
implementation differ only by a candidate-independent pre-link constant; tests
verify that their orderings and order-statistic-calibrated regions agree.

The unconditional transport cost is

$$
C(\omega)=E\int_0^1
\left|F_X^{-1}(u)-F_\omega^{-1}(u)\right|\,du.
$$

Its source orientation is lower-is-more-central. The implementation uses
$-C(\omega)$ only to expose a common higher-is-more-central calibration
interface.

Virta's main metric spatial depth is

$$
D_S(\mu)=1-\frac12E\{h(X_1,X_2,\mu)\},
$$

where

$$
h(x_1,x_2,\mu)=
\mathbf{1}_{\{\mu\notin\{x_1,x_2\}\}}
\frac{
d^2(x_1,\mu)+d^2(x_2,\mu)-d^2(x_1,x_2)
}{
d(x_1,\mu)d(x_2,\mu)
}.
$$

For $d(x,y)=|x-y|$,

$$
D_S(x)=1-\{P(X>x)-P(X<x)\}^2,
$$

and the exact empirical shortcut is

$$
D_{S,n}(x)=
1-\left\{
\frac{\#(X_i>x)-\#(X_i<x)}{n}
\right\}^2.
$$

The direct implementation is the ordered-pair V-statistic. Query ties make
$h$ zero, and training-point depths use self-inclusion in the full training
reference. A leave-one-out helper is tested but is not used in production.

Metric halfspaces are $H_{a,b}=\{y:d(y,a)\leq d(y,b)\}$. Metric halfspace
depth is

$$
D_H(x)=\inf_{H\ni x}P(H).
$$

On the real line,

$$
D_H(x)=\min\{F(x),1-F(x^-)\},
$$

with exact empirical formula

$$
D_{H,n}(x)=\frac1n
\min\{\#(X_i\leq x),\#(X_i\geq x)\}.
$$

Equality is included in both closed tails. Production uses this exact
one-dimensional formula; finite anchor enumeration is only a unit-test check.

For a continuous real-valued distribution,

$$
D_S(x)=4F(x)\{1-F(x)\}
=4D_H(x)\{1-D_H(x)\}.
$$

Because $D_H\in[0,1/2]$, the two depth rows have identical population
orderings. Automated checks verify the direct reductions, continuous
rankings, finite-sample grid and test memberships, ties, and the equivalence
between nonconformity $-D$ and high-depth acceptance. The two literature
definitions are both retained but are not described as independent orderings.

## Prediction-set geometry

Every score is calibrated with the manuscript's inclusive type-1 empirical
quantile. The primary candidate grid has spacing $0.01$. For the C-score,
transport rank, and transport cost, a score crossing is linearly interpolated
between its adjacent grid points. Empirical real-line depths are step
functions, so their Lebesgue geometry is instead computed exactly from the
atoms and gaps of the ordered training sample. Total length is

$$
\lambda(\widehat C)=
\sum_{k=1}^{K}(b_k-a_k)
$$

for the disjoint interval union
$\widehat C=\bigcup_{k=1}^{K}[a_k,b_k]$.

The domain expands by 25% until at least ten consecutive excluded grid points
lie beyond every bounded set on both sides. Selected replications are repeated
at spacing $0.005$ to check length and component stability.

If a depth calibration threshold is zero, the inclusive empirical depth set is
the whole real line. The implementation stores $(-\infty,\infty)$ and reports
infinite unconditional length and HDR error, an unbounded-set rate, and
conditional finite-replication summaries. It never truncates such a set at the
candidate domain.

The HDR diagnostic is the Lebesgue symmetric-difference length

$$
\lambda(\widehat C\mathbin{\triangle}H_{0.95}),
$$

with no denominator. It is not probability-mass capture, a fraction of HDR
length, or a Jaccard index. Exact true modes and, for genuinely bimodal
mixtures, the exact antimode at zero are scored directly rather than replaced
by nearby grid points.

## Local tests

Run from this module:

```sh
Rscript --vanilla tests/test_metric_depths.R
Rscript --vanilla tests/test_r1_equivalence.R
Rscript --vanilla tests/test_set_length.R
Rscript --vanilla tests/test_exact_depth_geometry.R
Rscript --vanilla tests/smoke_test.R
```

The smoke test includes fixed score identities, scale behavior, fast-versus-
direct transport cost, calibration boundary cases, unbounded-set propagation,
manifest checks, and a reduced end-to-end pipeline.

## Simulation commands

The diagnostic pilot is:

```sh
Rscript --vanilla run_all.R \
  --config=config/pilot_config.R \
  --output-dir=results/pilot \
  --workers=4
```

The original-seed baseline extension is:

```sh
Rscript --vanilla run_all.R \
  --config=config/paper_config.R \
  --output-dir=results/baseline \
  --workers=4
```

The fresh-seed confirmatory run is:

```sh
Rscript --vanilla run_all.R \
  --config=config/confirmatory_config.R \
  --output-dir=results/confirmatory \
  --workers=4
```

The same workflow can be run as separate population, simulation, aggregation,
and table stages:

```sh
Rscript --vanilla scripts/01_population_check.R \
  --config=config/confirmatory_config.R \
  --output-dir=results/confirmatory \
  --workers=4

Rscript --vanilla scripts/02_run_simulation.R \
  --config=config/confirmatory_config.R \
  --output-dir=results/confirmatory \
  --workers=4

Rscript --vanilla scripts/03_aggregate_results.R \
  --config=config/confirmatory_config.R \
  --output-dir=results/confirmatory

Rscript --vanilla scripts/04_make_markdown_tables.R \
  --config=config/confirmatory_config.R \
  --output-dir=results/confirmatory
```

For a batch scheduler, `scripts/run_batch.sh` records code checksums and the R
environment, runs the tests first, and stops at the first failure. Replication
seeds depend only on scenario and replication identifiers, so worker
scheduling cannot change simulated data. A cached replication is reused only
when both its configuration and implementation fingerprints match; `--force`
requests recomputation.

## Expected outputs

The simulation output directory contains:

- `replicates/`: cacheable replication RDS objects;
- `raw_replicates.csv`: one row per scenario, replication, and method;
- `tau_selection.csv`: the full midpoint criterion over the $\tau$ grid;
- `equivalence_checks.csv`: rank and depth identities and region agreement;
- `domain_diagnostics.csv`, `grid_refinement.csv`, and
  `depth_geometry_diagnostics.csv`;
- `population_curves.csv`, `population_summary.csv`, and population
  equivalence and tuning diagnostics;
- `modality_checks.csv` and `scale_equivariance.csv`;
- `aggregate_summary.csv`, `tau_summary.csv`, `domain_summary.csv`,
  `grid_refinement_summary.csv`, and `depth_geometry_summary.csv`;
- four independent recalculation CSVs for coverage, length, peaks, and HDR;
- `validation_checks.csv`;
- `results_table.md`, generated from the saved results; and
- configuration-bound simulation and population manifests.

The full configuration uses 200 independent replications in each of six
scenarios, with 200 training, 200 calibration, and 2,000 test observations.
The proposed score selects $\tau\in\{0,1,\ldots,20\}$ by the midpoint
proportion plus $0.0001\tau$. Because the same calibration sample selects
$\tau$ and sets the final threshold, fixed-score split-conformal validity does
not automatically carry over to the tuned score; coverage is reported
empirically.

## File map

```text
.
├── README.md
├── run_all.R
├── R/
│   ├── aggregation.R
│   ├── distributions.R
│   ├── markdown_tables.R
│   ├── metrics.R
│   ├── population.R
│   ├── prediction_sets.R
│   ├── profiles.R
│   ├── score_metric_halfspace_depth.R
│   ├── score_metric_spatial_depth.R
│   ├── score_proposed.R
│   ├── score_transport_cost.R
│   ├── score_transport_rank.R
│   ├── set_geometry.R
│   ├── set_length.R
│   ├── simulation.R
│   └── utilities.R
├── config/
│   ├── confirmatory_config.R
│   ├── paper_config.R
│   └── pilot_config.R
├── scripts/
│   ├── 01_population_check.R
│   ├── 02_run_simulation.R
│   ├── 03_aggregate_results.R
│   ├── 04_make_markdown_tables.R
│   └── run_batch.sh
└── tests/
    ├── smoke_test.R
    ├── test_exact_depth_geometry.R
    ├── test_metric_depths.R
    ├── test_r1_equivalence.R
    └── test_set_length.R
```
