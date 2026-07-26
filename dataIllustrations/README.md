# Real-data illustrations

This directory contains four applications of distance-profile centrality
scores:

- `colon.R` treats gene-specific normal-to-tumor optimal transport maps as
  functional objects.
- `temp.R` compares annual summer distributions of daily minimum temperature.
- `mnist.R` identifies representative and atypical handwriting within each
  digit class.
- `taxi.R` compares daily Manhattan transport networks and clusters their
  distance profiles.

All application scripts source the single implementation in `tolReg.R`, locate
inputs relative to this directory, and write numerical and figure output below
`results/`. They do not read from the simulation folders or the parent
project.

## Requirements

The applications were tested with R 4.5.1.

| Analysis | Direct packages |
|---|---|
| Colon | `ggplot2` 4.0.1, `patchwork` 1.3.2, `pracma` 2.4.4 |
| Temperature | `ggplot2` 4.0.1, `patchwork` 1.3.2 |
| MNIST | base R only |
| Taxi | `ggplot2` 4.0.1, `dbscan` 1.2.5 |

Install the non-base packages from an R session with:

```r
install.packages(c("ggplot2", "patchwork", "pracma", "dbscan"))
```

`tolReg.R` has an optional RSpectra branch for a general feature-MDS helper,
but none of these four scripts calls that helper; RSpectra is not required.

## Input-data tree

The scripts expect this exact layout:

```text
dataIllustrations/
├── data/
│   ├── colon/                         # user supplied; not distributed
│   │   ├── I2000.txt
│   │   ├── names2.txt
│   │   └── tissues.txt
│   ├── mnist/                         # user supplied; not distributed
│   │   ├── train-images.idx3-ubyte
│   │   └── train-labels.idx1-ubyte
│   ├── temperature/
│   │   └── tampa.csv                  # included
│   └── taxi/
│       └── taxi2019gl.rds             # included processed input
├── results/
│   ├── colon/
│   ├── mnist/
│   ├── temperature/
│   └── taxi/
├── colon.R
├── mnist.R
├── run_all.R
├── taxi.R
├── temp.R
└── tolReg.R
```

These SHA-256 checksums identify the exact inputs used for the reference
outputs:

| File | SHA-256 |
|---|---|
| `colon/I2000.txt` | `9fafa4d403a837085c48b70e20f150e4075d21bdcec5dd39c2672bbb6d68e7e8` |
| `colon/names2.txt` | `218f40a1be966b4ec88c6a2a330e0e4634bc92d914730ea4729c2c87b142af9d` |
| `colon/tissues.txt` | `d1993a2debb2b2046b5431f4909048766d2bad951526bbed8be7d6383d02f378` |
| `mnist/train-images.idx3-ubyte` | `ba891046e6505d7aadcbbe25680a0738ad16aec93bde7f9b65e87a2fc25776db` |
| `mnist/train-labels.idx1-ubyte` | `65a50cbbf4e906d70832878ad85ccda5333a97f0f4c3dd2ef09a8a9eef7101c5` |
| `temperature/tampa.csv` | `8b822d602ac260526e8227ea0b3c720432fb2ed5909bab4df5cac40dd965590d` |
| `taxi/taxi2019gl.rds` | `1d4150728ae5809cf0aa9ed1ab3feb14f3127c39411cbc24c8cd89abd527bea5` |

The colon and MNIST directories contain tracked placeholders so that their
expected locations are visible even though the source files are omitted.

## Commands

Run the following commands from the repository root. The scripts are
path-independent, but using one documented working directory makes output
locations unambiguous.

Run the two analyses whose inputs are included:

```sh
Rscript --vanilla dataIllustrations/run_all.R \
  --analyses=temperature,taxi
```

After placing the colon and MNIST inputs, run all four in the order colon,
temperature, MNIST, and taxi:

```sh
Rscript --vanilla dataIllustrations/run_all.R
```

Run applications individually with:

```sh
Rscript --vanilla dataIllustrations/colon.R
Rscript --vanilla dataIllustrations/temp.R
Rscript --vanilla dataIllustrations/mnist.R
Rscript --vanilla dataIllustrations/taxi.R
```

The runner checks all selected inputs before starting, so an all-analysis run
cannot stop halfway merely because a later input file is absent.

## Colon gene expression

### Purpose and preprocessing

`colon.R` reads a 2,000 by 62 expression matrix, tissue identifiers, and gene
annotations. It removes the normal-tissue sample labelled 12, recovers 21
matched normal/tumor pairs, and constructs one empirical optimal transport map
per gene on a 0.01 grid. The original analysis has one tied interpolation knot
at gene index 1,358; the script applies the recorded `1e-5` correction there.
It then fits functional C-scores with `tau = 0` and marks the lowest 5% of
scores as peripheral.

There is no random sampling in this analysis.

### Data source, citation, and access

The data accompany:

> Alon, U., Barkai, N., Notterman, D. A., Gish, K., Ybarra, S., Mack, D.,
> and Levine, A. J. (1999). Broad patterns of gene expression revealed by
> clustering analysis of tumor and normal colon tissues probed by
> oligonucleotide arrays. *PNAS*, 96(12), 6745--6750.
> [doi:10.1073/pnas.96.12.6745](https://doi.org/10.1073/pnas.96.12.6745)

The original study page is
`http://genomics-pubs.princeton.edu/oncology/affydata/index.html`; its legacy
server is intermittently unavailable. An
[archived copy of the source page](https://web.archive.org/web/20230602123200/http://genomics-pubs.princeton.edu/oncology/affydata/index.html)
links the `I2000`, `names`, and `tissues` pages. Save their plain-text data as
`I2000.txt`, `names2.txt`, and `tissues.txt`, respectively, and place them as
shown in the input tree. Verify the checksums before running the analysis. The
source does not state an explicit redistribution license, so the files are not
included here. Users are responsible for following the source terms.

### Outputs and success check

| Output | Content | Manuscript mapping |
|---|---|---|
| `results/colon/colon_c_scores.csv` | 2,000 scores, peripheral ranks, flags | numerical values for Figure `fig:colon` |
| `results/colon/colon.pdf` | transport maps, profile/object MDS, highlighted densities | Figure `fig:colon` |

A successful run prints `Colon analysis complete`, produces both nonempty
files, flags 100 genes, and includes T57619, T48804, and T61609 among the
highlighted peripheral genes.

## Tampa minimum temperatures

### Purpose and preprocessing

`temp.R` constructs one empirical distribution of June--August daily minimum
temperatures for each year from 1939 through 2023. It applies the two recorded
2019 source corrections used in the analysis: summer position 62 is set to 90,
and positions 26 and 27 are set to 78. Each year's 92 values are then sorted
before distributional C-scores and profile/object MDS are computed. The
implemented automatic criterion over candidates 0, 1, 2, 3, 5, 10, and 20
selects `tau = 20`.

There is no random sampling in this analysis.

### Data source, citation, and access

`data/temperature/tampa.csv` is included. It is a daily-summary export for
station `USW00012842`, Tampa International Airport, from
[NOAA/NCEI Climate Data Online](https://www.ncei.noaa.gov/cdo-web/).
NOAA-produced federal environmental data are in the public domain in the
United States; the [NCEI archive policy](https://www.ncei.noaa.gov/archive)
describes its open-data and CC0 practice. Cite NOAA/NCEI, the Daily Summaries
dataset, the station identifier, and the access date in downstream work.

### Outputs and success check

| Output | Content | Manuscript mapping |
|---|---|---|
| `results/temperature/temperature_c_scores.csv` | 85 annual scores and peripheral ranks | numerical values for Figure `fig:temp` |
| `results/temperature/temperature.pdf` | densities and profile/object MDS | Figure `fig:temp` |

A successful run prints `Temperature analysis complete`, produces both
nonempty files, and ranks 2023 and 2022 as the two most peripheral years. The
bundled scores match the automatic `tau = 20` fit. The manuscript caption for
Figure `fig:temp` reports `tau = 10`; that caption value does not reproduce the
bundled numerical output.

## MNIST handwriting

### Purpose, preprocessing, and seed

`mnist.R` validates the standard IDX headers, rotates each 28 by 28 image to
the figure orientation, and analyzes each digit class separately. For every
digit from 0 through 9 it resets `set.seed(1)`, samples 200 training images,
scales pixel values by 1/50, and computes Euclidean C-scores. Automatic tuning
also resets the midpoint-subsampling seed to 1 for each digit. The figure
shows the four highest- and four lowest-scoring images per digit. Automatic
tuning selects `tau = 20` for digits 0, 1, and 4, and `tau = 0` for digits 2,
3, 5, 6, 7, 8, and 9.

### Data source, citation, and access

Download
[`train-images-idx3-ubyte.gz`](https://storage.googleapis.com/cvdf-datasets/mnist/train-images-idx3-ubyte.gz)
and
[`train-labels-idx1-ubyte.gz`](https://storage.googleapis.com/cvdf-datasets/mnist/train-labels-idx1-ubyte.gz)
from the mirror listed by
[TensorFlow Datasets](https://www.tensorflow.org/datasets/catalog/mnist).
Decompress them and place the two uncompressed files under `data/mnist/` with
the names shown in the input tree.

Cite:

> LeCun, Y., Bottou, L., Bengio, Y., and Haffner, P. (1998).
> Gradient-based learning applied to document recognition.
> *Proceedings of the IEEE*, 86(11), 2278--2324.
> [doi:10.1109/5.726791](https://doi.org/10.1109/5.726791)

The distribution page does not state an explicit license for redistribution
of the IDX files, so they are not included here. Users are responsible for
checking and following the source terms.

### Outputs and success check

| Output | Content | Manuscript mapping |
|---|---|---|
| `results/mnist/mnist_c_scores.csv` | sampled source indices, 2,000 scores, ranks | numerical values for Figure `fig:mnist` |
| `results/mnist/mnist.png` | representative and atypical digits | Figure `fig:mnist` |

A successful run prints `MNIST analysis complete`, writes a 2,000-row score
table, and produces a nonempty 8,000 by 1,600 pixel image. The class-specific
tuning values above reproduce the bundled scores. The manuscript caption for
Figure `fig:mnist` reports a single `tau = 10`; that caption value does not
reproduce the bundled numerical output for any digit.

## Manhattan taxi networks

### Purpose and preprocessing

`taxi.R` reads 365 daily 13 by 13 graph Laplacians for 2019. It divides each
matrix by its maximum absolute upper-triangular entry, computes network
C-scores, and uses Wasserstein distances between profiles for profile MDS and
DBSCAN. Profile clustering uses `eps = 0.08`, while clustering of the original
network distances uses `eps = 0.5`; both use `minPts = 3`.
The implemented automatic tuning criterion selects `tau = 20`.

Automatic tuning samples 10,000 candidate midpoints after resetting the
midpoint-subsampling seed to 1. Weekday labels are computed from numeric
calendar weekdays, so output does not depend on the machine's language
locale.

### Data source, citation, and access

The analysis-ready `data/taxi/taxi2019gl.rds` is included. It contains only
365 aggregate graph Laplacians and date labels, not raw trip records or
traveller identifiers. It was derived from the publicly downloadable
[NYC Taxi and Limousine Commission 2019 trip records](https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page).
That page states that the records were supplied by technology providers,
disclaims their accuracy, and is governed by the linked
[NYC terms of use](https://www.nyc.gov/main/terms-of-use); it does not provide
a separate dataset-license grant.

The construction groups 66 Manhattan taxi zones into 13 regions following:

> Zhou, Y. and Müller, H.-G. (2022). Network regression with graph
> Laplacians. *Journal of Machine Learning Research*, 23, 1--41.
> [Article page](https://jmlr.org/papers/v23/22-0681.html)

The raw monthly trip files and raw-to-Laplacian preprocessing code are not
included. Consequently, the source trip records alone are not sufficient for
bitwise regeneration of `taxi2019gl.rds`; the included aggregate object is the
authoritative input for this analysis.

### Outputs and success check

| Output | Content | Manuscript mapping |
|---|---|---|
| `results/taxi/taxi_c_scores.csv` | daily scores, ranks, weekday/weekend, 5% flags | Figure `fig:taxi` |
| `results/taxi/taxi_dbscan_clusters.csv` | weekday counts by clustering method and cluster | Table `tab:taxi` |
| `results/taxi/taxi_cluster_assignments.csv` | daily profile/object cluster assignments | Table `tab:taxi` diagnostics |
| `results/taxi/taxi_profile_mds.pdf` | profile MDS colored by score | Figure `fig:taxi` |

A successful run prints `Taxi analysis complete` and produces all four
nonempty files. The profile clustering separates the main weekday and weekend
groups, places July 4 with the weekend-like group, and recovers September
23--25 as a separate group. The bundled scores match the automatic
`tau = 20` fit. The manuscript caption for Figure `fig:taxi` reports
`tau = 10`; that caption value does not reproduce the bundled numerical
output.

## Runtime and memory

Approximate single-process requirements on a recent laptop are:

| Analysis | Typical wall time | Peak memory guidance |
|---|---:|---:|
| Colon | 1--3 minutes | allow 2 GB |
| Temperature | under 30 seconds | below 1 GB |
| MNIST | 1--5 minutes | allow 2 GB |
| Taxi | 1--2 minutes | allow 2 GB |

Package versions, BLAS implementation, and hardware can change these values.
The scripts run sequentially; `run_all.R` does not hold one application's fit
in memory while starting the next.

## Raw, processed, and generated files

- Source inputs: the colon text files, MNIST IDX files, and temperature CSV.
- Analysis-ready processed input: `taxi2019gl.rds`.
- Intermediate objects: transport maps, profiles, pairwise distances, MDS
  coordinates, and cluster fits exist only in memory and are not written.
- Generated outputs: everything below `results/`.

The small generated outputs in `results/` are included as reference artifacts.
Rerunning an application replaces the corresponding reference CSV and figure;
compare numerical CSV content rather than PDF metadata or byte-level hashes
when R or graphics-library versions differ.

## Limitations

- Colon and MNIST cannot run until the omitted source files are supplied.
- The legacy colon source page and MNIST distribution page may be intermittently
  unavailable; the checksums above allow files obtained from an authorized
  mirror to be verified.
- Exact taxi raw-data preprocessing is outside this repository; only the
  aggregate analysis-ready object and the downstream analysis are reproducible
  here.
- The scripts reproduce the specified analyses; they are not general data
  downloaders or data-cleaning pipelines.
