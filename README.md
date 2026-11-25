# CroCoDeELR

CroCoDeELR is an **unofficial R port** of
[CroCoDeEL](https://github.com/metagenopolis/CroCoDeEL) 1.3.1. CroCoDeEL
detects cross-sample (well-to-well) contamination in shotgun metagenomic
data. It needs only a species abundance table: no negative controls and no
plate map.

This port is not affiliated with or endorsed by the authors of CroCoDeEL.
Please cite the original paper when you use it (see [Citation](#citation)).
For questions on the method, refer to the
[original repository](https://github.com/metagenopolis/CroCoDeEL).

## Installation

CroCoDeELR is written in pure R. Installing it does not require compilation
(no Rtools on Windows).

```r
# From GitHub
remotes::install_github("jun-yeong-heo/CroCoDeELR")

# From a local copy of this folder
install.packages("path/to/CroCoDeELR", repos = NULL, type = "source")

# Optional: needed only for train_model() and import_sklearn_model()
install.packages(c("ranger", "jsonlite"))
```

Check the installation on the toy dataset of CroCoDeEL:

```r
CroCoDeELR::self_test()
```

## Usage

```r
library(CroCoDeELR)

# Species in rows, samples in columns. Counts or relative abundances.
ab <- read_ab_table("species_abundance.tsv")

# Search contamination among all sample pairs.
events <- search_conta(ab, nproc = 4)

# Scatterplots of each event (PDF).
plot_conta(ab, events, file = "contamination_events.pdf")

# Save events in the CroCoDeEL format (readable by the CroCoDeEL
# interpretation interface).
write_conta_events(events, "contamination_events.tsv")
```

`search_conta()` also accepts a numeric matrix or a data.frame built in R,
with species names as row names.

All steps on files, as `crocodeel easy_wf`:

```r
easy_wf("species_abundance.tsv",
        events_file = "contamination_events.tsv",
        pdf_file = "contamination_events.pdf")
```

### Correspondence with the CroCoDeEL command line

| CroCoDeEL | CroCoDeELR |
|---|---|
| `crocodeel easy_wf` | `easy_wf()` |
| `crocodeel search_conta` | `search_conta()` + `write_conta_events()` |
| `crocodeel plot_conta` | `read_conta_events()` + `plot_conta()` |
| `crocodeel train_model` | `train_model()` |
| `crocodeel self_test` | `self_test()` |
| `-s2` (second table) | `ab_table_2` |
| `--filter-low-ab` | `filter_low_ab` |
| `-m` (model file) | `model = load_model(file)` |
| `--probability-cutoff`, `--rate-cutoff` | `probability_cutoff`, `rate_cutoff` |
| `--nproc` | `nproc` (default 1) |
| `--color-conta-species`, `--no-conta-line` | `color_conta_species`, `conta_line` |

### Input tables

The recommendations of CroCoDeEL apply. Meteor is the recommended profiler.
MetaPhlAn4 can be used at the SGB level (`t__` rows) with
`filter_low_ab = 20`. sylph can be used with GTDB representative genomes
(species-level `t__` rows). The authors advise against other profilers,
which often estimate subdominant species poorly.

### Models

The default model is the one shipped with CroCoDeEL 1.3.1
(`crocodeel_rf_Feb2026_2`), converted to R without retraining.

A model trained with the Python `crocodeel train_model` can be used after
conversion:

```sh
python <path>/export_sklearn_model.py model.joblib model.json
```

```r
system.file("python", package = "CroCoDeELR")   # <path> above
model <- import_sklearn_model("model.json", output = "model.rds")
events <- search_conta(ab, model = model)
```

A new model can be trained in R with 'ranger'. Sample names must follow
the CroCoDeEL training scheme (`conta_source_case_1`,
`conta_target_case_1`, `non_conta_source_case_1`, ...). The training
dataset of CroCoDeEL is available at <https://doi.org/10.57745/N6JSHQ>.

```r
res <- train_model(read_ab_table("training_dataset.meteor.tsv"),
                   model_file = "my_model.rds", report_file = "report.json",
                   nproc = 4)
```

## Differences from CroCoDeEL

| Component | CroCoDeEL 1.3.1 | CroCoDeELR |
|---|---|---|
| RANSAC random draws | NumPy generator (sklearn), seed 42 | R Mersenne-Twister (`sample.int`), seed 42 |
| Default model | scikit-learn Random Forest | Same trees, converted; same probabilities |
| Model training | scikit-learn | 'ranger' with equivalent settings |
| Train/test split | sklearn `train_test_split` | R generator, same test fraction |
| Plots | matplotlib | base R graphics |
| Parallelism | `multiprocessing` | PSOCK cluster (works on Windows) |

Everything else follows the original code: table preprocessing, candidate
species, RANSAC model selection and stopping rule, the six features, line
refinement, cutoffs, sorting and the output format. The only source of
differences in detection is the RANSAC random draws. The RANSAC seed is
fixed, so results are reproducible across runs and machines. The caller's
random number stream is not changed.

## Validation against CroCoDeEL 1.3.1

All scripts are in `validation/`. Results are in `validation/results/`.
Reference runs used CroCoDeEL 1.3.1 (Python 3.13, scikit-learn 1.9.1) and
CroCoDeELR 0.1.0 (R 4.3.3) on Linux, on 2026-11-25. Datasets are the
toy dataset of CroCoDeEL and public Meteor tables from the CroCoDeEL data
repository (<https://doi.org/10.57745/N6JSHQ>).

### 1. Default model

On 20,000 random feature vectors, the converted model gives the same
probabilities as scikit-learn: all 20,000 values are bitwise identical.
The test suite checks 300 more vectors, 100 of which sit exactly on split
thresholds.

### 2. Port without the random draws

The RANSAC subsets drawn by the Python tool were injected into the R port.
The R port then gave the same result as CroCoDeEL for every sample pair
with a contamination line: same pairs, same species lists, same formatted
rates and the same probabilities (maximum absolute difference 0). On the toy
dataset, feature values differed by at most 1e-14 (floating-point summation
order). The
test suite repeats this check on the toy dataset, in one-table and
two-table modes, by comparing the output files byte by byte.

| Dataset | Samples | Pairs with a line | Identical |
|---|---:|---:|:---:|
| toy | 14 | 49 | yes |
| PRJEB83730 | 17 | 110 | yes |
| PRJNA698986_P2 | 78 | 104 | yes |
| PRJNA698986_P3 | 91 | 294 | yes |
| PRJEB12449 | 110 | 1484 | yes |
| PRJEB10878 | 128 | 1984 | yes |
| PRJNA352475 | 182 | 1762 | yes |
| PRJEB6337 | 237 | 8935 | yes |

### 3. Port as distributed (R random draws)

**Criterion.** The R draws change which candidate line RANSAC picks for some
pairs. The question is whether this changes results more than the original
tool's own dependence on its draws. As a yardstick, CroCoDeEL was rerun with
RANSAC seeds 1, 2 and 3 instead of 42. The R port is considered within
tolerance when its difference from CroCoDeEL (seed 42) is in the range of
these seed-to-seed differences.

Events are pairs with probability ≥ 0.5. Jaccard = shared events / events
found by either run. Probability and rate differences are computed on shared
events. Rate differences are relative.

| Dataset | Events Py | Events R | Shared | Jaccard R vs Py | Jaccard Py seeds (min / mean) | Max abs. diff. prob. R / Py seeds | 95th pct rate diff. R / Py seeds |
|---|---:|---:|---:|---:|---:|---:|---:|
| toy | 7 | 7 | 7 | 1.00 | 1.00 / 1.00 | 0.00 / 0.00 | 1.7% / 1.6% |
| PRJEB83730 | 26 | 26 | 26 | 1.00 | 0.96 / 0.99 | 0.19 / 0.17 | 5.1% / 3.6% |
| PRJNA698986_P2 | 26 | 30 | 25 | 0.81 | 0.78 / 0.80 | 0.18 / 0.18 | 8.3% / 5.6% |
| PRJNA698986_P3 | 24 | 24 | 24 | 1.00 | 0.88 / 0.93 | 0.31 / 0.42 | 5.4% / 7.7% |
| PRJEB12449 | 77 | 73 | 72 | 0.92 | 0.88 / 0.89 | 0.28 / 0.28 | 6.3% / 6.8% |
| PRJEB10878 | 88 | 85 | 84 | 0.94 | 0.94 / 0.96 | 0.16 / 0.22 | 8.1% / 8.0% |
| PRJNA352475 | 96 | 93 | 86 | 0.83 | 0.78 / 0.82 | 0.24 / 0.29 | 10.0% / 11.3% |
| PRJEB6337 | 197 | 198 | 181 | 0.85 | 0.87 / 0.88 | 0.35 / 0.39 | 12.2% / 12.5% |

MetaPhlAn4 dataset PRJNA763023/PRJDB4176 (900 samples, `filter_low_ab = 20`,
one extra Python seed):

| Comparison | Events | Shared | Jaccard | Max abs. diff. prob. | 95th pct rate diff. |
|---|---:|---:|---:|---:|---:|
| R vs Py (seed 42) | 14,395 vs 14,329 | 14,295 | 0.991 | 0.29 | 4.1% |
| Py seed 1 vs Py seed 42 | 14,358 vs 14,329 | 14,298 | 0.994 | 0.26 | 3.9% |

**Result.**

- In 7 of 8 Meteor datasets, the event Jaccard index of R vs Python is
  within the range of Python seed changes. In PRJEB6337 it is 0.85, below
  the lowest of three Python seeds (0.87). Three seeds give only a rough
  estimate of that range.
- The median relative rate difference of shared events is 0 to 1.5%
  (Python seeds: 0 to 1.2%). The 95th percentile is 1.7 to 12.2% (Python
  seeds: 1.6 to 12.5%).
- In the Meteor datasets, events found by only one of the two runs had a
  probability below 0.5 in the other run (at most 0.49). They are pairs near
  the probability cutoff.
- With 900 samples, the R port finds 99% of the Python events. Python seed
  1 also finds 99%.

**Interpretation.** CroCoDeEL results depend on the RANSAC draws for pairs
near the probability cutoff. This holds for the original tool as well. The
R port differs from CroCoDeEL by about as much as CroCoDeEL differs from
itself with another seed. Exact Python results can only be reproduced with
NumPy's draws (section 2). The authors' advice to review the scatterplots
of each event applies here too.

### 4. Model training

Both tools were trained on the official CroCoDeEL training dataset
(`training_dataset.meteor.tsv`, 14,470 labeled pairs; 5,850 contaminated)
with default settings (`validation/r_train_eval.R`,
`validation/python_train_eval.py`, `validation/compare_training.R`).

**Pairs kept for training.** Python kept 8,762 pairs and R kept 8,767. The
difference comes from the RANSAC draws for non-contaminated pairs. 8,645
pairs were kept by both.

**Classification.** Each tool split its pairs 70/30 with its own random
generator. 768 pairs ended up in both test sets.

| Evaluation | Pairs | Accuracy | Precision | Recall | F1 | AUC |
|---|---:|---:|---:|---:|---:|---:|
| Python, own test set | 2,629 | 0.9985 | 0.9989 | 0.9989 | 0.9989 | 0.9997 |
| R, own test set | 2,631 | 0.9985 | 0.9994 | 0.9983 | 0.9989 | 1.0000 |
| Python, pairs held out by both | 768 | 0.9987 | 1.0000 | 0.9981 | 0.9991 | 0.9990 |
| R, pairs held out by both | 768 | 0.9987 | 1.0000 | 0.9981 | 0.9991 | 0.9999 |

On the 768 pairs held out by both, the two tools made the same call for
every pair. Precision, recall and F1 refer to the contaminated class.

**Rate estimation.** For the 5,840 contaminated pairs with a line in both
tools, the median absolute difference between the R and Python rates was
0.0006 log10 units (0.14%). The maximum was 0.09 log10 units. Both tools had
the same error against the simulated rate (median 0.055 log10 units).

**Real data.** The retrained models were applied to the seven Meteor
datasets of section 3, with identical features (`compare_trained_models.R`).

| Dataset | Events sklearn | Events ranger | Jaccard ranger vs sklearn | Jaccard sklearn vs default model | Probability correlation |
|---|---:|---:|---:|---:|---:|
| PRJEB83730 | 28 | 29 | 0.97 | 0.86 | 0.999 |
| PRJNA698986_P2 | 32 | 29 | 0.85 | 0.77 | 0.974 |
| PRJNA698986_P3 | 23 | 24 | 0.88 | 0.96 | 0.989 |
| PRJEB12449 | 80 | 69 | 0.82 | 0.87 | 0.979 |
| PRJEB10878 | 86 | 89 | 0.90 | 0.92 | 0.990 |
| PRJNA352475 | 95 | 85 | 0.77 | 0.84 | 0.973 |
| PRJEB6337 | 200 | 186 | 0.80 | 0.83 | 0.977 |

- The ranger and scikit-learn models agree on 77 to 97% of events.
- Retraining with scikit-learn alone gives a model that also differs from the
  shipped default model (Jaccard 0.77 to 0.96). With scikit-learn 1.9.1 the
  shipped model could not be reproduced exactly (largest probability
  difference 0.68 on the training pairs). The training data or library
  version of the shipped model may differ.
- In 4 of 7 datasets the ranger model called fewer events than the
  scikit-learn model.

For contamination detection, the default model (converted from
scikit-learn, section 1) remains the recommended choice. `train_model()` is
meant for custom models.

### 5. Run time

On 2 cores (Linux), `train_model(nproc = 2)` took 25 s on the official
training dataset (Python: 35 s), plus 24 s to read the table.
`search_conta(nproc = 2)` took 36 s for 237 samples
(Python: 18 s) and 144 s for 900 samples (Python: 250 s, while other jobs
were running). Run times are rough.

## License

CroCoDeELR is free software under the GNU General Public License, version 3
or later (GPL-3.0-or-later), the license of CroCoDeEL. See `COPYING`.

The R code is a translation of the CroCoDeEL Python code, copyright the
CroCoDeEL authors (Lindsay Goulet, Florian Plaza Oñate, Edi Prifti, Eugeni
Belda, Emmanuelle Le Chatelier, Guillaume Gautreau). The default model and
the toy dataset come from CroCoDeEL 1.3.1. The modifications are listed in
`NEWS.md`.

## Citation

Please cite the original method:

> Goulet, L. et al. CroCoDeEL: accurate control-free detection of
> cross-sample contamination in metagenomic data. *Nature Communications*
> (2026). <https://doi.org/10.1038/s41467-026-72637-9>

`citation("CroCoDeELR")` prints the same reference.
