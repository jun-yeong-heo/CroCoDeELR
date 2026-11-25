# CroCoDeELR 0.1.0 (2026-11-25)

First release. R translation of CroCoDeEL 1.3.1
(<https://github.com/metagenopolis/CroCoDeEL>, commit `f5bc209`).

Modifications relative to the original Python code (notice required by
GPL-3.0, section 5a):

* All modules were translated from Python to R:
  `ab_table_utils.py`, `conta_features.py`, `ransac.py`, `search_conta.py`,
  `conta_event.py`, `plot_conta.py`, `train_model.py`, `self_test.py`,
  `execution_description.py` and the `easy_wf` command.
* RANSAC candidate subsets are drawn with R's Mersenne-Twister generator
  (`sample.int()`, seed 42) instead of NumPy's generator. This is the only
  source of differences in contamination detection; see the README.
* The default Random Forest model `crocodeel_rf_Feb2026_2` was converted
  from scikit-learn to R without retraining
  (`inst/python/export_sklearn_model.py`, `data-raw/build_default_model.R`).
  Its predictions are identical to those of scikit-learn.
* `train_model()` trains the forest with 'ranger' instead of scikit-learn,
  with equivalent settings, and splits train/test sets with R's generator.
* Scatterplots are drawn with base R graphics instead of matplotlib.
* `self_test()` compares results with those expected for this port
  (`inst/extdata/contamination_events_R.tsv`).
* Parallel processing uses PSOCK clusters (`parallel` package) so that it
  also works on Windows. The default is one process.
* Abundance tables can be given as R matrices or data.frames, in addition
  to TSV files.
