# Validation scripts

Scripts used to compare CroCoDeELR with CroCoDeEL 1.3.1. They are not part
of the installed package (listed in `.Rbuildignore`). Results are in
`results/` and summarized in the main README.

Paths in the scripts point to the environment where the validation was run
(`/home/claude/...`). Adjust them before running again.

| Script | Purpose |
|---|---|
| `dump_numpy_subsets.py` | Write the RANSAC subsets drawn by the Python tool |
| `dump_python_features.py` | Features and probabilities of every pair (Python) |
| `helpers.R` | Inject NumPy subsets into the R port; R features of every pair |
| `run_r_searches.R` | R runs in two modes: R draws, NumPy draws |
| `python_seed_runs.py` | Python runs with another RANSAC seed (yardstick) |
| `compare_results.R` | Comparison tables for the Meteor datasets |
| `run_r_metaphlan.R`, `compare_metaphlan.R` | MetaPhlAn4 dataset (900 samples) |
| `r_train_eval.R`, `python_train_eval.py` | Training on the official dataset; per-pair results |
| `compare_training.R` | Classification and rate comparison on the training dataset |
| `compare_trained_models.R` | ranger vs scikit-learn models on real data |
| `make_figures.R` | Figures 1-3 of the main README (`man/figures/`) |
| `simulate_training_table.R` | Simulated training table (first check, before the official dataset was available) |

Python reference outputs were made with `crocodeel search_conta
--probability-cutoff 0` (all pairs with a contamination line).
