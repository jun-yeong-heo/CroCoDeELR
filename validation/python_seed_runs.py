"""Run the Python CroCoDeEL with another RANSAC seed (default is 42).

Measures how much the original tool's results depend on the RANSAC draws,
as a yardstick for the differences of the R port.
Usage: python python_seed_runs.py seed ab_table.tsv out.tsv [filter_low_ab] [prob_cutoff]
"""
import sys

from crocodeel import ab_table_utils, conta_features
from crocodeel.conta_event import ContaminationEventIO
from crocodeel.search_conta import Defaults, run_search_conta

if __name__ == "__main__":
    seed, ab_file, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    filter_low_ab = float(sys.argv[4]) if len(sys.argv) > 4 else None
    prob_cutoff = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    conta_features.ContaminationFeatureExtractor.RANSAC_RANDOM_STATE = seed
    with open(ab_file) as fh:
        tab = ab_table_utils.read_filter_normalize(fh, filter_low_ab)
    with open(Defaults.MODEL_FILE, "rb") as mfh:
        events = run_search_conta(tab, None, mfh, prob_cutoff, 0.0, nproc=2)
    with open(out, "w") as ofh:
        ContaminationEventIO.write_tsv(events, ofh)
