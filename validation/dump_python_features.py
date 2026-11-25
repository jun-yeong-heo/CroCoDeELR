"""Compute CroCoDeEL features and probabilities for every sample pair (Python).

Usage: python dump_python_features.py ab_table.tsv out.tsv
Output columns: source, target, 6 features, offset, probability, species.
"""
import sys
import warnings
from itertools import product

import joblib
import numpy as np

from crocodeel import ab_table_utils
from crocodeel.conta_features import ContaminationFeatureExtractor
from crocodeel.search_conta import Defaults

with open(sys.argv[1]) as fh:
    tab = ab_table_utils.read_filter_normalize(fh, None)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = joblib.load(Defaults.MODEL_FILE)

ext = ContaminationFeatureExtractor(tab)
rows = []
for s, t in product(tab.columns, tab.columns):
    if s == t:
        continue
    f = ext.extract(s, t)
    if f is None:
        continue
    rows.append((s, t, f))

X = np.vstack([f.values for _, _, f in rows]) if rows else np.empty((0, 6))
p = model.predict_proba(X)[:, 1] if rows else []
with open(sys.argv[2], "w") as fh:
    for (s, t, f), prob in zip(rows, p):
        vals = "\t".join(repr(float(v)) for v in f.values)
        fh.write(f"{s}\t{t}\t{vals}\t{repr(float(f.conta_line_offset))}\t{repr(float(prob))}\t{','.join(f.conta_line_species)}\n")
