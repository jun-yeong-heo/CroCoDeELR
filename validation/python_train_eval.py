"""Reproduce `crocodeel train_model` and keep per-pair results.

Same steps and seeds as crocodeel/train_model.py (features with NumPy RANSAC
draws, invalid pairs removed, train_test_split(test_size=0.3,
random_state=0), RandomForestClassifier(n_estimators=100, random_state=0)).
Writes one row per kept pair: source, target, label, set, probability.
Also writes the probabilities of the default CroCoDeEL model.
Usage: python python_train_eval.py training_table.tsv out.tsv model.joblib
"""
import sys
import warnings

import joblib
import numpy as np
from sklearn.model_selection import train_test_split

from crocodeel import ab_table_utils
from crocodeel.search_conta import Defaults
from crocodeel.train_model import (_compute_features, _reconstruct_sample_pairs,
                                   _train_model)

if __name__ == "__main__":
    ab_file, out, model_out = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(ab_file) as fh:
        tab = ab_table_utils.read_filter_normalize(fh, None)
    pairs = _reconstruct_sample_pairs(tab)
    labels = np.array([s.startswith("conta") for s, _ in pairs])
    feats = _compute_features(tab, pairs, 2)
    valid = ~np.any(np.isnan(feats), axis=1)
    feats, labels = feats[valid], labels[valid]
    kept = [p for p, v in zip(pairs, valid) if v]
    idx = np.arange(len(kept))
    f_tr, f_te, y_tr, y_te, i_tr, i_te = train_test_split(
        feats, labels, idx, test_size=0.3, random_state=0)
    model = _train_model(f_tr, y_tr, 100, 0, 2)
    model.set_params(n_jobs=1)
    joblib.dump(model, model_out, compress=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = joblib.load(Defaults.MODEL_FILE)
    p_new = model.predict_proba(feats)[:, 1]
    p_def = default.predict_proba(feats)[:, 1]
    test = np.zeros(len(kept), dtype=bool)
    test[i_te] = True
    with open(out, "w") as fh:
        fh.write("source\ttarget\tis_contaminated\tset\tprob_retrained\tprob_default\n")
        for k, (s, t) in enumerate(kept):
            fh.write(f"{s}\t{t}\t{labels[k]}\t{'test' if test[k] else 'train'}\t"
                     f"{float(p_new[k])!r}\t{float(p_def[k])!r}\n")
