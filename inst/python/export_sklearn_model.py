"""Export a CroCoDeEL scikit-learn Random Forest model to JSON.

This file is part of CroCoDeELR, an unofficial R port of CroCoDeEL
(https://github.com/metagenopolis/CroCoDeEL).
Copyright (C) the CroCoDeEL authors (original Python code) and
Jun-yeong Heo (R port).
License: GPL-3.0-or-later. See the COPYING file distributed with CroCoDeELR.

The JSON file stores, for every tree, the arrays that sklearn uses for
prediction (children, split feature, split threshold, node class values).
Read it in R with CroCoDeELR::import_sklearn_model().

Usage:
    python export_sklearn_model.py model.joblib model.json

Requirements: Python with joblib and the scikit-learn version that can
unpickle the model.
"""

import json
import sys
import warnings

import joblib
from sklearn.ensemble import RandomForestClassifier


def export_model(model_path: str, json_path: str) -> None:
    """Write the trees of a RandomForestClassifier to a JSON file."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = joblib.load(model_path)

    if not isinstance(model, RandomForestClassifier):
        raise SystemExit(f"{model_path} is not an sklearn RandomForestClassifier")

    if model.n_outputs_ != 1 or len(model.classes_) != 2:
        raise SystemExit("Only binary, single-output classifiers are supported")

    trees = []
    for estimator in model.estimators_:
        tree = estimator.tree_
        trees.append(
            {
                "children_left": tree.children_left.tolist(),
                "children_right": tree.children_right.tolist(),
                "feature": tree.feature.tolist(),
                # Python floats are written with repr(), which round-trips
                # float64 values exactly.
                "threshold": tree.threshold.tolist(),
                "value": tree.value[:, 0, :].tolist(),
            }
        )

    payload = {
        "format": "crocodeel_rf_sklearn_v1",
        "classes": [bool(c) for c in model.classes_],
        "n_features": int(model.n_features_in_),
        "sklearn_version": getattr(model, "_sklearn_version", "unknown"),
        "params": {
            key: value
            for key, value in model.get_params().items()
            if isinstance(value, (int, float, str, bool, type(None)))
        },
        "trees": trees,
    }

    with open(json_path, "w", encoding="utf8") as fh:
        json.dump(payload, fh)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    export_model(sys.argv[1], sys.argv[2])
