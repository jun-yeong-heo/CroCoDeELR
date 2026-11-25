"""Dump the RANSAC candidate subsets drawn by the Python CroCoDeEL (NumPy RNG).

Used only for validation: injecting these subsets into CroCoDeELR checks that
the rest of the R port reproduces the Python results exactly.
Usage: python dump_numpy_subsets.py max_n out.tsv
"""
import sys
from crocodeel.ransac import _candidate_subsets

max_n, out = int(sys.argv[1]), sys.argv[2]
with open(out, "w") as fh:
    for n in range(2, max_n + 1):
        s = _candidate_subsets(n, 30, 2, 42)
        for t in range(s.shape[0]):
            fh.write(f"{n}\t{s[t, 0] + 1}\t{s[t, 1] + 1}\n")
