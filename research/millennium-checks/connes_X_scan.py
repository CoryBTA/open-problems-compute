#!/usr/bin/env python3
"""Scan connes_weil_form.Form over X (36-function Legendre basis); record |Mellin zero - true zero| for the first 5 zeros.
Usage: connes_X_scan.py [n_basis]   -> writes connes_X_scan.json (also n_basis=24 if passed as a second run).
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from connes_weil_form import Form, mellin_zeros  # noqa: E402


def scan(Xs, n_basis):
    true = [float(mp.zetazero(k).imag) for k in range(1, 6)]
    out = []
    for X in Xs:
        form = Form(float(X), n_basis)
        w, V = np.linalg.eigh(0.5 * (form.A + form.A.T))
        zs = np.array(mellin_zeros(form, V[:, 0]))
        diffs = [float(abs(zs[np.argmin(abs(zs - t))] - t)) for t in true]
        out.append({"X": X, "n_basis": n_basis, "lowest_eig": float(w[0]), "n_mellin_zeros_below_130": int(len(zs)),
                    "true": true, "abs_diff_first5": diffs})
        print(X, n_basis, [f"{d:.2e}" for d in diffs], flush=True)
    return out


if __name__ == "__main__":
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    res = {"n_basis": nb, "runs": scan((13, 17, 23, 29), nb), "second_basis_24": scan((13, 17, 23, 29), 24)}
    Path(__file__).with_name("connes_X_scan.json").write_text(json.dumps(res, indent=1))
