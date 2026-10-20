#!/usr/bin/env python3
"""BSD leading-term check with a large coefficient window, aimed at curves whose predicted |Sha| exceeds 1.

Same model family and formula as bsd_leading_term.py (y^2 + y = x^3 + a2 x^2 + a4 x + a6, |Delta| = p prime, so
conductor p, k = 1, c_p = 1 split / 1 nonsplit). The window is a2 in -1..1, a4 in -1000..1000, a6 in -12000..12000,
|Delta| prime in 11..1500. Delta is computed vectorised; analyse() is the existing function. Writes bsd_sha_search.json.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("bsd_leading_term", HERE / "bsd_leading_term.py")
bl = importlib.util.module_from_spec(spec)
sys.modules["bsd_leading_term"] = bl
spec.loader.exec_module(bl)

A2, A4, A6, PMAX = (-1, 1), 1000, 12000, 1500


def main():
    a6 = np.arange(-A6, A6 + 1, dtype=np.int64)
    cands = []
    for a2 in range(A2[0], A2[1] + 1):
        for a4 in range(-A4, A4 + 1):
            b2, b4 = 4 * a2, 2 * a4
            b6 = 1 + 4 * a6
            b8 = 4 * a2 * a6 + a2 - a4 * a4
            delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
            hit = np.nonzero((np.abs(delta) >= 11) & (np.abs(delta) <= PMAX))[0]
            for i in hit:
                cands.append((a2, a4, int(a6[i]), abs(int(delta[i]))))
    found, seen = [], set()
    for a2, a4, a6v, d in cands:
        pk = bl.prime_power(d)
        if not pk or pk[1] != 1:
            continue
        p, k = pk
        if bl.invariants(a2, a4, a6v)[4] % p == 0:
            continue
        key = (p, tuple(bl.bsd.ap((0, a2, 1, a4, a6v), q) for q in (2, 3, 5, 7, 13)))
        if key in seen:
            continue
        seen.add(key)
        r = bl.analyse(a2, a4, a6v, p, k)
        if r:
            found.append(r)
    found.sort(key=lambda r: (r["N"], r["model"]))
    for r in found:
        print(f"N={r['N']:>5} model={r['model']} c_p={r['c_p']} tors<={r['torsion_bound']} L/Omega={r['L_over_Omega']:.6f} Sha_pred={r['sha_pred']:.6f} square={r['sha_is_square']}")
    out = {"window": {"a2": list(A2), "a4": [-A4, A4], "a6": [-A6, A6], "p_max": PMAX}, "n_candidates": len(cands), "n_curves": len(found),
           "n_not_square": sum(not r["sha_is_square"] for r in found), "curves": found}
    (HERE / "bsd_sha_search.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
