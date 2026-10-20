#!/usr/bin/env python3
"""Random 3-SAT satisfiable fraction vs clause ratio alpha = m/n (python-sat, Glucose4), n = 50/100/150.
Calibration only: finite-size crossing of fraction 0.5 vs the known asymptotic threshold ~4.267 (not proved; known numerically).
Writes sat_threshold.json.
"""
import json
import random
import time
from pathlib import Path

from pysat.solvers import Glucose4

ALPHAS = [3.8, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.7]
NS, TRIALS, SEED = (50, 100, 150), 200, 12345


def instance(n, m, rng):
    cl = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cl.append([v if rng.random() < 0.5 else -v for v in vs])
    return cl


def main():
    rng = random.Random(SEED)
    res = {}
    for n in NS:
        row = []
        for a in ALPHAS:
            m, sat, t0 = round(a * n), 0, time.time()
            for _ in range(TRIALS):
                with Glucose4(bootstrap_with=instance(n, m, rng)) as s:
                    sat += s.solve()
            row.append({"alpha": a, "m": m, "sat_fraction": sat / TRIALS, "seconds": round(time.time() - t0, 1)})
            print(n, a, sat / TRIALS, flush=True)
        # linear interpolation of the 0.5 crossing
        cross = None
        for r0, r1 in zip(row, row[1:], strict=False):
            if r0["sat_fraction"] >= 0.5 > r1["sat_fraction"]:
                cross = r0["alpha"] + (r0["sat_fraction"] - 0.5) / (r0["sat_fraction"] - r1["sat_fraction"]) * (r1["alpha"] - r0["alpha"])
        res[str(n)] = {"rows": row, "alpha_half": cross}
    out = {"solver": "Glucose4", "trials_per_point": TRIALS, "seed": SEED, "results": res}
    Path(__file__).with_name("sat_threshold.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
