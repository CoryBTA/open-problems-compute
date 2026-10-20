#!/usr/bin/env python3
"""Initial-data (t = 0) test of the corrected Step 7 hypothesis hG: G <= A*F with F > 1, where G = ||grad w||^2/Omega = R*F/4.
Family: Roberts-type Beltrami field at wavenumber k (|w| constant) plus eps times a Taylor-Green-like field (wavenumber 1) at N=32 and N=48 (same field).
Reports R, F, G for each (k, eps); writes step7_family_scan.json. Algebraic check on data, not on evolved NS solutions.
"""
import json
from pathlib import Path

import ns_ic_study as s
import numpy as np

HERE = Path(__file__).resolve().parent


def row(k, eps, N):
    sol = s.base.NS3D(N, 0.01)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ux = np.sin(k * Z) + eps * np.sin(X) * np.cos(Y) * np.cos(Z)
    uy = np.cos(k * Z) - eps * np.cos(X) * np.sin(Y) * np.cos(Z)
    uz = np.zeros_like(X)
    w = tuple(np.asarray(c) for c in sol.curl(np.fft.fftn(ux), np.fft.fftn(uy), np.fft.fftn(uz)))
    m = s.measure(sol, w)
    return {"k": k, "eps": eps, "N": N, "R": m["R"], "F": m["F"], "G": m["R"] * m["F"] / 4, "Omega": m["Omega"]}


if __name__ == "__main__":
    rows = [row(k, eps, N) for N in (32, 48) for k in (1, 2, 3, 4) for eps in (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0)]
    (HERE / "step7_family_scan.json").write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(f'N={r["N"]} k={r["k"]} eps={r["eps"]:<5g} F={r["F"]:<10.4g} R={r["R"]:<10.4g} G={r["G"]:.4g}')


def at_F(k, target, N=48):
    """Bisect eps (log scale) so that F is about `target`; returns the row (F is increasing in eps in this family for eps >= 0.01)."""
    lo, hi = 0.01, 200.0
    for _ in range(40):
        mid = (lo * hi) ** 0.5
        if row(k, mid, N)["F"] < target:
            lo = mid
        else:
            hi = mid
    return row(k, hi, N)
