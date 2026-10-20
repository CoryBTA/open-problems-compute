#!/usr/bin/env python3
"""Viscosity scan at N=64 (IFRK4, CFL 0.5, T=3, same run() as ns_ic_study.py) for one IC; usage: ns_scan_n64.py kida|rand1.
Writes ns_scan_n64_<ic>.json (temporary; `--merge` combines them into ns_scan_n64.json). Checkpoints per run (ckpt_scan_*.npz, git-ignored).
"""
import json
import sys
from pathlib import Path

import ns_ic_study as s

HERE = Path(__file__).resolve().parent
if __name__ == "__main__":
    if sys.argv[1] == "--merge":
        rows = []
        for ic in ("kida", "rand1"):
            rows += json.loads((HERE / f"ns_scan_n64_{ic}.json").read_text())
        (HERE / "ns_scan_n64.json").write_text(json.dumps(rows, indent=1))
        sys.exit()
    ic, rows = sys.argv[1], []
    for nu in (0.02, 0.01, 0.005):
        s.NU = nu
        r = s.run_robust(name=ic, N=64, T=3.0, ckpt=f"ckpt_scan_{ic}_{nu}.npz")
        r["N"], r["nu"] = 64, nu
        rows.append(r)
        print(json.dumps(r), flush=True)
        (HERE / f"ns_scan_n64_{ic}.json").write_text(json.dumps(rows, indent=1))
