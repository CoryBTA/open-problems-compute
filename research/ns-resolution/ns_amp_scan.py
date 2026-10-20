#!/usr/bin/env python3
"""Amplitude, viscosity and C(t) scan (IFRK4, CFL 0.5, T=3, run() of ns_ic_study.py with series=True).
usage: ns_amp_scan.py LABEL ic:amp:nu:N [ic:amp:nu:N ...]  -> writes ns_amp_scan_LABEL.json (rows incl. per-time R, F, C series).
Checkpoints per run (ckpt_amp_*.npz, git-ignored); rerun the same command to resume. `--merge` combines all ns_amp_scan_*.json into ns_amp_scan.json.
"""
import glob
import json
import sys
from pathlib import Path

import ns_ic_study as s

HERE = Path(__file__).resolve().parent
if __name__ == "__main__":
    if sys.argv[1] == "--merge":
        rows = []
        for f in sorted(glob.glob(str(HERE / "ns_amp_scan_*.json"))):
            rows += json.loads(Path(f).read_text())
        (HERE / "ns_amp_scan.json").write_text(json.dumps(rows))
        sys.exit()
    label, rows = sys.argv[1], []
    out = HERE / f"ns_amp_scan_{label}.json"
    if out.exists():
        rows = json.loads(out.read_text())
    for spec in sys.argv[2:]:
        ic, amp, nu, N = spec.split(":")
        if any(r["ic"] == ic and r["amp"] == float(amp) and r["nu"] == float(nu) and r["N"] == int(N) for r in rows):
            continue
        s.AMP, s.NU = float(amp), float(nu)
        r = s.run_robust(name=ic, N=int(N), T=3.0, ckpt=f"ckpt_amp_{ic}_{amp}_{nu}_{N}.npz", series=True)
        r["amp"], r["nu"], r["N"] = float(amp), float(nu), int(N)
        rows.append(r)
        print(json.dumps({k: v for k, v in r.items() if k != "series"}), flush=True)
        out.write_text(json.dumps(rows))
