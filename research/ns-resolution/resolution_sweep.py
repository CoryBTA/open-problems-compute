#!/usr/bin/env python3
"""Does R = ||grad omega||^2 / ||grad|omega|||^2 depend on grid size N at fixed physics?

Reuses the solver and diagnostic from ns_ratio_R.py (restored from git history, unchanged).
Config: Taylor-Green, amp=5, nu=0.01, T=2 (the 'TG5_Re500' run whose R at 48^3 was questioned).
CPU only (numpy). Writes resolution_sweep.json next to this file.
"""
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ns_ratio_R", HERE / "ns_ratio_R.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

sizes = [int(a) for a in sys.argv[1:]] or [24, 32, 48, 64, 96]
out = []
for N in sizes:
    r = mod.run_R_test(f"TG5_nu0.01_N{N}", 5, 0.01, 2.0, N=N)
    r["N"] = N
    out.append(r)
    (HERE / "resolution_sweep.json").write_text(json.dumps(out, indent=2, default=str))
print("done", [(r["N"], round(r["R_max"], 3)) for r in out])
