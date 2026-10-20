"""Dispatch for CI/VM compute jobs. Each job checkpoints and stops itself inside its wall-time budget.

Exit codes: 0 finished (result JSON written), 75 partial (checkpoint saved; run again to resume), 2 unknown job.
A partial run writes `<ckpt-dir>/<job>.partial.json` and no result file, so a half-finished job is never committed as a result.

Usage: python scripts/ci_run_job.py JOB --ckpt-dir DIR --budget-min N
Jobs: ns_kida_n128, ns_kida_n192, ns_tg5_n128 (checkpointed, resumable), connes_flint, bsd_flint, sat_threshold (short; run as is).
"""
import argparse
import json
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EX_PARTIAL = 75
SAFETY_MARGIN_MIN = 10  # the budget is wall time for the solver loop; leave room for install, cache save and push


def _ns(name, N, T, nu, result, ckpt_name):
    def job(ckpt_dir, budget_s):
        sys.path.insert(0, str(ROOT / "research" / "ns-resolution"))
        import ns_ic_study as m

        done = ROOT / "research" / "ns-resolution" / result
        if done.exists():  # the weekly schedule re-fires a finished job: do nothing
            return {"partial": False, "skipped": f"{result} exists"}
        m.NU = nu
        out = m.run_robust(name=name, N=N, T=T, ckpt=str(Path(ckpt_dir) / ckpt_name), budget_s=budget_s)
        if out.get("partial"):
            return out
        out.update({"N": N, "nu": nu, "source": "ci_run_job.py"})
        (ROOT / "research" / "ns-resolution" / result).write_text(json.dumps(out, indent=1))
        return out

    return job


def _script(path):
    """A short job that is one existing script run as __main__ (its own JSON output)."""

    def job(ckpt_dir, budget_s):
        saved, sys.argv = sys.argv, [str(ROOT / path)]  # these scripts read their own argv (basis size, ...)
        try:
            runpy.run_path(str(ROOT / path), run_name="__main__")
        finally:
            sys.argv = saved
        return {"partial": False}

    return job


JOBS = {
    "ns_kida_n128": _ns("kida", 128, 3.0, 0.01, "ns_kida_n128_ci.json", "ckpt_kida128.npz"),
    "ns_kida_n192": _ns("kida", 192, 3.0, 0.01, "ns_kida_n192.json", "ckpt_kida192.npz"),
    "ns_tg5_n128": _ns("tg5", 128, 2.0, 0.01, "ns_tg5_n128.json", "ckpt_tg5_128.npz"),
    "connes_flint": _script("research/millennium-checks/connes_flint_zeros.py"),
    "bsd_flint": _script("research/millennium-checks/bsd_flint_period.py"),
    "sat_threshold": _script("research/millennium-checks/sat_threshold.py"),
}


RESULTS = {  # repo-relative result file per job: a job is finished when this file exists on origin/results or origin/main
    "ns_kida_n128": "research/ns-resolution/ns_kida_n128_ci.json",
    "ns_kida_n192": "research/ns-resolution/ns_kida_n192.json",
    "ns_tg5_n128": "research/ns-resolution/ns_tg5_n128.json",
    "connes_flint": "research/millennium-checks/connes_flint_zeros.json",
    "bsd_flint": "research/millennium-checks/bsd_flint_period.json",
    "sat_threshold": "research/millennium-checks/sat_threshold.json",
}


def run_job(job, ckpt_dir, budget_min):
    """Run one job inside budget_min minutes; returns (exit code, result dict)."""
    if job not in JOBS:
        print(f"job {job!r} not implemented; implemented: {sorted(JOBS)}", file=sys.stderr)
        return 2, {}
    Path(ckpt_dir).mkdir(parents=True, exist_ok=True)
    budget_s = max(60, (budget_min - SAFETY_MARGIN_MIN) * 60)
    out = JOBS[job](ckpt_dir, budget_s)
    if out.get("partial"):
        (Path(ckpt_dir) / f"{job}.partial.json").write_text(json.dumps(out))
        return EX_PARTIAL, out
    return 0, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job")
    ap.add_argument("--ckpt-dir", default="ckpt")
    ap.add_argument("--budget-min", type=int, default=300)
    a = ap.parse_args()
    code, out = run_job(a.job, a.ckpt_dir, a.budget_min)
    print(json.dumps(out))
    sys.exit(code)


if __name__ == "__main__":
    main()
