"""Self-healing for the compute workflows: cancel stuck runs, re-run a failed run once, stop after repeated failures.

Called by .github/workflows/continuous-research.yml every 6 hours (and by Grove's daily watchdog from outside GitHub).
For each workflow it reads the recent runs and applies, in order:
  1. a run queued or in progress for longer than the workflow's limit is cancelled (the next tick starts it again; the
     Millennium jobs resume from their cached checkpoint);
  2. when the newest finished run failed and has had no re-run yet, re-run its failed jobs (transient runner, network or
     install errors); this happens at most once per run;
  3. after MAX_FAILURES consecutive failed runs it stops dispatching that workflow ("blocked") and opens one issue; the
     block lifts by itself when a run succeeds;
  4. an open "ended failure" issue is closed once the newest finished run succeeded.
Prints `block_<workflow>=1` lines for the dispatch step (also appended to $GITHUB_OUTPUT when set). Exit 0 always: healing
must not turn a healthy driver run red.

Usage: python scripts/ci_heal.py --repo OWNER/NAME [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

WORKFLOWS = {  # workflow file -> seconds a run may sit queued or in progress before it counts as stuck
    "schur-sweep.yml": 4 * 3600,
    "millennium-compute.yml": 7 * 3600,
}
MAX_FAILURES = 3
ACTIVE = {"queued", "in_progress", "waiting", "pending", "requested"}
FAILED = {"failure", "timed_out", "startup_failure"}


def _age(created: str, now: datetime) -> float:
    return (now - datetime.fromisoformat(created.replace("Z", "+00:00"))).total_seconds()


def decide(runs: list[dict], now: datetime, stuck_after: int, max_failures: int = MAX_FAILURES) -> dict:
    """Pure policy: `runs` newest first (gh run list JSON). Returns {"cancel": [ids], "rerun": id|None, "blocked": bool, "failures": n}."""
    cancel = [r["databaseId"] for r in runs if r["status"] in ACTIVE and _age(r["createdAt"], now) > stuck_after]
    live = [r for r in runs if r["status"] in ACTIVE and r["databaseId"] not in cancel]
    finished = [r for r in runs if r["status"] == "completed" and r.get("conclusion") not in ("cancelled", "skipped", None)]
    failures = 0
    for r in finished:
        if r["conclusion"] in FAILED:
            failures += 1
        else:
            break
    newest = finished[0] if finished else None
    rerun = None
    if newest and newest["conclusion"] in FAILED and newest.get("attempt", 1) == 1 and failures < max_failures and not live:
        rerun = newest["databaseId"]
    return {"cancel": cancel, "rerun": rerun, "blocked": failures >= max_failures, "failures": failures}


def _gh(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["gh", *args], capture_output=True, text=True)


def _runs(repo: str, workflow: str) -> list[dict]:
    out = _gh("run", "list", "-R", repo, "--workflow", workflow, "--limit", "20", "--json",
              "databaseId,status,conclusion,createdAt,attempt")
    return json.loads(out.stdout) if out.returncode == 0 and out.stdout.strip() else []


def _issue(repo: str, title: str, body: str) -> None:
    found = _gh("issue", "list", "-R", repo, "--state", "open", "--search", f'"{title}" in:title', "--json", "number")
    if found.returncode == 0 and not json.loads(found.stdout or "[]"):
        _gh("issue", "create", "-R", repo, "--title", title, "--body", body)


def _close_issues(repo: str, title_prefix: str, why: str) -> None:
    found = _gh("issue", "list", "-R", repo, "--state", "open", "--search", f'"{title_prefix}" in:title', "--json", "number,title")
    for it in json.loads(found.stdout or "[]") if found.returncode == 0 else []:
        if it["title"].startswith(title_prefix):
            _gh("issue", "close", str(it["number"]), "-R", repo, "--comment", why)


def heal(repo: str, dry_run: bool = False, now: datetime | None = None) -> dict[str, dict]:
    """Apply the policy to every workflow; returns the decision per workflow."""
    now = now or datetime.now(timezone.utc)
    report = {}
    for wf, limit in WORKFLOWS.items():
        name = wf.removesuffix(".yml")
        runs = _runs(repo, wf)
        d = decide(runs, now, limit)
        report[name] = d
        for rid in d["cancel"]:
            print(f"{name}: run {rid} stuck > {limit // 3600} h: cancelling")
            if not dry_run:
                _gh("run", "cancel", str(rid), "-R", repo)
        if d["rerun"]:
            print(f"{name}: run {d['rerun']} failed: re-running its failed jobs once")
            if not dry_run:
                _gh("run", "rerun", str(d["rerun"]), "-R", repo, "--failed")
        if d["blocked"]:
            print(f"{name}: {d['failures']} failed runs in a row: dispatch blocked until one succeeds")
            if not dry_run:
                _issue(repo, f"{name} failed {d['failures']} times in a row",
                       f"The driver stopped starting {name}. Read the latest runs at https://github.com/{repo}/actions/workflows/{wf}; "
                       "the block lifts by itself after one successful run.")
        elif d["failures"] == 0 and not dry_run:
            _close_issues(repo, f"{name} failed", "The latest run succeeded; closing.")
            _close_issues(repo, f"{name} ended", "The latest run succeeded; closing.")
        print(f"block_{name.replace('-', '_')}={1 if d['blocked'] else 0}")
    out = os.environ.get("GITHUB_OUTPUT")
    if out:
        with open(out, "a", encoding="utf-8") as fh:
            for name, d in report.items():
                fh.write(f"block_{name.replace('-', '_')}={1 if d['blocked'] else 0}\n")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    try:
        heal(a.repo, a.dry_run)
    except Exception as err:  # noqa: BLE001 - healing must never fail the driver
        print(f"::warning::ci_heal skipped: {err}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
