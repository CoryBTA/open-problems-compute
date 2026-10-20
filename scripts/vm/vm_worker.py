"""Always-on compute worker for a free VM (Oracle Cloud always-free ARM, 4 cores / 24 GB, no job time limit).

Loop forever: fetch origin, take the first job in compute_queue.json whose result file is not yet on origin/results or
origin/main, run it through scripts/ci_run_job.py (checkpoint in CKPT_DIR on disk, so a reboot or kill resumes it),
publish the result file to branch `results` (fast-forward push, never main or cloud), repeat. With nothing to do it sleeps
and re-checks, so adding a job to compute_queue.json on `cloud` is enough to start it.

Run under systemd (scripts/vm/millennium-worker.service). Needs: a clone of the repo on branch `main`, a git SSH deploy key
with write access, python deps from requirements.txt + requirements-extra.txt.
Env: MILLENNIUM_CKPT_DIR (default ~/millennium-ckpt), MILLENNIUM_IDLE_S (default 3600), MILLENNIUM_SEGMENT_MIN (default 600).
"""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
EX_PARTIAL = 75


def next_job(queue, finished_files, results):
    """First job of `queue` whose result file (results[job]) is not in `finished_files`; None when all are done."""
    for job in queue:
        if job not in results:
            raise KeyError(f"{job!r} is in compute_queue.json but has no entry in ci_run_job.RESULTS")
        if results[job] not in finished_files:
            return job
    return None


def git(*args, cwd=None, check=True):
    """Run a git command (in ROOT unless cwd is given; ROOT is read at call time) and return its stdout."""
    return subprocess.run(["git", *args], cwd=cwd or ROOT, check=check, capture_output=True, text=True).stdout


def files_on(ref):
    """Set of repo-relative paths in `ref` (empty if the ref does not exist)."""
    out = git("ls-tree", "-r", "--name-only", ref, check=False)
    return set(out.split())


def publish(result_relpath, job, results_wt):
    """Commit the result file onto branch `results` (own worktree) and push it, retrying when the remote moved."""
    if not results_wt.exists():
        git("worktree", "add", "--detach", str(results_wt), "origin/results" if files_on("origin/results") else "origin/main")
    for attempt in range(4):
        git("fetch", "origin")
        git("checkout", "-f", "--detach", "origin/results" if files_on("origin/results") else "origin/main", cwd=results_wt)
        dest = results_wt / result_relpath
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / result_relpath, dest)
        git("add", "--", result_relpath, cwd=results_wt)
        git("-c", "user.name=vm-worker", "-c", "user.email=vm-worker@users.noreply.github.com",
            "commit", "-m", f"compute results: {job} (vm worker)", cwd=results_wt, check=False)
        pushed = subprocess.run(["git", "push", "origin", "HEAD:refs/heads/results"], cwd=results_wt, capture_output=True, text=True)
        if pushed.returncode == 0:
            return True
        print(f"push rejected (attempt {attempt + 1}): {pushed.stderr.strip()[-200:]}", flush=True)
        time.sleep(10)
    return False


def main():
    import ci_run_job

    ckpt = Path(os.environ.get("MILLENNIUM_CKPT_DIR", Path.home() / "millennium-ckpt"))
    idle = int(os.environ.get("MILLENNIUM_IDLE_S", "3600"))
    segment = int(os.environ.get("MILLENNIUM_SEGMENT_MIN", "600"))
    results_wt = ROOT.parent / "millennium-results-wt"
    while True:
        git("fetch", "origin")
        git("merge", "--ff-only", "origin/main", check=False)  # pick up queue and code changes; a local edit blocks it, never loses work
        queue = json.loads((ROOT / "compute_queue.json").read_text())["jobs"]
        finished = files_on("origin/results") | files_on("origin/main")
        job = next_job(queue, finished, ci_run_job.RESULTS)
        if job is None:
            print(f"{time.strftime('%F %T')} queue empty; sleeping {idle}s", flush=True)
            time.sleep(idle)
            continue
        print(f"{time.strftime('%F %T')} running {job}", flush=True)
        code, out = ci_run_job.run_job(job, ckpt, segment)
        if code == EX_PARTIAL:
            print(f"{job}: segment done at t={out.get('t')}; continuing from the checkpoint", flush=True)
            continue
        rel = ci_run_job.RESULTS[job]
        if code == 0 and (ROOT / rel).exists():
            ok = publish(rel, job, results_wt)
            print(f"{job}: finished, published={ok}", flush=True)
            if git("ls-files", "--", rel):  # keep the working copy clean so the next fast-forward works
                git("checkout", "--", rel)
            else:
                (ROOT / rel).unlink()
        else:
            print(f"{job}: failed with code {code}; retrying in {idle}s", flush=True)
            time.sleep(idle)


if __name__ == "__main__":
    main()
