"""Print the next unfinished Millennium compute job (first entry of compute_queue.json whose result file is not on origin/results).

Used by .github/workflows/continuous-research.yml. Prints nothing when every queued job has a result. Exit 0 either way;
exit 2 when the queue names a job that ci_run_job.RESULTS does not know (a typo would otherwise stall silently).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "vm"))

from ci_run_job import RESULTS  # noqa: E402
from vm_worker import next_job  # noqa: E402


def finished_files(ref: str = "origin/results") -> set[str]:
    """Repo-relative paths on `ref`; empty when the branch does not exist yet."""
    out = subprocess.run(["git", "ls-tree", "-r", "--name-only", ref], cwd=ROOT, capture_output=True, text=True).stdout
    return set(out.split())


def main() -> int:
    queue = json.loads((ROOT / "compute_queue.json").read_text())["jobs"]
    try:
        job = next_job(queue, finished_files(), RESULTS)
    except KeyError as err:
        print(f"::error::{err}", file=sys.stderr)
        return 2
    if job:
        print(job)
    return 0


if __name__ == "__main__":
    sys.exit(main())
