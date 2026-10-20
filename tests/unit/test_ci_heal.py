"""scripts/ci_heal.py: the policy that cancels stuck runs, re-runs a failure once and stops after repeated failures."""
from __future__ import annotations

import importlib.util
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_spec = importlib.util.spec_from_file_location("ci_heal", Path(__file__).resolve().parents[2] / "scripts" / "ci_heal.py")
ch = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ch)

NOW = datetime(2026, 10, 20, 12, 0, tzinfo=timezone.utc)


def run(rid: int, status: str, conclusion: str | None = None, hours_ago: float = 1, attempt: int = 1) -> dict:
    created = (NOW - timedelta(hours=hours_ago)).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {"databaseId": rid, "status": status, "conclusion": conclusion, "createdAt": created, "attempt": attempt}


def test_healthy_history_needs_nothing():
    d = ch.decide([run(3, "completed", "success"), run(2, "completed", "failure", 5)], NOW, 4 * 3600)
    assert d == {"cancel": [], "rerun": None, "blocked": False, "failures": 0}


def test_a_stuck_run_is_cancelled_but_a_recent_one_is_not():
    d = ch.decide([run(2, "in_progress", hours_ago=6), run(1, "queued", hours_ago=0.5)], NOW, 4 * 3600)
    assert d["cancel"] == [2]


def test_first_failure_is_rerun_once_only():
    assert ch.decide([run(5, "completed", "failure")], NOW, 3600)["rerun"] == 5
    assert ch.decide([run(5, "completed", "failure", attempt=2)], NOW, 3600)["rerun"] is None  # already re-run
    assert ch.decide([run(6, "in_progress"), run(5, "completed", "failure")], NOW, 3600)["rerun"] is None  # something is running


def test_three_failures_in_a_row_block_and_a_success_lifts_it():
    fails = [run(i, "completed", "failure", hours_ago=10 - i, attempt=2) for i in (9, 8, 7)]
    assert ch.decide(fails, NOW, 3600) | {} == {"cancel": [], "rerun": None, "blocked": True, "failures": 3}
    assert ch.decide([run(10, "completed", "success"), *fails], NOW, 3600)["blocked"] is False


def test_cancelled_and_skipped_runs_do_not_count_as_failures():
    runs = [run(4, "completed", "cancelled"), run(3, "completed", "skipped"), run(2, "completed", "success")]
    assert ch.decide(runs, NOW, 3600)["failures"] == 0


def test_timed_out_counts_as_a_failure():
    assert ch.decide([run(2, "completed", "timed_out")], NOW, 3600)["failures"] == 1


def test_heal_applies_the_policy_with_gh_calls(monkeypatch):
    calls = []

    def fake_gh(*args):
        calls.append(args)

        class R:
            returncode = 0
            stdout = "[]"
        return R()

    monkeypatch.setattr(ch, "_gh", fake_gh)
    monkeypatch.setattr(ch, "_runs", lambda repo, wf: [run(7, "in_progress", hours_ago=9)] if wf == "millennium-compute.yml"
                        else [run(8, "completed", "failure")])
    report = ch.heal("o/r", now=NOW)
    assert ("run", "cancel", "7", "-R", "o/r") in calls
    assert ("run", "rerun", "8", "-R", "o/r", "--failed") in calls
    assert report["schur-sweep"]["rerun"] == 8 and report["millennium-compute"]["cancel"] == [7]
