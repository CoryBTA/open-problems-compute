"""research/schur-structured: the symmetric (Fredricksen-Sweet) class reaches the known Schur numbers and verifies."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_spec = importlib.util.spec_from_file_location("schur_structured", Path(__file__).resolve().parents[2] / "research" / "schur-structured" / "search.py")
ss = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ss)


def test_symmetric_class_holds_an_s4_witness():
    """n+1 = 45 is divisible by 3, so 15 and 30 stay unpaired; without that exception the class is empty (L#41235)."""
    w = ss.solve({"n": 44, "k": 4, "symmetric": True})
    assert isinstance(w, list) and ss.verify(w, 4) and ss.is_symmetric(w)
    assert ss.solve({"n": 45, "k": 4, "symmetric": True}) is None


@pytest.mark.timeout(60)
@pytest.mark.slow
def test_s5_160_found_in_symmetric_depth_class():
    w = ss.solve({"n": 160, "k": 5, "symmetric": True, "min_depth": 44})
    assert isinstance(w, list) and ss.verify(w, 5) and ss.is_symmetric(w) and ss.depth(w, 5) >= 44


def test_verifier_rejects_a_sum():
    assert not ss.verify([0, 0, 1, 1], 2)  # 1 + 1 = 2 inside set 0
    assert ss.verify([0, 1, 1, 0], 2)


def test_five_set_partitions_of_160_are_valid_and_distinct():
    """run6b enumerates symmetric depth-44 partitions of [1,160]; every one verifies, none repeats (the paper counts 5,840)."""
    spec = importlib.util.spec_from_file_location("run6b", Path(ss.__file__).with_name("run6b.py"))
    import sys
    sys.path.insert(0, str(Path(ss.__file__).parent))
    try:
        run6b = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(run6b)
    finally:
        sys.path.pop(0)
    parts = run6b.five_set_parts(20)
    assert len(parts) == 20 and len({tuple(p) for p in parts}) == 20
    assert all(ss.verify(p, 5) and ss.is_symmetric(p) and ss.depth(p, 5) >= 44 for p in parts)


@pytest.mark.timeout(60)
@pytest.mark.slow
def test_shards_cover_every_start_exactly_once(tmp_path):
    """The Actions matrix splits run6b's starts with --shard I/S; the shard sizes must add up to the unsharded count."""
    import json
    import subprocess
    import sys
    script = Path(ss.__file__).with_name("run6b.py")

    def starts(*extra: str) -> int:
        out = tmp_path / f"s{len(extra)}{abs(hash(extra))}.json"
        subprocess.run([sys.executable, str(script), "536", "--parts", "4", "--multiply", "--workers", "2", "--budget", "1000",
                        "--summary", str(out), *extra], check=True, cwd=script.parent, capture_output=True, timeout=50)
        return json.loads(out.read_text())["starts"]

    whole = starts()
    assert whole > 0
    assert starts("--shard", "0/3") + starts("--shard", "1/3") + starts("--shard", "2/3") == whole


def test_merge_labels_negative_only_when_every_shard_finished_without_timeout(tmp_path):
    spec = importlib.util.spec_from_file_location("merge_shards", Path(ss.__file__).with_name("merge_shards.py"))
    ms = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ms)

    def shard(i: int, **stats: int) -> None:
        (tmp_path / f"shard-{i}.json").write_text(__import__("json").dumps(
            {"n": 537, "starts": 10, "stats": {"sat": 0, "unsat": 10, "timeout": 0, **stats}}))

    shard(0)
    shard(1)
    assert ms.merge(tmp_path, 2, 537)["label"] == "negative (this class only)"
    assert ms.merge(tmp_path, 3, 537)["label"] == "inconclusive"  # a shard never reported
    shard(1, unsat=9, timeout=1)
    assert ms.merge(tmp_path, 2, 537)["label"] == "inconclusive"
    shard(1, unsat=9, sat=1)
    assert ms.merge(tmp_path, 2, 537)["label"].startswith("witness-found")


def test_finish_makes_a_missing_result_inconclusive_and_records_prefix(tmp_path):
    import json
    spec = importlib.util.spec_from_file_location("merge_shards", Path(ss.__file__).with_name("merge_shards.py"))
    ms = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ms)
    out = tmp_path / "result.json"
    ms.finish(out, 537, 100, "42")
    row = json.loads(out.read_text())
    assert row["label"] == "inconclusive" and row["prefix"] == 100 and row["run"] == "42"
    out.write_text(json.dumps({"n": 537, "label": "negative (this class only)"}))
    ms.finish(out, 537, 100, "43")
    assert json.loads(out.read_text())["label"] == "negative (this class only)"
