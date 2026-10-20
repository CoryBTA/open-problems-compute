# open-problems-compute

Reproducible computations for open problems, run on free CI. Nothing here is a new theorem: every claim carries a label
and a scope, and every search runs its known-answer cases first.

| Directory | What it is |
|---|---|
| `research/schur-structured` | Schur-number searches in the Fredricksen-Sweet symmetric class (S(6) >= 536 by extension). Sharded sweeps with a calibration gate. |
| `research/ns-resolution` | Navier-Stokes resolution studies (checkpointed, resumable). |
| `research/millennium-checks` | Small numerical checks (BSD, Riemann-type, SAT thresholds). |
| `research/lean-millennium` | Lean 4 statements and lemmas; CI fails on `sorry`. |
| `research/open-problems/schur.py` | Independent Schur-partition checker used to re-verify witnesses. |

## Workflows
- `continuous-research` (every 6 h): starts the next unfinished job from `compute_queue.json`, else the next item of
  `research/schur-structured/sweep_queue.json`; one job at a time; opens a deduplicated issue on failure or an empty queue.
- `schur-sweep`: calibration (n = 536 must be found and verified) then 16 shards; the label is
  `negative (this class only)` only if every shard finished with no timeout, otherwise `inconclusive`.
- `millennium-compute`: long numerical jobs that stop themselves inside the time limit, save a checkpoint and re-dispatch.
- `lean`: builds the Lean package and rejects `sorry`.

Results are appended to the `results` branch (never forced, never `main`). Labels: `reproduced`, `negative (class only)`,
`inconclusive`. A witness found by a search is not a record until `research/open-problems/schur.py` re-checks it and it is
compared with the literature.

License: MIT.
