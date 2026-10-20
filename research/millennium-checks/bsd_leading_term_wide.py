#!/usr/bin/env python3
"""Wider coefficient search for bsd_leading_term.py: a2 in -2..2, a4 in -80..80, a6 in -160..160, prime conductor 11..700.
Same analysis (analyse) as bsd_leading_term.py; writes bsd_leading_term_wide.json and lists any curve whose predicted |Sha| is not a perfect square.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bsd_leading_term as bl  # noqa: E402

bsd = bl.bsd


def main():
    found, seen = [], set()
    for a2 in range(-2, 3):
        for a4 in range(-80, 81):
            for a6 in range(-160, 161):
                _, _, _, delta, c4 = bl.invariants(a2, a4, a6)
                if delta == 0:
                    continue
                pk = bl.prime_power(abs(delta))
                if not pk or c4 % pk[0] == 0 or pk[0] > 700 or pk[0] < 11:
                    continue
                p, k = pk
                key = (p, tuple(bsd.ap((0, a2, 1, a4, a6), q) for q in (2, 3, 5, 7, 13)))
                if key in seen:
                    continue
                seen.add(key)
                r = bl.analyse(a2, a4, a6, p, k)
                if r:
                    found.append(r)
                    print(r["N"], r["model"], round(r["sha_pred"], 6), r["sha_is_square"], flush=True)
    found.sort(key=lambda r: (r["N"], r["model"]))
    bad = [r for r in found if not r["sha_is_square"]]
    print(f"{len(found)} rank-0 prime-conductor curves; {len(bad)} with non-square-integer Sha_pred")
    for r in bad:
        print("NOT SQUARE:", r["N"], r["model"], r["sha_pred"], "tors_bound", r["torsion_bound"], "c_p", r["c_p"])
    Path(__file__).with_name("bsd_leading_term_wide.json").write_text(json.dumps({"n_curves": len(found), "n_not_square": len(bad), "curves": found}, indent=1))


if __name__ == "__main__":
    main()
