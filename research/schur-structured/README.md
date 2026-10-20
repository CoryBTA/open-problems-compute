# schur-structured

Schur partitions in the symmetric class of Fredricksen and Sweet (EJC 7, 2000, #R32): i and (n+1)-i share a set, except
that when 3 divides n+1 the points (n+1)/3 and 2(n+1)/3 stay unpaired. "Depth >= m" means one set has no element below m.
Every witness is re-checked by `verify()` (independent of the SAT encoding). Nothing here is a new theorem.

| Problem | Script | Result | Label |
|---|---|---|---|
| S(2)=4, S(3)=13, S(4)=44, unrestricted | `search.py` | n and n+1 as expected (n=45 UNSAT, 13 s) | reproduced |
| Symmetric class, k=4 | `search.py` | symmetric 44 exists and verifies; symmetric 45 UNSAT | reproduced |
| S(5)=160, symmetric + depth >= 44 | `search.py` | verified witness in 0.7 s. The unstructured search in `research/open-problems` timed out at 1550 s | reproduced |
| S(6) >= 536, symmetric + depth >= 161, whole instance | `run6.py` | CaDiCaL timed out at n = 499 .. 544 (2800 s each) | timeout |
| S(6) >= 536 by extension | `run6b.py` | fixing the first 160 entries to a unit multiple of a depth-44 symmetric 5-set partition (5,840 base partitions, the exact count in Fredricksen-Sweet; 385,440 with multiples): verified witness for n=536 in 17 s, re-checked by `research/open-problems/schur.py` | reproduced |
| n = 537..600 by extension | `run6b.py`, `sweep.sh` | all 385,440 starts UNSAT for each of the 64 values of n, no timeouts | negative (this class only) |
| Prefix 120 instead of 160 (n = 536 gate, then 537..544) | `run6b.py --prefix 120` | n=536 SAT (verified, 51 s); n=537..544 all starts UNSAT, no timeouts | negative (this class only) |

## Corrections
- L#41211 claimed a palindromic constraint makes S(4)=44 UNSAT. That test omitted the 3 | n+1 exception; with it the class holds an S(4) witness (L#41235 supersedes it).
- S(5)=160 is exact (Heule 2017), so S(5) is a reproduction, not an open problem.

## Not claimed
- The negative covers only partitions whose first 160 entries are a unit multiple of a depth-44 symmetric 5-set partition of 160. It says nothing about other 5-set partitions of 160 or non-symmetric 6-set partitions.
- Prefix 120 is the first 120 entries of a 160-partition, which is a 5-set partition of [1, 120] but not every one. Bases not coming from a 160-partition are untried: enumerating all 5-set partitions of [1, m] for m < 160 is large.
- Anything about S(6) beyond the rows above. The published record is 536 (Fredricksen and Sweet 2000); a larger n is a new
  record only if it is verified, re-checked by an independent script, and compared against the current literature (Ageron et al. 2021 and later).
