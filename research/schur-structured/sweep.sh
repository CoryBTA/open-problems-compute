#!/bin/sh
# Sweep run6b.py over n (class: symmetric, depth >= 161, first 160 entries = unit multiple of a depth-44 symmetric 5-set partition of 160).
for n in $(seq "$1" "$2"); do
  echo "n=$n"; python research/schur-structured/run6b.py "$n" --parts 100000 --multiply --workers 24 --budget 100000 | tail -2
done
