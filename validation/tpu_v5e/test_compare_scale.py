#!/usr/bin/env python3
"""compare_v5e.scale_time on mocked run times (no simulator, no TPU): row tiles, strip interpolation and
extrapolation, the once-per-kernel layer setup, and the chain overlap."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_v5e as c

SETUP = 0.22                                  # a tile on the full config = inner tile + layer setup
sim = {("a", "s2"): (20.0 + SETUP, 3.0, 5.0), ("a", "s4"): (40.0 + SETUP, 3.0, 9.0),
       ("i", "s2"): (20.0, 3.0, 5.0), ("i", "s4"): (40.0, 3.0, 9.0),
       ("a", "s1"): (11.0, 2.0, 1.5), ("a", "one"): (7.0, 4.0, 6.0)}
near = lambda x, y: abs(x - y) < 1e-9

# one run, one tile
assert near(c.scale_time(sim, 1, ([("a", "one")], None), None), 7.0)
# strips: exact at the simulated counts, linear beyond (10 per strip) and between 1 and 2
first = ([("a", "s2"), ("a", "s4")], ((2, 4), 4))
assert near(c.scale_time(sim, 1, first, None), 40.0 + SETUP)
assert near(c.scale_time(sim, 1, ([("a", "s2"), ("a", "s4")], ((2, 4), 56)), None), 560.0 + SETUP)
assert near(c.scale_time(sim, 1, ([("a", "s1"), ("a", "s2")], ((1, 2), 1.5)), None), (11.0 + 20.0 + SETUP) / 2)
# row tiles: the layer setup once, not mt times
inner = ([("i", "s2"), ("i", "s4")], ((2, 4), 4))
assert near(c.scale_time(sim, 4, first, inner), 4 * 40.0 + SETUP)
assert not near(c.scale_time(sim, 4, first, inner), 4 * (40.0 + SETUP))
# chain: minus min(prologue, tail) of the first run, once
assert near(c.scale_time(sim, 1, ([("a", "one")], None), None, chain=True), 7.0 - 4.0)
assert near(c.scale_time(sim, 4, first, inner, chain=True), 4 * 40.0 + SETUP - 3.0)
# timing provenance of a GEMM csv: a guess is marked, an unrecognizable file is unknown
tk = c.timing_kind
assert tk([dict(wall_us="58.2", device_us="13.6", ops="-|fusion=0.81")]) == "device_trace(inferred)"
assert tk([dict(wall_us="58.2", device_us="13.6", ops="")]) == "unknown"              # differs, but no trace ops
assert tk([dict(wall_us="68.44", device_us="68.44"), dict(wall_us="45.65", device_us="45.65")]) == "loop_slope(inferred)"
assert tk([dict(wall_us="58.2", device_us="")]) == "wall"
assert tk([dict(device_us="36.31")]) == "unknown"                                      # device_us alone (gemm_int8.csv)
assert tk([dict(wall_us="1", device_us="1"), dict(wall_us="2", device_us="")]) == "unknown"   # partly missing
# the stored files
import csv
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
stored = lambda f: tk(list(csv.DictReader(open(os.path.join(R, f)))))
assert stored("2026-11-03_colab_v5e1_followup2/gemm.csv") == "loop_slope(inferred)"
assert stored("2026-11-03_colab_v5e1_followup2/gemm_int8.csv") == "unknown"
assert stored("2026-11-04_colab_v5e1_smallk2/gemm_oddchain.csv") == "unknown"
assert stored("2026-11-03_colab_v5e1_rot/gemm.csv") == "device_trace(inferred)"
# an explicit --timing wins and is one of the stated kinds
assert set(c.TIMING_KINDS) == {"device_trace", "device_slope", "loop_slope", "wall"}
print("ok")
