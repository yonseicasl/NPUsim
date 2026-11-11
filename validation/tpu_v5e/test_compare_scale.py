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
print("ok")
