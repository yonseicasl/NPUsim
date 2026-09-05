#!/usr/bin/env python3
"""Generate a nebula .wgh for a single connected-layer GEMM net.
Format (per nebula connected_layer::init_weight): bias[N] floats then weight[K*N]
floats, concatenated per layer.

Default: zero-filled (values are irrelevant for timing runs). With a seed argument the
weights are deterministic pseudo-random in [-0.5, 0.5] and the biases stay zero -- the
FUNCTIONAL verification fixture needs nonzero weights, or the accelerator-vs-reference
comparison is vacuously all-zero on both sides."""
import random
import struct
import sys, os

def main():
    if len(sys.argv) not in (4, 5):
        print("usage: gen_weights.py K N out.wgh [random_seed]"); sys.exit(1)
    K, N, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    seed = int(sys.argv[4]) if len(sys.argv) == 5 else None
    n_floats = N + K*N
    os.makedirs(os.path.dirname(out), exist_ok=True) if os.path.dirname(out) else None
    with open(out, "wb") as f:
        if seed is None:
            f.write(b"\x00" * (4*n_floats))
        else:
            rng = random.Random(seed)
            f.write(b"\x00" * (4*N))    # biases: zero (the accelerator models raw MACs)
            f.write(struct.pack("<%df" % (K*N),
                                *(rng.uniform(-0.5, 0.5) for _ in range(K*N))))
    kind = "zero" if seed is None else f"random(seed={seed})"
    print(f"{out}: {n_floats} floats ({4*n_floats} bytes) for K={K} N={N}, {kind}")

if __name__ == "__main__":
    main()
