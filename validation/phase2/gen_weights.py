#!/usr/bin/env python3
"""Generate a zero-filled nebula .wgh for a single connected-layer GEMM net.
Format (per nebula connected_layer::init_weight): bias[N] floats then weight[K*N]
floats, concatenated per layer. Values are irrelevant for timing runs."""
import struct, sys, os

def main():
    if len(sys.argv) != 4:
        print("usage: gen_weights.py K N out.wgh"); sys.exit(1)
    K, N, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    n_floats = N + K*N
    os.makedirs(os.path.dirname(out), exist_ok=True) if os.path.dirname(out) else None
    with open(out, "wb") as f:
        f.write(b"\x00" * (4*n_floats))
    print(f"{out}: {n_floats} floats ({4*n_floats} bytes) for K={K} N={N}")

if __name__ == "__main__":
    main()
