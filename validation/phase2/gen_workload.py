#!/usr/bin/env python3
"""Generate NPUsim net/map/weights for a GEMM point C[MxN] = A[MxK] x B[KxN]
mirroring Gemmini's 16x16 WS array.

- net: batch=M, input h*w*1 = K (grayscale spatial encoding), connected output=N.
- mapping: WS, 16x16 filled with PE_X=K16 / PE_Y=C16, B kept inner (PE temporal)
  so weights are reused across the batch like Gemmini's WS schedule.
- weights: zero-filled bias[N] + weight[K*N] floats.
Usage: gen_workload.py M K N  ->  net gemm_MxKxN, mapping ws.map
"""
import os, sys

DIM = 16
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."

def factor_hw(k):
    h = 1
    for cand in (16, 8, 32, 4, 64, 2, 128):
        if k % cand == 0 and k // cand <= 4096:
            return cand, k // cand
    return 1, k

def emit(m, k, n):
    name = f"gemm_{m}x{k}x{n}"
    assert m % 1 == 0 and k % DIM == 0 and n % DIM == 0, "K,N must be multiples of 16"
    h, w = factor_hw(k)

    # ---- net config ----
    net = f"""# GEMM {m}x{k}x{n} validation net (Gemmini-matched).
[net]
height={h}
width={w}
channels=1
batch={m}

[data]
dataset=imagenet
test=datasets/gemm/test.lst
labels=datasets/gemm/labels.lst
weight=weights/{name}.wgh
top=1

[connected]
output={n}
activation=linear

[softmax]
groups=1

[cost]
type=l2
"""
    open(f"{ROOT}/configs/networks/{name}.cfg", "w").write(net)

    # ---- mapping: totals K(col)=n, B=m, C=k ----
    # spatial: PE_X K=16, PE_Y C=16. B fully inner at PE (V3 edge_accumulation:
    # outputs accumulate at the array edge, so weights stay resident across the
    # whole batch like Gemmini's WS schedule). Per-PE input = b_pe bytes <= 512B LB.
    # per-PE K tile: bounded by n/DIM and by the accumulator slab m*k_pe*DIM <= 64K
    k_pe = max(1, min(DIM, n // DIM, max(1, 65536 // (m * DIM))))
    b_pe = m                           # full batch per PE (edge accumulation)
    b_glb = m // b_pe
    k_glb = n // (k_pe * DIM)        # remaining K after PE(16) x PE_X(16)
    c_glb = k // DIM                 # remaining C after PE_Y(16)
    assert k_pe * DIM * max(1, k_glb) * 1 or True
    # The GLB output partition (64KB) models Gemmini's edge accumulator; when the
    # output slab exceeds it, the outer K loop spills to DRAM (RTL likewise mvouts C
    # per N-strip). Fold count (N*K/256) is invariant to where the K loop lives.
    k_glb_glb, k_dram = max(1, k_glb), 1
    while m * k_pe * DIM * k_glb_glb > 65536 and k_glb_glb % 2 == 0:
        k_glb_glb //= 2
        k_dram *= 2
    rows = {
        "MAC":     [1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        "PE":      [k_pe, b_pe, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        "PE_X":    [DIM, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        "PE_Y":    [1, 1, 1, 1, DIM, 1, 1, 0, 0, 1, 1],
        "GLB":     [k_glb_glb, b_glb, 1, 1, c_glb, 1, 1, 0, 0, 1, 1],
        "CHIPS_X": [1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        "CHIPS_Y": [1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        "DRAM":    [k_dram, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
    }
    # sanity: totals
    tot_k = k_pe * DIM * max(1, k_glb)
    tot_b = b_pe * b_glb
    tot_c = DIM * c_glb
    assert tot_k == n and tot_b == m and tot_c == k, (tot_k, tot_b, tot_c)

    os.makedirs(f"{ROOT}/configs/mappings/gemmini/{name}", exist_ok=True)
    with open(f"{ROOT}/configs/mappings/gemmini/{name}/ws.map", "w") as f:
        f.write("# K, B, P, Q, C, R, S, H, W, GROUP, STRIDE\n\n[connected]\n")
        f.write(f"# {n}, {m}, 1, 1, {k}, 1, 1, 0, 0, 1, 1,\n")
        for lv, r in rows.items():
            f.write(f"{lv:7s} = " + ", ".join(str(x) for x in r) + ",\n")

    # ---- weights ----
    n_floats = n + k * n
    with open(f"{ROOT}/models/weights/{name}.wgh", "wb") as f:
        f.write(b"\x00" * (4 * n_floats))
    print(f"{name}: h*w={h}x{w}, PE(B{b_pe},K{k_pe}), GLB(K{max(1,k_glb)},B{b_glb},C{c_glb})")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        emit(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
    else:
        for m, k, n in [(64,64,64),(128,128,128),(256,256,256),(512,512,512),
                        (16,512,512),(512,512,64)]:
            emit(m, k, n)
