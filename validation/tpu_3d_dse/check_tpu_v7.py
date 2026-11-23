#!/usr/bin/env python3
"""configs/accelerators/tpu_v7.cfg against the published TPU v7 rooflines (2307 TFLOPS bf16, 7380 GB/s HBM): three
bf16 GEMMs on the analytical DRAM, isolated-program scope (input from DRAM, output written back),
program_setup_cycle 0. A sanity check of the parameterisation -- peak rate, HBM rate, MAC count, DRAM payload --
NOT a fidelity result: no TPU v7 measurement exists. The layer reports are copied to OUT_DIR.
Usage: python3 validation/tpu_3d_dse/check_tpu_v7.py OUT_DIR      (from the repository root; needs models/model)"""
import sys, os, re, subprocess, shutil
sys.path.insert(0, "validation/palm"); import palm_decode as palm
palm.CLOCK_HZ, palm.DIM, palm.MXUS, palm.MOSAIC_DRAM = 2.2e9, 256, 8, None
PEAK, HBM = 2307e12, 7380e9
OUT = os.path.abspath(sys.argv[1]); os.makedirs(OUT, exist_ok=True)
for M, K, N in ((512, 4096, 4096), (32, 8192, 8192), (128, 4096, 4096)):
    tmp, cfg, rund = palm.workdir(("tpu_v7",))
    c = os.path.join(cfg, "accelerators", "tpu_v7.cfg"); t = open(c).read()
    open(c, "w").write(re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", t))
    palm.gemm_files(cfg, rund, "g", "tpu_v7", M, K, N)
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfg); env.pop("NPUSIM_MOSAIC_DRAM", None)
    p = subprocess.run([os.path.join(palm.ROOT, "models", "model"), "run", "tpu_v7", "g", "ws"], cwd=rund, env=env, capture_output=True, text=True)
    if p.returncode:
        print("FAILED", M, K, N, (p.stdout + p.stderr)[-600:]); continue
    L = open(os.path.join(rund, "tpu_v7_g_layer_0.txt")).read()
    f = lambda rx: float(re.search(rx, L).group(1))
    crit, mac = f(r"Critical-path latency\s*:\s*([\d.]+)"), f(r"Busy scalar-MAC cycles\s*:\s*([\d.]+)")
    B = L[L.index("Busy cycles (ratio"):]
    dram = re.search(r"\* DRAM\s*:\s*[\d.]+ cycles \(\s*([\d.]+) %\)", B).group(1)
    comp_busy = re.search(r"\* PE \(compute\+LB\)\s*:\s*[\d.]+ cycles \(\s*([\d.]+) %\)", B).group(1)
    bott = re.search(r"Bottleneck level\s*:\s*(.*)", L).group(1).strip()
    blk = L[L.index("Multi-chip <-> DRAM transactions"):][:400]
    tx = [int(re.search(rf"\* {k}\s*:\s*(\d+)", blk).group(1)) * 4096 for k in ("Input data", "Weight", "Output data")]
    comp, mem = 2 * M * K * N / PEAK * 1e6, (2 * M * K + 2 * K * N + 2 * M * N) / HBM * 1e6
    us = crit / 2.2e9 * 1e6
    print(f"{M}x{K}x{N}: NPUsim {us:8.2f} us | roofline compute {comp:7.2f} us, memory {mem:7.2f} us -> x{us / max(comp, mem):.3f} of the larger | "
          f"MACs ok {mac == M * K * N} | DRAM payload {tx} expect {[2 * M * K, 2 * K * N, 2 * M * N]} | busy DRAM {dram}%, compute {comp_busy}% | bottleneck {bott}")
    shutil.copy(os.path.join(rund, "tpu_v7_g_layer_0.txt"), os.path.join(OUT, f"{M}x{K}x{N}_layer_0.txt"))
    shutil.rmtree(tmp)
