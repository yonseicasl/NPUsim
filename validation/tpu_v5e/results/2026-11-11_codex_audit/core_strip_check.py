import sys, os, re, subprocess
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, "validation/palm"); import palm_decode as palm
palm.CLOCK_HZ = 1.5e9; palm.MXUS = 4
tmp, cfg, run = palm.workdir(("tpu_v5e",))
p = os.path.join(cfg, "accelerators", "tpu_v5e.cfg"); txt = open(p).read().replace("program_setup_cycle = 1035", "program_setup_cycle = 0")
for knob in (0, 1):
    open(os.path.join(cfg, "accelerators", f"v5e_t{knob}.cfg"), "w").write(txt.replace("interlayer_overlap = 1", f"interlayer_overlap = 1\ntimeline_dram_tiles = {knob}"))
jobs = [(knob, M, K, s) for knob in (0, 1) for M, K in ((512, 2048),) for s in (1, 2, 4)]
def go(j):
    knob, M, K, s = j; name = f"t{M}_{K}_{s}"; acc = f"v5e_t{knob}"
    palm.gemm_files(cfg, run, name, acc, M, K, s * 512)
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfg); env.pop("NPUSIM_MOSAIC_DRAM", None)
    subprocess.run([os.path.join(palm.ROOT, "models", "model"), "run", acc, name, "ws"], cwd=run, env=env, capture_output=True)
    t = open(os.path.join(run, f"{acc}_{name}_layer_0.txt")).read()
    g = lambda rx: float(re.search(rx, t).group(1))
    h, tl = map(float, re.search(r"Overlap prologue/tail\s*:\s*([\d.]+)\s*/\s*([\d.]+)", t).groups())
    return j, g(r"Critical-path latency\s*:\s*([\d.]+)"), g(r"Compute-schedule latency\s*:\s*([\d.]+)"), h, tl
with ThreadPoolExecutor(6) as ex: res = list(ex.map(go, jobs))
for (knob, M, K, s), cp, cs, h, tl in res:
    roof = M * K * s * 512 / 65536
    print(f"knob{knob} M{M} K{K} strips{s}: critical {cp:8.0f} compute-sched {cs:8.0f} roofline {roof:8.0f} crit/roof {(cp/roof-1)*100:+5.1f}% prologue {h:6.0f} tail {tl:6.0f}", flush=True)
