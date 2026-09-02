#!/usr/bin/env python3
"""PyTorch-frontend acceptance gate: the executable-IR path must reproduce the
validated nebula path, and its identities must hold end to end.

WHY THIS EXISTS. The frontend's own unit tests validate lowering and the C++
loader in isolation. What they did not pin was the WHOLE-PATH claim: that a
model entering through executable IR produces the same timing/traffic the
RTL-validated nebula path produces for the equivalent workload. The 2026-08-28
analysis established that equivalence by measurement; this gate keeps it.

Checks (asserted; non-zero exit on failure):
  FE1  COMPILE DETERMINISM: recompiling each checked-in graph fixture yields an
       executable identical (modulo metadata) to the checked-in executable
       fixture, so the artifacts cannot drift from their source silently.
  FE2  NEBULA PARITY: the linear_relu fixture IS gemm_64x64x64 (64x64 Linear =
       the RTL-validated GEMM). Run through run-ir (gemmini_sfu.cfg) and through
       the nebula path (gemmini.cfg, same mapping rows), compute-schedule and
       every DRAM traffic counter must be IDENTICAL, and the critical-path
       delta must equal the SFU busy cycles exactly -- the only physical
       difference between the two configurations.
  FE3  CONV IDENTITIES: the conv_relu fixture's computation count equals
       P*Q*K*R*S*C from its own geometry, the fused SFU processes exactly
       K*P*Q elements, and the output-commit identity line is present.
  FE4  FAIL-LOUD: an activation without an [sfu] section and a mapping with a
       wrong op_id are both refused with the specific error, not simulated.
  FE6  POOL (REAL CAPTURE): the pool_chain fixture -- captured with REAL
       torch.export (Conv2d -> MaxPool2d -> AvgPool2d), not synthesized -- must
       reproduce its geometry identities end to end: conv MACs 8*8*8*3*3*3 =
       13,824; max pool = outputs x (window-1) VMAX = 128*3 = 384 scalar ops;
       average pool = outputs x (window-1) VADD + outputs VMUL = 96+32 = 128.
       The layer reports must name each pool kind and mode.
  FE7  POOL SFU CONTRACT: a pool with no [sfu] section is refused naming the
       operation, and an [sfu] whose supported_ops allowlist lacks vmax is
       refused naming the missing primitive -- a missing execution unit is an
       architecture fact, not a zero-cost default.
  FE5  DAG/LIFETIME/OPS: a residual branch executes all six operations, pins its
       skip tensor in the GLB, frees tensors after last use, and accounts exact
       SFU scalar-operation identities for BN/add/pooling/concat.
  FE8  CLASSIFIER E2E / ALIAS ELISION: the lenet fixture -- a REAL torch.export
       capture of Conv-ReLU-Pool x2 -> flatten -> Linear-ReLU -> Linear ->
       Softmax -- runs end to end as exactly 7 operations. The flatten is elided
       as an alias of the pool's storage: no eighth layer exists, the first
       Linear reads its input GLB-resident THROUGH the alias (no fictitious
       DRAM round trip), and every MAC/scalar identity is exact
       (86400/153600/30720/1200 MACs; 2592/768 pool, 49 softmax scalar ops).
"""
import json, os, re, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODELS = os.path.join(ROOT, "models")
FIX = os.path.join(ROOT, "frontend", "fixtures")
TMP = os.environ.get("TMPDIR", "/tmp")

def sh(args, cwd=ROOT):
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True)

def run_ir(config, exec_json, mapping, stem):
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=os.path.join(ROOT, "configs"))
    r = subprocess.run(["./model", "run-ir", config, exec_json, mapping, stem],
                       cwd=MODELS, capture_output=True, text=True, env=env)
    return r

def main():
    failures = []

    # FE1
    for graph, fixture in (
        ("linear_relu_graph.json", "linear_relu.exec.json"),
        ("conv_relu_graph.json", "conv_relu.exec.json"),
        ("residual_dag_graph.json", "residual_dag.exec.json"),
        ("pool_chain_graph.json", "pool_chain.exec.json"),
        ("lenet_graph.json", "lenet.exec.json"),
    ):
        out = os.path.join(TMP, "fe1_" + fixture)
        r = sh([sys.executable, "-m", "frontend.pytorch.cli", "compile",
                "--graph", os.path.join(FIX, graph), "--output", out])
        if r.returncode != 0:
            failures.append(f"FE1 {graph}: compile failed: {r.stderr.strip()[:200]}")
            continue
        strip = lambda p: {k: v for k, v in json.load(open(p)).items() if k != "metadata"}
        if strip(out) != strip(os.path.join(FIX, fixture)):
            failures.append(f"FE1 {fixture}: recompiled executable differs from the checked-in "
                            f"fixture -- the artifacts drifted from their source graph")

    # FE2
    r = run_ir("../configs/accelerators/gemmini_sfu.cfg",
               os.path.join(FIX, "linear_relu.exec.json"),
               os.path.join(FIX, "linear_relu.map"), "fe_gate_linear")
    if r.returncode != 0:
        failures.append(f"FE2 run-ir failed: {r.stderr.strip()[:200]}")
    nebula = sh(["./npusim.sh", "run", "gemmini", "gemm_64x64x64", "ws"])
    ir_path = os.path.join(MODELS, "gemmini_sfu_fe_gate_linear_layer_0.txt")
    nb_path = os.path.join(ROOT, "result", "gemmini", "gemm_64x64x64", "ws", "layer_0.txt")
    if not (os.path.exists(ir_path) and os.path.exists(nb_path)):
        failures.append("FE2 missing one of the comparison results")
    else:
        ir, nb = open(ir_path).read(), open(nb_path).read()
        def metric(text, name):
            return float(re.search(rf"{re.escape(name)}\s*:\s*([\d.]+)", text).group(1))
        for name in ("Compute-schedule latency",):
            if metric(ir, name) != metric(nb, name):
                failures.append(f"FE2 {name}: IR {metric(ir, name)} != nebula {metric(nb, name)}")
        def dram_rows(text):
            sec = text.split("Multi-chip <-> DRAM transactions", 1)[1]
            return re.findall(r"(Input data|Weight|Output data)\s*:\s*(\d+/\d+/\d+)", sec)[:3]
        if dram_rows(ir) != dram_rows(nb):
            failures.append(f"FE2 DRAM traffic differs: IR {dram_rows(ir)} vs nebula {dram_rows(nb)}")
        sfu_busy = metric(ir, "SFU busy cycles")
        crit_delta = metric(ir, "Critical-path latency") - metric(nb, "Critical-path latency")
        if abs(crit_delta - sfu_busy) > 0.5:
            failures.append(f"FE2 critical-path delta {crit_delta} != SFU busy {sfu_busy}; "
                            f"something besides the SFU moved between the two paths")

    # FE3
    r = run_ir("../configs/accelerators/gemmini_sfu.cfg",
               os.path.join(FIX, "conv_relu.exec.json"),
               os.path.join(FIX, "conv_relu.map"), "fe_gate_conv")
    conv_path = os.path.join(MODELS, "gemmini_sfu_fe_gate_conv_layer_0.txt")
    if r.returncode != 0 or not os.path.exists(conv_path):
        failures.append(f"FE3 conv run-ir failed: {r.stderr.strip()[:200]}")
    else:
        text = open(conv_path).read()
        geo = json.load(open(os.path.join(FIX, "conv_relu.exec.json")))["operations"][0]["geometry"]
        expected_macs = (geo["output_height"]*geo["output_width"]*geo["output_channels"]
                         *geo["filter_height"]*geo["filter_width"]*geo["input_channels"]
                         //geo["groups"])
        macs = int(re.search(r"# of computations\s*:\s*(\d+)", text).group(1))
        if macs != expected_macs:
            failures.append(f"FE3 computations {macs} != geometry-derived {expected_macs}")
        scalar = int(re.search(r"Scalar operations\s*:\s*(\d+)", text).group(1))
        outs = geo["output_height"]*geo["output_width"]*geo["output_channels"]
        if scalar != outs:
            failures.append(f"FE3 SFU scalar ops {scalar} != output volume {outs}")
        if "identity OK" not in text:
            failures.append("FE3 the output-commit identity line is missing")

    # FE4
    r = run_ir("../configs/accelerators/gemmini.cfg",
               os.path.join(FIX, "linear_relu.exec.json"),
               os.path.join(FIX, "linear_relu.map"), "fe_gate_nosfu")
    if r.returncode == 0 or "requires an SFU" not in r.stderr:
        failures.append("FE4 an activation without [sfu] was not refused with the specific error")
    bad_map = os.path.join(TMP, "fe4_bad.map")
    open(bad_map, "w").write(open(os.path.join(FIX, "linear_relu.map")).read()
                             .replace("op_id  = linear", "op_id  = wrong_name"))
    r = run_ir("../configs/accelerators/gemmini_sfu.cfg",
               os.path.join(FIX, "linear_relu.exec.json"), bad_map, "fe_gate_badid")
    if r.returncode == 0 or "no mapping section has op_id" not in r.stderr:
        failures.append("FE4 a wrong op_id binding was not refused with the specific error")

    # FE5
    r = run_ir(
        "../configs/accelerators/gemmini_sfu.cfg",
        os.path.join(FIX, "residual_dag.exec.json"),
        os.path.join(FIX, "residual_dag.map"),
        "fe_gate_dag",
    )
    dag_paths = [
        os.path.join(MODELS, f"gemmini_sfu_fe_gate_dag_layer_{index}.txt")
        for index in range(6)
    ]
    dag_network_path = os.path.join(MODELS, "gemmini_sfu_fe_gate_dag.txt")
    if r.returncode != 0 or not all(os.path.exists(path) for path in dag_paths):
        failures.append(f"FE5 residual DAG run-ir failed: {r.stderr.strip()[:200]}")
    elif not os.path.exists(dag_network_path):
        failures.append("FE5 residual DAG network result is missing")
    else:
        layers = [open(path).read() for path in dag_paths]
        network = open(dag_network_path).read()
        expected_kinds = ("conv2d", "batch_norm", "elementwise", "pool2d", "pool2d", "concat")
        for index, kind in enumerate(expected_kinds):
            if index == 0:
                continue
            if f"Executable DAG operation: {kind}" not in layers[index]:
                failures.append(f"FE5 layer {index} does not report operation kind {kind}")
        expected_scalar_ops = {1: 512, 2: 256, 3: 192, 4: 256, 5: 0}
        for index, expected in expected_scalar_ops.items():
            match = re.search(r"Scalar operations\s*:\s*(\d+)", layers[index])
            if match is None or int(match.group(1)) != expected:
                actual = "missing" if match is None else match.group(1)
                failures.append(
                    f"FE5 layer {index} scalar operations {actual} != expected {expected}"
                )
            latency = re.search(r"Critical-path latency\s*:\s*([\d.]+)", layers[index])
            energy = re.search(r"Layer total energy\s*:\s*([\d.]+)", layers[index])
            if latency is None or float(latency.group(1)) <= 0:
                failures.append(f"FE5 layer {index} has no positive critical-path latency")
            if energy is None or float(energy.group(1)) <= 0:
                failures.append(f"FE5 layer {index} has no positive modeled energy")
        residency_contracts = {
            0: ("output retained in GLB", "future-use input(s) pinned"),
            1: ("1/5 input tensor(s) in GLB", "output retained in GLB"),
            2: ("2/2 input tensor(s) in GLB", "output retained in GLB"),
            5: ("2/2 input tensor(s) in GLB", "output materialized to DRAM"),
        }
        for index, fragments in residency_contracts.items():
            for fragment in fragments:
                if fragment not in layers[index]:
                    failures.append(f"FE5 layer {index} is missing residency contract: {fragment}")
        if not re.search(r"Timing scope\s*:\s*6 of 6 layers\s+\(complete", network):
            failures.append("FE5 network timing scope is not complete for all 6 DAG operations")
        if not re.search(r"Energy scope\s*:\s*6 of 6 layers\s+\(complete", network):
            failures.append("FE5 network energy scope is not complete for all 6 DAG operations")

    # FE6 -- pool chain (real torch.export capture)
    r = run_ir("../configs/accelerators/gemmini_sfu.cfg",
               os.path.join(FIX, "pool_chain.exec.json"),
               os.path.join(FIX, "pool_chain.map"), "fe_gate_pool")
    pool_paths = [os.path.join(MODELS, f"gemmini_sfu_fe_gate_pool_layer_{i}.txt")
                  for i in range(3)]
    if r.returncode != 0 or not all(os.path.exists(p) for p in pool_paths):
        failures.append(f"FE6 pool chain run-ir failed: {r.stderr.strip()[:200]}")
    else:
        texts = [open(p).read() for p in pool_paths]
        macs = int(re.search(r"# of computations\s*:\s*(\d+)", texts[0]).group(1))
        if macs != 13824:
            failures.append(f"FE6 conv computations {macs} != 13824")
        for index, (marker, expected) in enumerate(
                ((r"pool2d \(max_pool2d\)", 384), (r"pool2d \(avg_pool2d\)", 128)), start=1):
            if not re.search(marker, texts[index]):
                failures.append(f"FE6 layer {index} does not identify itself via {marker}")
            scalar = re.search(r"Scalar operations\s*:\s*(\d+)", texts[index])
            if scalar is None or int(scalar.group(1)) != expected:
                actual = "missing" if scalar is None else scalar.group(1)
                failures.append(f"FE6 layer {index} scalar operations {actual} != {expected}")

    # FE7 -- pool SFU contract
    r = run_ir("../configs/accelerators/gemmini.cfg",
               os.path.join(FIX, "pool_chain.exec.json"),
               os.path.join(FIX, "pool_chain.map"), "fe_gate_pool_nosfu")
    if r.returncode == 0 or "requires an [sfu] section" not in r.stderr:
        failures.append("FE7 a pool without [sfu] was not refused with the specific error")
    novmax = os.path.join(TMP, "fe7_sfu_novmax.cfg")
    base = open(os.path.join(ROOT, "configs/accelerators/gemmini_sfu.cfg")).read()
    open(novmax, "w").write(base.replace("[sfu]", "[sfu]\nsupported_ops = relu", 1))
    r = run_ir(novmax, os.path.join(FIX, "pool_chain.exec.json"),
               os.path.join(FIX, "pool_chain.map"), "fe_gate_pool_novmax")
    if r.returncode == 0 or "supported_ops contract" not in r.stderr or "vmax" not in r.stderr:
        failures.append("FE7 an [sfu] without vmax was not refused naming the missing primitive")

    # FE8 -- classifier end to end across an elided flatten (real torch.export capture)
    r = run_ir("../configs/accelerators/gemmini_sfu.cfg",
               os.path.join(FIX, "lenet.exec.json"),
               os.path.join(FIX, "lenet.map"), "fe_gate_lenet")
    lenet_paths = [os.path.join(MODELS, f"gemmini_sfu_fe_gate_lenet_layer_{i}.txt")
                   for i in range(7)]
    if r.returncode != 0 or not all(os.path.exists(p) for p in lenet_paths):
        failures.append(f"FE8 lenet run-ir failed: {r.stderr.strip()[:200]}")
    else:
        texts = [open(p).read() for p in lenet_paths]
        if os.path.exists(os.path.join(MODELS, "gemmini_sfu_fe_gate_lenet_layer_7.txt")):
            failures.append("FE8 an eighth layer exists -- the flatten was NOT elided")
        if any("flatten" in text for text in texts):
            failures.append("FE8 a layer report mentions flatten -- the alias leaked into timing")
        for index, expected in {0: 86400, 2: 153600, 4: 30720, 5: 1200}.items():
            macs = int(re.search(r"# of computations\s*:\s*(\d+)", texts[index]).group(1))
            if macs != expected:
                failures.append(f"FE8 layer {index} computations {macs} != {expected}")
        for index, expected in {1: 2592, 3: 768, 6: 49}.items():
            scalar = re.search(r"Scalar operations\s*:\s*(\d+)", texts[index])
            if scalar is None or int(scalar.group(1)) != expected:
                actual = "missing" if scalar is None else scalar.group(1)
                failures.append(f"FE8 layer {index} scalar operations {actual} != {expected}")
        # The money check: the first Linear consumes the pool's storage THROUGH the
        # flatten alias -- GLB-resident input, occupancy handing 256 view bytes to 120.
        if "input GLB-resident" not in texts[4]:
            failures.append("FE8 the Linear after the elided flatten read its input from DRAM "
                            "-- the alias did not keep the pool storage resident")
        if not re.search(r"GLB occupancy 256 -> 120 /", texts[4]):
            failures.append("FE8 the alias boundary does not hand exactly 256 -> 120 GLB bytes")

    if failures:
        for f in failures:
            print(" - " + f)
        print(f"{len(failures)} check(s) FAILED")
        return 1
    print("FE1 all five fixtures recompile bit-identically (modulo metadata) ok")
    print("FE2 IR path == nebula path on the RTL-validated GEMM: compute-schedule and DRAM "
          "traffic identical; critical-path delta == SFU busy exactly ok")
    print("FE3 conv computations, SFU scalar ops, and the commit identity match the geometry ok")
    print("FE4 missing-SFU and wrong-op_id runs fail loudly with the specific errors ok")
    print("FE5 residual DAG, tensor lifetime, and BN/add/pool/concat identities all match ok")
    print("FE6 real-capture pool chain: conv MACs and max/avg scalar-op identities exact ok")
    print("FE7 pool without [sfu] and [sfu] without vmax both refused with specific errors ok")
    print("FE8 lenet classifier runs e2e as 7 ops; flatten elided as an alias with exact "
          "MAC/scalar identities and a GLB-resident handoff across the boundary ok")
    print("ALL FRONTEND CHECKS PASSED")
    return 0

if __name__ == "__main__":
    sys.exit(main())
