#!/usr/bin/env python3
"""Traffic-counter validation: analytic identities vs simulator counters.

Checks, per validated workload:
  T1  output DRAM bytes == output volume (single-write convention, DR6)
  T2  weight DRAM bytes >= weight volume; ratio reported (1.0 = fetch-once)
  T3  input  DRAM bytes == exact mapped input union (including convolution halo)
  T4  GLB->array bytes >= DRAM bytes per type (on-chip streams at least off-chip)
  T5  GLB->array weight bytes == weight events x array weight tile bytes (PA1)
  T6  computation cycles == MACs / active PEs
  T7  per-stage busy == the printed axes combined by that stage's documented rule
      (the PE stage also folds in the compute schedule and the format-IP axis, and
      combines them by the printed double-/single-buffered rule -- P1-A/P4-1).
      LAYER SCOPE ONLY: network.txt sums busy and axes independently over layers, and
      sum_layer(max(axes)) != max(sum_layer(axes)), so network axes do not reduce to
      network busy. That rollup contract is pinned by check_network_rollup() in
      validation/check_timing.py instead (L3).
  T8  critical path >= the validated compute-schedule metric
  T9  a GLB-bypassed datatype costs ZERO GLB SRAM access but still carries non-zero
      GLB<->PE-array fabric traffic (P1-B: bypass skips the SRAM, not the wires)
  T10 the GLB access axis never exceeds its own datatype combination of the reported
      per-datatype GLB access cycles, using the printed rule (sum for a shared port, max
      for separate partitions) -- L1. Equality on a single-chip config; the bound also
      holds with several chips, because the axis takes the max ACROSS chips of each chip's
      own datatype total while the per-datatype rows are already maxed across chips.
  T11 convolution results must report an applied halo contract whose unique element and
      post-coalescing transaction counts equal T3; non-overlapping GEMMs must report that
      halo reuse is not needed.
  T12 the psum round trip is exact: stores - loads == the number of output tiles. Every partial
      sum the array stores must be loaded back to be accumulated, and each tile takes one extra
      store when it is finished. A store leg without its load leg is not a cheaper schedule, it
      is an unaccumulated partial sum -- and it understates the GLB traffic bound while leaving
      the compute-schedule latency untouched, so no other check here notices (E20-3b; Eyeriss
      conv3 once read 19656:312, and the 512x512x512 GEMM 96:0).

Units are derived from each config's link bitwidths and data formats:
  eyeriss: fp16 (2 B/elem), GLB link 256 b (32 B/txn), DRAM link 64 b (8 B/txn)
  gemmini: int8 (1 B/elem), GLB link 128 b (16 B/txn), DRAM link 64 b (8 B/txn)
"""
import os, re, sys

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."

# name -> (result file, elem_bytes, glb_txn_B, dram_txn_B,
#          in_vol, wt_vol, out_vol (elements, incl. padding/halo), macs, active)
def alexnet_cases():
    b = 4
    cases = {}
    # (name, layer, K, P, Q, C, R, S, stride, groups, pad_P)
    shapes = [("conv1", 0, 96, 55, 55, 3, 11, 11, 4, 1, 56),
              ("conv2", 2, 256, 27, 27, 96, 5, 5, 1, 2, 28),
              ("conv3", 4, 384, 13, 13, 256, 3, 3, 1, 1, 13),
              ("conv4", 5, 384, 13, 13, 384, 3, 3, 1, 2, 13),
              ("conv5", 6, 256, 13, 13, 384, 3, 3, 1, 2, 13)]
    active = {"conv1": 154, "conv2": 140, "conv3": 156, "conv4": 156, "conv5": 156}
    for name, L, K, P, Q, C, R, S, st, g, Pp in shapes:
        ih = (Pp - 1)*st + R
        iw = (Q - 1)*st + S
        in_vol = b*C*ih*iw
        wt_vol = K*C*R*S//g
        out_vol = b*K*Pp*Q
        macs = b*K*Pp*Q*(C//g)*R*S
        cases[name] = (f"result/eyeriss/alexnet/silicon/layer_{L}.txt", 2, 32, 8,
                       in_vol, wt_vol, out_vol, macs, active[name])
    return cases

def gemm_cases():
    cases = {}
    for m, k, n in [(64,64,64),(256,256,256),(512,512,512)]:
        # NPUsim models neither bias MACs nor a bias stream, so the reference
        # volumes exclude the bias[n] appended in the .wgh file.
        cases[f"gemm{m}x{k}x{n}"] = (
            f"result/gemmini/gemm_{m}x{k}x{n}/ws/layer_0.txt", 1, 16, 8,
            m*k, k*n, m*n, m*k*n, 256)
    return cases

def timeline_invariants(txt, name, failures):
    """T7/T8 (CE3/CE7/P1-A/P4-1): each stage's busy value must equal its printed axes
    combined by that stage's documented rule, and the critical path must be at least the
    validated compute-schedule metric. Applies to LAYER results only -- see the module
    docstring for why the network rollup has a different contract (L3).

    Memory stages: busy == max(access, link, overlap).
    PE stage: the compute schedule and the format-IP axis are axes too, and the
    combination rule depends on local-buffer double buffering (printed):
      double-buffered -> max(compute, access, link, overlap, format)
      single-buffered -> compute + max(access, link, overlap) + format
    The single-buffered rule is NOT a sum of all axes: the overlap axis is already the
    pipelined makespan of the access/link work, so summing the three would charge the
    same LB<->MAC transactions two or three times (P1-A)."""
    tl = txt.split("Layer timeline")[1].split("MAC result")[0]
    if "[layer: busy = max of these axes]" not in tl:
        failures.append(f"{name}: T7 expects a layer-scope result; the report does not "
                        f"declare the layer axis contract")
        return False
    sched = float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", tl).group(1))
    crit = float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", tl).group(1))
    fmt = float(re.search(r"PE format-IP axis\s*:\s*([\d.]+)", tl).group(1))
    double_buffered = "double-buffered" in re.search(r"PE local buffer\s*:\s*(\S+)", tl).group(1)
    busy = [float(v) for v in re.findall(r"\* [^:\n]+:\s*([\d.]+) cycles \(", tl)]
    axes = [[float(a), float(b), float(c)] for a, b, c in
            re.findall(r":\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+) cycles", tl)]
    ok = True
    for stage in range(5):
        if stage < 4:
            expect = max(axes[stage])
        elif double_buffered:
            expect = max(axes[stage] + [sched, fmt])
        else:
            expect = sched + max(axes[stage]) + fmt
        if abs(busy[stage] - expect) > 0.5:
            failures.append(f"{name}: T7 stage{stage} busy {busy[stage]} != expected {expect}")
            ok = False
    if crit + 0.5 < sched:
        failures.append(f"{name}: T8 critical {crit} < compute-schedule {sched}")
        ok = False
    return ok

def parse(path):
    txt = open(os.path.join(ROOT, path)).read()
    comp = float(re.search(r"Computation cycle\s*:\s*([\d.]+)", txt).group(1))
    glb_sec = txt.split("Global buffer result")[1]
    tx = re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)",
                    glb_sec.split("serialized")[1].split("# of request")[0])
    glb = {k: int(v) for k, v in tx}
    # P1-B: the simulator names the GLB-bypassed (direct-stream) datatypes explicitly.
    bypassed = set(re.search(r"GLB-bypassed \(direct stream\)\s*:(.*)", glb_sec).group(1).split())
    glb_shared = "sum (shared port)" in re.search(r"GLB datatype rule\s*:(.*)", txt).group(1)
    glb_access_axis = [float(a) for a, _, _ in
                       re.findall(r":\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+) cycles",
                                  txt.split("Busy-cycle axes")[1])][2]
    glb_access = {k: float(v) for k, v in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.]+) cycles",
        glb_sec.split("Access cycle")[1].split("Energy")[0])}
    ev = re.findall(r"(Input data|Weight|Output data)\s*:\s*(\d+)\n",
                    glb_sec.split("# of data transfer")[1].split("PE-array")[0])
    events = {k: int(v) for k, v in ev}
    dram_sec = txt.split("Multi-chip <-> DRAM transactions")[1]
    dx = re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", dram_sec)[:3]
    dram = {k: int(v) for k, v in dx}
    tiles = txt.split("======= PE array =======")[1].split("========================")[0]
    array_tile = {k: int(v) for k, v in
                  re.findall(r"(Input data|Weight|Output data)\s*:\s*(\d+)", tiles)}
    halo = re.search(r"Input halo reuse\s*:\s*(.*)", txt).group(1).strip()
    psum = re.search(r"Psum round trip\s*:\s*(\d+) loads / (\d+) stores", txt)
    psum_legs = (int(psum.group(1)), int(psum.group(2))) if psum else None
    return (comp, glb, events, dram, array_tile, bypassed, glb_access, glb_shared,
            glb_access_axis, halo, psum_legs)

def main():
    cases = {}
    cases.update(alexnet_cases())
    cases.update(gemm_cases())
    failures = 0
    invariant_failures = []
    print(f"{'case':>16s} {'T1 out':>7s} {'T2 wt':>6s} {'T3 in':>6s} {'T4':>4s} {'T5 glb-wt':>9s} {'T6 comp':>8s} {'T7/8':>5s} {'T9':>4s} {'T10':>4s} {'T11':>4s} {'T12':>4s}")
    for name, (path, eb, gtxn, dtxn, in_v, wt_v, out_v, macs, active) in cases.items():
        (comp, glb, events, dram, tile, bypassed, glb_access, glb_shared,
         glb_axis, halo, psum_legs) = parse(path)
        rows = {}
        # T1: output DRAM == volume
        out_b, out_ref = dram["Output data"]*dtxn, out_v*eb
        rows["T1"] = out_b/out_ref
        # T2: >= volume; T3: the input is the exact dense union, including halo.
        rows["T2"] = dram["Weight"]*dtxn/(wt_v*eb)
        rows["T3"] = dram["Input data"]*dtxn/(in_v*eb)
        # T4: GLB >= DRAM per type (bytes). P1-B: this now applies to a GLB-bypassed
        # datatype too -- bypass skips the SRAM, not the fabric, so the stream is still
        # counted and the on-chip volume must still cover the off-chip volume.
        # 2026-08-25: OUTPUT is now checked too. It was silently absent from this tuple while
        # the docstring said "per type", so the one datatype whose on-chip volume can legitimately
        # fall BELOW the off-chip volume if a write-back is missed was the one not covered. Every
        # output byte reaches DRAM through this boundary, so the inequality must hold for it.
        t4 = all(glb[t]*gtxn + 1e-9 >= dram[t]*dtxn
                 for t in ("Input data", "Weight", "Output data"))
        # T5: GLB weight bytes == events x array tile bytes (holds for bypassed weight too)
        t5 = glb["Weight"]*gtxn/(events["Weight"]*tile["Weight"]*eb) if events["Weight"] else float("nan")
        t5_ok = abs(t5 - 1) < 0.01
        # T6
        t6 = comp/(macs/active)
        t78 = timeline_invariants(open(os.path.join(ROOT, path)).read(), name, invariant_failures)
        # T9 (P1-B): the bypass contract, checked instead of exempted. A bypassed
        # datatype must show zero GLB SRAM access cycles (no fill write, no read) AND
        # non-zero GLB<->PE-array fabric transactions; a non-bypassed one must show
        # non-zero access cycles whenever it carries traffic.
        label_of = {"Input data": "input", "Weight": "weight", "Output data": "output"}
        t9_ok = True
        for t, label in label_of.items():
            if label in bypassed:
                if glb_access[t] != 0.0 or glb[t] == 0:
                    t9_ok = False
                    invariant_failures.append(
                        f"{name}: T9 {label} is GLB-bypassed but access={glb_access[t]} "
                        f"transactions={glb[t]} (expected access 0 with non-zero traffic)")
            elif glb[t] > 0 and glb_access[t] <= 0.0:
                t9_ok = False
                invariant_failures.append(
                    f"{name}: T9 {label} is GLB-resident with {glb[t]} transactions but "
                    f"zero access cycles")
        # T10 (L1): the GLB access axis must not exceed its own datatype combination of the
        # reported per-datatype access cycles. Summing the type-combined read and fill sides
        # separately (the pre-fix form) breaks this on a separate buffer: it adds the peaks
        # of two different partitions, which actually run in parallel.
        glb_combined = sum(glb_access.values()) if glb_shared else max(glb_access.values())
        t10_ok = glb_axis <= glb_combined + 0.5
        if not t10_ok:
            invariant_failures.append(
                f"{name}: T10 GLB access axis {glb_axis} exceeds its datatype combination "
                f"{glb_combined} of the per-type access cycles "
                f"({'sum' if glb_shared else 'max'} rule)")
        # T11 (E20-3): pin the report-level halo provenance to the exact T3 counter. This catches
        # a future change that preserves total input bytes while silently disabling the capacity
        # gate or falling back to a fractional repetition approximation.
        if name.startswith("conv"):
            match = re.fullmatch(
                r"applied; (\d+) -> (\d+) input elements, working set (\d+) B fits GLB, "
                r"DRAM serialized (\d+) -> (\d+)", halo)
            t11_ok = bool(match)
            if match:
                replicated, unique, working, pre, post = map(int, match.groups())
                t11_ok = (replicated > unique == in_v and working > 0 and pre > post and
                          post == dram["Input data"])
            if not t11_ok:
                invariant_failures.append(
                    f"{name}: T11 invalid applied halo contract '{halo}' for input volume "
                    f"{in_v} and DRAM transactions {dram['Input data']}")
        else:
            t11_ok = halo == "not needed; per-repetition streaming already moves exactly the union"
            if not t11_ok:
                invariant_failures.append(f"{name}: T11 unexpected halo contract '{halo}'")
        # T12 (E20-3b): pair the psum legs. The upper slack is one final store per output
        # tile -- that psum is never reloaded, it leaves over the off-chip store path.
        if psum_legs is None:
            t12_ok = False
            invariant_failures.append(f"{name}: T12 result reports no psum round trip line")
        else:
            loads, stores = psum_legs
            output_tiles = out_v/tile["Output data"]
            if stores == 0:
                t12_ok = loads == 0
                if not t12_ok:
                    invariant_failures.append(
                        f"{name}: T12 {loads} psum loads against 0 stores -- nothing was written "
                        f"out, so there is nothing to read back")
            else:
                # Exact, not a band. Each output tile is stored once every time it is swapped out
                # and read back once every time it returns, plus ONE final store when it is
                # finished (that last one leaves over the off-chip path instead of being reloaded).
                # So stores - loads is precisely the number of output tiles. Every validated case
                # meets this exactly; an inequality here would have let the 512x512x512 GEMM's
                # reduction-depth over-count hide inside the slack.
                t12_ok = abs(stores - loads - output_tiles) < 0.5
                if not t12_ok:
                    invariant_failures.append(
                        f"{name}: T12 psum legs unpaired: {loads} loads, {stores} stores, "
                        f"expected stores - loads == {output_tiles:.0f} output tiles")
        ok = (abs(rows["T1"]-1) < 0.01 and rows["T2"] > 0.999 and
              abs(rows["T3"]-1) < 0.01 and t4 and t5_ok and abs(t6-1) < 0.01 and
              t78 and t9_ok and t10_ok and t11_ok and t12_ok)
        if not ok:
            failures += 1
        print(f"{name:>16s} {rows['T1']:7.3f} {rows['T2']:6.2f} {rows['T3']:6.2f} "
              f"{'ok' if t4 else 'FAIL':>4s} {t5:9.3f} {t6:8.3f} {'ok' if t78 else 'FAIL':>5s} "
              f"{'ok' if t9_ok else 'FAIL':>4s} {'ok' if t10_ok else 'FAIL':>4s} "
              f"{'ok' if t11_ok else 'FAIL':>4s} {'ok' if t12_ok else 'FAIL':>4s}"
              f"{'' if ok else ' <-- FAIL'}")
    for f in invariant_failures:
        print("  " + f)
    print("\nT1/T3/T5/T6 must be 1.000; T2 >= 1 (ratio = refetch factor); T4 GLB>=DRAM; "
          "T7 busy==combined axes; T8 critical>=compute-schedule; T9 bypass = no SRAM access, "
          "fabric traffic still charged; T10 GLB axis <= its datatype combination; "
          "T11 exact halo contract reported; T12 psum loads and stores paired")
    print(f"{'ALL CHECKS PASSED' if failures == 0 else f'{failures} case(s) FAILED'}")
    return 1 if failures else 0

if __name__ == "__main__":
    sys.exit(main())
