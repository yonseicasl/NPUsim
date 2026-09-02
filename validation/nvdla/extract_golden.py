#!/usr/bin/env python3
"""Extract cycle + DBB traffic goldens from CRC-passing nv_small Verilator runs.

Cycle boundary (same definition the first data point used, README "Cycle boundary"):
  start = the LAST engine D_OP_ENABLE_0 CSB write accepted (the six enables are
          programmed SDP, CACC, CMAC_A, CMAC_B, CSC, CDMA; CDMA's is last)
  end   = the last output DBB write beat ("DBB: write, last tick")
  cycles = (end_tick - start_tick) / 2   (the harness runs 2 ticks per clock)

Traffic: every DBB beat is 8 B (64-bit, burst length 1 -- nv_small
PRIMARY_MEMIF_MAX_BURST_LENGTH_1). Beats are classified into feature/weight/dest
by the address regions the trace's own .cfg declares (mem_load with the *_feature/
*_dt .dat file = feature, *_wt.dat = weight; the mem_init region without a load =
dest). physical_bytes = beats x 8 as measured; distinct_bytes = distinct beat
addresses x 8 (re-reads collapse). The logical tensor volume is computed by the
comparison script from the trace's register-programmed dims, not here.

Only CRC-passing runs are emitted: a failing run's cycle count is not golden.
(2026-08-25 finding: the RAMDP *_E2 model's value corruption had ZERO timing
effect -- the fixed EMULATION-model binary reproduces failing runs' tick counts
exactly -- but golden status still requires the functional PASS.)
"""
import re, os, sys, glob, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
EMU = os.path.join(ROOT, "ext/nvdla/outdir/nv_small/verilator_emu")
TRACES = os.path.join(ROOT, "ext/nvdla/verif/tests/trace_tests/nv_small")

OP_ENABLE_WORDS = {0x3010//4: "CDMA", 0x4008//4: "CSC", 0x5008//4: "CMAC_A",
                   0x6008//4: "CMAC_B", 0x7008//4: "CACC", 0x9038//4: "SDP"}

def regions(trace):
    cfg = glob.glob(os.path.join(TRACES, trace, "*.cfg"))[0]
    text = open(cfg).read()
    inits = re.findall(r"mem_init\(pri_mem, (0x[0-9a-fA-F]+), (0x[0-9a-fA-F]+)", text)
    loads = re.findall(r'mem_load\(pri_mem, (0x[0-9a-fA-F]+), "([^"]+)"', text)
    out = {}
    loaded = set()
    for addr, fname in loads:
        kind = "weight" if re.search(r"_wt\.dat$|_in_wt|_in_weight", fname) else "feature"
        size = next(int(s,16) for a,s in inits if a == addr)
        out[kind] = (int(addr,16), size); loaded.add(addr)
    for a, s in inits:
        if a not in loaded:
            out["dest"] = (int(a,16), int(s,16))
    return out

def extract(trace):
    log = os.path.join(EMU, f"run_{trace}", "run.log")
    if not os.path.exists(log): return None
    text = open(log, errors="replace").read()
    if "*** PASS" not in text: return ("NOT_PASS", None)
    reg = regions(trace)
    start = None
    for m in re.finditer(r"\((\d+)\) write to nvdla: addr ([0-9a-f]+), data 0*1\b", text):
        if int(m.group(2), 16) in OP_ENABLE_WORDS:
            start = int(m.group(1))
    reads, writes = {}, {}
    for m in re.finditer(r"\(\d+\) DBB: read request from dla, addr ([0-9a-f]+) burst (\d+)", text):
        a = int(m.group(1), 16) & ~7
        for b in range(int(m.group(2)) + 1):
            reads[a + 8*b] = reads.get(a + 8*b, 0) + 1
    last_write = None
    for m in re.finditer(r"\((\d+)\) DBB: write request from dla, addr ([0-9a-f]+)", text):
        a = int(m.group(2), 16) & ~7
        writes[a] = writes.get(a, 0) + 1
    for m in re.finditer(r"\((\d+)\) DBB: write, last tick", text):
        last_write = int(m.group(1))
    crc = re.search(r"CRC matched [0-9a-f]+ \+ [0-9a-f]+ -> ([0-9a-f]+)", text)
    def classify(table):
        out = {k: [0, 0] for k in ("feature", "weight", "dest", "other")}
        for a, n in table.items():
            for k, (base, size) in reg.items():
                if base <= a < base + size:
                    out[k][0] += 8*n; out[k][1] += 8
                    break
            else:
                out["other"][0] += 8*n; out["other"][1] += 8
        return out
    r, w = classify(reads), classify(writes)
    cycles = (last_write - start) // 2 if start and last_write else None
    return ("OK", {
        "cycles": cycles, "start_tick": start, "end_tick": last_write,
        "crc": crc.group(1) if crc else "",
        "in_phys": r["feature"][0], "in_distinct": r["feature"][1],
        "wt_phys": r["weight"][0], "wt_distinct": r["weight"][1],
        "out_phys": w["dest"][0], "out_distinct": w["dest"][1],
        "read_other": r["other"][0], "write_other": w["other"][0] + w["feature"][0] + w["weight"][0],
    })

def main():
    rows = []
    for d in sorted(os.listdir(EMU)):
        if not d.startswith("run_"): continue
        t = d[4:]
        res = extract(t)
        if res is None: continue
        st, g = res
        if st != "OK":
            print(f"{t:<40} {st}"); continue
        print(f"{t:<40} PASS cycles={g['cycles']} in={g['in_phys']}B wt={g['wt_phys']}B out={g['out_phys']}B"
              + (f" [unclassified r={g['read_other']} w={g['write_other']}]" if g['read_other'] or g['write_other'] else ""))
        rows.append({"trace": t, **g})
    if "--write" in sys.argv and rows:
        path = os.path.join(ROOT, "validation/nvdla/golden_rtl_full.csv")
        with open(path, "w", newline="") as f:
            wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wcsv.writeheader(); wcsv.writerows(rows)
        print(f"\nwrote {path} ({len(rows)} traces)")

if __name__ == "__main__":
    main()
