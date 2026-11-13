#!/usr/bin/env python3
"""Per-op byte accounting of profiled XLA programs, read from the XProf xplane.pb of each kept trace.

What the numbers are: STATIC COMPILER ACCOUNTING. `bytes_accessed` and `memory_access_breakdown` are stats of the
op's event METADATA (one record per HLO op, shared by all its executions), written by the compiler's cost
analysis -- not per-execution hardware counters, and the traces hold no per-DMA events (the device lines
"Async XLA Ops", "TC Overlay", "Pallas Primitives", "_counters_" are empty). memory_space is printed as the raw
id: 1 for operands with the default layout (no S(n) tag in the HLO), 3 for the destination of a copy tagged S(1);
XProf's own naming of these (HBM, VMEM) is not verifiable from the artifacts.

Usage: python3 trace_bytes.py TRACES_DIR OUT.csv      (TRACES_DIR/<program>/**/*.xplane.pb; no schema library)
"""
import csv
import glob
import os
import re
import struct
import sys

def fields(buf):
    i, n = 0, len(buf)
    while i < n:
        key = 0; s = 0
        while True:
            b = buf[i]; i += 1; key |= (b & 0x7F) << s; s += 7
            if b < 0x80: break
        f, wt = key >> 3, key & 7
        if wt == 0:
            v = 0; s = 0
            while True:
                b = buf[i]; i += 1; v |= (b & 0x7F) << s; s += 7
                if b < 0x80: break
        elif wt == 2:
            ln = 0; s = 0
            while True:
                b = buf[i]; i += 1; ln |= (b & 0x7F) << s; s += 7
                if b < 0x80: break
            v = buf[i:i + ln]; i += ln
        elif wt == 1: v = buf[i:i + 8]; i += 8
        elif wt == 5: v = buf[i:i + 4]; i += 4
        else: raise ValueError(wt)
        yield f, wt, v
def stat(buf, names):
    d = dict((f, v) for f, _, v in fields(buf))
    name = names.get(d.get(1), d.get(1))
    val = next((d[k] for k in (2, 3, 4, 5, 6, 7) if k in d), None)
    if 2 in d: val = struct.unpack("<d", d[2])[0]
    if 7 in d: val = ("ref", names.get(d[7], d[7]))
    return name, val
def breakdown(b):
    out = []
    for f, _, v in fields(b):
        if f == 1:
            m = dict((ff, vv) for ff, _, vv in fields(v))
            out.append((m.get(1, 0), m.get(2, 0), m.get(3, 0)))        # (operation_type, memory_space, bytes)
    return out
def ops(path):
    """{op name: {stat: value}} for the XLA op event metadata of the TPU plane, plus line names."""
    res, lines = {}, []
    for f, _, plane in fields(open(path, "rb").read()):
        if f != 1: continue
        p = list(fields(plane))
        pname = next((v.decode() for ff, _, v in p if ff == 2), "")
        if "TPU" not in pname: continue
        names = {}
        for ff, _, v in p:
            if ff == 5:                                   # map<int64, XStatMetadata>
                e = dict((a, b) for a, _, b in fields(v)); m = dict((a, b) for a, _, b in fields(e[2]))
                names[e[1]] = m.get(2, b"").decode()
        for ff, _, v in p:
            if ff == 3:
                lines.append(next((b.decode() for a, _, b in fields(v) if a == 2), "?"))
            if ff == 4:                                   # map<int64, XEventMetadata>
                e = dict((a, b) for a, _, b in fields(v))
                md = list(fields(e[2])); nm = next((b.decode() for a, _, b in md if a == 2), "")
                st = dict(stat(b, names) for a, _, b in md if a == 5)
                if st: res[nm] = st
    return res, lines

def main(traces, out):
    rows = []
    for d in sorted(os.listdir(traces)):
        pb = glob.glob(os.path.join(traces, d, "**", "*.xplane.pb"), recursive=True)
        if not pb:
            continue
        for long_name, st in ops(pb[0])[0].items():
            bd = st.get("memory_access_breakdown")
            bd = breakdown(bd) if isinstance(bd, (bytes, bytearray)) else []
            acc = lambda op, space: sum(b for o, s, b in bd if o == op and s == space)
            rows.append(dict(program=d, op=re.match(r"%([\w.\-]+)", long_name).group(1),
                             hlo_category=(st.get("hlo_category") or b"").decode() if isinstance(st.get("hlo_category"), bytes) else st.get("hlo_category"),
                             bytes_accessed=st.get("bytes_accessed"), model_flops=st.get("model_flops"),
                             read_space1=acc(1, 1), write_space1=acc(2, 1), read_space3=acc(1, 3), write_space3=acc(2, 3),
                             other=";".join(f"{o}:{s}:{b}" for o, s, b in bd if (o, s) not in ((1, 1), (2, 1), (1, 3), (2, 3))),
                             operands=" | ".join(re.findall(r"(\w+\[[\d,]*\]\{[^}]*\})", long_name.split("=", 1)[1])[:6])))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    # the balance loop: one fused op per iteration -- what does it read and write?
    bad = n = 0
    for r in rows:
        m = re.fullmatch(r"balance_(bf16|int8)_(\d+)x(\d+)_loop", r["program"])
        if m and "fusion" in r["op"]:
            M, K = int(m.group(2)), int(m.group(3))
            want = K * K * (1 if m.group(1) == "int8" else 2) + M * K * 2
            n += 1
            bad += not (want <= r["read_space1"] <= want + 16 and r["write_space1"] <= 8 and not r["read_space3"])
    print(f"wrote {len(rows)} op records of {len({r['program'] for r in rows})} programs -> {out}")
    if n:
        print(f"balance loop fused op: {n - bad} of {n} read weights + M*K*2 input bytes (+ <= 16 B of scalars) from "
              f"space 1, write <= 8 B, touch no space 3")


if __name__ == "__main__":
    main(*sys.argv[1:3])
