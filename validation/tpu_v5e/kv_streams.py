#!/usr/bin/env python3
"""Model head-major contiguous KV reads with MOSAIC dram_bw_probe (no TPU required).

This is a DRAM supply model, not a QK/softmax/AV attention simulator. Up to
--max-probe-mib bytes are simulated exactly. Larger streams use the slope
between half that cap and the cap, retaining the first-access intercept once.
Use a cap >= the largest stream for exact-only simulation. No measured TPU
latencies enter the model. The JSON sidecar records every command and hash.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

MIB = 1 << 20


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def kv_bytes(batch, context):
    if batch <= 0 or context <= 0:
        raise ValueError('batch and context must be positive')
    return batch * context * 8 * 128 * 2 * 2  # heads, head_dim, bf16, K+V


def parse_probe(output, expected_bytes):
    clock = re.search(r'tCK=([\d.]+)ns', output)
    timing = re.search(r'ticks=(\d+).*read_bytes=(\d+) write_bytes=(\d+)', output)
    if not clock or not timing:
        raise ValueError('unrecognized dram_bw_probe output')
    ticks, reads, writes = map(int, timing.groups())
    if reads != expected_bytes or writes != 0:
        raise ValueError(f'probe traffic mismatch: {reads} reads, {writes} writes; expected {expected_bytes}')
    us = ticks * float(clock.group(1)) / 1000
    if not math.isfinite(us) or us <= 0:
        raise ValueError('probe returned non-positive/non-finite latency')
    return us


def stream_time(size, cap, probe):
    if size <= 0 or size % MIB:
        raise ValueError('probe accepts only positive whole MiB streams; refusing to round KV bytes')
    if cap < 2 * MIB or cap % (2 * MIB):
        raise ValueError('probe cap must be an even number of MiB, >= 2')
    if size <= cap:
        return probe(size), 'exact'
    lo, hi = probe(cap // 2), probe(cap)
    if hi <= lo:
        raise ValueError('non-positive steady-state DRAM slope')
    return hi + (size - cap) * (hi - lo) / (cap // 2), 'extrapolated'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('results', type=Path)
    ap.add_argument('--dram', type=Path, required=True)
    ap.add_argument('--probe', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--max-probe-mib', type=int, default=32)
    a = ap.parse_args()
    source = a.results / 'results.json'
    data = json.loads(source.read_text())
    cases = sorted({(r['B'], r['L']) for r in data['decode_stack'] if r['variant'] == 'full'})
    if not cases:
        ap.error('no full decode stack cases')
    cap = a.max_probe_mib * MIB
    if cap < 2 * MIB or cap % (2 * MIB):
        ap.error('--max-probe-mib must be positive and even')
    for b, length in cases:
        if kv_bytes(b, length) % MIB:
            ap.error('KV size is not a whole MiB; this probe cannot represent it exactly')
    samples = {}
    def probe(size):
        if size not in samples:
            command = [str(a.probe.resolve()), str(a.dram.resolve()), str(size // MIB), '0']
            print(f'Probing {size // MIB} MiB', flush=True)
            p = subprocess.run(command, text=True, capture_output=True, check=True, timeout=3600)
            samples[size] = {'us': parse_probe(p.stdout, size), 'command': command, 'stdout': p.stdout}
        return samples[size]['us']
    rows = []
    for b, length in cases:
        size = kv_bytes(b, length)
        us, method = stream_time(size, cap, probe)
        rows.append(dict(B=b, L=length, kv_bytes=size, kv_us=us, method=method))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    meta = dict(model='contiguous head-major K+V DRAM supply, no attention compute',
                geometry=dict(kv_heads=8, head_dim=128, element_bytes=2),
                max_probe_bytes=cap, samples=samples,
                inputs={str(p.resolve()): digest(p) for p in (source, a.dram, a.probe)},
                output_sha256=digest(a.out),
                limitation='Fresh single sequential stream; no layout, competing traffic, QK/softmax/AV or vector costs. Large streams may be extrapolated; see method.')
    a.out.with_suffix(a.out.suffix + '.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(f'Wrote {len(rows)} cases to {a.out}')


if __name__ == '__main__':
    main()
