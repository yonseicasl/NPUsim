#!/usr/bin/env python3
"""Reproduce the fixed 102-point TPU-v5e cohort; do not count gemm_only twice.

The cohort is retrospective validation, not a held-out test set. Excludes warm
buffer reuse, isolated copy-dominated small kernels and while-loop artifacts.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
from kv_streams import digest, kv_bytes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('results', type=Path)
    ap.add_argument('--tag', default='hbm2e_tpu_v5e_refsb')
    ap.add_argument('--stack-json', type=Path, help='override full decode with compare_stack JSON')
    ap.add_argument('--compute-csv', type=Path, help='override compute cohort with a new comparison')
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    root, tag = a.results, '_' + a.tag if a.tag != 'analytic' else ''
    files = set()
    def read(path):
        files.add(path)
        with path.open() as f:
            rows = list(csv.DictReader(f))
        keys = [(r['M'], r['K'], r['N']) for r in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f'duplicate shape in {path}')
        return rows
    groups = {}
    def errors(rows, kernel=False):
        result = []
        for r in rows:
            if r['src'] != 'device':
                raise ValueError('wall timings must not enter the device validation cohort')
            measured = float(r['kernel_us' if kernel else 'measured_us'])
            predicted = float(r['npusim_us'])  # primary metric retains once-per-program cost
            if not all(math.isfinite(x) and x > 0 for x in (measured, predicted)):
                raise ValueError('invalid comparison time')
            result.append(100 * (predicted / measured - 1))
        return result
    f2 = root / '2026-11-03_colab_v5e1_followup2'
    f3 = root / '2026-11-03_colab_v5e1_followup3'
    streamed = [r for r in read(f2 / f'compare_gemm{tag}.csv') if (int(r['K']), int(r['N'])) in ((4096,6144), (4096,4096), (4096,28672), (14336,4096))]
    groups['streamed_bf16'] = errors(streamed)
    groups['streamed_int8'] = errors(read(f2 / f'compare_gemm_int8_tpu_v5e_int8w{tag}.csv'))
    groups['chained_odd'] = errors(read(root / '2026-11-04_colab_v5e1_smallk2' / f'compare_gemm_oddchain{tag}_chain.csv'))
    compute_path = a.compute_csv or root / '2026-11-03_colab_v5e1_rot' / f'compare_gemm{tag}.csv'
    compute = [r for r in read(compute_path) if int(r['M']) >= 512]
    # The baseline cohort contains exactly these shapes, not merely 15 arbitrary rows.
    baseline_compute = read(root / '2026-11-03_colab_v5e1_rot' / f'compare_gemm{tag}.csv')
    shape = lambda r: (r['M'], r['K'], r['N'])
    if {shape(r) for r in compute} != {shape(r) for r in baseline_compute if int(r['M']) >= 512}:
        raise ValueError('compute cohort differs from baseline')
    groups['compute_isolated_kernel'] = errors(compute, kernel=True)
    projections = [(4096,6144), (4096,4096), (4096,28672), (14336,4096)]
    if a.stack_json:
        files.add(a.stack_json)
        data = json.loads(a.stack_json.read_text())
        stack_rows = [r for r in data['predictions'] if r['variant']=='full' and r['model']==a.tag]
        source = f3 / 'results.json'
        files.add(source)
        expected_stack = {(r['B'], r['L']): r['us_per_layer'] for r in json.loads(source.read_text())['decode_stack'] if r['variant']=='full'}
        if len(stack_rows) != len(expected_stack) or {(r['B'],r['L']) for r in stack_rows} != expected_stack.keys():
            raise ValueError('decode cohort differs from baseline')
        for r in stack_rows:
            if r['measured_us'] != expected_stack[r['B'],r['L']]:
                raise ValueError('decode measurement differs from baseline')
            if not math.isfinite(r['predicted_us']) or r['predicted_us'] <= 0:
                raise ValueError('invalid decode prediction')
        groups['decode_full'] = [100*(r['predicted_us']/r['measured_us']-1) for r in stack_rows]
    else:
        source = f3 / 'results.json'
        files.add(source)
        raw = json.loads(source.read_text())
        sim = {(int(r['M']),int(r['K']),int(r['N'])):float(r['npusim_us']) for r in read(f3 / f'compare_gemm{tag}.csv')}
        groups['decode_full'] = [100*((sum(sim[r['B'],k,n] for k,n in projections)+kv_bytes(r['B'],r['L'])/819e3)/r['us_per_layer']-1)
                                 for r in raw['decode_stack'] if r['variant']=='full']
    expected = dict(streamed_bf16=20, streamed_int8=20, chained_odd=35, compute_isolated_kernel=15, decode_full=12)
    if {k:len(v) for k,v in groups.items()} != expected:
        raise ValueError(f'cohort count mismatch: { {k: len(v) for k,v in groups.items()} }')
    groups['pooled'] = sum(groups.values(), [])
    metrics = {k:dict(n=len(v), mean_pct=statistics.mean(v), mape_pct=statistics.mean(map(abs,v)),
                      max_abs_pct=max(map(abs,v))) for k,v in groups.items()}
    result = dict(metrics=metrics, inputs={str(p.resolve()):digest(p) for p in sorted(files)},
                  caveat='Fixed retrospective cohort, not held-out. gemm_only diagnostic excluded to avoid double counting full decode.')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    for k,v in metrics.items():
        print(f"{k:26} n={v['n']:3d} mean={v['mean_pct']:+.3f}% MAPE={v['mape_pct']:.3f}% worst={v['max_abs_pct']:.3f}%")


if __name__ == '__main__':
    main()
