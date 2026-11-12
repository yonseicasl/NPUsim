#!/usr/bin/env python3
"""Audit measured int8/bf16 pairs and optimized HLO; never fit simulator constants.

--source points to the validation/tpu_v5e directory containing the original
followup2, followup3 and Nov4 profile datasets. Output is separate from inputs.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import statistics
from pathlib import Path

RUNS = ('2026-11-03_colab_v5e1_followup2',
        '2026-11-03_colab_v5e1_followup3',
        '2026-11-04_colab_v5e1_profile')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--prediction-bf16', type=Path)
    ap.add_argument('--prediction-int8', type=Path)
    args = ap.parse_args()
    if bool(args.prediction_bf16) != bool(args.prediction_int8):
        ap.error('provide both prediction CSVs or neither')
    args.output.mkdir(parents=True, exist_ok=True)
    hashes, pairs, hlo, profiles = {}, [], [], []

    def read(path):
        data = path.read_bytes()
        hashes[str(path.relative_to(args.source))] = hashlib.sha256(data).hexdigest()
        return data.decode()

    for run in RUNS:
        folder = args.source / 'results' / run
        data = json.loads(read(folder / 'results.json'))
        indexed = {}
        for r in data['stream_gemm']:
            key = (r['M'], r['K'], r['N'], r['weights'], r['cold'])
            if key in indexed:
                raise ValueError(f'duplicate {run}: {key}')
            indexed[key] = r
        if len(indexed) != 80:
            raise ValueError(f'{run}: expected 80 cold/warm dtype rows')
        for key, bf in indexed.items():
            m, k, n, dtype, cold = key
            if dtype != 'bf16':
                continue
            it = indexed[(m, k, n, 'int8', cold)]
            if min(it['us'], bf['us']) <= 0:
                raise ValueError('nonpositive time')
            # An algebraic diagnostic, not conversion latency. Halving bf16
            # also halves shared work and is invalid in the compute regime.
            pairs.append(dict(run=run, M=m, K=k, N=n, cold=cold,
                bf16_us=bf['us'], int8_us=it['us'],
                int8_minus_half_bf16_us=it['us']-.5*bf['us'],
                bf16_weight_GBps=2*k*n/bf['us']/1000,
                int8_weight_GBps=k*n/it['us']/1000,
                logical_activation_plus_output_bf16_bytes=2*m*(k+n),
                bf16_copies=bf['copies'], int8_copies=it['copies']))
        for path in sorted((folder/'hlo').glob('*gemm*.txt')):
            body = read(path)
            for line in body.splitlines():
                if '"window_config"' not in line:
                    continue
                config = json.loads(line.split('backend_config=', 1)[1])
                win = config['window_config']
                hlo.append(dict(run=run, case=path.stem,
                    emitter=config.get('convolution_algorithm_config', {}).get('emitter'),
                    window=win,
                    # Scope is evidence only: s8 operand survives in the module.
                    contains_s8=bool(re.search(r'\bs8\[', body)),
                    convolution_lines=[x.strip() for x in body.splitlines()
                                       if ' convolution(' in x]))
        stats = folder/'xprof/hlo_stats.csv'
        if stats.exists():
            for row in csv.DictReader(io.StringIO(read(stats))):
                if not re.match(r'(gemm_|stream_gemm_)', row['case']):
                    continue
                if 'fusion' in row['category'] and 'dot_general' in row['tf_op_name']:
                    profiles.append(dict(run=run, **{k: row[k] for k in ('case','category',
                        'occurrences','avg_time','avg_self_time','dma_stall_percent',
                        'hbm_bw','bound_by')}))
    with (args.output/'paired.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]));w.writeheader();w.writerows(pairs)
    summaries=[]
    for run in RUNS:
        for m in (1,8,32,128,512):
            rows=[r for r in pairs if r['run']==run and r['M']==m and r['cold']]
            values=[r['int8_minus_half_bf16_us'] for r in rows]
            summaries.append(dict(run=run,M=m,count=len(rows),
                residual_mean_us=statistics.mean(values),
                residual_min_us=min(values),residual_max_us=max(values)))
    comparisons, comparison_summary = [], []
    if args.prediction_bf16:
        for dtype, path in [('bf16', args.prediction_bf16), ('int8', args.prediction_int8)]:
            raw = path.read_bytes()
            hashes['prediction:' + path.name] = hashlib.sha256(raw).hexdigest()
            predictions = {}
            for row in csv.DictReader(io.StringIO(raw.decode())):
                key = tuple(int(row[k]) for k in ('M', 'K', 'N'))
                if key in predictions:
                    raise ValueError(f'duplicate prediction: {key}')
                predictions[key] = row
            for pair in pairs:
                if not pair['cold']:
                    continue
                key = tuple(pair[k] for k in ('M','K','N'))
                row = predictions[key]
                predicted = float(row['npusim_us'])
                # Current saved predictions already omit program setup. Never
                # silently subtract a constant from arbitrary cached results.
                if abs(predicted-float(row['npusim_noprog_us'])) > .011:
                    raise ValueError('cached predictions must be in-program')
                observed = [('loop_wall_slope', pair[dtype+'_us'])]
                prefix = 'gemm_' + 'x'.join(map(str, key)) + '_' + dtype + '_cold'
                special = 'stream_gemm_8x4096x4096_bf16'
                matches = [r for r in profiles if r['run']==pair['run'] and
                           (r['case']==prefix or (key==(8,4096,4096) and
                            dtype=='bf16' and r['case']==special))]
                if len(matches)>1:
                    raise ValueError('ambiguous fusion profile')
                if matches:
                    observed.append(('profile_convolution_fusion', float(matches[0]['avg_self_time'])))
                for scope, measured in observed:
                    comparisons.append(dict(run=pair['run'], dtype=dtype,
                        M=key[0], K=key[1], N=key[2], measurement_scope=scope,
                        measured_us=measured, npusim_inprogram_us=predicted,
                        error_pct=100*(predicted/measured-1)))
        with (args.output/'scope_comparison.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
        for key in sorted({(r['run'],r['dtype'],r['measurement_scope']) for r in comparisons}):
            errs=[r['error_pct'] for r in comparisons if (r['run'],r['dtype'],r['measurement_scope'])==key]
            comparison_summary.append(dict(run=key[0], dtype=key[1],measurement_scope=key[2],
                count=len(errs),mean_error_pct=statistics.mean(errs),mape_pct=statistics.mean(map(abs,errs))))
    result=dict(timing_boundary='median wall-time slope of dependent GEMM loop; CSV device_us is not direct trace duration',
        interpretation='int8-minus-half-bf16 is not isolated conversion latency; contains shared work and changes in compiler scheduling. Warm means one weight copy, not proven on-chip residency.',
        input_sha256=hashes,paired_rows=len(pairs),summaries=summaries,
        comparison_summary=comparison_summary,
        optimized_hlo=hlo,profile_fusions=profiles)
    (args.output/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(paired_rows=len(pairs),hlo_windows=len(hlo),
                         profile_fusions=len(profiles),hashed_files=len(hashes))))


if __name__ == '__main__':
    main()
