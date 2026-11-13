import csv, hashlib, json, math, re, statistics as st
from pathlib import Path
P=Path('/home/bogil/MOSAIC/.claude/worktrees/npusim-engine/simulator/3rd/npusim/.claude/worktrees/tpu-measure')
V=P/'validation/tpu_v5e'
O=V/'results/2026-11-12_window_candidate'
S=V/'results/2026-11-12_colab_v5e1_balance'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
measure={(r['dtype'],int(r['M']),int(r['K'])):r for r in csv.DictReader((S/'balance.csv').open()) if r['schedule']=='loop'}
for line in (O/'FROZEN.sha256').read_text().splitlines():
    digest,rel=line.split(maxsplit=1)
    assert sha((V/rel).resolve())==digest, rel
out={}
for tag in ('analytic','hbm2e_tpu_v5e_refsb'):
    f=O/f'balance_npusim_{tag}.csv'
    if not f.exists(): print(tag,'PENDING'); continue
    rows=list(csv.DictReader(f.open()))
    assert len(rows)==28
    meta=json.loads(f.with_suffix('.json').read_text())
    assert meta['balance_csv_sha256']==sha(S/'balance.csv')
    assert meta['frozen']['binary_sha256']==sha(P/'models/model')
    assert meta['frozen']['library_sha256']==sha(P/'library/libnpusim.so')
    for n,h in meta['frozen']['configs'].items(): assert h==sha(P/'configs/accelerators'/f'{n}.cfg')
    for n,h in (meta['dram'] or {}).items(): assert h==sha(Path(n))
    assert meta['window_switch_cycle']['value']==590
    gen=O/f'generated_{tag}'
    assert len(list(gen.iterdir()))==84
    assert len(list(gen.rglob('*.txt')))==168
    seen=set()
    for r in rows:
        dt,M,K=r['dtype'],int(r['M']),int(r['K'])
        key=(dt,M,K); assert key not in seen; seen.add(key)
        ref=measure[key]
        assert float(r['measured_us'])==float(ref['fused_ops_sum_per_iter_us'])
        assert r['hlo_sha256']==sha(S/'hlo'/f'balance_{dt}_{M}x{K}_loop.txt')
        assert r['window_model']=='count-based equal-share'
        for v in 'ABC':
            name=f'bal_{dt}_{M}x{K}_{v}'.lower(); d=gen/name
            net=(d/'network.cfg').read_text()
            val=lambda k:int(re.search(rf'(?m)^{k}=(\d+)',net).group(1))
            assert (val('batch'),val('height')*val('width')*val('channels'),val('output'))==(M,K,K)
            assert f'weight=weights/{name}.wgh' in net
            mp=[[int(t) for t in line.split('=')[1].split(',') if t.strip()] for line in next(d.rglob('ws.map')).read_text().splitlines() if '=' in line]
            assert tuple(math.prod(x[c] for x in mp) for c in (1,4,0))==(M,K,K)
            config=(d/f'{name}.cfg').read_text()
            kv=lambda k:re.findall(rf'(?m)^{k}\s*=\s*(\d+)',config)
            assert kv('program_setup_cycle')==['0'] and kv('activation_residency')==['2']
            assert kv('offchip_window_count')==([] if v=='A' else [r['windows']])
            assert kv('window_switch_cycle')==(['590'] if v=='C' else [])
            layer=(d/f'{name}_{name}_layer_0.txt').read_text()
            network=(d/f'{name}_{name}.txt').read_text()
            ns=float(re.search(r'Critical-path latency\s*:\s*([\d.]+) cycles',network).group(1))
            ls=float(re.search(r'Critical-path latency\s*:\s*([\d.]+) cycles',layer).group(1))
            assert ns==ls
            assert abs(float(r[f'npusim_{v}_us'])-ls/1500)<0.006
            assert int(re.search(r'# of computations\s*:\s*(\d+)',layer).group(1))==M*K*K
            resid=re.search(r'Tensor residency\s*:\s*(.*)',layer).group(1).strip()
            assert r[f'residency_{v}']==resid
            e=100*(ls/1500/float(r['measured_us'])-1)
            assert abs(e-float(r[f'err_{v}_pct']))<0.06
            r[f'exact_error_{v}']=e
            if v=='A':
                compute=float(re.search(r'Computation cycle\s*:\s*([\d.]+)',layer).group(1))
                dram=float(re.search(r'\* DRAM\s*:\s*([\d.]+) cycles',layer).group(1))
                r['compute_dram_ratio']=compute/dram
        assert r['residency_A']==r['residency_B']==r['residency_C']
        delta=float(r['npusim_C_us'])-float(r['npusim_B_us'])
        assert -0.011<=delta<=int(r['windows'])*590/1500+0.011
        r['regime']=('compute/DRAM<0.8' if r['compute_dram_ratio']<.8 else 'compute/DRAM>1.25' if r['compute_dram_ratio']>1.25 else 'compute/DRAM[0.8,1.25]')
        r['edge']='partial' if float(r['window_coverage'])>1 else 'full'
    assert seen==set(measure)
    groups={'all':rows}
    for field in ('dtype','edge','regime','residency_C'):
        for value in sorted(set(r[field] for r in rows)):
            groups[f'{field}={value}']=[r for r in rows if r[field]==value]
    metrics={}
    for key,rs in groups.items():
        metrics[key]={'n':len(rs), **{v:{'mape':st.mean(abs(r[f'exact_error_{v}']) for r in rs),'bias':st.mean(r[f'exact_error_{v}'] for r in rs)} for v in 'ABC'}}
    worst=sorted(rows,key=lambda r:abs(r['exact_error_C'])-abs(r['exact_error_A']),reverse=True)
    regressions=[{'case':f"{r['dtype']}:{r['M']}:{r['K']}",'delta_abs_error_pp':abs(r['exact_error_C'])-abs(r['exact_error_A']),'errors':{v:r[f'exact_error_{v}'] for v in 'ABC'},'residency':r['residency_C']} for r in worst[:8]]
    out[tag]={'artifact_checks':'PASS 84 inputs, 168 reports, measurements and hashes','groups':metrics,'worst_regressions':regressions,'improved_C_vs_A':sum(abs(r['exact_error_C'])<abs(r['exact_error_A']) for r in rows)}
    print(tag,json.dumps(out[tag],indent=2))
Path('/tmp/npusim_codex_full_batch_review.json').write_text(json.dumps(out,indent=2))

