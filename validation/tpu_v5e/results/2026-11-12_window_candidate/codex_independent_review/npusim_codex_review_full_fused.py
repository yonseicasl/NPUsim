import csv, hashlib, json, math, re, statistics as st
from pathlib import Path
P=Path('/home/bogil/MOSAIC/.claude/worktrees/npusim-engine/simulator/3rd/npusim/.claude/worktrees/tpu-measure')
V=P/'validation/tpu_v5e'; O=V/'results/2026-11-12_window_candidate'; G=O/'full_fused_output'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for line in (G/'FROZEN.sha256').read_text().splitlines():
    digest,rel=line.split(maxsplit=1)
    assert sha(Path('/tmp/npusim_codex_full_fused_binary')/Path(rel).name)==digest
cases={(dt,M,K) for dt in ('bf16','int8') for M in (32,64,128,192,256,384,512) for K in (4096,8192)}
find=lambda s,p:re.search(p,s).group(1)
def tx(s,title):
    block=s[s.index(title):].splitlines()[1:4]
    return tuple(int(find(l,r':\s*(\d+)/')) for l in block)
out={}
for tag in ('analytic','hbm2e_tpu_v5e_refsb'):
    key=tag+'_fusedout'; f=G/f'balance_npusim_{key}.csv'
    if not f.exists(): print(key,'PENDING'); continue
    rows=list(csv.DictReader(f.open())); assert len(rows)==28
    meta=json.loads(f.with_suffix('.json').read_text())
    assert meta['activations']=='fused_output'
    assert meta['frozen']['binary_sha256']==sha(Path('/tmp/npusim_codex_full_fused_binary/model'))
    assert meta['frozen']['library_sha256']==sha(Path('/tmp/npusim_codex_full_fused_binary/libnpusim.so'))
    assert meta['window_switch_cycle']['value']==590
    for n,h in meta['frozen']['configs'].items(): assert h==sha(P/'configs/accelerators'/f'{n}.cfg')
    source=Path(meta['measurement']['run'])
    assert meta['balance_csv_sha256']==sha(source/'balance.csv')
    for n,h in (meta['dram'] or {}).items(): assert h==sha(Path(n))
    gen=G/f'generated_{key}'; oldgen=O/f'generated_{tag}'
    assert len(list(gen.iterdir()))==84 and len(list(gen.rglob('*.txt')))==168
    oldrows={(r['dtype'],int(r['M']),int(r['K'])):r for r in csv.DictReader((O/f'balance_npusim_{tag}.csv').open())}
    seen=set(); checked=[]
    for r in rows:
        dt,M,K=r['dtype'],int(r['M']),int(r['K']); shape=(dt,M,K)
        assert shape not in seen; seen.add(shape); oldrow=oldrows[shape]
        assert r['hlo_sha256']==oldrow['hlo_sha256']
        assert r['measured_us']==oldrow['measured_us']
        rs={'case':shape,'old_errors':{},'new_errors':{},'bytes':None}
        for v in 'ABC':
            name=f'bal_{dt}_{M}x{K}_{v}'.lower(); d=gen/name
            net=(d/'network.cfg').read_text(); val=lambda k:int(find(net,rf'(?m)^{k}=(\d+)'))
            assert (val('batch'),val('height')*val('width')*val('channels'),val('output'))==(M,K,K)
            assert f'weight=weights/{name}.wgh' in net
            mp=[[int(t) for t in l.split('=')[1].split(',') if t.strip()] for l in next(d.rglob('ws.map')).read_text().splitlines() if '=' in l]
            assert tuple(math.prod(x[c] for x in mp) for c in (1,4,0))==(M,K,K)
            cfg=(d/f'{name}.cfg').read_text()
            assert find(cfg,r'(?m)^activation_residency\s*=\s*(\d+)')=='0'
            assert find(cfg,r'(?m)^fused_output_consumer\s*=\s*(\d+)')=='1'
            assert find(cfg,r'(?m)^program_setup_cycle\s*=\s*(\d+)')=='0'
            if v!='A': assert find(cfg,r'(?m)^offchip_window_count\s*=\s*(\d+)')==r['windows']
            if v=='C': assert find(cfg,r'(?m)^window_switch_cycle\s*=\s*(\d+)')=='590'
            report=(d/f'{name}_{name}_layer_0.txt').read_text()
            old=(oldgen/name/f'{name}_{name}_layer_0.txt').read_text()
            assert int(find(report,r'# of computations\s*:\s*(\d+)'))==M*K*K
            assert find(report,r'Tensor residency\s*:\s*(.*)').startswith('fused_output_consumer:')
            cycles=float(find(report,r'Critical-path latency\s*:\s*([\d.]+)'))
            assert abs(float(r[f'npusim_{v}_us'])-cycles/1500)<.006
            bits=tuple(map(int,find(cfg[cfg.index('[dram]'):],r'line_size\s*=\s*([\d:]+)').split(':')))
            counts=tx(report,'Multi-chip <-> DRAM transactions')
            actual=tuple(a*b//8 for a,b in zip(counts,bits))
            assert tuple(map(int,r['dram_line_bits'].split(':')))==bits
            assert tuple(int(r[k]) for k in ('dram_input_bytes','dram_weight_bytes','dram_output_bytes'))==actual
            assert actual==(M*K*2,K*K*(1 if dt=='int8' else 2),0),actual
            assert tx(report,'PE-array <-> GLB transactions')==tx(old,'PE-array <-> GLB transactions')
            assert find(report,r'Psum round trip\s*:\s*(.*)')==find(old,r'Psum round trip\s*:\s*(.*)')
            assert find(report,r'Computation cycle\s*:\s*([\d.]+)')==find(old,r'Computation cycle\s*:\s*([\d.]+)')
            assert find(report,r'Fold fill cycle\s*:\s*([\d.]+)')==find(old,r'Fold fill cycle\s*:\s*([\d.]+)')
            rs['old_errors'][v]=100*(float(find(old,r'Critical-path latency\s*:\s*([\d.]+)'))/1500/float(r['measured_us'])-1)
            rs['new_errors'][v]=100*(cycles/1500/float(r['measured_us'])-1)
            rs['bytes']=actual
        rs['dtype']=dt
        rs['old_resident']='resident' if 'input GLB-resident' in oldrow['residency_C'] else 'fallback'
        rs['edge']='partial' if float(r['window_coverage'])>1 else 'full'
        ratio=(2*M*K*K/196.6e12)/(K*K*(1 if dt=='int8' else 2)/750e9)
        rs['nominal_regime']='C/S<0.6' if ratio<.6 else 'C/S>1.6' if ratio>1.6 else 'C/S[0.6,1.6]'
        checked.append(rs)
    assert seen==cases
    groups={'all':checked}
    for field in ('dtype','edge','old_resident','nominal_regime'):
        for value in sorted({r[field] for r in checked}):
            groups[field+'='+value]=[r for r in checked if r[field]==value]
    metrics={}
    for group,rs in groups.items():
        metrics[group]={'n':len(rs)}
        for mode in ('old_errors','new_errors'):
            metrics[group][mode]={v:{'mape':st.mean(abs(r[mode][v]) for r in rs),'bias':st.mean(r[mode][v] for r in rs)} for v in 'ABC'}
    regressions=[{'case':r['case'],'A':r['new_errors']['A'],'C':r['new_errors']['C'],'delta_abs_pp':abs(r['new_errors']['C'])-abs(r['new_errors']['A'])} for r in checked]
    regressions.sort(key=lambda r:r['delta_abs_pp'],reverse=True)
    out[key]={'checks':'PASS 84 input sets, 168 reports, hashes, exact I/W/O bytes, preserved on-chip transactions/psum/compute/fill','rows':checked,'metrics':metrics,'worst_regressions':regressions[:6]}
    print(key,json.dumps({k:v for k,v in out[key].items() if k!='rows'},indent=2))
Path('/tmp/npusim_codex_full_fused_review.json').write_text(json.dumps(out,indent=2))

