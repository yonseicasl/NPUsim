import csv,gzip,hashlib,json,re,statistics as st
from pathlib import Path
from zipfile import ZipFile
V=Path('/home/bogil/MOSAIC/.claude/worktrees/npusim-engine/simulator/3rd/npusim/.claude/worktrees/tpu-measure/validation/tpu_v5e')
src=Path('/home/bogil/measure_tpu_v5e1.zip')
rows={(r['dtype'],int(r['M']),int(r['K'])):r for r in csv.DictReader((V/'results/2026-11-12_colab_v5e1_balance/balance.csv').open()) if r['schedule']=='loop'}
result=[]
with ZipFile(src) as z:
    names=[n for n in z.namelist() if 'balance_' in n and '_loop/' in n and n.endswith('perfetto_trace.json.gz')]
    assert len(names)==28
    for n in names:
        m=re.search(r'balance_(bf16|int8)_(\d+)x(\d+)_loop/',n); key=(m[1],int(m[2]),int(m[3]))
        raw=z.read(n); contents=json.loads(gzip.decompress(raw)); ev=contents.get('traceEvents',[]) if isinstance(contents,dict) else contents
        proc={e['pid']:e['args']['name'] for e in ev if e.get('ph')=='M' and e.get('name')=='process_name'}
        threads={(e['pid'],e.get('tid')):e['args']['name'] for e in ev if e.get('ph')=='M' and e.get('name')=='thread_name'}
        candidates=sorted((p for p,name in proc.items() if 'TPU' in name),key=lambda p:proc[p]);assert candidates
        pid=candidates[0]
        events=[e for e in ev if e.get('ph')=='X' and e.get('pid')==pid]
        mods=[e for e in events if 'Module' in threads.get((pid,e.get('tid')),'')]
        conv=[e for e in events if 'Ops' in threads.get((pid,e.get('tid')),'') and 'convolution' in str(e.get('args',{}).get('hlo_category',''))]
        assert mods and conv
        counts=[0]*len(mods);bad=[]
        for e in conv:
            containers=[i for i,mod in enumerate(mods) if float(e['ts'])>=float(mod['ts'])-1e-6 and float(e['ts'])+float(e['dur'])<=float(mod['ts'])+float(mod['dur'])+1e-6]
            if len(containers)!=1:bad.append((e['name'],len(containers)))
            else:counts[containers[0]]+=1
        assert not bad,(key,bad[:3])
        assert set(counts)=={8},(key,{c:counts.count(c) for c in set(counts)})
        observed=sum(float(e['dur']) for e in conv)/len(mods)/8
        recorded=float(rows[key]['fused_ops_sum_per_iter_us'])
        assert abs(observed-recorded)<=.00051,(key,observed,recorded)
        result.append(dict(case=key,modules=len(mods),conv_events=len(conv),events_each_module=8,reconstructed_us=observed,csv_us=recorded,trace_member=n,trace_sha256=hashlib.sha256(raw).hexdigest()))
assert len({tuple(r['case']) for r in result})==28
out=dict(zip_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),checks='PASS all 28 original TPU traces: exactly8 convolution events contained in every module, each event uniquely assigned, reconstructed time matches CSV rounding',cases=result)
Path('/tmp/npusim_codex_balance_trace_audit.json').write_text(json.dumps(out,indent=2))
print(out['checks'])
print('Total modules:',sum(r['modules'] for r in result),'total convolution events:',sum(r['conv_events'] for r in result))

