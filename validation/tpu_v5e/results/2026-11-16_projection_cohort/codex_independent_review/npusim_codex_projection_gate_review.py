import csv,json,re,hashlib,math,sys
from pathlib import Path
P=Path('/home/bogil/MOSAIC/.claude/worktrees/npusim-engine/simulator/3rd/npusim/.claude/worktrees/tpu-measure')
V=P/'validation/tpu_v5e'; D=V/'results/2026-11-16_projection_cohort'; G=D/'gate'
sys.path.insert(0,str(V))
from trace_bytes import fields,stat,breakdown
# Share only protobuf wire/stat decoding; independently reconstruct intervals, counts and means.
def raw_audit(pb):
 planes=[list(fields(v)) for f,w,v in fields(pb.read_bytes()) if f==1]
 planes=[p for p in planes if b'TPU' in next((v for f,w,v in p if f==2),b'')]
 assert len(planes)==1
 p=planes[0]; names={}; meta={}
 for f,w,v in p:
  if f==5:
   e=dict((a,b) for a,_,b in fields(v)); m=dict((a,b) for a,_,b in fields(e[2]))
   names[e[1]]=m.get(2,b'').decode()
 for f,w,v in p:
  if f==4:
   e=dict((a,b) for a,_,b in fields(v)); md=list(fields(e[2]))
   meta[e[1]]=dict(stat(b,names) for a,_,b in md if a==5)
 ids={k for k,s in meta.items() if b'convolution' in s.get('hlo_category',b'')}
 mods=[]; conv=[]
 for f,w,v in p:
  if f!=3: continue
  ln=list(fields(v)); name=next((b.decode() for a,_,b in ln if a==2),'')
  if 'Module' not in name and not ('Ops' in name and 'Async' not in name): continue
  ts=next((b for a,_,b in ln if a==3),None)
  if not any(a==4 for a,_,b in ln): continue
  ts = 0 if ts is None else ts  # protobuf scalar default; captures use offsets relative to zero
  for a,_,b in ln:
   if a!=4: continue
   e=dict((x,y) for x,_,y in fields(b))
   if 'Module' not in name and e.get(1) not in ids: continue
   assert 3 in e and e[3]>0
   # Proto3 zero offsets may be omitted; duration must be positive.
   start=ts*1000+e.get(2,0); interval=(start,start+e[3])
   (mods if 'Module' in name else conv).append(interval)
 if not ids: return dict(ids=0,conv=len(conv),modules=len(mods))
 assert len(ids)==1 and mods
 counts=[0]*len(mods)
 for a,b in conv:
  owners=[i for i,(c,d) in enumerate(mods) if c<=a and b<=d]
  assert len(owners)==1
  counts[owners[0]]+=1
 assert set(counts)=={8}
 bd=breakdown(meta[next(iter(ids))]['memory_access_breakdown'])
 return dict(ids=1,conv=len(conv),modules=len(mods),us=sum(b-a for a,b in conv)/len(conv)/1e6,
             read=sum(n for op,space,n in bd if op==1 and space==1),
             write=sum(n for op,space,n in bd if op==2 and space==1),
             other=sum(n for op,space,n in bd if space!=1))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
mj=json.loads((D/'manifest.json').read_text()); src=Path(mj['run'])
for f,h in mj['source_files'].items(): assert sha(src/f)==h
for f,h in mj['raw_archive'].items(): assert sha(src.parent/f)==h
rows=list(csv.DictReader((D/'manifest.csv').open())); assert len(rows)==40
audits=[]; support=0
for r in rows:
 assert sha(src/r['hlo'])==r['hlo_sha256']
 for f,h in json.loads(r['trace_files']).items(): assert sha(src/'traces'/r['tag']/f)==h
 pbs=list((src/'traces'/r['tag']).rglob('*.xplane.pb')); assert len(pbs)==1
 x=raw_audit(pbs[0]); x['tag']=r['tag']; audits.append(x)
 if r['supported']=='True':
  support+=1; M,K,N=map(int,(r['M'],r['K'],r['N']))
  assert x['modules']==int(r['modules'])==1 and x['conv']==8
  assert abs(x['us']-float(r['conv_event_us']))<1e-10
  assert x['read']==M*K*2+K*N*(1 if r['dtype']=='int8' else 2)
  assert x['write']==4 and x['other']==0
 else: assert int(r['M'])==1 and x['ids']==0
assert support==32
results=[]; run_count=0
for cp in G.glob('balance_npusim_*.csv'):
 tag=cp.stem.removeprefix('balance_npusim_'); meta=json.loads(cp.with_suffix('.json').read_text())
 assert meta['manifest_sha256']==sha(D/'manifest.csv')
 for f,h in (meta.get('dram') or {}).items(): assert sha(f)==h
 fr=meta['frozen']
 assert sha(P/'models/model')==fr['binary_sha256']
 assert sha(P/'library/libnpusim.so')==fr['library_sha256']
 assert sha(G/'balance_npusim.py.gate_snapshot')==fr['driver_sha256']
 for n,h in fr['configs'].items(): assert sha(P/'configs/accelerators'/f'{n}.cfg')==h
 errs=[]
 for r in csv.DictReader(cp.open()):
  M,K,N=map(int,(r['M'],r['K'],r['N'])); group=[]
  for v in 'ABC':
   name=f"prj_{r['dtype']}_{M}x{K}x{N}_{v.lower()}"; d=G/f'generated_{tag}'/name
   net=(d/'network.cfg').read_text()
   val=lambda k:int(re.search(rf'(?m)^{k}=(\d+)',net)[1])
   assert (val('batch'),val('channels')*val('height')*val('width'),val('output'))==(M,K,N)
   mapping=list((d/'mappings').rglob('ws.map')); assert len(mapping)==1
   ms=[[int(x) for x in l.split('=')[1].split(',') if x.strip()] for l in mapping[0].read_text().splitlines() if '=' in l]
   assert tuple(math.prod(x[c] for x in ms) for c in (1,4,0))==(M,K,N)
   t=(d/f'{name}_{name}_layer_0.txt').read_text(); cfg=(d/f'{name}.cfg').read_text()
   number=lambda key:float(re.search(re.escape(key)+r'\s*:\s*([\d.]+)',t)[1])
   crit=number('Critical-path latency'); us=crit/1500
   assert abs(us-float(r[f'npusim_{v}_us']))<=0.0051
   assert number('# of computations')==M*K*N
   assert 'fused_output_consumer = 1' in cfg
   bits=list(map(int,re.search(r'(?m)^line_size\s*=\s*([\d:]+)',cfg.split('[dram]')[1])[1].split(':')))
   block=t.split('Multi-chip <-> DRAM transactions')[1]
   tx=[int(re.search(r'\* '+key+r'\s*:\s*(\d+)',block)[1]) for key in ('Input data','Weight','Output data')]
   assert [a*b//8 for a,b in zip(tx,bits)]==[M*K*2,K*N*(1 if r['dtype']=='int8' else 2),0]
   busy=t.split('Busy cycles (ratio of critical path)')[1].split('Bottleneck level')[0]
   busyvals=[float(x) for x in re.findall(r':\s*([\d.]+) cycles',busy)]; assert len(busyvals)==5
   group.append((v,crit,100*(us/float(r['measured_us'])-1),busyvals))
   run_count+=1
  for level in range(4): assert len({g[3][level] for g in group})==1
  assert abs(group[2][3][4]-group[1][3][4]-int(r['windows'])*590)<0.11
  errs.append(dict(case=f"{r['dtype']}:{M}:{K}:{N}",calibration=r['w_calibration_case']=='True',errors={g[0]:g[2] for g in group}))
 non=[x for x in errs if not x['calibration']]
 results.append(dict(dram=tag,cases=errs,noncal_mape={v:sum(abs(x['errors'][v]) for x in non)/len(non) for v in 'ABC'}))
out=dict(status='PASS',limits='Reused trace_bytes protobuf wire/stat decoding only; reconstructed event intervals/aggregation separately. Report times have 0.1-cycle print precision; no TPU rerun or new simulator runs.',manifest_cases=40,supported=32,unsupported=8,conv_events=sum(x['conv'] for x in audits),runs=run_count,metrics=results,trace_audit=audits)
Path('/tmp/npusim_codex_projection_gate_review.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='trace_audit'},indent=2))
