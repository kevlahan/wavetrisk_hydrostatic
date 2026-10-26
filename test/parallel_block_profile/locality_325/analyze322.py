"""Whole-step workload sweep; integer census is not a phase-time attribution."""
from pathlib import Path
import json,math,statistics as st,sys

def census(directory,n,steps,block):
 root=Path(directory);times=[];work=[];startup=[];states=[]
 env=json.loads((root/'environment.json').read_text())
 bridge=env.get('WAVETRISK_OWNER_MAP322','off')=='B256'
 if bridge and (n!=256 or not block or steps!=4):raise ValueError('Invalid empty-owner scope')
 expected={f'campaign312-rank-{r}.txt' for r in range(n)}
 if {p.name for p in root.glob('campaign312-rank-*.txt')}!=expected:raise ValueError('Missing rank census')
 for rank in range(n):
  lines=[s.split() for s in (root/f'campaign312-rank-{rank}.txt').read_text().splitlines()]
  h=lines[0]
  if len(h)!=5 or h[0]!='H' or list(map(int,h[1:4]))!=[rank,steps,4*steps if block else 0]:raise ValueError('Rank census header')
  startup.append(float(h[4]));ts=[];ss=[];ws=[[0]*5 for _ in range(steps)];last_key=None;last=[0]*5
  rr=[]
  for i in range(steps):
   f=lines[1+i]
   if len(f)!=6 or f[:2]!=['T',str(i+1)]:raise ValueError('Timestep census ordering')
   values=list(map(float,f[2:5]));active=int(f[5])
   if any(not math.isfinite(v) or v<=0 for v in values) or active<=0:raise ValueError('Invalid timestep metric')
   ts.append(values[0]);ss.append([values[1],values[2],active])
  for f in lines[1+steps:]:
   if len(f)!=13 or f[0]!='R':raise ValueError('Malformed RHS census')
   v=list(map(int,f[1:]));step,gen,mask,nd,groups,prim,scalar,vel,rv,rm,nnode,call=v
   if any(x<0 for x in v) or not 1<=step<=steps:raise ValueError('Invalid RHS census')
   if bridge and rank>=160:
    if any(x!=0 for x in (nd,groups,prim,scalar,vel,rv,rm,nnode)):raise ValueError('Extra rank has numerical domains/groups/work/payload')
   elif groups<1 or nnode<1:raise ValueError('Unexpected empty numerical owner')
   key=(gen,mask);cur=[prim,scalar,vel,rv,rm]
   delta=[a-b for a,b in zip(cur,last)] if key==last_key else cur
   if any(x<0 for x in delta):raise ValueError('Counters decreased within preparation epoch')
   ws[step-1]=[a+b for a,b in zip(ws[step-1],delta)];last_key=key;last=cur;rr.append(v)
  if block and [r[0] for r in rr]!=[i for i in range(1,steps+1) for _ in range(4)]:raise ValueError('Four RK stages per step required')
  if states and ss!=states[0]:raise ValueError('Rank physical state metadata differs')
  times.append(ts);states.append(ss);work.append(ws)
 data=[]
 for i in range(steps):
  row=dict(step=i+1,rank0_seconds=times[0][i],max_rank_seconds=max(t[i] for t in times),physical_time=states[0][i][0],dt=states[0][i][1],active_wavelets=states[0][i][2])
  if block:
   totals=[sum(w[i][j] for w in work) for j in range(5)];peak=[max(w[i][j] for w in work) for j in range(5)]
   row.update(global_work=totals[:3],peak_rank_work=peak[:3],mean_rank_work=[t/n for t in totals[:3]],peak_to_mean=[p*n/t if t else None for p,t in zip(peak[:3],totals[:3])],exported_values=totals[3],send_messages=totals[4])
  data.append(row)
 groups={}
 if block:
  for label in ('initial','final'):
   rows=[line.split() for line in (root/f'groups312-{label}.txt').read_text().splitlines()];h=rows[0]
   if h[0]!='H' or int(h[1])!=n or len(rows)!=1+int(h[2])+int(h[3]):raise ValueError('Execution-group snapshot extent')
   fine=[list(map(int,f[1:])) for f in rows if f[0]=='L'];physical=[list(map(int,f[1:])) for f in rows if f[0]=='E']
   if len(fine)!=int(h[2]) or len(physical)!=int(h[3]):raise ValueError('Execution group membership')
   expected_owners=set(range(160 if bridge else n))
   if {v[4] for v in physical}!=expected_owners:raise ValueError('Execution owner set differs')
   if bridge and any(v[4]!=(i//4) for i,v in enumerate(physical)):raise ValueError('B owner map is not exact160-owner map')
   if any(v[4]!=physical[v[6]-1][4] for v in fine):raise ValueError('Coalescing changed ownership')
   roots={};current=None
   for v in fine:
    if v[3]==5:current=(v[1],v[2]);roots[current]=set()
    if current is None or current[0]!=v[1]:raise ValueError('Original-root traversal incomplete')
    roots[current].add(v[4])
   groups[label]=dict(fine_units=len(fine),execution_groups=len(physical),original_roots=len(roots),mixed_original_roots=sum(len(x)>1 for x in roots.values()),rank_group_counts=[sum(v[4]==r for v in physical) for r in range(n)])
 return dict(steps=data,groups=groups,initialization_max_seconds=max(startup),first_step_seconds=data[0]['max_rank_seconds'] if steps else None,
             ordinary_mean_seconds=st.mean(x['max_rank_seconds'] for x in data[1:]) if steps>1 else None)

def snapshots(directory,n):
 root=Path(directory)
 expected={f'snapshot312-rank-{r}.txt' for r in range(n)}
 if {p.name for p in root.glob('snapshot312-rank-*.txt')}!=expected:raise ValueError('Missing snapshot preservation record')
 totals={9900:[0,0],9901:[0,0]}
 for rank in range(n):
  rows=[list(map(int,line.split())) for line in (root/f'snapshot312-rank-{rank}.txt').read_text().splitlines()]
  if len(rows)!=2 or any(len(v)!=3 for v in rows) or [v[0] for v in rows]!=[9900,9901] or any(x<0 for v in rows for x in v):raise ValueError('Snapshot preservation ordering/extent')
  for index,a,b in rows:totals[index][0]+=a;totals[index][1]+=b
 return {str(k):dict(restored_deletion_flags=v[0],restored_level_entries_or_lengths=v[1]) for k,v in totals.items()}
