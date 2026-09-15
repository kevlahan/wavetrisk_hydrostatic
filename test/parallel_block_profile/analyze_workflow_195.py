"""Validate common workflow envelopes and report paired excess and observer effects."""
from pathlib import Path
import json,math,statistics
PHASES=('setup','advance','physics','adapt_reconstruct','remap','output_restart')

def collect(root,ranks,mode):
 files=list(Path(root).glob('workflow195-rank-*.txt'))
 if len(files)!=ranks:raise ValueError('Missing workflow rank files')
 data={}
 for file in files:
  rank=int(file.stem.rsplit('-',1)[1]);steps={};seen=set();events=[];header=None
  if rank in data or not 0<=rank<ranks:raise ValueError('Invalid workflow rank')
  for line in file.read_text().splitlines():
   f=line.split()
   if not f:continue
   if f[0]=='schema':
    if header is not None:raise ValueError('Duplicate schema')
    header=list(map(int,f[1:]));continue
   if f[0]=='step':
    seq,restart,remap=map(int,f[1:4])
    if seq in steps:raise ValueError('Duplicate step')
    steps[seq]={'restart':restart,'remap':remap,'wall':float(f[4]),'phases':{},'transactions':{}}
   elif f[0]=='phase':
    seq,p=map(int,f[1:3]);key=('p',seq,p)
    if key in seen or not 0<=p<6:raise ValueError('Duplicate/invalid phase')
    seen.add(key);steps[seq]['phases'][p]=float(f[3])
   elif f[0]=='transaction':
    seq,p,tag,calls=map(int,f[1:5]);key=('t',seq,p,tag)
    if key in seen or not 0<=p<6 or tag not in (1,2) or calls<0:raise ValueError('Invalid transaction')
    seen.add(key);steps[seq]['transactions'][(p,tag)]={'calls':calls,'seconds':float(f[5])}
   elif f[0]=='event':events.append((list(map(int,f[1:5])),list(map(float,f[5:7]))))
   else:raise ValueError('Unknown workflow record')
  if header!=[1,mode,rank] or not steps:raise ValueError('Wrong workflow mode/schema')
  if set(steps)!=set(range(1,len(steps)+1)):raise ValueError('Noncontiguous steps')
  for seq,s in steps.items():
   if s['restart'] not in (0,1) or s['remap'] not in (0,1):raise ValueError('Invalid workflow flags')
   if set(s['phases'])!=set(range(6)) or len(s['transactions'])!=12:raise ValueError('Incomplete phase/transaction coverage')
   times=[s['wall'],*s['phases'].values(),*(v['seconds'] for v in s['transactions'].values())]
   if not all(math.isfinite(v) and v>=0 for v in times):raise ValueError('Invalid workflow time')
   if abs(sum(s['phases'].values())-s['wall'])>1e-6:raise ValueError('Workflow conservation failure')
   for (p,tag),v in s['transactions'].items():
    if v['seconds']>s['phases'][p]+1e-6:raise ValueError('Transaction exceeds containing phase')
   if mode==0 and (any(s['phases'][p] for p in range(1,6)) or any(v['calls'] for v in s['transactions'].values())):
    raise ValueError('Minimal mode contains phase instrumentation')
  if mode!=2 and events:raise ValueError('Unexpected event trace')
  if mode==2:
   totals={};counts={};phase_events={}
   for expected,(info,times) in enumerate(events,1):
    seq,tag,p,ordinal=info;a,b=times
    if ordinal!=expected or seq not in steps or tag not in (0,1,2) or p not in range(6):raise ValueError('Invalid event identity')
    if not all(math.isfinite(t) for t in times) or a<0 or b<a or b>steps[seq]['wall']+1e-6:raise ValueError('Invalid event interval')
    key=seq,p,tag;totals[key]=totals.get(key,0.)+b-a;counts[key]=counts.get(key,0)+1
    if tag==0:phase_events.setdefault(seq,[]).append((p,a,b))
   for seq,s in steps.items():
    ev=sorted(phase_events.get(seq,[]))
    if len(ev)!=6 or [e[0] for e in ev]!=list(range(6)):raise ValueError('Missing outer events')
    if abs(ev[0][1])>1e-6 or any(abs(ev[i][2]-ev[i+1][1])>1e-6 for i in range(5)):raise ValueError('Noncontiguous phase events')
    for info,times in events:
     if info[0]==seq and info[1]>0:
      p=info[2]
      if times[0]<ev[p][1]-1e-6 or times[1]>ev[p][2]+1e-6:raise ValueError('Transaction outside phase interval')
    for p in range(6):
     if abs(totals[(seq,p,0)]-s['phases'][p])>1e-6:raise ValueError('Event/phase times disagree')
     for tag in (1,2):
      v=s['transactions'][(p,tag)]
      if counts.get((seq,p,tag),0)!=v['calls'] or abs(totals.get((seq,p,tag),0)-v['seconds'])>1e-6:raise ValueError('Event/transaction totals disagree')
  data[rank]=steps
 if any(set(v)!=set(data[0]) for v in data.values()):raise ValueError('Unmatched workflow steps')
 rows=[]
 for seq in data[0]:
  values=[data[r][seq] for r in range(ranks)]
  if len({(v['restart'],v['remap']) for v in values})!=1:raise ValueError('Rank event flags disagree')
  slow=max(range(ranks),key=lambda r:values[r]['wall'])
  rows.append({'sequence':seq,'restart':values[0]['restart'],'remap':values[0]['remap'],
   'wall_mean':statistics.mean(v['wall'] for v in values),'wall_max':values[slow]['wall'],'slow_rank':slow,
   'phase_mean':[statistics.mean(v['phases'][p] for v in values) for p in range(6)],
   'phase_on_slow_rank':[values[slow]['phases'][p] for p in range(6)],
   'transaction_mean':{str(tag):{'calls_global':sum(v['transactions'][(p,tag)]['calls'] for v in values for p in range(6)),
     'seconds':statistics.mean(sum(v['transactions'][(p,tag)]['seconds'] for p in range(6)) for v in values)} for tag in (1,2)}})
 ordinary=[s for s in rows if not s['restart']]
 if not ordinary:raise ValueError('No ordinary workflow steps')
 return {'mode':mode,'ranks':ranks,'steps':rows,
 'ordinary_wall_max_median':statistics.median(s['wall_max'] for s in ordinary),
 'ordinary_rank_mean_wall':statistics.mean(s['wall_mean'] for s in ordinary),
 'ordinary_phase_mean':{n:statistics.mean(s['phase_mean'][p] for s in ordinary) for p,n in enumerate(PHASES)},
 'scope':'Same observer and outer workflow boundaries. Max is selected once per step, not independently per phase. Transaction 1 is WT_after_step; transaction 2 is implementation-specific tendency preparation and must not be assumed equivalent. Local event clocks are not synchronized across ranks.'}

def summarize(runs):
 out={};pairs={}
 for mode in ('bare','minimal','coarse','detail'):
  selected={}
  for side in ('block','legacy'):
   rows=[v for v in runs.values() if v['phase']=='timing-'+mode and v['side']==side]
   if not rows:raise ValueError('Missing timing mode '+mode)
   selected[side]=rows
  if len(selected['block'])!=len(selected['legacy']):raise ValueError('Unequal paired repeats')
  out[mode]={}
  for side,rows in selected.items():
   out[mode][side]={'elapsed_median':statistics.median(r['elapsed_seconds'] for r in rows),
    'printed_step_median':statistics.median(r['ordinary_step_median'] for r in rows),'repeats':len(rows)}
   if mode!='bare':out[mode][side]['wall_step_median']=statistics.median(r['workflow']['ordinary_wall_max_median'] for r in rows)
  pairs[mode]=[]
  for b,l in zip(selected['block'],selected['legacy']):
   if b['round']!=l['round']:raise ValueError('Unmatched round')
   pair={'round':b['round'],'elapsed_gap':b['elapsed_seconds']-l['elapsed_seconds'],
         'printed_step_gap':b['ordinary_step_median']-l['ordinary_step_median']}
   if mode!='bare':pair['wall_step_gap']=b['workflow']['ordinary_wall_max_median']-l['workflow']['ordinary_wall_max_median']
   if mode in ('coarse','detail'):
    bsteps=b['workflow']['steps'];lsteps=l['workflow']['steps']
    if [(s['sequence'],s['restart'],s['remap']) for s in bsteps]!=[(s['sequence'],s['restart'],s['remap']) for s in lsteps]:raise ValueError('Different phase workload')
    for bs,ls in zip(bsteps,lsteps):
     if bs['transaction_mean']['1']['calls_global']!=ls['transaction_mean']['1']['calls_global']:raise ValueError('Unmatched complete wavelet calls')
    gap=b['workflow']['ordinary_rank_mean_wall']-l['workflow']['ordinary_rank_mean_wall']
    delta={n:b['workflow']['ordinary_phase_mean'][n]-l['workflow']['ordinary_phase_mean'][n] for n in PHASES}
    residual=gap-sum(delta.values())
    if abs(residual)>2e-6:raise ValueError('Excess does not reconcile')
    pair.update(rank_mean_gap=gap,phase_excess=delta,residual=residual,
      whole_window_transaction_calls={side:{str(tag):sum(s['transaction_mean'][str(tag)]['calls_global'] for s in steps) for tag in (1,2)} for side,steps in [('block',bsteps),('legacy',lsteps)]})
   pairs[mode].append(pair)
 calibration={side:{mode:100*(out[mode][side]['wall_step_median']/out['minimal'][side]['wall_step_median']-1)
   for mode in ('coarse','detail')} for side in ('block','legacy')}
 return {'modes':out,'paired_differences':pairs,'observer_increase_percent':calibration,
  'status':'measured_phase_accounting_requires_causal_review',
  'scope':'Bare executables are authoritative timing references; minimal/coarse/detail share the same diagnostic binary per side. Positive signed differences mean block excess. Phase accounting is not an estimate of removable cost. Inspect phase calls, rank-local events, detailed counters and production source before claiming a cause.'}
