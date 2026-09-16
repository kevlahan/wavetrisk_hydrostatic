"""Full-step disjoint operation ledger, calibration, and explicit unresolved costs."""
from pathlib import Path
import json,math,statistics as st
NR,NC=100,42
CONTEXTS={0:'dynamics remainder',1:'dynamics WT',2:'tendency tag',3:'adaptation WT',4:'adaptation remainder',5:'step setup',6:'physics',7:'remap',8:'output/restart'}
PHASES={'dynamics WT':[1],'dynamics outside WT':[0,2],'adaptation/reconstruction':[3,4],'step setup':[5],'physics':[6],'remap':[7],'output/restart':[8]}

def finite(x):return math.isfinite(x) and x>=0

def common(path):
 rows={}
 for line in path.read_text().splitlines():
  f=line.split()
  if f[0]=='step':rows[int(f[1])]={'restart':int(f[2]),'remap':int(f[3]),'wall':float(f[4]),'phase':{},'tx':{},'tc':{}}
  elif f[0]=='phase':rows[int(f[1])]['phase'][int(f[2])]=float(f[3])
  elif f[0]=='transaction':
   r=rows[int(f[1])];key=int(f[2]),int(f[3]);r['tx'][key]=float(f[5]);r['tc'][key]=int(f[4])
 for r in rows.values():
  p,t=r['phase'],r['tx'];r['context']={0:p[1]-t[1,1]-t[1,2],1:t[1,1],2:t[1,2],3:t[3,1],4:p[3]-t[3,1],5:p[0],6:p[2],7:p[4],8:p[5]}
  if any(not math.isfinite(x) or x< -1e-6 for x in r['context'].values()):raise ValueError('Overlapping common contexts')
  if abs(sum(r['context'].values())-r['wall'])>1e-6:raise ValueError('Common full-step conservation')
 return rows

def collect(root,ranks,mode,side=None):
 root=Path(root);files=list(root.glob('operations198-rank-*.txt'));data={};samples={};dependencies=[]
 if len(files)!=ranks:raise ValueError('Missing operations rank files')
 for file in files:
  rank=int(file.stem.rsplit('-',1)[1]);rows={};schema=None
  if rank in data or rank not in range(ranks):raise ValueError('Invalid/duplicate rank')
  ref=common(root/f'workflow195-rank-{rank}.txt')
  for line in file.read_text().splitlines():
   f=line.split();kind=f[0]
   if kind=='schema':
    if schema is not None:raise ValueError('Duplicate schema')
    schema=list(map(int,f[1:]));continue
   if kind=='step':
    n=int(f[1])
    if n in rows:raise ValueError('Duplicate step')
    rows[n]={'restart':int(f[2]),'remap':int(f[3]),'wall':float(f[4]),'self':{},'calls':{},'count':{},'volume':{}}
    continue
   if kind not in ('self','count','volume'):raise ValueError('Unknown operation record')
   n,c,p=map(int,f[1:4])
   if n not in rows or c not in range(9):raise ValueError('Invalid context/step')
   r=rows[n];key=c,p
   if key in r[kind]:raise ValueError('Duplicate operation record')
   if kind=='self':
    if len(f)!=6 or p not in range(NR+1) or int(f[4])<0 or not finite(float(f[5])):raise ValueError('Invalid operation value')
    r['self'][key]=float(f[5]);r['calls'][key]=int(f[4])
   elif kind=='count':
    if len(f)!=5 or p not in range(1,NC+1) or int(f[4])<0:raise ValueError('Invalid work count')
    r['count'][key]=int(f[4])
   else:
    if len(f)!=6 or p not in range(1,57) or min(map(int,f[4:]))<0:raise ValueError('Invalid volume')
    r['volume'][key]=list(map(int,f[4:]))
  if schema!=[2,mode,rank] or set(rows)!=set(ref) or set(rows)!=set(range(1,len(rows)+1)):raise ValueError('Incomplete operations schema/steps')
  for n,r in rows.items():
   q=ref[n]
   if r['restart']!=q['restart'] or r['remap']!=q['remap'] or not finite(r['wall']) or abs(r['wall']-q['wall'])>1e-6:raise ValueError('Full-step flags/wall mismatch')
   if len(r['self'])!=9*(NR+1) or len(r['count'])!=9*NC:raise ValueError('Missing dense operation records')
   if mode==0:
    if any(r['self'].values()) or any(r['calls'].values()) or any(r['count'].values()) or r['volume']:raise ValueError('Disabled observer recorded operations')
   else:
    for c in range(9):
     if abs(sum(r['self'][c,p] for p in range(NR+1))-q['context'][c])>1e-6:raise ValueError('Per-context conservation failure')
    if r['calls'][1,65]!=q['tc'][1,1] or r['calls'][3,65]!=q['tc'][3,1]:raise ValueError('WT coverage failure')
   if mode<2 and any(r['count'].values()):raise ValueError('Timing-only run recorded work counters')
   r['envelopes']=q['context']
  data[rank]=rows
  sf=root/f'storage198-rank-{rank}.txt';lines=sf.read_text().splitlines();head=lines[0].split()
  if len(head)!=5 or head[:3]!=['schema','1',str(rank)]:raise ValueError('Invalid storage schema')
  stride,offset=map(int,head[3:]);sr={}
  if stride!=257 or offset not in range(stride):raise ValueError('Invalid sample schedule')
  for line in lines[1:]:
   f=line.split()
   if len(f)!=8 or f[0]!='sample':raise ValueError('Invalid storage row')
   n,c,p,visits,takes=map(int,f[1:6]);seconds,square=map(float,f[6:]);key=n,c,p
   if key in sr or n not in rows or c not in range(9) or p not in range(1,5) or visits<0 or takes!=(visits+offset)//stride or not finite(seconds) or not finite(square):raise ValueError('Invalid storage sample values')
   if mode!=3 and (visits or takes or seconds or square):raise ValueError('Unexpected storage sampling')
   if takes==0 and (seconds or square):raise ValueError('Storage times without samples')
   if takes and square+1e-12<seconds*seconds/takes:raise ValueError('Storage moments inconsistent')
   sr[key]=(visits,takes,seconds,square)
  if len(sr)!=len(rows)*9*4:raise ValueError('Missing dense storage records')
  samples[rank]=sr
  lines=(root/f'dependencies198-rank-{rank}.txt').read_text().splitlines()
  if lines[0].split()!=['schema','1',str(rank)]:raise ValueError('Invalid dependency schema')
  for sequence,line in enumerate(lines[1:],1):
   f=line.split()
   if len(f)!=14:raise ValueError('Invalid dependency row')
   v=list(map(int,f[:12]));t=list(map(float,f[12:]))
   if v[11]!=sequence or v[0] not in rows or v[1] not in range(9) or v[2] not in (1,2,3) or min(v[7:10])<0 or not all(finite(x) for x in t):raise ValueError('Invalid dependency values')
   dependencies.append([rank,*v,*t])
  if mode<2 and len(lines)>1:raise ValueError('Unexpected dependencies in timing-only pass')
 if len({tuple(r) for r in data.values()})!=1:raise ValueError('Rank step coverage differs')
 for n in data[0]:
  if len({(data[k][n]['restart'],data[k][n]['remap']) for k in data})!=1:raise ValueError('Rank flags differ')
 if mode and side in ('block','legacy'):
  required=range(83,101) if side=='block' else (88,89,90,97,100)
  missing=[p for p in required if not sum(v['calls'][c,p] for rows in data.values() for v in rows.values() for c in range(9))]
  if missing:raise ValueError('Missing full-step scope coverage: '+repr(missing))
 if mode==3 and side=='block' and not sum(v[1] for rows in samples.values() for v in rows.values()):raise ValueError('Storage probe binary recorded no samples')
 ordinary=[n for n,r in data[0].items() if not r['restart']]
 if not ordinary:raise ValueError('No ordinary steps')
 if side=='block' and mode>=2:
  waits=sum(v['calls'][c,19] for rows in data.values() for v in rows.values() for c in range(9))
  if waits!=len(dependencies):raise ValueError('Incomplete inverse dependency coverage')
 if side=='legacy' and dependencies:raise ValueError('Unexpected legacy inverse implementation records')
 steps=[]
 for n in data[0]:
  slow=max(data,key=lambda k:data[k][n]['wall']);q=data[slow][n]
  steps.append({'sequence':n,'restart':q['restart'],'remap':q['remap'],'slow_rank':slow,'wall_max':q['wall'],'wall_mean':st.mean(data[k][n]['wall'] for k in data),'contexts':{c:{'mean':st.mean(data[k][n]['envelopes'][c] for k in data),'slow_rank':q['envelopes'][c],'calls_global':[sum(data[k][n]['calls'][c,p] for k in data) for p in range(NR+1)],'counters_global':[sum(data[k][n]['count'][c,p] for k in data) for p in range(1,NC+1)]} for c in range(9)}})
 ledger=[]
 for c in range(9):
  for p in range(NR+1):
   rankmeans=[st.mean(data[k][n]['self'][c,p] for n in ordinary) for k in data]
   ledger.append({'context':c,'operation':p,'mean':st.mean(rankmeans),'min_rank_mean':min(rankmeans),'median_rank_mean':st.median(rankmeans),'max_rank_mean':max(rankmeans),'slow_rank_path_mean':st.mean(data[s['slow_rank']][s['sequence']]['self'][c,p] for s in steps if not s['restart']),'calls_per_rank_step':sum(data[k][n]['calls'][c,p] for k in data for n in ordinary)/(ranks*len(ordinary))})
 storage=[]
 for c in range(9):
  for p in range(1,5):
   vals=[samples[k][n,c,p] for k in data for n in ordinary];tot=[sum(v[i] for v in vals) for i in range(4)]
   estimates=[v[0]*v[2]/v[1] for v in vals if v[1]]
   storage.append({'context':c,'operation':p,'calls':tot[0],'samples':tot[1],'sample_seconds':tot[2],'sample_square_seconds':tot[3],'cells_without_samples':sum(v[0]>0 and v[1]==0 for v in vals),'estimated_inclusive_rank_step_seconds':sum(estimates)/(ranks*len(ordinary)) if estimates else None})
 return {'schema':2,'mode':mode,'ranks':ranks,'side':side,'steps':steps,'ordinary_wall_mean':st.mean(s['wall_mean'] for s in steps if not s['restart']),'ordinary_slow_rank_wall_mean':st.mean(s['wall_max'] for s in steps if not s['restart']),'ledger':ledger,'work_per_rank_step':{c:[sum(data[k][n]['count'][c,p] for k in data for n in ordinary)/(ranks*len(ordinary)) for p in range(1,NC+1)] for c in range(9)},'dependency_records':len(dependencies),'storage':storage,'scope':'Self times are disjoint; actual slowest full-step rank is selected once per step. Storage samples are inclusive, nested and systematic; never add them to this ledger. Operation names do not by themselves prove matched mathematics or recoverable cost.'}

def summarize(runs):
 result={'timings':{},'calibration':{},'ranking':[],'scope':'Signed full-step location differences, not estimates of removable time. Unclassified and unmatched scopes remain explicit.'}
 for side in ('baseline','candidate','legacy'):
  values=[v for k,v in runs.items() if k.startswith(side+'-') and v['phase']=='timing-bare']
  result['timings'][side]={'ordinary_step_medians':[v['ordinary_step_median'] for v in values],'median':st.median(v['ordinary_step_median'] for v in values)}
 for side in ('candidate','legacy'):
  bare=result['timings'][side]['median'];cal={}
  for mode in ('minimal','timers','work'):
   probes=[v for k,v in runs.items() if k.startswith(side+'-'+mode+'-')]
   cal[mode]={'printed_step_medians':[v['ordinary_step_median'] for v in probes],'observer_wall_means':[v['operations']['ordinary_wall_mean'] for v in probes],'envelope_means':{name:[sum(s['contexts'].get(c,s['contexts'].get(str(c)))['mean'] for s in v['operations']['steps'] if not s['restart'] for c in ids)/sum(not s['restart'] for s in v['operations']['steps']) for v in probes] for name,ids in PHASES.items()}}
   cal[mode]['printed_change_from_bare_percent']=100*(st.median(cal[mode]['printed_step_medians'])/bare-1)
  for name,reference in (('timers','minimal'),('work','timers')):
   changes={key:st.mean(cal[name]['envelope_means'][key])-st.mean(cal[reference]['envelope_means'][key]) for key in PHASES}
   cal[name]['parent_change_seconds']=changes
   cal[name]['review_flags']=[key for key,delta in changes.items() if abs(delta)>max(0.001,0.02*st.mean(cal[reference]['envelope_means'][key]))]
  cal['timers']['bare_step_overhead_flag']=abs(cal['timers']['printed_change_from_bare_percent'])>2.0
  result['calibration'][side]=cal
 def rows(side):return [v['operations'] for k,v in runs.items() if k.startswith(side+'-timers-')]
 left,right=rows('legacy'),rows('candidate')
 for i in range(9*(NR+1)):
  a=[v['ledger'][i] for v in left];b=[v['ledger'][i] for v in right];c,p=a[0]['context'],a[0]['operation']
  result['ranking'].append({'context':c,'operation':p,'legacy_seconds':st.mean(v['mean'] for v in a),'block_seconds':st.mean(v['mean'] for v in b),'excess_seconds':st.mean(v['mean'] for v in b)-st.mean(v['mean'] for v in a),'paired_excess_seconds':[y['mean']-x['mean'] for x,y in zip(a,b)],'legacy_calls':st.mean(v['calls_per_rank_step'] for v in a),'block_calls':st.mean(v['calls_per_rank_step'] for v in b),'actual_slow_rank_path_excess':st.mean(v['slow_rank_path_mean'] for v in b)-st.mean(v['slow_rank_path_mean'] for v in a),'status':'unclassified' if p==0 else 'location only; semantic matching and causality require source audit'})
 result['ranking'].sort(key=lambda x:x['excess_seconds'],reverse=True)
 expected=st.mean(v['ordinary_wall_mean'] for v in right)-st.mean(v['ordinary_wall_mean'] for v in left)
 if abs(sum(v['excess_seconds'] for v in result['ranking'])-expected)>1e-6:raise ValueError('Ranked ledger does not conserve full-step excess')
 result['instrumented_mean_excess']=expected
 result['accounting']={'bare_median_excess':result['timings']['candidate']['median']-result['timings']['legacy']['median'],
  'instrumented_rank_mean_excess':expected,'signed_ledger_sum':sum(v['excess_seconds'] for v in result['ranking']),
  'positive_excess_sum':sum(max(0.,v['excess_seconds']) for v in result['ranking']),
  'negative_offset_sum':sum(min(0.,v['excess_seconds']) for v in result['ranking']),
  'actual_slow_rank_path_excess':sum(v['actual_slow_rank_path_excess'] for v in result['ranking']),
  'scope':'Bare median, additive rank mean and actual slow-rank path are distinct statistics; none is rescaled to force agreement.'}
 result['parent_ranking']=sorted([{'task':name,'legacy_seconds':sum(v['legacy_seconds'] for v in result['ranking'] if v['context'] in ids),
  'block_seconds':sum(v['block_seconds'] for v in result['ranking'] if v['context'] in ids),
  'excess_seconds':sum(v['excess_seconds'] for v in result['ranking'] if v['context'] in ids)} for name,ids in PHASES.items()],key=lambda v:v['excess_seconds'],reverse=True)

 result['storage_probes']={k:{'printed_step_median':v['ordinary_step_median'],'mode':v['operations']['mode'],'wall_mean':v['operations']['ordinary_wall_mean'],'storage':v['operations']['storage'],'work_per_rank_step':v['operations']['work_per_rank_step']} for k,v in runs.items() if v['kind']=='candidate_storage'}
 controls=[v for k,v in result['storage_probes'].items() if '-control-' in k]
 sampled=[v for k,v in result['storage_probes'].items() if '-control-' not in k]
 if controls and sampled:
  control=st.median(v['printed_step_median'] for v in controls);sample=st.median(v['printed_step_median'] for v in sampled)
  minimal=st.median(result['calibration']['candidate']['minimal']['printed_step_medians'])
  result['storage_calibration']={'control_median':control,'sampled_median':sample,
   'disabled_probe_binary_change_from_minimal_percent':100*(control/minimal-1),
   'sampled_change_from_own_control_percent':100*(sample/control-1),
   'scope':'Binary layout/disabled-hook effect and combined clocks/counting/sampling effect; not an overhead correction.'}
 result['calibration_scope']='2% or 1ms parent-change flags are review thresholds, not confidence bounds. Storage-control uses identical probe binary with sampling disabled. Sampling estimates are inclusive, may overlap, and are excluded from ranking.'
 result['optimization_ready']=False
 result['unresolved']='No automatically inferred recoverable-time estimate. Observer effects, fused scopes, dependency causes and hardware behavior require review.'
 return result
