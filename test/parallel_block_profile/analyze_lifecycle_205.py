"""Disjoint operation accounting with explicit remainder, calibration and MPI owners."""
import json,math,statistics as st
from pathlib import Path
from analyze_operations_199 import common,CONTEXTS,PHASES
from instrument_lifecycle_205 import REGIONS,EVENTS
NR=len(REGIONS)

def timing(root,ranks,mode):
    root=Path(root);data={};files=list(root.glob('operations205-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing Stage205 rank files')
    for file in files:
        rank=int(file.stem.rsplit('-',1)[1]);rows={};schema=None
        if rank not in range(ranks) or rank in data:raise ValueError('Invalid rank')
        ref=common(root/f'workflow195-rank-{rank}.txt')
        for line in file.read_text().splitlines():
            f=line.split()
            if not f:raise ValueError('Empty record')
            if f[0]=='schema':
                if schema is not None:raise ValueError('Duplicate schema')
                schema=list(map(int,f[1:]));continue
            if f[0]=='step':
                if len(f)!=5:raise ValueError('Invalid step')
                n=int(f[1])
                if n in rows:raise ValueError('Duplicate step')
                rows[n]={'restart':int(f[2]),'remap':int(f[3]),'wall':float(f[4]),'self':{},'calls':{},'mpi':{},'inclusive':{}}
                continue
            if f[0] not in ('self','mpi') or len(f)!=(7 if f[0]=='self' else 6):raise ValueError('Unknown/invalid operation row')
            n,c,p=map(int,f[1:4]);v=float(f[5])
            if n not in rows or c not in range(9) or p not in range(NR+1) or not math.isfinite(v) or v<0:raise ValueError('Invalid operation value')
            r=rows[n]
            if f[0]=='self':
                key=c,p;count=int(f[4])
                if key in r['self'] or count<0:raise ValueError('Duplicate/invalid self row')
                inc=float(f[6])
                if not math.isfinite(inc) or inc<0 or (p and inc+1e-6<v):raise ValueError('Invalid inclusive scope')
                r['self'][key]=v;r['calls'][key]=count;r['inclusive'][key]=inc
            else:
                parent=int(f[4]);key=c,p,parent
                if p not in (26,27,28,29) or parent not in range(NR+1) or key in r['mpi']:raise ValueError('Invalid MPI owner')
                r['mpi'][key]=v
        if schema!=[1,mode,rank] or set(rows)!=set(ref) or set(rows)!=set(range(1,len(rows)+1)):raise ValueError('Schema/step mismatch')
        for n,r in rows.items():
            q=ref[n]
            if not math.isfinite(r['wall']) or r['wall']<0 or abs(r['wall']-q['wall'])>1e-6 or (r['restart'],r['remap'])!=(q['restart'],q['remap']):raise ValueError('Workflow mismatch')
            if len(r['self'])!=9*(NR+1):raise ValueError('Incomplete operation ledger')
            for c in range(9):
                if abs(sum(r['self'][c,p] for p in range(NR+1))-q['context'][c])>1e-6:raise ValueError('Context conservation failure')
                for p in (26,27,28,29):
                    if abs(sum(v for (cc,pp,_),v in r['mpi'].items() if (cc,pp)==(c,p))-r['self'][c,p])>1e-6:raise ValueError('MPI owner conservation failure')
            if mode==0 and (any(v for (c,p),v in r['self'].items() if p) or any(r['calls'].values()) or any(r['inclusive'].values()) or r['mpi']):raise ValueError('Disabled observer did work')
        data[rank]=rows
    if len({tuple(v) for v in data.values()})!=1:raise ValueError('Rank step mismatch')
    for n in data[0]:
        if len({(data[k][n]['restart'],data[k][n]['remap']) for k in data})!=1:raise ValueError('Rank flag mismatch')
    ordinary=[n for n,r in data[0].items() if not r['restart']]
    if not ordinary:raise ValueError('No ordinary steps')
    slow={n:min(k for k in data if data[k][n]['wall']==max(r[n]['wall'] for r in data.values())) for n in ordinary}
    ledger=[];mpi=[]
    for c in range(9):
        for p in range(NR+1):
            means=[st.mean(data[k][n]['self'][c,p] for n in ordinary) for k in range(ranks)]
            counts=[st.mean(data[k][n]['calls'][c,p] for n in ordinary) for k in range(ranks)]
            ledger.append({'context':c,'operation':p,'name':REGIONS.get(p,'unselected'),
                'seconds':st.mean(means),'inclusive_seconds':st.mean(data[k][n]['inclusive'][c,p] for k in range(ranks) for n in ordinary),
                'rank_seconds':means,'rank_calls':counts,'calls':st.mean(counts),
                'actual_slow_rank_path_seconds':st.mean(data[slow[n]][n]['self'][c,p] for n in ordinary)})
        keys={(p,parent) for rows in data.values() for n in ordinary for cc,p,parent in rows[n]['mpi'] if cc==c}
        for p,parent in sorted(keys):
            vals=[st.mean(data[k][n]['mpi'].get((c,p,parent),0) for n in ordinary) for k in range(ranks)]
            mpi.append({'context':c,'operation':p,'owner':parent,'owner_name':REGIONS.get(parent,'unassigned'),'seconds':st.mean(vals),'rank_seconds':vals})
    wall=st.mean(data[k][n]['wall'] for k in range(ranks) for n in ordinary)
    if abs(sum(v['seconds'] for v in ledger)-wall)>1e-6:raise ValueError('Full-step conservation failure')
    return {'mode':mode,'ranks':ranks,'ordinary_steps':ordinary,'wall_rank_mean':wall,
            'slow_rank_wall_mean':st.mean(data[slow[n]][n]['wall'] for n in ordinary),
            'slow_rank_ties':{str(n):[k for k in sorted(data) if data[k][n]['wall']==data[slow[n]][n]['wall']] for n in ordinary},
            'slow_rank_rule':'Lowest tied rank; a descriptive path, not a proven critical path. Do not rank optimizations by this secondary statistic.','ledger':ledger,'mpi_owners':mpi,
            'scope':'Exclusive wall time (inclusive_seconds is a nested union, never additive);  MPI time is also broken down by immediate enclosing scope, not extra time. Calls execute inside the original conditional. Selected MPI excludes posting and unselected scopes. Rank counts are not weighted mathematical work.'}

def lifecycle(root,ranks,mode,ordinary):
    root=Path(root);files=list(root.glob('lifecycle205-rank-*.txt'));data={};sites={};plans={};events={};transitions={};duplicates=[]
    if len(files)!=ranks:raise ValueError('Missing lifecycle ranks')
    for f in files:
        rank=int(f.stem.rsplit('-',1)[1]);lines=f.read_text().splitlines();head=lines.pop(0).split()
        if rank not in range(ranks) or rank in data or head[:4]!=['schema','1',str(mode),str(rank)] or len(head)!=5 or int(head[4])!=len(lines):raise ValueError('Invalid lifecycle header')
        data[rank]=[];seen_geometry=set();last_step=0
        known_steps=common(root/f'workflow195-rank-{rank}.txt')
        for i,line in enumerate(lines,1):
            v=line.split()
            if len(v)!=14 or v[0]!='event' or int(v[1])!=i:raise ValueError('Lifecycle sequence error')
            step,context,owner,rk,kind,obj,gen,a,b,c=map(int,v[2:12])
            start,duration=map(float,v[12:]);
            if not math.isfinite(start) or not math.isfinite(duration) or duration<0 or (kind==40 and start<0):raise ValueError('Invalid collective clock')
            if step<last_step or step not in known_steps or context not in range(9) or owner not in range(NR+1) or rk not in range(5) or kind not in EVENTS or gen<0:raise ValueError('Invalid lifecycle event')
            last_step=step
            if not mode:raise ValueError('Disabled observer recorded lifecycle events')
            data[rank].append([step,context,owner,rk,kind,obj,gen,a,b,c,start,duration])
            if step not in ordinary:continue
            events[str(kind)]=events.get(str(kind),0)+1
            if kind==1:
                if obj not in (1,2,3,4) or b not in (0,1):raise ValueError('Invalid plan reuse event')
                key=str(obj);q=plans.setdefault(key,{'requests':0,'reused':0,'rebuilt':0,'generations':set()})
                q['requests']+=1;q['reused' if b else 'rebuilt']+=1;q['generations'].add((rank,gen))
                if obj<4 and b!=int(a==gen):raise ValueError('Reuse predicate contradicts generation')
            if kind==5:
                key=rank,gen
                if key in seen_geometry:duplicates.append({'rank':rank,'step':step,'generation':gen,'kind':'geometry reseed within generation; audit inputs before declaring redundancy'})
                seen_geometry.add(key)
            if kind in (6,7):
                key=f'{kind}:{a}:{b}';q=transitions.setdefault(key,{'events':0,'stage_generation_visits':{}});q['events']+=1
                visit=f'{rank}:{step}:{rk}:{gen}';q['stage_generation_visits'][visit]=q['stage_generation_visits'].get(visit,0)+1
            if kind==40:
                if obj<=0 or min(a,b)<-1:raise ValueError('Invalid collective counts')
                q=sites.setdefault(str(obj),{'calls':0,'send_elements':0,'recv_elements':0,'unknown_count_calls':0,'owners':set(),'seconds':0.0,'rank_seconds':{}})
                q['calls']+=1;q['owners'].add(owner);q['seconds']+=duration;q['rank_seconds'][str(rank)]=q['rank_seconds'].get(str(rank),0.0)+duration
                if min(a,b)<0:q['unknown_count_calls']+=1
                else:q['send_elements']+=a;q['recv_elements']+=b
            if kind==13 and min(a,b)<0:raise ValueError('Negative producer increments')
        # Legacy has selected timers but no block lifecycle events; an empty tape is valid.
    divisor=ranks*len(ordinary)
    for q in plans.values():q['rank_generation_count']=len(q.pop('generations'))
    for q in sites.values():q['owners']=sorted(q['owners']);q['calls_per_rank_step']=q['calls']/divisor;q['seconds_per_rank_step']=q['seconds']/divisor
    for q in transitions.values():
        vals=list(q.pop('stage_generation_visits').values());q.update(max_visits_in_same_rank_step_stage_generation=max(vals),repeat_visits=sum(max(0,n-1) for n in vals))
    return {'schema':1,'ranks':ranks,'mode':mode,'event_counts':events,'plans':plans,'collective_sites':sites,'field_transitions':transitions,
      'geometry_reseeds':duplicates,'rank_events':data,'scope':'Counts and lifecycle predicates, not proof that live fields are unchanged. Site volumes use the MPI datatype at that site, not inferred DRAM bytes.'}

def collect(root,ranks,mode):
    q=timing(root,ranks,mode);q['lifecycle']=lifecycle(root,ranks,mode,q['ordinary_steps'])
    expected=sum(v['seconds'] for v in q['ledger'] if v['operation'] in (26,27))
    actual=sum(v['seconds_per_rank_step'] for v in q['lifecycle']['collective_sites'].values())
    if abs(expected-actual)>1e-6:raise ValueError('Collective site clock conservation failed')
    return q

def summarize(result):
    runs=list(result['runs'].values());out={'bare':{},'calibration':{},'signed_locations':[],'lifecycle':{},'scope':'Selected batch accounting only. Unselected includes required computation and unobserved MPI. Inclusive parents overlap; do not add. Timing/trace evidence does not prove removability.'}
    for side in ('baseline','candidate','legacy'):
        bare=[q['ordinary_step_median'] for q in runs if q['side']==side and q['phase']=='bare']
        if len(bare)!=2:raise ValueError('Two bare controls required')
        out['bare'][side]={'samples':bare,'median_seconds':st.median(bare)}
    for side in ('baseline','legacy'):
        records=[];out['lifecycle'][side]=[]
        for rep in (1,2):
            selected={p:[q for q in runs if q['side']==side and q['phase']==p and q['round']==rep] for p in ('coarse','operations')}
            if any(len(v)!=1 for v in selected.values()):raise ValueError('Missing/duplicate adjacent calibration pair')
            control,detail=(selected[p][0] for p in ('coarse','operations'))
            for q in (control,detail):
                flags=[];delta=100*(q['ordinary_step_median']/out['bare'][side]['median_seconds']-1);changes={}
                if abs(delta)>3:flags.append('printed median differs from bare by >3%')
                for phase in ('whole step',*PHASES):
                    base=control['hardware']['phases'][phase]['wall_seconds_rank_mean'];value=q['hardware']['phases'][phase]['wall_seconds_rank_mean']
                    changes[phase]=value-base
                    if abs(value-base)>max(.001,.02*base):flags.append(phase)
                if q['hardware']['quality_flags']:flags.append('hardware clock quality')
                records.append({'label':q['label'],'control':control['label'],'printed_change_percent':delta,'phase_changes_seconds':changes,'flags':flags})
            life=detail['operations']['lifecycle'];out['lifecycle'][side].append({k:v for k,v in life.items() if k!='rank_events'})
        out['calibration'][side]=records
    for c in range(9):
        for p in range(NR+1):
            pairs=[];means={};incs={}
            for side in ('baseline','legacy'):
                rows=[next(x for x in q['operations']['ledger'] if (x['context'],x['operation'])==(c,p)) for q in runs if q['side']==side and q['phase']=='operations']
                means[side]=st.mean(x['seconds'] for x in rows);incs[side]=st.mean(x['inclusive_seconds'] for x in rows)
            for rep in (1,2):
                vals={q['side']:next(x['seconds'] for x in q['operations']['ledger'] if (x['context'],x['operation'])==(c,p)) for q in runs if q['side'] in ('baseline','legacy') and q['phase']=='operations' and q['round']==rep}
                pairs.append(vals['baseline']-vals['legacy'])
            out['signed_locations'].append({'context':c,'operation':p,'name':REGIONS.get(p,'unselected'),'exclusive':means,'inclusive_nonadditive':incs,'difference':means['baseline']-means['legacy'],'paired_difference':pairs})
    out['signed_locations'].sort(key=lambda r:r['difference'],reverse=True)
    total=st.mean(q['operations']['wall_rank_mean'] for q in runs if q['side']=='baseline' and q['phase']=='operations')-st.mean(q['operations']['wall_rank_mean'] for q in runs if q['side']=='legacy' and q['phase']=='operations')
    if abs(total-sum(x['difference'] for x in out['signed_locations']))>1e-6:raise ValueError('Signed ledger does not close')
    out['instrumented_rank_mean_excess']=total
    flags=[f for rows in out['calibration'].values() for q in rows for f in q['flags']]
    out['decision']='observer_effects_require_review' if flags else 'ready_for_consumer_lifetime_audit'
    out['optimization_authorized']=False
    out['next_gate']='Calibrated, non-overlapping avoidable costs plus explicit consumer/dependency proof must support >=0.10 s/step before selecting a prototype. Bare runs must then validate actual savings.'
    return out

def write_report(root,result):
    s=result['summary'];lines=['# Stage205 lifecycle confirmation','',s['decision'],'',s['scope'],'',s['next_gate'],'',
      '| Context | Selected location | Baseline − legacy (s) | Reversed differences (s) |','|---|---|---:|---|']
    for q in s['signed_locations'][:25]:lines.append(f"| {CONTEXTS[q['context']]} | {q['name']} | {q['difference']:.6f} | {q['paired_difference']} |")
    for side, repeats in s['lifecycle'].items():
        lines+=['',f'## {side}: lifecycle observations']
        for i,life in enumerate(repeats,1):
            lines+=['',f'Repeat {i}; totals across ordinary steps and ranks.']
            for key,q in life.get('plans',{}).items():
                lines.append(f"- Plan {key}: {q['requests']} requests, {q['reused']} reuses, {q['rebuilt']} rebuild predicates.")
            lines.append(f"- Repeated geometry seeds within a generation: {len(life.get('geometry_reseeds',[]))}; inspect storage validity before interpreting.")
            lines+=['','| Collective site | Source | Seconds/rank/step | Calls/rank/step |','|---|---|---:|---:|']
            for key,q in sorted(life.get('collective_sites',{}).items(),key=lambda kv:kv[1]['seconds_per_rank_step'],reverse=True)[:12]:
                src=q.get('source',{});label=str(src.get('file',''))+':'+str(src.get('line',''))
                lines.append(f"| {key} | {label} | {q['seconds_per_rank_step']:.6f} | {q['calls_per_rank_step']:.2f} |")
    lines+=['','Read lifecycle-summary.json and raw lifecycle205-rank-*.txt with consumer-map-205.json. Repeated generation or stage alone is not proof of unchanged live fields. Collective-site durations reconcile with selected MPI time; counts identify datatype elements, not independently removable latency. Start times are rank-local step offsets, not synchronized global arrival times.']
    (root/'REVIEW.md').write_text('\n'.join(lines)+'\n')


def attach_sites(q,inventory,side):
    sites={str(s['collective_site']):{'file':name,**s} for name,item in inventory['implementations'][side].items() for s in item['sites'] if s['collective_site'] is not None}
    for key,row in q['lifecycle']['collective_sites'].items():
        if key not in sites:raise ValueError('Unknown collective call site '+key)
        row['source']=sites[key]
    if side=='baseline' and q['mode']==1:
        life=q['lifecycle'];required={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,40,41}
        if not required.issubset(map(int,life['event_counts'])):raise ValueError('Missing lifecycle producer/consumer coverage')
        if set(life['plans'])!={'1','2','3','4'}:raise ValueError('Missing plan lifetime coverage')
        for rank,rows in life['rank_events'].items():
            for step in q['ordinary_steps']:
                stages=[v[7] for v in rows if v[0]==step and v[4]==15]
                if stages!=[1,2,3,4]:raise ValueError('RK producer stage coverage differs')
        if sum(v[7] for rows in life['rank_events'].values() for v in rows if v[4]==13)<=0:raise ValueError('No scalar physical production observed')
