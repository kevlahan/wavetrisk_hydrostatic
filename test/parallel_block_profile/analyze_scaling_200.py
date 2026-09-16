"""Matched scaling evidence; measured locations are never labelled recoverable savings."""
import csv,json,math,statistics as st
from pathlib import Path
from analyze_operations_199 import common,PHASES
from analyze_workflow_195 import collect as workflow
MPI={13,18,19,34,37,38,39,40,60,72,75}

def corr(x,y):
    if len(x)<2 or len(set(x))<2 or len(set(y))<2:return None
    return st.correlation(x,y)

def collect(root,ranks,mode):
    root=Path(root);wf=workflow(root,ranks,1);data={}
    files=list(root.glob('scaling200-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing scaling rank files')
    for f in files:
        rank=int(f.stem.rsplit('-',1)[1]);ref=common(f.with_name(f'workflow195-rank-{rank}.txt'));rows={}
        if rank not in range(ranks) or rank in data:raise ValueError('Duplicate/invalid rank')
        lines=f.read_text().splitlines()
        if lines[0].split()!=['schema','1',str(mode),str(rank)]:raise ValueError('Wrong scaling schema/mode')
        for line in lines[1:]:
            v=line.split()
            if not v:continue
            if len(v)<4:raise ValueError('Truncated scaling record')
            kind=v[0];n=int(v[1])
            if kind=='step':
                if n in rows:raise ValueError('Duplicate step')
                rows[n]={'restart':int(v[2]),'remap':int(v[3]),'wall':float(v[4]),'self':{},'calls':{},'count':{},'volume':{}}
                continue
            if kind not in ('self','count','volume') or n not in rows:raise ValueError('Invalid scaling record')
            c,p=map(int,v[2:4]);r=rows[n];key=c,p
            if c not in range(9) or key in r[kind]:raise ValueError('Invalid/duplicate context record')
            if kind=='self':
                if p not in MPI|{0} or len(v)!=6:raise ValueError('Unexpected detailed timer')
                r['calls'][key]=int(v[4]);r[kind][key]=float(v[5])
                if r['calls'][key]<0 or not math.isfinite(r[kind][key]) or r[kind][key]<0:raise ValueError('Invalid timer')
            elif kind=='count':
                if p not in range(1,43) or len(v)!=5 or int(v[4])<0:raise ValueError('Invalid counter')
                r[kind][key]=int(v[4])
            else:
                if p not in range(1,57) or len(v)!=6 or min(map(int,v[4:]))<0:raise ValueError('Invalid volume')
                r[kind][key]=list(map(int,v[4:]))
        if set(rows)!=set(ref):raise ValueError('Step coverage differs')
        for n,r in rows.items():
            q=ref[n]
            if not math.isfinite(r['wall']) or abs(r['wall']-q['wall'])>1e-6 or (r['restart'],r['remap'])!=(q['restart'],q['remap']):raise ValueError('Workflow mismatch')
            if len(r['self'])!=108 or len(r['count'])!=378:raise ValueError('Incomplete dense coverage')
            for c in range(9):
                if abs(sum(r['self'][c,p] for p in MPI|{0})-q['context'][c])>1e-6:raise ValueError('Context conservation failure')
            if mode==0 and (any(r['calls'].values()) or any(r['count'].values()) or r['volume'] or any(v for (c,p),v in r['self'].items() if p)):raise ValueError('Coarse control recorded MPI/work')
        data[rank]=rows
    if any(set(v)!=set(data[0]) for v in data.values()):raise ValueError('Rank step mismatch')
    ordinary=[n for n,r in data[0].items() if not r['restart']];cells=[]
    if not ordinary:raise ValueError('No ordinary steps')
    if mode==1 and not sum(sum(r['calls'].values()) for rows in data.values() for r in rows.values()):raise ValueError('Missing MPI coverage')
    for n in ordinary:
        slow=max(data,key=lambda rank:data[rank][n]['wall'])
        for name,contexts in {'whole step':list(range(9)),**PHASES}.items():
            rankrows=[]
            for rank in range(ranks):
                r=data[rank][n];local=sum(r['self'][c,0] for c in contexts);mpi=sum(r['self'][c,p] for c in contexts for p in MPI)
                counts=[sum(r['count'][c,p] for c in contexts) for p in range(1,43)]
                rankrows.append({'rank':rank,'outside_mpi_seconds':local,'mpi_seconds':mpi,'wall_seconds':local+mpi,'counts':counts,
                  'mpi_calls':sum(r['calls'][c,p] for c in contexts for p in MPI),
                  'volume':{str(p):[sum(r['volume'].get((c,p),[0,0])[i] for c in contexts) for i in (0,1)] for p in range(1,57) if any((c,p) in r['volume'] for c in contexts)}})
            local=[r['outside_mpi_seconds'] for r in rankrows];mpi=[r['mpi_seconds'] for r in rankrows];work=[r['counts'][34] for r in rankrows]
            cells.append({'step':n,'phase':name,'rank_rows':rankrows,'outside_mpi_mean':st.mean(local),'outside_mpi_max':max(local),
              'outside_mpi_max_over_mean':max(local)/st.mean(local) if st.mean(local) else None,
              'outside_mpi_median':st.median(local),'outside_mpi_min':min(local),
              'mpi_mean':st.mean(mpi),'wall_mean':st.mean(r['wall_seconds'] for r in rankrows),'wall_max':max(r['wall_seconds'] for r in rankrows),
              'phase_on_full_step_slowest_rank':rankrows[slow]['wall_seconds'],'full_step_slowest_rank':slow,
              'work35_global':sum(work),'outside_mpi_seconds_per_work35':sum(local)/sum(work) if sum(work) else None,
              'local_mpi_correlation':corr(local,mpi),'local_work35_correlation':corr(local,work)})
    return {'mode':mode,'ranks':ranks,'ordinary_steps':ordinary,'workflow':wf,'cells':cells,
      'scope':'MPI is a selected scope subset, including progress and producer wait. Outside MPI is elapsed time, not hardware CPU time. Phase maxima are not additive. Counter 35 means temperature direct-flux node evaluations; do not assume all work is represented by it.'}

def summarize(result):
    runs=result['runs'];cases={}
    for case in result['cases']:
        row={'bare':{},'observer':{}}
        for side in ('candidate','legacy'):
            selected=[r for r in runs.values() if r['case']==case and r['side']==side]
            vals=[r['ordinary_step_median'] for r in selected if r['phase']=='bare']
            if len(vals)!=2:raise ValueError('Two bare repeats required')
            row['bare'][side]={'samples':vals,'median':st.median(vals),'range':max(vals)-min(vals)}
            row['observer'][side]={}
            for phase in ('coarse','probe'):
                r=next(r for r in selected if r['phase']==phase);a=r['scaling'];parents={}
                for name in PHASES:
                    cells=[c for c in a['cells'] if c['phase']==name]
                    parents[name]={k:st.mean(c[k] for c in cells) for k in ('outside_mpi_mean','outside_mpi_max','mpi_mean','wall_mean','phase_on_full_step_slowest_rank','work35_global')}
                    if phase=='probe':
                        parents[name]['work_counters_global_per_step']=[st.mean(sum(q['counts'][i] for q in c['rank_rows']) for c in cells) for i in range(42)]
                        parents[name]['work_counters_rank_max_per_step']=[st.mean(max(q['counts'][i] for q in c['rank_rows']) for c in cells) for i in range(42)]
                        ratios=[c['outside_mpi_max_over_mean'] for c in cells if c['outside_mpi_max_over_mean'] is not None]
                        parents[name]['outside_mpi_rank_max_over_mean']=st.mean(ratios) if ratios else None
                        parents[name]['mpi_calls_per_rank_step']=st.mean(st.mean(q['mpi_calls'] for q in c['rank_rows']) for c in cells)
                    if phase=='coarse':
                        for key in ('outside_mpi_mean','outside_mpi_max','mpi_mean','work35_global'):parents[name][key]=None
                row['observer'][side][phase]={'printed_median':r['ordinary_step_median'],'parents':parents,
                   'printed_change_from_bare_percent':100*(r['ordinary_step_median']/st.median(vals)-1)}
            coarse=row['observer'][side]['coarse'];probe=row['observer'][side]['probe']
            probe['parent_change_from_coarse_seconds']={n:probe['parents'][n]['wall_mean']-coarse['parents'][n]['wall_mean'] for n in PHASES}
            probe['review_flags']=[n for n,v in probe['parent_change_from_coarse_seconds'].items() if abs(v)>max(.001,.02*coarse['parents'][n]['wall_mean'])]
        row['gap_seconds']=row['bare']['candidate']['median']-row['bare']['legacy']['median']
        row['ratio']=row['bare']['candidate']['median']/row['bare']['legacy']['median']
        row['probe_parent_excess']={n:row['observer']['candidate']['probe']['parents'][n]['wall_mean']-row['observer']['legacy']['probe']['parents'][n]['wall_mean'] for n in PHASES}
        cases[case]=row
    for name,row in cases.items():
        n=result['cases'][name]['ranks']
        row['speedup_vs_21_one_socket']={side:cases['n21-s1']['bare'][side]['median']/row['bare'][side]['median'] for side in ('candidate','legacy')}
        row['efficiency_vs_21_one_socket']={side:v/(n/21) for side,v in row['speedup_vs_21_one_socket'].items()}
    topology={side:cases['n83-s2']['bare'][side]['median']-cases['n83-s1']['bare'][side]['median'] for side in ('candidate','legacy')}
    return {'cases':cases,'two_socket_minus_one_socket_seconds':topology,
      'optimization_ready':False,'material_target_seconds':0.10,
      'decision_rule':'Seek a reproducible mechanism plausibly recovering at least 0.10 s of the 83-rank gap. This threshold screens proposals; measured timing locations or scaling differences are not removable-time estimates.',
      'limitations':['Two bare repeats detect gross consistency, not precise confidence intervals. One coarse/probe pair per side/case; repeat only ambiguous cases.',
       'Scaling changes cache footprint, NUMA traffic, ownership granularity and rank work; it cannot uniquely separate these causes.',
       'MPI delay is not added to local-work savings. Do not extrapolate a one-percent micro-optimization from this campaign.']}

def write_reports(out,result):
    s=result['summary'];lines=['# Stage 200 paired scaling evidence','','| Placement | Blocks s | Legacy s | Gap s | Ratio |','|---|---:|---:|---:|---:|']
    for case,v in s['cases'].items():lines.append(f"| {case} | {v['bare']['candidate']['median']:.4f} | {v['bare']['legacy']['median']:.4f} | {v['gap_seconds']:.4f} | {v['ratio']:.3f} |")
    lines+=['','Bare runs are authoritative. Phase/MPI/work measurements are separately calibrated, not a recovery estimate.','',s['decision_rule'],'','Inspect scaling.csv, observer calibration in results.json and raw per-rank records before selecting a code change.']
    (out/'scaling.md').write_text('\n'.join(lines)+'\n')
    with (out/'scaling.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['case','side','mode','step','phase','rank','outside_mpi_s','mpi_s','wall_s','mpi_calls','temperature_flux_nodes'])
        for r in result['runs'].values():
            if 'scaling' not in r:continue
            for c in r['scaling']['cells']:
                for q in c['rank_rows']:w.writerow([r['case'],r['side'],r['phase'],c['step'],c['phase'],q['rank'],q['outside_mpi_seconds'] if r['phase']=='probe' else '',q['mpi_seconds'] if r['phase']=='probe' else '',q['wall_seconds'],q['mpi_calls'] if r['phase']=='probe' else '',q['counts'][34] if r['phase']=='probe' else ''])
