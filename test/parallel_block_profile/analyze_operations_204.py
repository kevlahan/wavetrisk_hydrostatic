"""Disjoint operation accounting with explicit remainder, calibration and MPI owners."""
import json,math,statistics as st
from pathlib import Path
from analyze_operations_199 import common,CONTEXTS,PHASES
from instrument_operations_204 import REGIONS
NR=len(REGIONS)

def collect(root,ranks,mode):
    root=Path(root);data={};files=list(root.glob('operations204-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing Stage204 rank files')
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
                if p not in (46,47,48) or parent not in range(NR+1) or key in r['mpi']:raise ValueError('Invalid MPI owner')
                r['mpi'][key]=v
        if schema!=[1,mode,rank] or set(rows)!=set(ref) or set(rows)!=set(range(1,len(rows)+1)):raise ValueError('Schema/step mismatch')
        for n,r in rows.items():
            q=ref[n]
            if not math.isfinite(r['wall']) or r['wall']<0 or abs(r['wall']-q['wall'])>1e-6 or (r['restart'],r['remap'])!=(q['restart'],q['remap']):raise ValueError('Workflow mismatch')
            if len(r['self'])!=9*(NR+1):raise ValueError('Incomplete operation ledger')
            for c in range(9):
                if abs(sum(r['self'][c,p] for p in range(NR+1))-q['context'][c])>1e-6:raise ValueError('Context conservation failure')
                for p in (46,47,48):
                    if abs(sum(v for (cc,pp,_),v in r['mpi'].items() if (cc,pp)==(c,p))-r['self'][c,p])>1e-6:raise ValueError('MPI owner conservation failure')
            if mode==0 and (any(v for (c,p),v in r['self'].items() if p) or any(r['calls'].values()) or any(r['inclusive'].values()) or r['mpi']):raise ValueError('Disabled observer did work')
        data[rank]=rows
    if len({tuple(v) for v in data.values()})!=1:raise ValueError('Rank step mismatch')
    for n in data[0]:
        if len({(data[k][n]['restart'],data[k][n]['remap']) for k in data})!=1:raise ValueError('Rank flag mismatch')
    ordinary=[n for n,r in data[0].items() if not r['restart']]
    if not ordinary:raise ValueError('No ordinary steps')
    if mode and not sum(row['calls'].get((1,3),0) for rows in data.values() for row in rows.values()):
        raise ValueError('Missing scalar inverse kernel coverage')
    slow={n:max(data,key=lambda k:data[k][n]['wall']) for n in ordinary}
    ledger=[];mpi=[]
    for c in range(9):
        for p in range(NR+1):
            means=[st.mean(data[k][n]['self'][c,p] for n in ordinary) for k in range(ranks)]
            counts=[st.mean(data[k][n]['calls'][c,p] for n in ordinary) for k in range(ranks)]
            ledger.append({'context':c,'operation':p,'name':REGIONS.get(p,'unassigned'),
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
            'slow_rank_wall_mean':st.mean(data[slow[n]][n]['wall'] for n in ordinary),'ledger':ledger,'mpi_owners':mpi,
            'scope':'Exclusive wall time (inclusive_seconds is a nested union, never additive);  MPI time is also broken down by immediate enclosing scope, not extra time. Calls count wrapped sites (including false conditional calls). Rank counts are not weighted mathematical work.'}

def summarize(result):
    runs=list(result['runs'].values());out={'bare':{},'calibration':{},'comparisons':{},'inclusive_envelopes':{}}
    for side in ('baseline','candidate','legacy'):
        bare=[r['ordinary_step_median'] for r in runs if r['side']==side and r['phase']=='bare']
        if len(bare)!=2:raise ValueError('Two reversed bare controls required')
        out['bare'][side]={'samples':bare,'median_seconds':st.median(bare)}
        coarse=[r for r in runs if r['side']==side and r['phase']=='coarse']
        details=[r for r in runs if r['side']==side and r['phase']=='operations']
        if len(coarse)!=1 or len(details)!=2 or {r['round'] for r in details}!={1,2}:raise ValueError('Incomplete observer controls')
        ref=coarse[0]['hardware']['phases'];cal=[]
        for r in coarse+details:
            delta=100*(r['ordinary_step_median']/st.median(bare)-1);flags=[];changes={}
            if abs(delta)>3:flags.append('printed median differs from bare by >3%')
            for phase in ('whole step',*PHASES):
                value=r['hardware']['phases'][phase]['wall_seconds_rank_mean'];base=ref[phase]['wall_seconds_rank_mean']
                changes[phase]=value-base
                if abs(value-base)>max(.001,.02*base):flags.append(phase)
            if r['hardware']['quality_flags']:flags.append('hardware clock quality')
            cal.append({'label':r['label'],'printed_change_percent':delta,'phase_changes_seconds':changes,'flags':flags})
        out['calibration'][side]=cal
    # Complete routine envelopes answer whether a parent still costs time after
    # child transfers/MPI have been separated. These rows MUST NOT be added.
    for side in ('baseline','candidate','legacy'):
        detail=[r['operations'] for r in runs if r['side']==side and r['phase']=='operations']
        out['inclusive_envelopes'][side]=[{k:row[k] for k in ('context','operation','name')} |
            {'inclusive_seconds':st.mean(q['ledger'][i]['inclusive_seconds'] for q in detail)}
            for i,row in enumerate(detail[0]['ledger']) if row['operation']]
    for left,right,name in [('legacy','baseline','Stage202 minus legacy'),('legacy','candidate','Stage203 minus legacy'),('candidate','baseline','Stage202 minus Stage203')]:
        pairs=[]
        for rep in (1,2):
            pair={s:next(r for r in runs if r['side']==s and r['phase']=='operations' and r['round']==rep)['operations'] for s in (left,right)}
            pairs.append(pair)
        ranking=[]
        for i in range(9*(NR+1)):
            a=[v[left]['ledger'][i] for v in pairs];b=[v[right]['ledger'][i] for v in pairs]
            c,p=a[0]['context'],a[0]['operation'];diff=[y['seconds']-x['seconds'] for x,y in zip(a,b)]
            ranking.append({'context':c,'context_name':CONTEXTS[c],'operation':p,'name':a[0]['name'],
                'left_seconds':st.mean(v['seconds'] for v in a),'right_seconds':st.mean(v['seconds'] for v in b),
                'excess_seconds':st.mean(diff),'paired_excess_seconds':diff,'positive_in_both_pairs':all(x>0 for x in diff),
                'actual_slow_rank_path_excess':st.mean(y['actual_slow_rank_path_seconds']-x['actual_slow_rank_path_seconds'] for x,y in zip(a,b)),
                'interpretation':'unassigned' if p==0 else 'exclusive location difference; audit work/consumers before claiming removable cost'})
        total=st.mean(v[right]['wall_rank_mean']-v[left]['wall_rank_mean'] for v in pairs)
        if abs(sum(r['excess_seconds'] for r in ranking)-total)>1e-6:raise ValueError('Signed ranking does not close')
        ranking.sort(key=lambda r:r['excess_seconds'],reverse=True)
        unassigned={s:st.mean(sum(v[s]['ledger'][i]['seconds'] for i in range(9*(NR+1)) if v[s]['ledger'][i]['operation']==0) for v in pairs) for s in (left,right)}
        mpi_keys={(v['context'],v['operation'],v['owner']) for pair in pairs for side in (left,right) for v in pair[side]['mpi_owners']}
        mpi_rows=[]
        for key in sorted(mpi_keys):
            values={side:[sum(v['seconds'] for v in pair[side]['mpi_owners'] if (v['context'],v['operation'],v['owner'])==key) for pair in pairs] for side in (left,right)}
            mpi_rows.append({'context':key[0],'operation':key[1],'owner':key[2],
                'owner_name':REGIONS.get(key[2],'unassigned'),'left_seconds':st.mean(values[left]),'right_seconds':st.mean(values[right]),
                'excess_seconds':st.mean(values[right])-st.mean(values[left])})
        out['comparisons'][name]={'mpi_owner_breakdown':mpi_rows,'left':left,'right':right,'instrumented_excess_seconds':total,'unassigned_seconds':unassigned,'ranking':ranking}
    clean=not any(row['flags'] for rows in out['calibration'].values() for row in rows)
    out['decision']='review_operation_evidence' if clean else 'observer_effects_require_review'
    out['scope']='Diagnostic only; no promotion. Signed exclusive rows conserve instrumented excess. Calibration failure blocks causal ranking. Residual/driver time remains explicit; MPI wait is not independently removable. Repeated scopes show stability, not statistical confidence.'
    return out

def write_report(root,result):
    s=result['summary'];lines=['# Stage204 operation evidence','',s['decision'],'',s['scope'],'']
    for label,q in s['comparisons'].items():
        lines.extend([f'## {label}','',f"Instrumented rank-mean excess: {q['instrumented_excess_seconds']:.6f} s/ordinary step.",f"Unassigned absolute times: {q['unassigned_seconds']}",'',
                      '| Context | Operation | Difference (s) | Reversed pairs (s) |','|---|---|---:|---|'])
        for r in q['ranking'][:20]:lines.append(f"| {r['context_name']} | {r['name']} | {r['excess_seconds']:.6f} | {r['paired_excess_seconds']} |")
        lines.append('')
    (root/'REVIEW.md').write_text('\n'.join(lines)+'\n')
