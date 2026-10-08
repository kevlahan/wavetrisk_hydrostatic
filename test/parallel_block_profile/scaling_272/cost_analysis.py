"""Validate complete bounded records; retain overhead failures and signed costs."""
import math,re,statistics
from pathlib import Path
from phases import NAMES, IDS, RETIRED
WORK={8:('plan rebuilds','unused','unused'),9:('published interior nodes','published interior values','publication calls'),10:('installed interior slots','installed boundary slots','actual installations'),12:('rebuilt columns','rebuilt physical layers','rebuilt blocks'),27:('posted remote requests','sent remote doubles','received remote doubles')}
WORK.update({IDS['OWNERS']:('summed Domain ownership','summed block ownership','stage samples'), IDS['PATCHES']:('summed Domain patch inventory','summed block patch inventory','stage samples'),IDS['BASIC']:('producer patch-layer visits','Domain-level-layer visits','unused'),IDS['SCALAR_LOCAL']:('restricted parent patches','restricted edges','replayed divergence values'), IDS['CAPTURE_MPI']:('capture exchanges','sent remote doubles','received remote doubles'),IDS['RESTRICT_MPI']:('posted remote requests','sent remote doubles','received remote doubles'),IDS['TEMP_EDGES']:('posted remote requests','sent remote doubles','received remote doubles'),IDS['MIG_PAYLOAD']:('sent logical bytes','received logical bytes','exchanges'),IDS['GHOST_PLAN']:('rebuilds','reuses','unused'),IDS['WRITEBACK_PLAN']:('rebuilds','reuses','unused'),IDS['EDGE_METRIC']:('rebuilds','reuses','unused')})
INCLUDED=(1,2,3,4,5,6,9,10,11)

def read_costs(directory,ranks,block,oracle=False,domains=160):
    directory=Path(directory)
    files=sorted(directory.glob('cost-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing/extra observer rank files')
    all_rows=[]
    for rank in range(ranks):
        p=directory/f'cost-rank-{rank:06d}.txt';lines=p.read_text().splitlines()
        if lines[0].split()!=['STAGE272_COST_V1',str(rank),'12',str(len(NAMES))]:raise ValueError('Observer header differs')
        if len(lines)!=1+12*len(NAMES):raise ValueError('Observer rectangle incomplete')
        rank_rows=[]
        for s in range(12):
            rows=[]
            for phase in range(len(NAMES)):
                fields=lines[1+s*len(NAMES)+phase].split()
                if len(fields)!=10 or list(map(int,fields[:2]))!=[s+1,phase+1]:raise ValueError('Observer row identity differs')
                calls=int(fields[2]);v=list(map(float,fields[3:7]));w=list(map(int,fields[7:]))
                if calls<0 or min(w)<0 or any(not math.isfinite(x) or x<0 for x in v):raise ValueError('Invalid observer values')
                if any(v[i]>v[i+2]+1e-7 for i in (0,1)):raise ValueError('Self exceeds inclusive time')
                if not calls and any(v):raise ValueError('Time without scope call')
                rows.append(dict(calls=calls,self_wall=v[0],self_cpu=v[1],inclusive_wall=v[2],inclusive_cpu=v[3],work=w))
            if rows[0]['calls']!=1 or rows[1]['calls']!=1 or rows[2]['calls']!=4 or rows[3]['calls']<4:raise ValueError('Incomplete RK4 transaction scopes')
            for metric in ('wall','cpu'):
                total=math.fsum(row['self_'+metric] for row in rows)
                elapsed=rows[0]['inclusive_'+metric]
                if elapsed<=0 or abs(total-elapsed)>max(1e-6,elapsed*1e-6):raise ValueError('Observer phase conservation failed')
            if block:
                expected=7 if oracle else 4
                if rows[9]['work'][2]!=expected or rows[10]['calls']!=3 or rows[25]['calls']!=1:raise ValueError('Installation/provisional/accepted transaction count differs')
            elif any(rows[p-1]['calls'] for p in (7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,26,27)+tuple(range(28,len(NAMES)+1))):
                raise ValueError('Native block scope in legacy main')
            if any(rows[p-1]['calls'] or any(rows[p-1]['work']) for p in RETIRED):raise ValueError('Retired fine-grained observer scope executed')
            if block:
                for symbol in ('COMPAT','CAPTURE_BEGIN','CAPTURE_PUBLISH','VELOCITY_REMAINDER','RK_ASSEMBLE','SCALAR_RESTRICT','EXNER_TRANSACTION','EXNER_KERNEL','ACCUM_RETAIN'):
                    if rows[IDS[symbol]-1]['calls']<4:raise ValueError('Missing decomposed stage transaction: '+symbol)
                for symbol in ('GHOST_PLAN','WRITEBACK_PLAN','EDGE_METRIC'):
                    row=rows[IDS[symbol]-1]
                    if row['work'][0]+row['work'][1]!=row['calls']:raise ValueError('Rebuild/reuse accounting differs: '+symbol)
            rank_rows.append(rows)
        all_rows.append(rank_rows)
    result=dict(rank_files=ranks,steps=12,ordinary_steps=[s+1 for s in INCLUDED],phases={},longest_local_step_rank_phases={})
    for p,name in enumerate(NAMES):
        rows=[all_rows[r][s][p] for r in range(ranks) for s in INCLUDED]
        result['phases'][name]={metric:statistics.mean(row[metric] for row in rows) for metric in ('self_wall','self_cpu','inclusive_wall','inclusive_cpu','calls')}
        result['phases'][name]['mean_work']=[statistics.mean(row['work'][i] for row in rows) for i in range(3)]
        if p+1 in WORK:result['phases'][name]['work_labels']=list(WORK[p+1])
        longest=[max(range(ranks),key=lambda r:all_rows[r][s][0]['inclusive_wall']) for s in INCLUDED]
        result['longest_local_step_rank_phases'][name]=statistics.mean(all_rows[rank][s][p]['self_wall'] for rank,s in zip(longest,INCLUDED))
    result['mean_rank_step_seconds']=result['phases']['step']['inclusive_wall']
    result['mean_longest_local_step_seconds']=math.fsum(result['longest_local_step_rank_phases'].values())
    if block:result['ownership_groups']=ownership_groups(all_rows,domains)
    result['limitations']='Self times add within a rank/step. Inclusive times overlap. Longest-local-step composition selects one rank per step, not independent region maxima; unsynchronized local durations do not prove the global critical path. CPU includes MPI spinning. Work counters cover named inverse/cache products, not all solver traffic; doubles are logical payload counts, not physical memory traffic.'
    result['mean_scope_entries_per_step']=math.fsum(v['calls'] for v in result['phases'].values())
    return result

def calibrate(batch):
    if len(batch)!=4 or [r['observer'] for r in batch]!=[False,True,True,False]:raise ValueError('Require bare-observer-observer-bare')
    a=[batch[i]['ordinary']['median'] for i in (0,3)];b=[batch[i]['ordinary']['median'] for i in (1,2)]
    control=statistics.mean(a);observed=statistics.mean(b)
    shifts=[b[0]-a[0],b[1]-a[1]]
    limits=dict(mean=max(.003,.01*control),pair=max(.005,.02*control),spread=.02*control)
    passed=abs(observed-control)<=limits['mean'] and max(map(abs,shifts))<=limits['pair'] and max(max(a)-min(a),max(b)-min(b))<=limits['spread']
    return dict(status='compatible_at_prespecified_resolution' if passed else 'observer_costs_unqualified',control_seconds=control,observer_seconds=observed,shift_seconds=observed-control,paired_shifts_seconds=shifts,limits_seconds=limits,control_spread_seconds=max(a)-min(a),observer_spread_seconds=max(b)-min(b),scope='Two independent runs per binary; compatibility does not prove each small scope unbiased.')

def compare_costs(rows,calibration):
    means={}
    for family in ('main','stage266'):
        selected=[r['costs'] for r in rows if r['family']==family and r['observer']]
        if len(selected)!=2:raise ValueError('Require two observers per control family')
        means[family]={name:{metric:statistics.mean(row['phases'][name][metric] for row in selected) for metric in ('self_wall','self_cpu','inclusive_wall','inclusive_cpu','calls')} for name in NAMES}
    signed={name:means['stage266'][name]['self_wall']-means['main'][name]['self_wall'] for name in NAMES}
    valid=all(row['status']=='compatible_at_prespecified_resolution' for row in calibration.values())
    return dict(status='calibrated_diagnostic_only' if valid else 'observer_costs_unqualified',mean_regions=means,signed_self_wall_seconds=signed,additive_mean_rank_difference_seconds=math.fsum(signed.values()),ranking=sorted(signed,key=signed.get,reverse=True),recoverable_budget_seconds=None,limitation='Current measured locations only if both calibrations pass; no omission or full replacement has been tested. Retain negative offsets, inclusive nesting and calibration failures.')


def distribution(values):
    ordered=sorted(values)
    if not ordered:return None
    return dict(mean=statistics.mean(ordered),minimum=ordered[0],median=statistics.median(ordered),
                p90=ordered[math.ceil(.9*len(ordered))-1],maximum=ordered[-1])


def ownership_groups(all_rows,domains):
    """Group actual rank/step ownership; never infer owners from rank indices."""
    groups={'Domain_owners':[], 'without_Domains':[]}
    membership={key:{} for key in groups}
    for s in INCLUDED:
        total=0
        for r,steps in enumerate(all_rows):
            rows=steps[s];own=rows[IDS['OWNERS']-1];patch=rows[IDS['PATCHES']-1]
            d,b,samples=own['work']
            if samples!=4 or own['calls']!=4 or patch['work'][2]!=4:
                raise ValueError('Missing four ownership snapshots per step')
            # Geometric ownership is stable within this RK transaction. Block
            # and patch inventories can change; retain sums and sample counts so
            # stage averages can be recovered without confusing them with visits.
            if d%4:raise ValueError('Mixed Domain ownership within a step requires finer grouping')
            total+=d
            key='Domain_owners' if d else 'without_Domains'
            if not d and rows[IDS['BASIC']-1]['work'][0]:
                raise ValueError('Domain-free rank reports Domain producer patch work')
            groups[key].append(rows);membership[key].setdefault(str(s+1),[]).append(r)
        if total!=4*domains:raise ValueError('Global geometric Domain census differs')
        for symbol in ('CAPTURE_MPI','RESTRICT_MPI','TEMP_EDGES','INVERSE_MPI'):
            sent=sum(steps[s][IDS[symbol]-1]['work'][1] for steps in all_rows)
            received=sum(steps[s][IDS[symbol]-1]['work'][2] for steps in all_rows)
            if sent!=received:raise ValueError('Remote payload conservation differs: '+symbol)
    result={}
    for key,samples in groups.items():
        result[key]=dict(rank_step_samples=len(samples),ranks_by_step=membership[key],phases={})
        for p,name in enumerate(NAMES):
            result[key]['phases'][name]={metric:distribution([row[p][metric] for row in samples])
                for metric in ('self_wall','inclusive_wall','calls')}
            result[key]['phases'][name]['work']=[distribution([row[p]['work'][i] for row in samples]) for i in range(3)]
    return result


def summarize_scaling(rows,calibration):
    costs={str(n):[r['costs'] for r in rows if r['ranks']==n and r['observer']] for n in (160,256)}
    if any(len(v)!=2 for v in costs.values()):raise ValueError('Require two observer runs at each count')
    means={n:{name:statistics.mean(c['phases'][name]['self_wall'] for c in samples) for name in NAMES}
           for n,samples in costs.items()}
    delta={name:means['256'][name]-means['160'][name] for name in NAMES}
    qualified=all(c['status']=='compatible_at_prespecified_resolution' for c in calibration.values())
    return dict(status='calibrated_diagnostic_only' if qualified else 'observer_costs_unqualified',
        mean_rank_self_seconds=means, change_256_minus_160_seconds=delta,
        ranking_by_increase=sorted(delta,key=delta.get,reverse=True),
        ownership_runs={n:[c['ownership_groups'] for c in samples] for n,samples in costs.items()},
        limitation='Distributions pool rank/ordinary-step samples within each run, not independent repeats. Inclusive scopes overlap. Local durations are unsynchronized; MPI scopes include late-arrival waiting and transfer, CPU includes spinning. Patch inventories differ from executed patch-layer visits. No recoverable speedup or causal dominance established automatically.')
