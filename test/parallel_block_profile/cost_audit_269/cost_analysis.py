"""Validate complete bounded records; retain overhead failures and signed costs."""
import math,re,statistics
from pathlib import Path
NAMES=('step','RK','tendency context','wavelet transaction','adapt algorithm','Domain inverse','native inverse','inverse setup','canonical publication','block installation','provisional refresh','block hydrostatic cache','native recomposition','catalogue','source blocks','migration','post-grid refresh','native compression','native scalar forward','native vector forward','block ghosts','native velocity source','native velocity gradient','Domain surface pressure','Domain boundary','accepted state refresh','inverse MPI wait/traffic')
WORK={8:('plan rebuilds','unused','unused'),9:('published interior nodes','published interior values','publication calls'),10:('installed interior slots','installed boundary slots','actual installations'),12:('rebuilt columns','rebuilt physical layers','rebuilt blocks'),27:('posted remote requests','sent remote doubles','received remote doubles')}
INCLUDED=(1,2,3,4,5,6,9,10,11)

def read_costs(directory,ranks,block,oracle=False):
    directory=Path(directory)
    files=sorted(directory.glob('cost-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing/extra observer rank files')
    all_rows=[]
    for rank in range(ranks):
        p=directory/f'cost-rank-{rank:06d}.txt';lines=p.read_text().splitlines()
        if lines[0].split()!=['STAGE269_COST_V1',str(rank),'12',str(len(NAMES))]:raise ValueError('Observer header differs')
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
                total=sum(row['self_'+metric] for row in rows)
                elapsed=rows[0]['inclusive_'+metric]
                if elapsed<=0 or abs(total-elapsed)>max(1e-6,elapsed*1e-6):raise ValueError('Observer phase conservation failed')
            if block:
                expected=7 if oracle else 4
                if rows[9]['work'][2]!=expected or rows[10]['calls']!=3 or rows[25]['calls']!=1:raise ValueError('Installation/provisional/accepted transaction count differs')
            elif any(rows[p-1]['calls'] for p in (7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,26,27)):
                raise ValueError('Native block scope in legacy main')
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
    result['mean_longest_local_step_seconds']=sum(result['longest_local_step_rank_phases'].values())
    result['limitations']='Self times add within a rank/step. Inclusive times overlap. Longest-local-step composition selects one rank per step, not independent region maxima; unsynchronized local durations do not prove the global critical path. CPU includes MPI spinning. Work counters cover named inverse/cache products, not all solver traffic; doubles are logical payload counts, not physical memory traffic.'
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
    return dict(status='calibrated_diagnostic_only' if valid else 'observer_costs_unqualified',mean_regions=means,signed_self_wall_seconds=signed,additive_mean_rank_difference_seconds=sum(signed.values()),ranking=sorted(signed,key=signed.get,reverse=True),recoverable_budget_seconds=None,limitation='Current measured locations only if both calibrations pass; no omission or full replacement has been tested. Retain negative offsets, inclusive nesting and calibration failures.')
