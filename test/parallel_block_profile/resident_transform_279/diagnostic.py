"""Predeclared resident-transform comparisons; no automatic promotion or fine runs."""
import statistics
from common import runtime
RANKS=(83,160,256)

def environment(kind, profile='off', oracle=False, catalogue=False):
    if profile not in ('off','coarse'):raise ValueError('Only off/coarse modes are admitted')
    if profile!='off' and (kind not in ('previous','candidate') or oracle or catalogue):
        raise ValueError('Coarse attribution requires optimized, oracle-free, audit-free repair/control')
    env=runtime(oracle,catalogue,candidate=kind in ('candidate','checked','previous'))
    env.update(WAVETRISK_PROFILE_PARALLEL_BLOCKS=str(int(profile!='off')),
               WAVETRISK_PROFILE_BLOCK_DETAIL=str(int(profile!='off')),WAVETRISK_PROFILE_RHS_FINE='0')
    env['WAVETRISK_RESIDENT_TRANSFORM']=str(int(kind in ('candidate','checked')))
    env['WAVETRISK_RESIDENT_TRANSFORM_CHECK']=str(int(kind in ('candidate','checked') and oracle))
    return env

def schedule():
    rows=[]
    for n in RANKS:
        controls=['main','stage266'] if n<=160 else ['stage266']
        rows.extend((n,k,'off') for k in controls)
        rows.extend((n,k,'off') for k in ('previous','candidate','candidate','previous','candidate','previous','previous','candidate'))
        rows.extend((n,k,'off') for k in reversed(controls))
        # Alternate the paired coarse order; bare timing is always separate.
        rows.extend((n,k,'coarse') for k in (('candidate','previous') if n==160 else ('previous','candidate')))
    return rows

def selected(rows,n,kind,profile='off'):
    return [r for r in rows if r['phase']=='measurement' and (r['ranks'],r['kind'],r['profile'])==(n,kind,profile)]
def medians(rows):return [r['ordinary']['median'] for r in rows]
def mean(rows):return statistics.mean(medians(rows))
def spread(rows):return max(medians(rows))-min(medians(rows))
def region(rows,index):return statistics.mean(r['attribution']['regions'].get(str(index), {'mean_on_slowest_step_rank_seconds':0.0})['mean_on_slowest_step_rank_seconds'] for r in rows)
def check_measurements(rows):
    actual=[(r['ranks'],r['kind'],r['profile']) for r in rows if r['phase']=='measurement']
    if actual!=schedule():raise ValueError('Incomplete or reordered resident-transform schedule')
    if any(r.get('status')!='passed' or r['oracle'] or r['catalogue'] for r in rows if r['phase']=='measurement'):
        raise ValueError('Invalid, oracle or ownership-audit measurement')

def summary(rr):
    totals=[sum(x['seconds'] for x in r['ordinary']['included']) for r in rr]
    means=[total/len(r['ordinary']['included']) for total,r in zip(totals,rr)]
    return dict(run_medians_seconds=medians(rr),mean_run_median_seconds=mean(rr),run_spread_seconds=spread(rr),
        ordinary_run_totals_seconds=totals,ordinary_run_means_seconds=means,mean_ordinary_seconds=statistics.mean(means),
        resources=[r['resources'] for r in rr])

def summarize_diagnostic(rows):
    check_measurements(rows);by_rank={}
    for n in RANKS:
        variants={}
        for k,p in [('main','off'),('stage266','off'),('previous','off'),('candidate','off'),('previous','coarse'),('candidate','coarse')]:
            rr=selected(rows,n,k,p)
            if rr:variants[k+'-'+p]=summary(rr)
        off=mean(selected(rows,n,'candidate'));old=mean(selected(rows,n,'previous'))
        ca=selected(rows,n,'candidate','coarse');pa=selected(rows,n,'previous','coarse')
        remaps=ca[0]['attribution']['remap_groups']['remap']['step_indices']
        if remaps!=[4,9] or pa[0]['attribution']['remap_groups']['remap']['step_indices']!=remaps:raise ValueError('Remap membership differs')
        groups={}
        for k in ('previous','candidate','stage266','main'):
            rr=selected(rows,n,k)
            if rr:
                groups[k]={label:dict(run_medians_seconds=[statistics.median(x['seconds'] for x in r['ordinary']['included'] if (x['index'] in remaps)==event) for r in rr],
                   run_means_seconds=[statistics.mean(x['seconds'] for x in r['ordinary']['included'] if (x['index'] in remaps)==event) for r in rr]) for label,event in [('remap',True),('non_remap',False)]}
        regions={i:dict(name=ca[0]['attribution']['regions'][i]['name'],previous_seconds=region(pa,int(i)),candidate_seconds=region(ca,int(i)),
                    candidate_minus_previous_seconds=region(ca,int(i))-region(pa,int(i))) for i in ca[0]['attribution']['regions']}
        entry=dict(variants=variants,ordinary_subgroups=groups,candidate_minus_previous_percent=100*(off/old-1),
            ordinary_mean_candidate_minus_previous_percent=100*(variants['candidate-off']['mean_ordinary_seconds']/variants['previous-off']['mean_ordinary_seconds']-1),
            observed_run_spread_seconds=max(spread(selected(rows,n,k)) for k in ('candidate','previous','stage266')),
            profile_perturbation_percent={k:100*(mean(selected(rows,n,k,'coarse'))/mean(selected(rows,n,k))-1) for k in ('previous','candidate')},
            gap_to_stage266_seconds=off-mean(selected(rows,n,'stage266')),coarse_regions=regions)
        if n<=160:entry['gap_to_main_seconds']=off-mean(selected(rows,n,'main'))
        by_rank[str(n)]=entry
    scaling={}
    for k in ('previous','candidate','stage266'):
        speed=mean(selected(rows,160,k))/mean(selected(rows,256,k))
        mean_speed=by_rank['160']['variants'][k+'-off']['mean_ordinary_seconds']/by_rank['256']['variants'][k+'-off']['mean_ordinary_seconds']
        scaling[k]=dict(speedup_160_to_256=speed,relative_efficiency_percent=100*speed*160/256,
            ordinary_mean_speedup_160_to_256=mean_speed,ordinary_mean_relative_efficiency_percent=100*mean_speed*160/256)
    return dict(automatic_promotion=False,performance_acceptance='requires_review',by_ranks=by_rank,scaling=scaling,
        scope='Primary: mean of independent bare-run ordinary medians. Secondary: ordinary means/totals and remap/non-remap groups. Steps are not independent repetitions. Single paired coarse profiles are attribution, not speed measurements. Regions85-90 separate new transform preparation, execution and exchange; predecessor has none. Region78 is not a pure transport comparison. Review absolute time, run spread, RSS and 160-to-256 scaling separately. Preserve ownership/retirement gates; >=640 scaling and legacy parity remain open.')
