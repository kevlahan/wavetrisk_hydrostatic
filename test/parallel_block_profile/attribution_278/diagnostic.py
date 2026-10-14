"""Bounded diagnostic schedule and explicit, non-causal interpretation of profiles."""
import statistics
from common import runtime
RANKS=(83,160,256)

def environment(kind, profile='off', oracle=False, catalogue=False):
    if profile not in ('off','coarse','fine'):
        raise ValueError('Unknown profile mode')
    if profile!='off' and (kind!='candidate' or oracle or catalogue):
        raise ValueError('Attribution requires optimized, oracle-free, audit-free candidate')
    env=runtime(oracle,catalogue,candidate=kind in ('candidate','checked','v5'))
    env.update(WAVETRISK_PROFILE_PARALLEL_BLOCKS=str(int(profile!='off')),
               WAVETRISK_PROFILE_BLOCK_DETAIL=str(int(profile!='off')),
               WAVETRISK_PROFILE_RHS_FINE=str(int(profile=='fine')))
    return env

def schedule():
    """Six repeated measurements/rank; reverse the binary/profile order halfway."""
    rows=[]
    for n in RANKS:
        if n<=160:rows.append((n,'main','off'))
        rows.append((n,'stage266','off'))
        rows.extend((n,k,p) for k,p in [('v5','off'),('candidate','off'),('candidate','coarse'),
                                     ('candidate','coarse'),('candidate','off'),('v5','off')])
    return rows

def selected(rows,n,kind,profile='off'):
    return [r for r in rows if r['phase']=='measurement' and (r['ranks'],r['kind'],r['profile'])==(n,kind,profile)]

def medians(rows):return [r['ordinary']['median'] for r in rows]
def mean(rows):return statistics.mean(medians(rows))
def spread(rows):return max(medians(rows))-min(medians(rows))
def region(rows,index):
    return statistics.mean(r['attribution']['regions'][str(index)]['mean_on_slowest_step_rank_seconds'] for r in rows)

def check_measurements(rows):
    actual=[(r['ranks'],r['kind'],r['profile']) for r in rows if r['phase']=='measurement']
    if actual!=schedule():raise ValueError('Incomplete or reordered diagnostic schedule')
    if any(r.get('status')!='passed' or r['oracle'] or r['catalogue'] for r in rows if r['phase']=='measurement'):
        raise ValueError('Invalid, oracle or ownership-audit measurement')

def fine_decisions(rows):
    check_measurements(rows)
    result={}
    for n in RANKS:
        off=selected(rows,n,'candidate'); coarse=selected(rows,n,'candidate','coarse')
        gap=mean(off)-mean(selected(rows,n,'stage266'))
        noise=max(spread(off),spread(coarse),spread(selected(rows,n,'v5')))
        exchange=region(coarse,78)
        # This is a diagnostic admission rule, never a performance acceptance threshold.
        budget_material=gap>0 and exchange>=max(gap,noise)
        result[str(n)]=dict(run_fine=budget_material,exchange_seconds=exchange,
            gap_to_stage266_seconds=gap,observed_run_spread_seconds=noise,
            reasons=['Exchange budget can cover the positive 266 gap and exceeds observed run spread'] if budget_material else [])
    low=selected(rows,160,'candidate','coarse'); high=selected(rows,256,'candidate','coarse')
    whole_growth=mean(selected(rows,256,'candidate'))-mean(selected(rows,160,'candidate'))
    exchange_growth=region(high,78)-region(low,78)
    uncertainty=max(result['160']['observed_run_spread_seconds'],result['256']['observed_run_spread_seconds'])
    if whole_growth>0 and exchange_growth>uncertainty:
        for n in ('160','256'):
            result[n]['run_fine']=True
            result[n]['reasons'].append('Exchange growth accompanies the 160-to-256 slowdown and exceeds observed run spread')
    for row in result.values():row['scope']='Admission to a bounded fine/off pair; not proof of network dominance. Run spread is a confounder screen, not a confidence interval.'
    return result

def summarize_diagnostic(rows):
    check_measurements(rows)
    by_rank={}
    for n in RANKS:
        variants={}
        for k,p in [('main','off'),('stage266','off'),('v5','off'),('candidate','off'),('candidate','coarse')]:
            rr=selected(rows,n,k,p)
            if rr:variants[k+'-'+p]=dict(run_medians_seconds=medians(rr),mean_run_median_seconds=mean(rr),
                run_spread_seconds=spread(rr),resources=[r['resources'] for r in rr])
        off=mean(selected(rows,n,'candidate'));old=mean(selected(rows,n,'v5'));coarse=selected(rows,n,'candidate','coarse')
        remaps=coarse[0]['attribution']['remap_groups']['remap']['step_indices']
        if any(r['attribution']['remap_groups']['remap']['step_indices']!=remaps for r in coarse):
            raise ValueError('Remap profile membership changed')
        groups={}
        for k in ('v5','candidate','stage266','main'):
            rr=selected(rows,n,k)
            if rr:
                groups[k]={label:[statistics.median(x['seconds'] for x in r['ordinary']['included'] if (x['index'] in remaps)==event) for r in rr]
                    for label,event in [('remap',True),('non_remap',False)]}
        regions={i:dict(name=coarse[0]['attribution']['regions'][i]['name'],
                        mean_on_slowest_step_rank_seconds=region(coarse,int(i))) for i in coarse[0]['attribution']['regions']}
        entry=dict(variants=variants,ordinary_subgroup_run_medians_seconds=groups,
            profile_off_vs_original_v5_percent=100*(off/old-1),coarse_vs_off_percent=100*(mean(coarse)/off-1),
            gap_to_stage266_seconds=off-mean(selected(rows,n,'stage266')),coarse_regions=regions)
        if n<=160:entry['gap_to_main_seconds']=off-mean(selected(rows,n,'main'))
        fine=[r for r in rows if r['phase']=='fine' and r['ranks']==n]
        control=[r for r in rows if r['phase']=='fine-control' and r['ranks']==n]
        if fine:
            if len(fine)!=1 or len(control)!=1:raise ValueError('Fine run lacks its off control')
            entry['fine_pair']=dict(fine_median_seconds=mean(fine),off_median_seconds=mean(control),
                observed_fine_vs_off_percent=100*(mean(fine)/mean(control)-1),
                scope='Single paired runs; diagnostic perturbation and run variation are not separated.')
        by_rank[str(n)]=entry
    low=by_rank['160']['coarse_regions'];high=by_rank['256']['coarse_regions']
    deltas={i:dict(name=high[i]['name'],seconds_256_minus_160=high[i]['mean_on_slowest_step_rank_seconds']-low[i]['mean_on_slowest_step_rank_seconds']) for i in high}
    scaling={}
    for k in ('v5','candidate','stage266'):
        speed=mean(selected(rows,160,k))/mean(selected(rows,256,k))
        scaling[k]=dict(speedup_160_to_256=speed,relative_efficiency_percent=100*speed*160/256)
    return dict(automatic_promotion=False,repair_validated=False,by_ranks=by_rank,scaling=scaling,
        coarse_region_deltas_256_minus_160=deltas,fine_decisions=fine_decisions(rows),
        scope='Diagnostic source only. Run medians are repetitions; steps are not independent repetitions. Profile changes and MPI waits do not prove causation or isolate transport. Absolute time, memory, ownership retirement and useful >=640-rank scaling remain separate gates.')
