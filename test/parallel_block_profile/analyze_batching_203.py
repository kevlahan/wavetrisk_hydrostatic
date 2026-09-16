"""A substantial-gain screen for the Stage203 prototype; no confidence claims."""
import statistics as st
import re
from collections import defaultdict
from analyze_hardware_201 import PHASES

def summarize(result):
    runs=list(result['runs'].values());out={'bare':{},'paired_recovery_seconds':[], 'phases':{},'calibration':{}}
    for side in ('baseline','candidate','legacy'):
        values=[r['ordinary_step_median'] for r in runs if r['side']==side and r['phase']=='bare']
        if len(values)!=2:raise ValueError('Two bare repetitions required per implementation')
        out['bare'][side]={'samples':values,'median_seconds':st.median(values)}
    for repeat in (1,2):
        pair={s:next(r for r in runs if r['side']==s and r['phase']=='bare' and r['round']==repeat) for s in ('baseline','candidate')}
        out['paired_recovery_seconds'].append(pair['baseline']['ordinary_step_median']-pair['candidate']['ordinary_step_median'])
    b,c,l=(out['bare'][s]['median_seconds'] for s in ('baseline','candidate','legacy'))
    out.update(recovery_seconds=b-c,recovery_percent=100*(b-c)/b,baseline_legacy_ratio=b/l,candidate_legacy_ratio=c/l,
               remaining_gap_seconds=c-l,fraction_of_measured_gap_removed=(b-c)/(b-l) if b>l else None)
    for side in ('baseline','candidate'):
        coarse=next(r for r in runs if r['side']==side and r['phase']=='coarse')
        for r in (r for r in runs if r['side']==side and 'hardware' in r):
            changes={}
            for phase in PHASES:
                ref=coarse['hardware']['phases'][phase]['wall_seconds_rank_mean'];v=r['hardware']['phases'][phase]['wall_seconds_rank_mean']
                changes[phase]={'seconds':v-ref,'flag':abs(v-ref)>.001 and abs(v-ref)>.02*ref}
            bare=out['bare'][side]['median_seconds'];delta=100*(r['ordinary_step_median']/bare-1)
            out['calibration'][r['label']]={'printed_change_percent':delta,'phase_changes':changes,
                'review_required':abs(delta)>3 or any(x['flag'] for x in changes.values()) or bool(r['hardware']['quality_flags'])}
    for phase in ('whole step',*PHASES):
        rows=[]
        for repeat in (1,2):
            pair={s:next(r for r in runs if r['side']==s and r['phase']=='core' and r['round']==repeat) for s in ('baseline','candidate')}
            v={s:pair[s]['hardware']['phases'][phase] for s in pair}
            bi,ci=(v[s]['counts_rank_mean_per_step']['instructions'] for s in ('baseline','candidate'))
            bc,cc=(v[s]['counts_rank_mean_per_step']['cycles'] for s in ('baseline','candidate'))
            rows.append({'repeat':repeat,'instruction_reduction_percent':100*(bi-ci)/bi if bi else None,
                         'cycles_reduction_percent':100*(bc-cc)/bc if bc else None,'baseline':v['baseline'],'candidate':v['candidate'],
                         'observer_review_required':any(out['calibration'][pair[s]['label']]['review_required'] for s in pair)})
        out['phases'][phase]=rows
    legacy_coarse=next(r for r in runs if r['side']=='legacy' and r['phase']=='coarse')
    delta=100*(legacy_coarse['ordinary_step_median']/out['bare']['legacy']['median_seconds']-1)
    out['calibration'][legacy_coarse['label']]={'printed_change_percent':delta,
        'review_required':abs(delta)>3 or bool(legacy_coarse['hardware']['quality_flags']),
        'scope':'One fresh legacy coarse control; no independent repeated phase calibration.'}
    out['matched_coarse_phases']={}
    for phase in ('whole step',*PHASES):
        vals={side:next(r for r in runs if r['side']==side and r['phase']=='coarse')['hardware']['phases'][phase]['wall_seconds_rank_mean'] for side in ('baseline','candidate','legacy')}
        out['matched_coarse_phases'][phase]={**vals,'remaining_excess':vals['candidate']-vals['legacy'],'recovered':vals['baseline']-vals['candidate']}
    substantial=all(x>=.10 for x in out['paired_recovery_seconds'])
    clean=not any(v['review_required'] for v in out['calibration'].values())
    out['substantial_gain_screen_passed']=substantial
    out['decision']='review_for_promotion' if substantial and clean else 'review_observer_effects' if substantial else 'do_not_promote_insufficient_gain'
    out['scope']='Two reversed bare repeats are a prototype screen, not statistical confidence. Require >=0.10 s/step in BOTH pairs. PMU includes MPI progress/spinning; no per-symbol causal savings claim. All numerical/placement gates must pass. Memory counters and sampling are deliberately omitted.'
    return out

def write_report(root,result):
    s=result['summary'];lines=['# Stage203 prototype review','',f"Decision: **{s['decision']}** (subject to review).",'',
        '| Version | Bare repetitions (s/step) | Median |','|---|---|---:|']
    for side,row in s['bare'].items():lines.append(f"| {side} | {row['samples']} | {row['median_seconds']:.6f} |")
    lines.extend(['',f"Recovery: {s['recovery_seconds']:.6f} s/step; paired recoveries {s['paired_recovery_seconds']}.",
                  f"Candidate/legacy ratio: {s['candidate_legacy_ratio']:.4f}.",'',s['scope']])
    lines.extend(['','Same-job coarse phase comparison (rank-mean seconds/ordinary step):','',
                  '| Phase | Stage202 | Stage203 | Legacy | Remaining excess |','|---|---:|---:|---:|---:|'])
    for phase,row in s['matched_coarse_phases'].items():
        lines.append(f"| {phase} | {row['baseline']:.6f} | {row['candidate']:.6f} | {row['legacy']:.6f} | {row['remaining_excess']:.6f} |")
    lines.extend(['','Phase means are distinct from printed bare medians. The one legacy coarse control is checked against its bare timing; inspect observer calibration before interpreting residuals.'])
    (root/'REVIEW.md').write_text('\n'.join(lines)+'\n')


def inverse_ledger(log,ranks):
    rows=defaultdict(list)
    for rank,values in re.findall(r'^Stage203 inverse rank=(\d+) work=([ \d]+)$',log,re.M):
        values=[int(x) for x in values.split()]
        if len(values)!=6:raise ValueError('Invalid inverse ledger record')
        rows[int(rank)].append(values)
    if set(rows)!=set(range(ranks)):raise ValueError('Missing inverse ledger ranks')
    if len({len(v) for v in rows.values()})!=1:raise ValueError('Inverse transaction counts differ between ranks')
    totals={str(k):[sum(v[i] for v in values) for i in range(6)] for k,values in rows.items()}
    return {'columns':['replaced_gather_calls','replaced_scatter_calls','scalar_outputs','inner_cells',
                       'canonical_restored_values','exit_interior_scalar_values'],
            'per_rank_totals':totals,'transactions_per_rank':len(rows[0]),
            'scope':'Whole-run diagnostic work, including checkpoint steps. Counts are not seconds; do not add waiting savings to producer savings.'}


def phase_balance(hardware):
    """Keep absolute and relative imbalance separate; do not add phase maxima."""
    result={}
    for phase in ('whole step',*PHASES):
        steps=defaultdict(list)
        for cell in hardware['cells']:
            if cell['phase']==phase:steps[cell['step']].append(cell)
        values={}
        for field in ('wall_seconds','thread_cpu_seconds'):
            means=[st.mean(c[field] for c in cells) for cells in steps.values()]
            maxima=[max(c[field] for c in cells) for cells in steps.values()]
            values[field]={'rank_mean':st.mean(means),'mean_step_rank_max':st.mean(maxima),
                           'mean_step_max_minus_mean':st.mean(x-m for x,m in zip(maxima,means)),
                           'mean_step_max_over_mean':st.mean(x/m for x,m in zip(maxima,means) if m) if any(means) else None}
        result[phase]=values
    return {'phases':result,'scope':'Ordinary-step rank envelopes include MPI progress/waiting. Maxima in different phases need not belong to the same rank and must not be summed as a critical path.'}
