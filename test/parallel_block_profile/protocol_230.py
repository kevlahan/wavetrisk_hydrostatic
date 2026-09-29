"""Architectural cleanup: regression guard, no minimum performance gain."""
import math
import statistics
from protocol_229 import ordinary

POLICY = dict(clear_regression_fraction=0.03, near_threshold_fraction=0.005,
              maximum_confirmation_batches=2, minimum_gain=None)


def screen(a,b):
    if not all(math.isfinite(v) and v>0 for v in (a,b)):raise ValueError('Invalid run medians')
    change=b/a-1
    return dict(baseline=a,candidate=b,gain_seconds=a-b,regression_fraction=change,
                decision='review_regression' if change>=POLICY['clear_regression_fraction'] else 'advance',
                scope='Single exploratory pair; no confirmed speedup claim.')


def confirm(batches):
    a=[];b=[];paired=[]
    for batch in batches:
        if [r[0] for r in batch]!=['baseline','candidate','candidate','baseline']:
            raise ValueError('Require A-B-B-A')
        v=[r[1] for r in batch]
        if not all(math.isfinite(x) and x>0 for x in v):raise ValueError('Invalid run medians')
        a.extend((v[0],v[3]));b.extend(v[1:3]);paired.extend((v[1]/v[0]-1,v[2]/v[3]-1))
    change=statistics.mean(b)/statistics.mean(a)-1
    target=POLICY['clear_regression_fraction']
    near=abs(change-target)<=POLICY['near_threshold_fraction']
    disagree=min(paired)<target<=max(paired)
    if near or disagree:
        decision='repeat_abba' if len(batches)<POLICY['maximum_confirmation_batches'] else 'inconclusive'
    else:decision='review_regression' if change>=target else 'within_regression_guard'
    return dict(decision=decision,mean_baseline_run_median=statistics.mean(a),
                mean_candidate_run_median=statistics.mean(b),gain_seconds=statistics.mean(a)-statistics.mean(b),
                regression_fraction=change,paired_regression_fractions=paired,
                near_threshold=near,disagreement=disagree,independent_runs_per_version=len(a),
                scope='Separate run medians; screens excluded. No minimum gain; no timestep-level significance test.')
