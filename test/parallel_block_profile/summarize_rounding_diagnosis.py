"""Keep unmodified-reference failures distinct from controlled causal replays."""
import argparse
import json
from pathlib import Path
from checkpoint_gate import accept,POLICY
from compare_restart_checkpoints import compare
from compare_physics_probe import load
from experiment import digest


def checkpoints(reference,candidate,exact=False):
    result={}
    for cp in (5,6):
        name=f'test_checkpoint_{cp:04d}.bin.zst'
        c=compare(reference/name,candidate/name)
        result[str(cp)]={'reference_sha256':digest(reference/name),'candidate_sha256':digest(candidate/name),
                        'failures':accept(c,exact),'comparison':c}
    return result


def event(left,right,step,domain,node,level,field):
    pattern=f'physics-input-step-{step}-rank-*.bin'
    a=load(str(left/pattern))[(step,domain,node)]
    b=load(str(right/pattern))[(step,domain,node)]
    nz=a['header'][9];k=level-1
    offset=3+nz+7*k+{'U':0,'V':1,'W':2}[field]
    raw_a=float(a['raw'][offset]);raw_b=float(b['raw'][offset])
    x=float(a['fields'][field][k]);y=float(b['fields'][field][k])
    midpoint=(x+y)/2
    result={'step':step,'domain':domain,'node':node,'level':level,'field':field,
            'legacy_double':raw_a,'block_double':raw_b,'double_difference':abs(raw_a-raw_b),
            'legacy_float32':x,'block_float32':y,'float32_difference':abs(x-y),
            'legacy_bits':int(a['fields'][field].view('u4')[k]),
            'block_bits':int(b['fields'][field].view('u4')[k]),'midpoint':midpoint,
            'opposite_sides_of_midpoint':min(raw_a,raw_b)<midpoint<max(raw_a,raw_b)}
    if not result['opposite_sides_of_midpoint'] or abs(result['legacy_bits']-result['block_bits'])!=1:
        raise ValueError('Event is not a one-ULP rounding-midpoint crossing')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('reference','focused','global-trace','one-replay','two-replay','out'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    result={'policy_unchanged':POLICY,
            'scope':'Two observed float32 midpoint crossings; replay is causal evidence, not production acceptance',
            'original_legacy_vs_block':checkpoints(a.reference/'legacy',a.reference/'block'),
            'one_replay_vs_block':checkpoints(a.one_replay/'legacy',a.reference/'block'),
            'two_replay_vs_block':checkpoints(a.two_replay/'legacy',a.reference/'block'),
            'instrumentation_transparency':{
                'focused_legacy':checkpoints(a.reference/'legacy',a.focused/'legacy',True),
                'focused_block':checkpoints(a.reference/'block',a.focused/'block',True),
                'global_legacy_replay':checkpoints(a.one_replay/'legacy',a.global_trace/'legacy',True),
                'global_block':checkpoints(a.reference/'block',a.global_trace/'block',True)},
            'events':[
                event(a.focused/'legacy',a.focused/'block',260,77,123,26,'U'),
                event(a.global_trace/'legacy',a.global_trace/'block',262,8,123,28,'W')]}
    for reports in result['instrumentation_transparency'].values():
        if any(c['failures'] for c in reports.values()):raise ValueError('Instrumentation transparency failed')
    log=(a.two_replay/'legacy/run.log').read_text()
    if log.count('PHYSICS REPLAY U:')!=1 or log.count('PHYSICS REPLAY W:')!=1:
        raise ValueError('Expected exactly the two declared replay interventions')
    if 'Total cpu time' not in log:raise ValueError('Causal run incomplete')
    result['two_replay_run_identity']=json.loads((a.two_replay/'legacy/identity.json').read_text())
    result['two_replay_screen_passed']=not any(c['failures'] for c in result['two_replay_vs_block'].values())
    with a.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
