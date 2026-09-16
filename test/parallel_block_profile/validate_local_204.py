#!/usr/bin/env python3
"""Check observer-only J4 restart equivalence and disjoint clock accounting."""
import argparse,json
from pathlib import Path
from local_j4 import compare_windows
from analyze_operations_204 import collect
from analyze_hardware_201 import collect as hardware
from analyze_scaling_200 import collect as scaling
from cluster_199 import workload_gate
from compare_restart_checkpoints import compare
from checkpoint_gate import accept
from build_restriction_204 import digest,verify

def validate(a):
    manifest=json.loads((a.kit/'stage204-sources.json').read_text())
    result={'stage':204,'status':'passed','scope':'Four-rank J4/J6 numerical and observer checks only; no performance claim. Cluster calibration remains mandatory.','runs':{}}
    for side,key in [('baseline','baseline_profile_sources'),('candidate','profiled_sources'),('legacy','legacy_profile_sources')]:
        root=getattr(a,side);ref=getattr(a,'reference_'+side);r=json.loads((root/'results.json').read_text());rr=json.loads((ref/'results.json').read_text())
        if r['status']!='passed' or r['ranks']!=4:raise ValueError('Incomplete local run '+side)
        if r['environment'].get('WAVETRISK_OPERATIONS204')!='1':raise ValueError('Operation observer disabled')
        verify(Path(r['build']),manifest[key])
        if digest(root/'climateJ4')!=r['binary_sha256']:raise ValueError('Local binary identity differs')
        exact=compare_windows(ref,root,'localJ4J6_checkpoint_0002.bin.zst',expected_new=[3])
        if not exact['passed'] or not workload_gate(rr['analysis'],r['analysis'])['passed']:raise ValueError('Observer changed numerics/workload '+side)
        op=collect(root,4,1);h=hardware(root,4,0);sc=scaling(root,4,0)
        if h['quality_flags'] or len(sc['workflow']['steps'])!=12:raise ValueError('Invalid clock/step coverage')
        result['runs'][side]={'path':str(root),'status':'passed','binary_sha256':r['binary_sha256'],'source_inventory':manifest[key],
            'steps':len(r['analysis']['step_records']),'exact_to_unmodified_reference':str(ref),'reloads':r['reloads'],
            'unassigned_by_context':{row['context']:row['seconds'] for row in op['ledger'] if row['operation']==0},
            'observed_operations':sorted({row['operation'] for row in op['ledger'] if row['calls']>0}),
            'observer_conservation':'passed (exclusive contexts, MPI owners, inclusive bounds and full step)',
            'build_log_sha256':digest(Path('/private/tmp/stage204-build-'+side+'-final.log'))}
    cross=compare(a.candidate/'localJ4J6_checkpoint_0003.bin.zst',a.legacy/'localJ4J6_checkpoint_0003.bin.zst',expected_domains=40)
    failures=accept(cross,exact=False)
    if failures:raise ValueError(failures)
    result['legacy_screen']={'passed':True,'scope':'Unchanged fixed tolerance, not bitwise equality.'}
    result['tests']=['Actual Fortran observer at O0/fcheck=all: nested/context-crossing scopes, disabled mode, checkpoint exclusion, MPI ownership, corrupted row rejection.',
      'Synthetic signed-ranking conservation; calibration flags and missing repeats rejected.',
      'Real seven-build patch/source/hash/submit sequence with compiler stand-in; modified executable rejected.']
    result['excluded_pilots']=['Initial sandboxed MPI launch failed before solver initialization.','Earlier 56-scope candidate/baseline pilots passed exact restart but are superseded by the final 65-scope observers.']
    a.out.write_text(json.dumps(result,indent=2)+'\n');print('PASS Stage204 local observer and numerical validation')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('kit','baseline','candidate','legacy','reference-baseline','reference-candidate','reference-legacy','out'):p.add_argument('--'+n,type=Path,required=True)
    validate(p.parse_args())
