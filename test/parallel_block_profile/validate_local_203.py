#!/usr/bin/env python3
"""Record local restart evidence without making a performance claim."""
import argparse,json
from pathlib import Path
from local_j4 import compare_windows
from cluster_199 import workload_gate
from compare_restart_checkpoints import compare
from checkpoint_gate import accept
from analyze_hardware_201 import collect
from analyze_scaling_200 import collect as scaling
from prepare_cluster_194 import file_digest
from analyze_batching_203 import inverse_ledger

def validate(a):
    report={'stage':203,'status':'passed','scope':'J4/J6, four ranks. Local correctness only; tests ran under variable/concurrent load. No timing claim. Combined checked-plus-oracle acceptance is deferred to the 83-rank cluster run.', 'runs':{}}
    ref=json.loads((a.reference/'results.json').read_text())
    repo=Path(__file__).resolve().parents[2]
    for name in ('candidate','runtime','oracle','observer'):
        root=getattr(a,name);r=json.loads((root/'results.json').read_text())
        if r['status']!='passed' or r['ranks']!=4:raise ValueError('Incomplete/incorrect local run '+name)
        if file_digest(root/'climateJ4')!=r['binary_sha256']:raise ValueError('Binary changed '+name)
        for n in ('parallel_block.f90','parallel_block_inverse.f90','parallel_block_mpi.f90','parallel_block_scalar_storage.f90'):
            if r['source_sha256'].get('src/'+n)!=file_digest(repo/'src'/n):raise ValueError('Prototype source differs '+name+': '+n)

        exact=compare_windows(a.reference,root,'localJ4J6_checkpoint_0002.bin.zst',expected_new=[3])
        if not exact['passed']:raise ValueError(exact['failures'])
        if not workload_gate(ref['analysis'],r['analysis'])['passed']:raise ValueError('Workload differs '+name)
        log=(root/'run.log').read_text()
        if any('Stage203 '+k+' = on' not in log for k in ('scalar fields','inverse workspace')):raise ValueError('Prototype disabled '+name)
        if name=='oracle' and any(r['environment'].get('WAVETRISK_VALIDATE_BLOCK_'+k)!='1' for k in ('DYNAMICS','ADAPTATION','REMAP')):raise ValueError('Missing oracles')
        row=dict(path=str(root),status='passed',binary_sha256=r['binary_sha256'],exact_to_stage202=True,reloads=r['reloads'],post_restart_remap=True,steps=len(r['analysis']['step_records']))
        if name=='observer':
            h=collect(root,4,0);s=scaling(root,4,0)
            if h['quality_flags']:raise ValueError(h['quality_flags'])
            if len(s['workflow']['steps'])!=row['steps']:raise ValueError('Observer step coverage differs')
            row.update(hardware_rank_coverage=h['ranks'],coarse_clock_consistency='passed',inverse_ledger=inverse_ledger(log,4))
        report['runs'][name]=row
    legacy=compare(a.candidate/'localJ4J6_checkpoint_0003.bin.zst',a.legacy/'localJ4J6_checkpoint_0003.bin.zst',expected_domains=40)
    failures=accept(legacy,exact=False)
    if failures:raise ValueError(failures)
    report['legacy_screen']={'passed':True,'reference':str(a.legacy),'scope':'Existing fixed numerical screen; not bitwise equivalence.'}
    repo=Path(__file__).resolve().parents[2]
    report['prototype_sources']={n:file_digest(repo/n) for n in ('Makefile','src/parallel_block.f90','src/parallel_block_inverse.f90','src/parallel_block_mpi.f90','src/parallel_block_scalar_storage.f90')}
    a.out.write_text(json.dumps(report,indent=2)+'\n');print('PASS local Stage203 restart, workload, numerical and observer gates')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('reference','candidate','runtime','oracle','observer','legacy','out'):p.add_argument('--'+n,type=Path,required=True)
    validate(p.parse_args())
