#!/usr/bin/env python3
"""Freeze unchanged Stage197 numerics with small scaling observers and protocol."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from prepare_cluster_199 import HELPERS as PREVIOUS_HELPERS
HELPERS=('build_restriction_200.py','build_200.sbatch','submit_restriction_200.py','cluster_200.py','cluster_200.sbatch',
 'analyze_scaling_200.py','placement_200.py','rank_launch_200.py','collect_cluster_logs_200.py')

def prepare(repo,out,previous):
    repo=repo.resolve();out=out.resolve();previous=previous.resolve();base=Path(__file__).parent
    if out.exists() or out.with_suffix('.tar.gz').exists() or out.is_relative_to(repo):raise ValueError('Fresh external kit/archive required')
    old=json.loads((previous/'stage199-sources.json').read_text())
    for n,h in old['candidate_sources'].items():
        if file_digest(repo/n)!=h:raise ValueError('Production numerics differ: '+n)
    out.mkdir();manifest=dict(old)
    for name in PREVIOUS_HELPERS:shutil.copy2(previous/name,out/name)
    for name in HELPERS:shutil.copy2(base/name,out/name)
    for name in ('stage196-baseline.bundle','physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch','operation-instrumentation.json'):
        shutil.copy2(previous/name,out/name)
    inventory=json.loads((out/'operation-instrumentation.json').read_text())
    mpi=[13,18,19,34,37,38,39,40,60,72,75]
    inventory={k:v for k,v in inventory.items() if k in ('counter_names','volume_names','region_names')}
    inventory['region_names']={k:v for k,v in inventory['region_names'].items() if int(k) in [0,*mpi]}
    inventory.update(stage=200,selected_mpi_operation_ids=mpi,
      modes={'0':'Common phases only; MPI and work are unavailable, not zero measured work.','1':'Common phases, selected MPI scopes and integer counters.'},
      scope='No operation/batch/dependency tape or per-access clocks. MPI includes producer waiting and progress; outside-MPI elapsed includes unclassified work and MPI calls not covered by designated scopes.')
    (out/'operation-instrumentation.json').write_text(json.dumps(inventory,indent=2)+'\n')
    shutil.copy2(previous/'stage199-candidate.patch',out/'stage200-candidate.patch')
    if file_digest(out/'stage200-candidate.patch')!=old['candidate_patch_sha256']:raise ValueError('Candidate patch changed')
    for side in ('block','legacy'):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            if side=='block':
                raw=subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']])
                with tarfile.open(fileobj=io.BytesIO(raw)) as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(out/'stage200-candidate.patch')],cwd=root,check=True)
                original=old['candidate_sources'];oldpatch='stage199-observer.patch';expected=old['observer_patch_sha256'];key='profiled_sources';target='stage200-observer.patch'
            else:
                with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
                subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
                original=old['legacy_io_sources'];oldpatch='legacy-profile-199.patch';expected=old['legacy_profile_patch_sha256'];key='legacy_profile_sources';target='legacy-profile-200.patch'
            before={n:(root/n).read_bytes() for n in original if not (root/n).is_symlink()}
            if file_digest(previous/oldpatch)!=expected:raise ValueError('Frozen observer changed')
            subprocess.run(['git','apply','--whitespace=nowarn',str(previous/oldpatch)],cwd=root,check=True)
            for n,h in old[key].items():
                if file_digest(root/n)!=h:raise ValueError('Frozen profile source differs: '+n)
            (root/'src/workflow_profile_195.f90').write_bytes((base/'workflow_profile_200.f90').read_bytes())
            (root/'src/workflow_clock_199.c').write_bytes((base/'workflow_clock_200.c').read_bytes())
            manifest[key]={n:file_digest(root/n) for n in old[key]}
            (out/target).write_text(patch(before,root,[n for n,h in manifest[key].items() if h!=original.get(n)]))
    manifest.update(observer_patch_sha256=file_digest(out/'stage200-observer.patch'),legacy_profile_patch_sha256=file_digest(out/'legacy-profile-200.patch'),
                    observer_sha256=file_digest(base/'workflow_profile_200.f90'),diagnosis_only=True,stage=200)
    (out/'stage200-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for source,target in [('parallel-blocks-cluster-200.md','README.md'),('parallel-blocks-stage200.md','EXPERIMENT.md')]:shutil.copy2(repo/'docs'/source,out/target)
    (out/'helpers-sha256.json').write_text(json.dumps({n:file_digest(out/n) for n in (*PREVIOUS_HELPERS,*HELPERS)},indent=2)+'\n')
    validation=repo/'docs/parallel-blocks-stage200-local-results.json'
    if validation.is_file():shutil.copy2(validation,out/'LOCAL_VALIDATION.json')
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    report={'stage':200,'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],
      'numerical_candidate_patch_sha256':manifest['candidate_patch_sha256'],'observer_patch_sha256':manifest['observer_patch_sha256'],
      'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'scope':'Five builds, unchanged numerical candidate/legacy; 42 sequential solver runs in four configurations.'}
    print(json.dumps(report,indent=2));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('repo','out','stage199-kit'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.stage199_kit)
