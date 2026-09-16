#!/usr/bin/env python3
"""Freeze production candidate, independent Stage196 reference and focused probes."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from instrument_workflow_195 import instrument as workflow
from instrument_operations_196 import instrument as operations
from instrument_operations_197 import upgrade,REGIONS,COUNTERS
BASELINE='9bc3a19556dab79b3cba3944cf5df53517389da3'
HELPERS=('build_restriction_197.py','build_197.sbatch','submit_restriction_197.py','cluster_197.py','cluster_197.sbatch',
 'rank_metrics_190.py','rank_launch_190.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
 'profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs.py','analyze_workflow_195.py','analyze_operations_197.py')

def prepare(repo,out,previous):
 repo=repo.resolve();out=out.resolve();previous=previous.resolve()
 if out.exists() or out.is_relative_to(repo):raise ValueError('Fresh kit outside repository required')
 old=json.loads((previous/'stage196-sources.json').read_text());tracked=old['baseline_sources']
 if subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()!=BASELINE:raise ValueError('Prepare from Stage196 HEAD')
 out.mkdir();manifest={k:old[k] for k in ('legacy_revision','legacy_original_sources','legacy_io_sources','legacy_io_patch_sha256','legacy_archive_sha256')}
 for name in HELPERS:shutil.copy2(Path(__file__).with_name(name),out/name)
 for name in ('physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch'):shutil.copy2(previous/name,out/name)
 subprocess.run(['git','-C',str(repo),'bundle','create',str(out/'stage196-baseline.bundle'),'HEAD','^a0c75adb'],check=True)
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',BASELINE]))) as t:t.extractall(root,filter='data')
  for n,h in tracked.items():
   if file_digest(root/n)!=h:raise ValueError('Stage196 baseline mismatch '+n)
  before={n:(root/n).read_bytes() for n in tracked if not (root/n).is_symlink()}
  for n in before:shutil.copy2(repo/n,root/n)
  candidate={n:file_digest(root/n) for n in tracked}
  changed=[n for n in tracked if candidate[n]!=tracked[n]]
  if changed!=['src/parallel_block_mpi.f90']:raise ValueError('Unexpected production scope '+str(changed))
  (out/'stage197-candidate.patch').write_text(patch(before,root,changed))
  before={n:(root/n).read_bytes() for n in tracked if not (root/n).is_symlink()}
  workflow(root,True);operations(root,True);block_inventory=upgrade(root,True)
  names=[*tracked,'src/workflow_profile_195.f90'];profiled={n:file_digest(root/n) for n in names}
  (out/'stage197-observer.patch').write_text(patch(before,root,[n for n in names if profiled[n]!=candidate.get(n)]))
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['legacy_io_sources'] if not (root/n).is_symlink()}
  p=previous/'legacy-profile-196.patch'
  if file_digest(p)!=old['legacy_profile_patch_sha256']:raise ValueError('Legacy observer input changed')
  subprocess.run(['git','apply','--whitespace=nowarn',str(p)],cwd=root,check=True)
  legacy_inventory=upgrade(root,False)
  legacy={n:file_digest(root/n) for n in old['legacy_profile_sources']}
  (out/'legacy-profile-197.patch').write_text(patch(before,root,[n for n in legacy if legacy[n]!=old['legacy_io_sources'].get(n)]))
 manifest.update(baseline_commit=BASELINE,baseline_stage=196,baseline_sources=tracked,candidate_sources=candidate,profiled_sources=profiled,
  makefile_sha256=tracked['Makefile'],candidate_makefile_sha256=candidate['Makefile'],
  candidate_patch_sha256=file_digest(out/'stage197-candidate.patch'),observer_patch_sha256=file_digest(out/'stage197-observer.patch'),
  legacy_profile_sources=legacy,legacy_profile_patch_sha256=file_digest(out/'legacy-profile-197.patch'),
  observer_sha256=file_digest(Path(__file__).with_name('workflow_profile_197.f90')))
 assert profiled['src/workflow_profile_195.f90']==legacy['src/workflow_profile_195.f90']==manifest['observer_sha256']
 (out/'stage197-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
 scope=json.loads((previous/'operation-instrumentation.json').read_text())
 scope['region_names'].update({str(k):v for k,v in REGIONS.items()});scope['counter_names'].update({str(k):v for k,v in COUNTERS.items()})
 scope['stage197_inventory']={'block':block_inventory,'legacy':legacy_inventory}
 scope['counter_scope']='Counters 4,5,12-15 repaired. Selected work counters only; legacy direct-flux scope includes shared step1 work and not all special boundary fluxes. Do not interpret scope ratios as arithmetic-only cost. Dependencies use rank-local clocks.'
 from analyze_operations_197 import GROUPS,BOUNDARY_GROUPS
 scope.update(groups=GROUPS,boundary_groups=BOUNDARY_GROUPS)
 (out/'operation-instrumentation.json').write_text(json.dumps(scope,indent=2)+'\n')
 shutil.copy2(repo/'docs/parallel-blocks-cluster-197.md',out/'README.md')
 archive=out.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
 report={'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':BASELINE,
  'candidate_patch_sha256':manifest['candidate_patch_sha256'],'observer_patch_sha256':manifest['observer_patch_sha256'],
  'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'observer_sha256':manifest['observer_sha256'],
  'helpers_sha256':{n:file_digest(out/n) for n in HELPERS}}
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('repo','out','stage196-kit'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();prepare(a.repo,a.out,a.stage196_kit)
