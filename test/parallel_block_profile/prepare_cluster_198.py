#!/usr/bin/env python3
"""Freeze diagnostic-only Stage198 observers around validated Stage197 numerics."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from instrument_operations_198 import upgrade,REGIONS,COUNTERS
HELPERS=('build_restriction_198.py','build_198.sbatch','submit_restriction_198.py','cluster_198.py','cluster_198.sbatch','analyze_operations_198.py',
 'rank_metrics_190.py','rank_launch_190.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py','profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs_198.py','analyze_workflow_195.py')

def prepare(repo,out,previous):
 repo=repo.resolve();out=out.resolve();previous=previous.resolve()
 if out.exists() or out.is_relative_to(repo):raise ValueError('Fresh kit outside repository required')
 old=json.loads((previous/'stage197-sources.json').read_text())
 for n,h in old['candidate_sources'].items():
  if file_digest(repo/n)!=h:raise ValueError('Numerical source changed since validated Stage197: '+n)
 out.mkdir();manifest=dict(old)
 for name in HELPERS:shutil.copy2(Path(__file__).with_name(name),out/name)
 for name in ('stage196-baseline.bundle','physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch'):shutil.copy2(previous/name,out/name)
 shutil.copy2(previous/'stage197-candidate.patch',out/'stage198-candidate.patch')
 if file_digest(out/'stage198-candidate.patch')!=old['candidate_patch_sha256']:raise ValueError('Candidate numerics changed')
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'stage198-candidate.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['candidate_sources'] if not (root/n).is_symlink()}
  if file_digest(previous/'stage197-observer.patch')!=old['observer_patch_sha256']:raise ValueError('Previous block observer changed')
  subprocess.run(['git','apply',str(previous/'stage197-observer.patch')],cwd=root,check=True)
  block_inventory=upgrade(root,True)
  profiled={n:file_digest(root/n) for n in [*old['profiled_sources'],'src/workflow_clock_198.c']}
  (out/'stage198-observer.patch').write_text(patch(before,root,[n for n in profiled if profiled[n]!=old['candidate_sources'].get(n)]))
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['legacy_io_sources'] if not (root/n).is_symlink()}
  if file_digest(previous/'legacy-profile-197.patch')!=old['legacy_profile_patch_sha256']:raise ValueError('Previous legacy observer changed')
  subprocess.run(['git','apply','--whitespace=nowarn',str(previous/'legacy-profile-197.patch')],cwd=root,check=True)
  legacy_inventory=upgrade(root,False)
  legacy={n:file_digest(root/n) for n in [*old['legacy_profile_sources'],'src/workflow_clock_198.c']}
  (out/'legacy-profile-198.patch').write_text(patch(before,root,[n for n in legacy if legacy[n]!=old['legacy_io_sources'].get(n)]))
 manifest.update(profiled_sources=profiled,legacy_profile_sources=legacy,observer_patch_sha256=file_digest(out/'stage198-observer.patch'),legacy_profile_patch_sha256=file_digest(out/'legacy-profile-198.patch'),observer_sha256=file_digest(Path(__file__).with_name('workflow_profile_198.f90')),numerical_candidate_stage=197,diagnosis_only=True)
 assert profiled['src/workflow_profile_195.f90']==legacy['src/workflow_profile_195.f90']==manifest['observer_sha256']
 (out/'stage198-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
 scope=json.loads((previous/'operation-instrumentation.json').read_text());scope['region_names'].update({str(k):v for k,v in REGIONS.items()});scope['counter_names'].update({str(k):v for k,v in COUNTERS.items()});scope['stage198_inventory']={'block':block_inventory,'legacy':legacy_inventory}
 scope['region_names']['0']='unclassified self work'
 scope['scope']='Full-step exclusive location ledger. Timed call-site envelopes include branch checks; legacy fused kernels are not individually matched arithmetic. Storage probes are inclusive systematic samples, not additive or hardware-byte measurements. Modes independently calibrated.'
 (out/'operation-instrumentation.json').write_text(json.dumps(scope,indent=2)+'\n')
 shutil.copy2(repo/'docs/parallel-blocks-cluster-198.md',out/'README.md')
 shutil.copy2(repo/'docs/parallel-blocks-stage198.md',out/'DIAGNOSIS.md')
 archive=out.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
 report={'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],'numerical_candidate_stage':197,'numerical_candidate_patch_sha256':manifest['candidate_patch_sha256'],'observer_patch_sha256':manifest['observer_patch_sha256'],'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'helpers_sha256':{n:file_digest(out/n) for n in HELPERS}}
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('repo','out','stage197-kit'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();prepare(a.repo,a.out,a.stage197_kit)
