#!/usr/bin/env python3
"""Freeze diagnostic-only Stage199 observers around validated Stage197 numerics."""
import argparse,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import file_digest
from prepare_cluster_195 import patch
from instrument_operations_199 import upgrade,REGIONS,COUNTERS
HELPERS=('build_restriction_199.py','build_199.sbatch','submit_restriction_199.py','cluster_199.py','cluster_199.sbatch','analyze_operations_199.py',
 'rank_metrics_190.py','rank_launch_190.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py','profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs_199.py','analyze_workflow_195.py','analyze_batches_199.py')

def prepare(repo,out,previous):
 repo=repo.resolve();out=out.resolve();previous=previous.resolve()
 if out.exists() or out.is_relative_to(repo):raise ValueError('Fresh kit outside repository required')
 old=json.loads((previous/'stage198-sources.json').read_text())
 for n,h in old['candidate_sources'].items():
  if file_digest(repo/n)!=h:raise ValueError('Numerical source changed since validated Stage197: '+n)
 out.mkdir();manifest=dict(old)
 for name in HELPERS:shutil.copy2(Path(__file__).with_name(name),out/name)
 for name in ('stage196-baseline.bundle','physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch'):shutil.copy2(previous/name,out/name)
 shutil.copy2(previous/'stage198-candidate.patch',out/'stage199-candidate.patch')
 if file_digest(out/'stage199-candidate.patch')!=old['candidate_patch_sha256']:raise ValueError('Candidate numerics changed')
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(fileobj=io.BytesIO(subprocess.check_output(['git','-C',str(repo),'archive',old['baseline_commit']]))) as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'stage199-candidate.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['candidate_sources'] if not (root/n).is_symlink()}
  if file_digest(previous/'stage198-observer.patch')!=old['observer_patch_sha256']:raise ValueError('Previous block observer changed')
  subprocess.run(['git','apply',str(previous/'stage198-observer.patch')],cwd=root,check=True)
  block_inventory=upgrade(root,True)
  profiled={n:file_digest(root/n) for n in [n for n in old['profiled_sources'] if n!='src/workflow_clock_198.c']+['src/workflow_clock_199.c']}
  (out/'stage199-observer.patch').write_text(patch(before,root,[n for n in profiled if profiled[n]!=old['candidate_sources'].get(n)]))
 with tempfile.TemporaryDirectory() as d:
  root=Path(d)
  with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['legacy_io_sources'] if not (root/n).is_symlink()}
  if file_digest(previous/'legacy-profile-198.patch')!=old['legacy_profile_patch_sha256']:raise ValueError('Previous legacy observer changed')
  subprocess.run(['git','apply','--whitespace=nowarn',str(previous/'legacy-profile-198.patch')],cwd=root,check=True)
  legacy_inventory=upgrade(root,False)
  legacy={n:file_digest(root/n) for n in [n for n in old['legacy_profile_sources'] if n!='src/workflow_clock_198.c']+['src/workflow_clock_199.c']}
  (out/'legacy-profile-199.patch').write_text(patch(before,root,[n for n in legacy if legacy[n]!=old['legacy_io_sources'].get(n)]))
 manifest.update(profiled_sources=profiled,legacy_profile_sources=legacy,observer_patch_sha256=file_digest(out/'stage199-observer.patch'),legacy_profile_patch_sha256=file_digest(out/'legacy-profile-199.patch'),observer_sha256=file_digest(Path(__file__).with_name('workflow_profile_199.f90')),numerical_candidate_stage=197,diagnosis_only=True)
 assert profiled['src/workflow_profile_195.f90']==legacy['src/workflow_profile_195.f90']==manifest['observer_sha256']
 (out/'stage199-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
 scope=json.loads((previous/'operation-instrumentation.json').read_text());scope['region_names'].update({str(k):v for k,v in REGIONS.items()});scope['counter_names'].update({str(k):v for k,v in COUNTERS.items()});scope['stage199_inventory']={'block':block_inventory,'legacy':legacy_inventory}
 scope['batch_schema']={'version':1,'columns':['step','context','sequence','parent','operation','layer','level','position','closed','work','range_requests','logical_slots','full_ranges','read_visits','address_resolutions','write_visits','column_requests','local_start','local_end','mpi_subset_seconds'],
  'all_or_unknown_layer':0,'all_or_unknown_level':[0,-1],
  'work_units':{'105':'local Domain level patches','106':'local Domain level patches','107':'local native mass Domains','108':'send-buffer doubles plus local-copy entries','109':'receive-buffer doubles','117':'catalog descriptors','119':'source block leaf-node ownership weight'},
  'unlisted_work':'zero means not recorded, not zero computational work',
  'mpi_subset_operation_ids':[13,18,19,34,37,38,39,40,60,72,75],
  'omitted_event_ids':[102,103,128],
  'scope':'Inclusive intervals and counts. Patch producer/allocation clocks remain in exclusive ledger. Do not sum nested batches or compare independent rank-clock origins.'}
 scope['region_names']['0']='unclassified self work'
 scope['scope']='Full-step exclusive location ledger. Timed call-site envelopes include branch checks; legacy fused kernels are not individually matched arithmetic. Ordered batch records are inclusive; storage access census has no per-access clocks. Rank clock origins are independent. Modes independently calibrated.'
 (out/'operation-instrumentation.json').write_text(json.dumps(scope,indent=2)+'\n')
 shutil.copy2(repo/'docs/parallel-blocks-cluster-199.md',out/'README.md')
 shutil.copy2(repo/'docs/parallel-blocks-stage199.md',out/'DIAGNOSIS.md')
 archive=out.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
 report={'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':manifest['baseline_commit'],'numerical_candidate_stage':197,'numerical_candidate_patch_sha256':manifest['candidate_patch_sha256'],'observer_patch_sha256':manifest['observer_patch_sha256'],'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'helpers_sha256':{n:file_digest(out/n) for n in HELPERS}}
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('repo','out','stage198-kit'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();prepare(a.repo,a.out,a.stage198_kit)
