#!/usr/bin/env python3
"""Freeze operation observers against committed Stage 195 and its exact legacy."""
import argparse,io,json,re,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import sources,file_digest,digest
from prepare_cluster_195 import patch
from instrument_operations_196 import instrument
from analyze_operations_196 import GROUPS,BOUNDARY_GROUPS
BASELINE='5dfe0db1'
HELPERS=('build_restriction_196.py','build_196.sbatch','submit_restriction_196.py','cluster_196.py','cluster_196.sbatch',
 'rank_metrics_190.py','rank_launch_190.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
 'profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs.py','analyze_workflow_195.py','analyze_operations_196.py')

def prepare(repo,out,previous):
 repo=repo.resolve();out=out.resolve();previous=previous.resolve()
 if out.exists() or out.is_relative_to(repo):raise ValueError('New kit outside repository required')
 base=subprocess.check_output(['git','-C',str(repo),'rev-parse',BASELINE],text=True).strip()
 old=json.loads((previous/'stage195-sources.json').read_text())
 entries=sources(repo,base)
 for n,h in old['baseline_sources'].items():
  if entries[n]!=h:raise ValueError('Committed numerical baseline differs from Stage 195: '+n)
 tracked=old['baseline_sources'];out.mkdir()
 for n in HELPERS:shutil.copy2(Path(__file__).with_name(n),out/n)
 for n in ('physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch'):
  shutil.copy2(previous/n,out/n)
 subprocess.run(['git','-C',str(repo),'bundle','create',str(out/'stage195-baseline.bundle'),'HEAD','^a0c75adb'],check=True)
 inventories={}
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp);raw=subprocess.check_output(['git','-C',str(repo),'archive',base])
  with tarfile.open(fileobj=io.BytesIO(raw)) as t:t.extractall(root,filter='data')
  before={n:(root/n).read_bytes() for n in tracked if not (root/n).is_symlink()}
  p=previous/'stage195-candidate.patch'
  if file_digest(p)!=old['candidate_patch_sha256']:raise ValueError('Previous candidate patch changed')
  subprocess.run(['git','apply',str(p)],cwd=root,check=True)
  inventories['block']=instrument(root,True)
  names=[*tracked,'src/workflow_profile_195.f90'];candidate={n:file_digest(root/n) for n in names}
  changed=[n for n in names if candidate[n]!=tracked.get(n)]
  (out/'stage196-candidate.patch').write_text(patch(before,root,changed))
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp)
  with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['legacy_io_sources'] if not (root/n).is_symlink()}
  p=previous/'legacy-profile-195.patch'
  if file_digest(p)!=old['legacy_profile_patch_sha256']:raise ValueError('Previous legacy patch changed')
  subprocess.run(['git','apply','--whitespace=nowarn',str(p)],cwd=root,check=True)
  inventories['legacy']=instrument(root,False)
  names=list(old['legacy_profile_sources']);profiled={n:file_digest(root/n) for n in names}
  changed=[n for n in names if profiled[n]!=old['legacy_io_sources'].get(n)]
  (out/'legacy-profile-196.patch').write_text(patch(before,root,changed))
 manifest={**{k:old[k] for k in ('legacy_revision','legacy_original_sources','legacy_io_sources','legacy_io_patch_sha256','legacy_archive_sha256')},
 'baseline_commit':base,'baseline_stage':195,'baseline_sources':tracked,'candidate_sources':candidate,
 'makefile_sha256':tracked['Makefile'],'candidate_makefile_sha256':candidate['Makefile'],
 'candidate_patch_sha256':file_digest(out/'stage196-candidate.patch'),
 'legacy_profile_sources':profiled,'legacy_profile_patch_sha256':file_digest(out/'legacy-profile-196.patch'),
 'observer_sha256':file_digest(Path(__file__).with_name('workflow_profile_196.f90'))}
 assert candidate['src/workflow_profile_195.f90']==profiled['src/workflow_profile_195.f90']==manifest['observer_sha256']
 (out/'stage196-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
 text=(repo/'src/parallel_block_profile.f90').read_text()
 region_names=re.findall(r"'([^']+)'",text[text.index('detail_name(DP_COUNT)'):text.index('logical, public, save :: detail_enabled')])
 counter_names=re.findall(r"'([^']+)'",text[text.index('detail_counter_name(DC_COUNT)'):text.index('detail_count(DC_COUNT)')])
 region_names=['unclassified dynamics self work',*region_names,'native mass begin exchange','native mass finish exchange','Domain full import preflight/install','Domain full import pack/zero/traversal','Domain full import MPI','Domain boundary import preflight','Domain boundary import zero/pack','Domain boundary import MPI','Domain boundary import install/local copies','Domain publication commit/residual']
 counter_names += ['native mass send buffer doubles','native mass receive buffer doubles','native mass local-copy entries','native mass posted requests','full import buffer doubles zeroed','full import patches installed','boundary import buffer doubles zeroed','boundary import local boundaries','boundary import remote boundaries']
 assert len(region_names)==78 and len(counter_names)==27
 text=(repo/'src/parallel_block_mpi.f90').read_text()
 a=text.index('BLOCK_PROFILE_PHASE_COUNT) = [character');b=text.index('real(dp), save :: block_profile_seconds',a)
 volume_names=re.findall(r'"([^"\n]+)"',text[a:b])
 assert len(volume_names)==56
 (out/'operation-instrumentation.json').write_text(json.dumps({'inventories':inventories,'groups':GROUPS,'boundary_groups':BOUNDARY_GROUPS,
 'region_names':dict(enumerate(region_names)),'counter_names':dict(enumerate(counter_names,1)),'block_volume_names':dict(enumerate(volume_names,1)),
 'counter_scope':'Logical work only; generic Domain sends/receives are doubles, not DRAM traffic. Native volume records cover instrumented families, not an exhaustive cross-implementation byte census. Boundary routine calls include start and finish scopes.'},indent=2)+'\n')
 shutil.copy2(repo/'docs/parallel-blocks-cluster-196.md',out/'README.md')
 archive=out.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
 report={'kit':str(out),'archive':str(archive),'sha256':file_digest(archive),'baseline_commit':base,
         'observer_sha256':manifest['observer_sha256'],'candidate_patch_sha256':manifest['candidate_patch_sha256'],
         'legacy_profile_patch_sha256':manifest['legacy_profile_patch_sha256'],'helpers_sha256':{n:file_digest(out/n) for n in HELPERS}}
 print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('repo','out','stage195-kit'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();prepare(a.repo,a.out,a.stage195_kit)
