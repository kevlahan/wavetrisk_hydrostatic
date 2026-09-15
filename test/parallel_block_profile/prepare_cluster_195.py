#!/usr/bin/env python3
"""Freeze paired diagnostic-only changes on committed Stage 194 and fixed legacy."""
import argparse,difflib,io,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from prepare_cluster_194 import sources,file_digest,digest
from instrument_workflow_195 import instrument
BASELINE='adb2e91e'
HELPERS=('build_restriction_195.py','build_195.sbatch','submit_restriction_195.py','cluster_195.py','cluster_195.sbatch',
 'rank_metrics_190.py','rank_launch_190.py','analyze.py','checkpoint_gate.py','compare_restart_checkpoints.py',
 'profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs.py','analyze_workflow_195.py')

def patch(before,root,names):
 return ''.join(''.join(difflib.unified_diff(before.get(n,b'').decode().splitlines(True),(root/n).read_text().splitlines(True),
  fromfile='a/'+n if n in before else '/dev/null',tofile='b/'+n)) for n in names if not (root/n).is_symlink())

def prepare(repo,out,previous):
 repo=repo.resolve();out=out.resolve();previous=previous.resolve()
 if out.exists() or out.is_relative_to(repo):raise ValueError('New kit outside repository required')
 base=subprocess.check_output(['git','-C',str(repo),'rev-parse',BASELINE],text=True).strip()
 old=json.loads((previous/'stage194-sources.json').read_text())
 entries=sources(repo,base)
 for n,h in old['candidate_sources'].items():
  if entries[n]!=h:raise ValueError('Committed Stage 194 differs from validated kit: '+n)
 tracked={n:h for n,h in entries.items() if (n.startswith('src/') and n.lower().endswith(('.f90','.f','.c','.h'))) or n in ('Makefile','test/climate/climate.f90')}
 out.mkdir()
 for n in HELPERS:shutil.copy2(Path(__file__).with_name(n),out/n)
 for n in ('physics-Makefile.inc','prepare_physics_dependencies.sh','legacy-source.tar.gz','legacy_checkpoint_io_190.patch'):
  shutil.copy2(previous/n,out/n)
 subprocess.run(['git','-C',str(repo),'bundle','create',str(out/'stage194-baseline.bundle'),'HEAD','^a0c75adb'],check=True)
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp)
  raw=subprocess.check_output(['git','-C',str(repo),'archive',base])
  with tarfile.open(fileobj=io.BytesIO(raw)) as t:t.extractall(root,filter='data')
  before={n:(root/n).read_bytes() for n in ('Makefile','test/climate/climate.f90','src/main.f90','src/adapt.f90','src/time_integr.f90')}
  instrument(root,True)
  changed=[*before,'src/workflow_profile_195.f90']
  text=patch(before,root,changed);(out/'stage195-candidate.patch').write_text(text)
  candidate={n:file_digest(root/n) for n in [*tracked,'src/workflow_profile_195.f90']}
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp)
  with tarfile.open(out/'legacy-source.tar.gz') as t:t.extractall(root,filter='data')
  subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
  before={n:(root/n).read_bytes() for n in old['legacy_io_sources'] if not (root/n).is_symlink()}
  # Existing legacy source inventory contains text files only except excluded figures.
  subprocess.run(['git','apply',str(previous/'legacy-profile-194.patch')],cwd=root,check=True)
  instrument(root,False)
  names=[*old['legacy_io_sources'],'src/parallel_block_profile.f90','src/legacy_profile_io.f90','src/workflow_profile_195.f90']
  profiled={n:file_digest(root/n) for n in names}
  changed=[n for n in names if profiled[n]!=old['legacy_io_sources'].get(n)]
  text=patch(before,root,changed);(out/'legacy-profile-195.patch').write_text(text)
 manifest={'baseline_commit':base,'baseline_stage':194,'baseline_sources':tracked,'candidate_sources':candidate,
 'makefile_sha256':tracked['Makefile'],'candidate_makefile_sha256':candidate['Makefile'],
 'candidate_patch_sha256':digest((out/'stage195-candidate.patch').read_bytes()),
 'legacy_revision':old['legacy_revision'],'legacy_original_sources':old['legacy_original_sources'],
 'legacy_io_sources':old['legacy_io_sources'],'legacy_io_patch_sha256':old['legacy_io_patch_sha256'],
 'legacy_archive_sha256':old['legacy_archive_sha256'],'legacy_profile_sources':profiled,
 'legacy_profile_patch_sha256':digest((out/'legacy-profile-195.patch').read_bytes()),
 'observer_sha256':file_digest(Path(__file__).with_name('workflow_profile_195.f90'))}
 (out/'stage195-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
 shutil.copy2(repo/'docs/parallel-blocks-cluster-195.md',out/'README.md')
 archive=out.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
 print(json.dumps({'kit':str(out),'archive':str(archive),'sha256':digest(archive.read_bytes())},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('repo','out','stage194-kit'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();prepare(a.repo,a.out,a.stage194_kit)
