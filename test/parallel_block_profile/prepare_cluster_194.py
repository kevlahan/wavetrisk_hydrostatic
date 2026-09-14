#!/usr/bin/env python3
"""Freeze the Stage 194 numerical patch and pinned legacy source into a new kit."""
import argparse
import hashlib
import difflib
from legacy_instrumentation_191 import instrument
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

BASELINE='30c83114'
LEGACY='b82647b36eb825aed0c0aae4e92adcf4ab46d608'
HELPERS=('build_restriction_194.py','build_194.sbatch','submit_restriction_194.py','cluster_194.py','cluster_194.sbatch',
         'rank_metrics_190.py','rank_launch_190.py','legacy_checkpoint_io_190.patch','analyze.py','checkpoint_gate.py',
         'compare_restart_checkpoints.py','profile_protocol.py','experiment.py','legacy_metrics_191.py','collect_cluster_logs.py','analyze_inverse_192.py','analyze_scalar_inverse_193.py','analyze_inner_inverse_194.py')


def digest(data): return hashlib.sha256(data).hexdigest()


def sources(repo,ref):
    raw=subprocess.check_output(['git','-C',str(repo),'archive','--format=tar',ref])
    with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
        return {m.name:digest(m.linkname.encode() if m.issym() else tar.extractfile(m).read())
                for m in tar.getmembers() if m.isfile() or m.issym()}


def file_digest(path):
    return digest(str(path.readlink()).encode() if path.is_symlink() else path.read_bytes())


def prepare(repo,out,physics_include,dependency_helper,stage193_kit):
    repo=repo.resolve();out=out.resolve()
    if out.exists() or out.is_relative_to(repo): raise ValueError('Use a new output outside the repository')
    baseline=subprocess.check_output(['git','-C',str(repo),'rev-parse',BASELINE],text=True).strip()
    original=sources(repo,LEGACY)
    base=sources(repo,baseline)
    bootstrap={n:h for n,h in base.items() if n.startswith('src/') and n.lower().endswith(('.f90','.f','.c','.h'))}
    previous=json.loads((stage193_kit/'stage193-sources.json').read_text())
    baseline_patch=(stage193_kit/'stage193-candidate.patch').read_bytes()
    if digest(baseline_patch)!='e00d6dbb6d010c34c1f0fcc56b0d127d80551cda7d4174e0078cc7527d0cfc33':
        raise ValueError('Require the validated Stage 193 kit')
    if previous['baseline_sources']!=bootstrap: raise ValueError('Stage 193 bootstrap differs')
    numerical=previous['candidate_sources']
    candidate={n:file_digest(repo/n) for n in numerical}
    changed=[n for n in candidate if candidate[n]!=numerical[n]]
    if set(changed)!={'src/parallel_block_mpi.f90','src/parallel_block_profile.f90'}: raise ValueError('Unexpected numerical patch scope: '+repr(changed))
    if subprocess.run(['git','-C',str(repo),'merge-base','--is-ancestor',baseline,'HEAD']).returncode:
        raise ValueError('Stage 191 baseline must be an ancestor of HEAD for the bundle')
    out.mkdir(parents=True)
    subprocess.run(['git','-C',str(repo),'bundle','create',str(out/'stage191-baseline.bundle'),
                    'HEAD','^a0c75adb'],check=True)
    for name in HELPERS: shutil.copy2(Path(__file__).with_name(name),out/name)
    shutil.copy2(physics_include,out/'physics-Makefile.inc')
    shutil.copy2(dependency_helper,out/'prepare_physics_dependencies.sh')
    shutil.copy2(repo/'docs/parallel-blocks-cluster-194.md',out/'README.md')
    # Reconstruct the previously tested uncommitted Stage 193 state from its
    # frozen patch, so neither the user nor a remote branch must recreate it.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        reference=Path(tmp)
        archive_bytes=subprocess.check_output(['git','-C',str(repo),'archive','--format=tar',baseline])
        with tarfile.open(fileobj=io.BytesIO(archive_bytes)) as tar: tar.extractall(reference,filter='data')
        subprocess.run(['git','apply','-'],input=baseline_patch,cwd=reference,check=True)
        for name,expected in numerical.items():
            if file_digest(reference/name)!=expected: raise ValueError('Stage 193 source identity differs: '+name)
        old_profile=(reference/'src/parallel_block_profile.f90').read_bytes()
        patch=''.join(''.join(difflib.unified_diff((reference/name).read_text().splitlines(True),
            (repo/name).read_text().splitlines(True),fromfile='a/'+name,tofile='b/'+name)) for name in changed).encode()
    (out/'stage193-baseline.patch').write_bytes(baseline_patch)
    (out/'stage194-candidate.patch').write_bytes(patch)
    archive=out/'legacy-source.tar.gz'
    subprocess.run(['git','-C',str(repo),'archive','--format=tar.gz','--output='+str(archive),LEGACY],check=True)
    # Verify application of the I/O-only delta before publishing its after hashes.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        subprocess.run(['tar','-xf',str(archive),'-C',str(root)],check=True)
        subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
        after={n:file_digest(root/n) for n in original}
        before_text={n:(root/n).read_bytes() for n in original if not (root/n).is_symlink()}
        instrument(repo,root)
        # Preserve the legacy 65-region schema; new split timers are block-only.
        (root/'src/parallel_block_profile.f90').write_bytes(old_profile)
        profile_names=[*original,'src/parallel_block_profile.f90','src/legacy_profile_io.f90']
        profiled={n:file_digest(root/n) for n in profile_names}
        timer_patch=[]
        for name in profile_names:
            if (root/name).is_symlink(): continue
            data=(root/name).read_bytes()
            before=before_text.get(name,b'')
            if data==before: continue
            timer_patch.extend(difflib.unified_diff(before.decode().splitlines(True),data.decode().splitlines(True),
                fromfile='a/'+name if name in before_text else '/dev/null',tofile='b/'+name))
        (out/'legacy-profile-194.patch').write_text(''.join(timer_patch))
    if [n for n in after if after[n]!=original[n]]!=['src/checkpoint.f90']:
        raise ValueError('Legacy numerical sources unexpectedly changed')
    manifest={'baseline_commit':baseline,'bootstrap_sources':bootstrap,'baseline_stage':193,
              'baseline_patch_sha256':digest(baseline_patch),'baseline_sources':numerical,'candidate_sources':candidate,
              'makefile_sha256':base['Makefile'],'candidate_patch_sha256':digest(patch),
              'legacy_revision':LEGACY,'legacy_original_sources':original,'legacy_io_sources':after,
              'legacy_archive_sha256':digest(archive.read_bytes()),
              'baseline_bundle_sha256':digest((out/'stage191-baseline.bundle').read_bytes()),
              'legacy_io_patch_sha256':digest((out/'legacy_checkpoint_io_190.patch').read_bytes()),
              'legacy_profile_sources':profiled,'legacy_profile_patch_sha256':digest((out/'legacy-profile-194.patch').read_bytes())}
    (out/'stage194-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as tar: tar.add(out,arcname=out.name)
    print(json.dumps({'kit':str(out),'archive':str(archive),'sha256':digest(archive.read_bytes())},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--physics-include',type=Path,required=True);p.add_argument('--dependency-helper',type=Path,required=True)
    p.add_argument('--stage193-kit',type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.physics_include,a.dependency_helper,a.stage193_kit)
