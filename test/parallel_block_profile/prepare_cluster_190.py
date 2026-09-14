#!/usr/bin/env python3
"""Freeze the Stage 190 numerical patch and pinned legacy source into a new kit."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

BASELINE='a0c75adb'
LEGACY='b82647b36eb825aed0c0aae4e92adcf4ab46d608'
HELPERS=('build_restriction_190.py','build_190.sbatch','submit_restriction_190.py','cluster_190.py','cluster_190.sbatch',
         'rank_metrics_190.py','rank_launch_190.py','legacy_checkpoint_io_190.patch','analyze.py','checkpoint_gate.py',
         'compare_restart_checkpoints.py','profile_protocol.py','experiment.py')


def digest(data): return hashlib.sha256(data).hexdigest()


def sources(repo,ref):
    raw=subprocess.check_output(['git','-C',str(repo),'archive','--format=tar',ref])
    with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
        return {m.name:digest(m.linkname.encode() if m.issym() else tar.extractfile(m).read())
                for m in tar.getmembers() if m.isfile() or m.issym()}


def file_digest(path):
    return digest(str(path.readlink()).encode() if path.is_symlink() else path.read_bytes())


def prepare(repo,out,physics_include,dependency_helper):
    repo=repo.resolve();out=out.resolve()
    if out.exists() or out.is_relative_to(repo): raise ValueError('Use a new output outside the repository')
    baseline=subprocess.check_output(['git','-C',str(repo),'rev-parse',BASELINE],text=True).strip()
    original=sources(repo,LEGACY)
    base=sources(repo,baseline)
    numerical={n:h for n,h in base.items() if n.startswith('src/') and n.lower().endswith(('.f90','.f','.c','.h'))}
    candidate={n:file_digest(repo/n) for n in numerical}
    changed=[n for n in candidate if candidate[n]!=numerical[n]]
    if changed!=['src/parallel_block_mpi.f90']: raise ValueError('Unexpected numerical patch scope: '+repr(changed))
    head=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    if head!=baseline: raise ValueError('Freeze this kit from the Stage 189 baseline checkout plus candidate edits')
    out.mkdir(parents=True)
    subprocess.run(['git','-C',str(repo),'bundle','create',str(out/'stage189-baseline.bundle'),
                    'HEAD','^61509429'],check=True)
    for name in HELPERS: shutil.copy2(Path(__file__).with_name(name),out/name)
    shutil.copy2(physics_include,out/'physics-Makefile.inc')
    shutil.copy2(dependency_helper,out/'prepare_physics_dependencies.sh')
    shutil.copy2(repo/'docs/parallel-blocks-cluster-190.md',out/'README.md')
    patch=subprocess.check_output(['git','-C',str(repo),'diff',baseline,'--','src/parallel_block_mpi.f90'])
    (out/'stage190-candidate.patch').write_bytes(patch)
    archive=out/'legacy-source.tar.gz'
    subprocess.run(['git','-C',str(repo),'archive','--format=tar.gz','--output='+str(archive),LEGACY],check=True)
    # Verify application of the I/O-only delta before publishing its after hashes.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        subprocess.run(['tar','-xf',str(archive),'-C',str(root)],check=True)
        subprocess.run(['git','apply',str(out/'legacy_checkpoint_io_190.patch')],cwd=root,check=True)
        after={n:file_digest(root/n) for n in original}
    if [n for n in after if after[n]!=original[n]]!=['src/checkpoint.f90']:
        raise ValueError('Legacy numerical sources unexpectedly changed')
    manifest={'baseline_commit':baseline,'baseline_sources':numerical,'candidate_sources':candidate,
              'makefile_sha256':base['Makefile'],'candidate_patch_sha256':digest(patch),
              'legacy_revision':LEGACY,'legacy_original_sources':original,'legacy_io_sources':after,
              'legacy_archive_sha256':digest(archive.read_bytes()),
              'baseline_bundle_sha256':digest((out/'stage189-baseline.bundle').read_bytes()),
              'legacy_io_patch_sha256':digest((out/'legacy_checkpoint_io_190.patch').read_bytes())}
    (out/'stage190-sources.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive=out.with_suffix('.tar.gz')
    with tarfile.open(archive,'w:gz') as tar: tar.add(out,arcname=out.name)
    print(json.dumps({'kit':str(out),'archive':str(archive),'sha256':digest(archive.read_bytes())},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--physics-include',type=Path,required=True);p.add_argument('--dependency-helper',type=Path,required=True)
    a=p.parse_args();prepare(a.repo,a.out,a.physics_include,a.dependency_helper)
