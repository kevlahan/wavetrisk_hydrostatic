"""Prepare immutable accepted266/main controls and an explicit connected RHS overlay."""
import argparse,hashlib,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path
from common import BASELINE,MAIN,save,sha,verify,source_scope,verify_packaged_grid
HERE=Path(__file__).resolve().parent

def inventory(root):
    return {str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*'))
            if (p.is_file() or p.is_symlink()) and '__pycache__' not in p.parts and p.name!='files-sha256.json'}

def prepare(repo,anchor,destination):
    repo,anchor=map(lambda p:Path(p).resolve(strict=True),(repo,anchor))
    verify(anchor);identity=json.loads((anchor/'IDENTITY.json').read_text())
    if identity['stage']!=266 or identity['revision']!='deferred-installation' or identity['main_commit']!=MAIN:
        raise ValueError('Require sealed accepted266 anchor')
    destination=Path(destination).resolve()
    if destination.exists():raise ValueError('Fresh kit directory required')
    overlay=['Makefile','src/multi_level.f90','src/time_integr.f90','src/parallel_block.f90',
             'src/parallel_block_mass.f90','src/parallel_block_mpi.f90','src/parallel_block_primitive_walk.f90']
    overlay+=sorted(str(p.relative_to(repo)) for p in (repo/'src').glob('parallel_block_rhs*.f90'))
    overlay=sorted(overlay)
    destination.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk278-source-') as tmp:
        root=Path(tmp)
        with tarfile.open(anchor/'sources.tar.gz') as stream:stream.extractall(root,filter='data')
        prior=json.loads((anchor/'source-files.json').read_text())
        for kind in ('candidate','main'):
            for name,expected in prior[kind].items():
                if sha(root/kind/name)!=expected:raise ValueError('Anchor source differs: '+kind+'/'+name)
        for name in identity['committed_inputs_checked']:
            actual=subprocess.check_output(['git','show',BASELINE+':'+name],cwd=repo)
            if hashlib.sha256(actual).hexdigest()!=prior['candidate'][name]:raise ValueError('Accepted commit differs: '+name)
        shutil.rmtree(root/'stage265t');(root/'candidate').rename(root/'stage266')
        shutil.copytree(root/'stage266',root/'candidate',symlinks=True)
        for name in overlay:shutil.copyfile(repo/name,root/'candidate'/name)
        manifests={kind:inventory(root/kind) for kind in ('stage266','main','candidate')}
        changed=sorted(n for n in set(manifests['stage266'])|set(manifests['candidate'])
                       if manifests['stage266'].get(n)!=manifests['candidate'].get(n))
        if changed!=overlay:raise ValueError('Explicit overlay contains unchanged/missing files')
        save(destination/'source-files.json',manifests)
        with tarfile.open(destination/'sources.tar.gz','w:gz',format=tarfile.PAX_FORMAT) as stream:
            for kind in manifests:stream.add(root/kind,arcname=kind)
    for p in anchor.iterdir():
        if p.is_file() and p.suffix in ('.py','.sh') and p.name!='repeat.sh':shutil.copy2(p,destination/p.name)
    # Use the established launch/control/build/analysis graph; development-only
    # component drivers are intentionally outside the cluster execution graph.
    names=('common.py','prepare.py','run.py','build.py','control.py','review.py','evidence.py','accounting.py',
           'launch.sh','status.sh','collect.sh','timing.sbatch','placement_278.py','contract.py','campaign_gate.py',
           'seal.py','campaign-spec.json','README.md','test_campaign.py')
    for name in names:
        if (HERE/name).is_file():shutil.copy2(HERE/name,destination/name)
    for name in ('grids','fixture-J4'):shutil.copytree(anchor/name,destination/name)
    identity.update(stage=278,revision='connected-rhs',analysis_schema_version=4,stage266_commit=BASELINE,
                    changed_source=overlay,anchor_manifest_sha256=sha(anchor/'files-sha256.json'),
                    scope='Resident connected RHS and real RK, paired retirement, bounded inverse/boundary adapters; partition unchanged.')
    identity.pop('stage265t_commit',None)
    identity.update(candidate_source_sha256={n:manifests['candidate'][n] for n in overlay},
                    baseline_source_sha256={n:manifests['stage266'].get(n) for n in overlay})
    save(destination/'overlay-files.json',overlay);save(destination/'IDENTITY.json',identity)
    source_scope(destination);verify_packaged_grid(destination)
    save(destination/'files-sha256.json',inventory(destination));verify(destination)
    print('Prepared unsealed kit: '+str(destination))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('repo','anchor','destination'):p.add_argument('--'+name,type=Path,required=True)
    prepare(**vars(p.parse_args()))
