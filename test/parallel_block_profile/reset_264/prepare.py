"""Create a clean self-contained reset anchor from recorded source/fixture seals."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import MAIN, J5_HASHES, save, sha, verify

HERE = Path(__file__).resolve().parent
HELPERS = ('analyze.py', 'capability_223.py', 'checkpoint_gate.py', 'compare_restart_checkpoints.py',
           'experiment.py', 'profile_protocol.py', 'protocol_229.py', 'run_local_230.py',
           'run_local_220.py', 'local_j4.py', 'validation_230.py', 'placement_230.py',
           'placement_200.py', 'rank_launch_223.py', 'rank_launch_190.py',
           'prepare_physics_dependencies.sh')


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() or p.is_symlink()}


def clean_physics(tree):
    root = tree / 'src/physics/simple_physics/phyparam'
    for name in ('obj', 'include', 'lib', 'xml', 'tmp', 'outputs'):
        path = root / name
        if path.exists():
            shutil.rmtree(path)
    for path in sorted((root / 'driver').rglob('*'), reverse=True):
        if path.is_dir() and path.name.endswith('.dSYM'):
            shutil.rmtree(path)
        elif path.is_file() and (path.suffix in ('.mod', '.o', '.a', '.so') or path.name in ('driver', 'driver.exe')):
            path.unlink()


def prepare(repo, oldkit, main_archive, fixture, destination):
    repo, oldkit, main_archive, fixture = [Path(p).resolve(strict=True) for p in (repo, oldkit, main_archive, fixture)]
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError('Fresh kit destination required')
    if sha(main_archive) != '23dbbd8abbcdf3d1fa8d39445de912d09e3806f9a381381a5135240db224f8e5':
        raise ValueError('Established direct-main archive differs')
    prior = json.loads((repo / 'docs/parallel-blocks-stage230-kit.json').read_text())
    for name, expected in json.loads((oldkit / 'helpers-sha256.json').read_text()).items():
        if sha(oldkit / name) != expected:
            raise ValueError('Original Stage230 kit file differs: ' + name)
    reference = subprocess.check_output(['git', 'rev-parse', '4ea64012'], cwd=repo, text=True).strip()
    main = subprocess.check_output(['git', 'rev-parse', 'main'], cwd=repo, text=True).strip()
    if main != MAIN:
        raise ValueError('Main advanced; audit before replacing the historical benchmark')
    destination.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix='rk264-source-') as tmp:
        root = Path(tmp)
        stage230 = root / 'stage230'
        stage230.mkdir()
        with tarfile.open(oldkit / 'stage230-source.tar.gz') as stream:
            stream.extractall(stage230, filter='data')
        subprocess.run(['git', 'apply', str(oldkit / 'stage230-candidate.patch')], cwd=stage230, check=True)
        candidates = json.loads((oldkit / 'stage230-sources.json').read_text())['candidate_sources']
        for name, expected in candidates.items():
            if sha(stage230 / name) != expected:
                raise ValueError('Frozen Stage230 source differs: ' + name)
        for path in stage230.iterdir():
            if path.name not in ('Makefile', 'src', 'test'):
                shutil.rmtree(path) if path.is_dir() else path.unlink()
        for path in (stage230 / 'test').iterdir():
            if path.name != 'climate':
                shutil.rmtree(path) if path.is_dir() else path.unlink()
        clean_physics(stage230)
        # Confirm exact committed numerical source, beyond the historical kit's seal.
        checked = []
        names = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', reference,
                                         'Makefile', 'src', 'test/climate'], cwd=repo, text=True).splitlines()
        for name in names:
            path = stage230 / name
            if (name == 'Makefile' or name.endswith(('.f90', '.F90', '.h'))):
                mode = subprocess.check_output(['git', 'ls-tree', reference, name], cwd=repo, text=True).split()[0]
                body = subprocess.check_output(['git', 'show', reference + ':' + name], cwd=repo)
                actual = str(path.readlink()).encode() if mode == '120000' else path.read_bytes()
                if actual != body:
                    raise ValueError('Stage230 commit/source differs: ' + name)
                checked.append(name)
        with tarfile.open(main_archive) as stream:
            stream.extractall(root, filter='data')
        main_tree = root / 'main'
        (root / 'rk249-main-legacy').rename(main_tree)
        for name, expected in json.loads((main_tree / 'files-sha256.json').read_text()).items():
            with (main_tree / name).open('rb') as stream:
                if hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
                    raise ValueError('Original main archive source differs: ' + name)
        main_identity = json.loads((main_tree / 'MAIN_IDENTITY.json').read_text())
        if main_identity['main_commit'] != MAIN or main_identity['changed_source'] != ['src/checkpoint.f90']:
            raise ValueError('Unexpected legacy solver adjustment')
        (main_tree / 'files-sha256.json').unlink()
        shutil.rmtree(main_tree / 'src/physics/simple_physics')
        shutil.copytree(stage230 / 'src/physics/simple_physics', main_tree / 'src/physics/simple_physics', symlinks=True)
        for path in main_tree.iterdir():
            if path.name not in ('Makefile', 'src', 'test', 'MAIN_IDENTITY.json'):
                shutil.rmtree(path) if path.is_dir() else path.unlink()
        manifests = {kind: inventory(root / kind) for kind in ('stage230', 'main')}
        save(destination / 'source-files.json', manifests)
        with tarfile.open(destination / 'sources.tar.gz', 'w:gz', format=tarfile.PAX_FORMAT) as stream:
            for kind in ('stage230', 'main'):
                stream.add(root / kind, arcname=kind)
    for path in HERE.iterdir():
        if path.is_file() and path.suffix in ('.py', '.sh', '.sbatch', '.md'):
            shutil.copy2(path, destination / path.name)
    helper_root = repo / 'test/parallel_block_profile'
    for name in HELPERS:
        shutil.copy2(helper_root / name, destination / name)
    (destination / 'grids').mkdir()
    for name in ('HR95JT_WT_003', 'HR95JT_WT_004'):
        shutil.copy2(fixture / 'grids' / name, destination / 'grids' / name)
    if sha(destination / 'grids/HR95JT_WT_004') != J5_HASHES['grids/HR95JT_WT_004']:
        raise ValueError('Established J5 coarse grid differs')
    local = destination / 'fixture-J4'
    local.mkdir()
    for name in ('simple.in', 'localJ4J6_checkpoint_0002.bin.zst'):
        shutil.copy2(fixture / name, local / name)
    shutil.copytree(destination / 'grids', local / 'grids')
    if sha(local / 'localJ4J6_checkpoint_0002.bin.zst') != '285df4e715661c2a74836f5781f9cb099c777f9fff4e456a0b3593061dddf34c':
        raise ValueError('Established local J4 seed differs')
    save(destination / 'IDENTITY.json', dict(stage=264, stage230_commit=reference, main_commit=MAIN,
         main_identity=main_identity, stage230_original_kit_sha256=prior['sha256'],
         main_original_archive_sha256=sha(main_archive), committed_inputs_checked=checked,
         j5_fixture_hashes=J5_HASHES, scope='Exact Stage230 versus legacy main; no owner-graph migration or new solver optimization.'))
    save(destination / 'files-sha256.json', inventory(destination))
    verify(destination)
    print('Prepared: ' + str(destination))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'oldkit', 'main_archive', 'fixture', 'destination'):
        parser.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    prepare(**vars(parser.parse_args()))
