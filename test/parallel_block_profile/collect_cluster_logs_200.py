#!/usr/bin/env python3
"""Collect one cluster job's review evidence, excluding checkpoints and binaries."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile


def collect(stage, job, home, results=None, tree=None):
    home = home.resolve()
    if results is None:
        matches = []
        for path in home.glob(f'stage{stage}*results-*/results.json'):
            record = json.loads(path.read_text())
            if str(record.get('slurm_environment', {}).get('SLURM_JOB_ID')) == job:
                matches.append(path.parent)
        if len(matches) != 1:
            raise ValueError(f'Expected one results directory for job {job}; found {matches}. Use --results PATH.')
        results = matches[0]
    results = results.resolve(strict=True)
    record_path = results / 'results.json'
    if record_path.exists():
        recorded_job = json.loads(record_path.read_text()).get('slurm_environment', {}).get('SLURM_JOB_ID')
        if recorded_job is not None and str(recorded_job) != job:
            raise ValueError('Results belong to a different Slurm job')
    tree = (tree or home / f'wt-stage{stage}').resolve()
    archive = home / f'stage{stage}-review-{job}.tar.gz'
    if archive.exists():
        raise FileExistsError(f'Review archive already exists: {archive}')
    files = {}
    for path in sorted(results.rglob('*')):
        if path.is_file() and not path.is_symlink() and path.suffix in {'.json', '.log', '.txt', '.resources', '.csv', '.md'}:
            files['results/' + str(path.relative_to(results))] = path
    roots = [tree, *sorted(tree.parent.glob(tree.name + '-legacy*'))]
    for root in roots:
        for path in sorted(root.glob('build-*.log')) + sorted(root.glob('*-builds.json')):
            if path.is_file() and not path.is_symlink():
                files[f'builds/{root.name}/{path.name}'] = path
    slurm = home / f'wavetrisk-{stage}-{job}.log'
    if slurm.is_file():
        files['slurm/' + slurm.name] = slurm
    if not files:
        raise ValueError('No review evidence found')
    inventory = {'stage': stage, 'job': job, 'results_directory': str(results),
                 'source_tree': str(tree), 'files': {},
                 'scope': 'Diagnostic evidence only; checkpoints and executables excluded. Collect after the job stops.'}
    with tarfile.open(archive, 'x:gz') as tar:
        for name, path in sorted(files.items()):
            data = path.read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mode = 0o644
            tar.addfile(info, io.BytesIO(data))
            inventory['files'][name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
        data = (json.dumps(inventory, indent=2) + '\n').encode()
        info = tarfile.TarInfo('collection.json')
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return archive


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', type=int, required=True)
    p.add_argument('--job', required=True)
    p.add_argument('--home', type=Path, default=Path.home())
    p.add_argument('--results', type=Path, help='Explicit directory for custom paths or setup failures')
    p.add_argument('--tree', type=Path, help='Custom build tree; default HOME/wt-stageSTAGE')
    a = p.parse_args()
    if a.stage < 1 or not a.job.isdigit(): p.error('Require a positive stage and numeric job ID')
    print(collect(a.stage, a.job, a.home, a.results, a.tree))


if __name__ == '__main__':
    main()
