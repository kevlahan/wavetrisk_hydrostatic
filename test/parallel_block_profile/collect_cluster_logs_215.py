#!/usr/bin/env python3
"""Collect review evidence atomically; verify the inventory before publication."""
import argparse,hashlib,io,json,os,tarfile,tempfile
from pathlib import Path


def verify_archive(archive):
    with tarfile.open(archive,'r:gz') as tar:
        members=tar.getmembers();names=[m.name for m in members]
        if len(set(names))!=len(names) or any(not m.isfile() for m in members):raise ValueError('Duplicate/non-file archive members')
        if 'collection.json' not in names or 'results/results.json' not in names:raise ValueError('Incomplete review: missing manifest or main results')
        inventory=json.load(tar.extractfile('collection.json'))
        if set(names)!={'collection.json',*inventory['files']}:raise ValueError('Inventory/member set differs')
        for name,expected in inventory['files'].items():
            data=tar.extractfile(name).read()
            if len(data)!=expected['bytes'] or hashlib.sha256(data).hexdigest()!=expected['sha256']:raise ValueError('Corrupt member: '+name)
        results=json.load(tar.extractfile('results/results.json'))
        job=results.get('slurm_environment',{}).get('SLURM_JOB_ID')
        if job is not None and str(job)!=inventory['job']:raise ValueError('Result job differs from inventory')
    return inventory


def collect(stage,job,home,results=None,tree=None):
    home=home.resolve()
    if results is None:
        matches=[p.parent for p in home.glob(f'stage{stage}*results-*/results.json')
                 if str(json.loads(p.read_text()).get('slurm_environment',{}).get('SLURM_JOB_ID'))==job]
        if len(matches)!=1:raise ValueError('Specify --results: expected one result directory')
        results=matches[0]
    results=results.resolve(strict=True);record=json.loads((results/'results.json').read_text())
    recorded_job=record.get('slurm_environment',{}).get('SLURM_JOB_ID')
    if recorded_job is not None and str(recorded_job)!=job:raise ValueError('Results belong to another job')
    tree=(tree or home/f'wt-stage{stage}').resolve();archive=home/f'stage{stage}-review-{job}.tar.gz'
    if archive.exists():raise FileExistsError('Archive already exists; use --verify: '+str(archive))
    files={}
    for path in sorted(results.rglob('*')):
        if path.is_file() and not path.is_symlink() and (path.suffix in {'.json','.log','.txt','.resources','.csv','.md','.data'} or path.name=='test_log' or path.name.endswith('.data.gz')):
            files['results/'+str(path.relative_to(results))]=path
    for root in [tree,*sorted(tree.parent.glob(tree.name+'-legacy*'))]:
        for path in sorted(root.glob('build-*.log'))+sorted(root.glob('*-builds.json')):
            if path.is_file() and not path.is_symlink():files[f'builds/{root.name}/{path.name}']=path
    slurm=home/f'wavetrisk-{stage}-{job}.log'
    if slurm.is_file():files['slurm/'+slurm.name]=slurm
    inventory=dict(stage=stage,job=job,results_directory=str(results),source_tree=str(tree),files={},
                   scope='Complete review inventory; binaries/checkpoints excluded. Final filename is published only after verification.')
    descriptor,name=tempfile.mkstemp(prefix=f'.stage{stage}-review-{job}.partial-',suffix='.tar.gz',dir=home)
    temporary=Path(name)
    try:
        with os.fdopen(descriptor,'wb') as stream:
            with tarfile.open(fileobj=stream,mode='w:gz') as tar:
                for name,path in sorted(files.items()):
                    data=path.read_bytes();info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644
                    tar.addfile(info,io.BytesIO(data));inventory['files'][name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
                data=(json.dumps(inventory,indent=2)+'\n').encode();info=tarfile.TarInfo('collection.json');info.size=len(data)
                tar.addfile(info,io.BytesIO(data))
            stream.flush();os.fsync(stream.fileno())
        verify_archive(temporary)
        # Atomic, same-filesystem publication; unlike rename, never overwrites.
        os.link(temporary,archive)
    finally:
        temporary.unlink(missing_ok=True)
    return archive


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',type=Path)
    p.add_argument('--stage',type=int);p.add_argument('--job');p.add_argument('--home',type=Path,default=Path.home())
    p.add_argument('--results',type=Path);p.add_argument('--tree',type=Path);a=p.parse_args()
    if a.verify:
        m=verify_archive(a.verify);print(f"PASS: {len(m['files'])} file hashes verified in {a.verify}");return
    if a.stage is None or a.stage<1 or not a.job or not a.job.isdigit():p.error('Require --stage and numeric --job')
    print(collect(a.stage,a.job,a.home,a.results,a.tree))
if __name__=='__main__':main()
