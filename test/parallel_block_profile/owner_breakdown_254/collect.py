"""Create or verify a Stage254 review archive after the Slurm job ends."""
import argparse
import hashlib
import io
import json
import shutil
import subprocess
from pathlib import Path
import tarfile


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify(archive):
    with tarfile.open(archive) as tar:
        members = {member.name: member for member in tar if member.isfile()}
        manifests = [name for name in members if name.endswith('/REVIEW_SHA256.json')]
        if len(manifests) != 1:
            raise ValueError('Missing or duplicate review manifest')
        manifest = json.load(tar.extractfile(members[manifests[0]]))
        prefix = manifests[0].removesuffix('REVIEW_SHA256.json')
        if set(members) != {prefix+name for name in manifest} | set(manifests):
            raise ValueError('Review archive file list differs')
        for name, expected in manifest.items():
            if hashlib.sha256(tar.extractfile(members[prefix+name]).read()).hexdigest() != expected:
                raise ValueError('Review archive hash differs: '+name)
    return dict(archive=str(archive.resolve()), sha256=digest(archive),
                verified_files=len(manifest))


def collect(results, archive, slurm_log=None, state=None):
    results, archive = results.resolve(), archive.resolve()
    if archive.exists():
        raise ValueError('Use a fresh archive: '+str(archive))
    has_results = (results/'results.json').is_file()
    if not has_results and (slurm_log is None or not slurm_log.is_file()):
        raise ValueError('No results.json or Slurm log exists; inspect the recorded job first')
    record = json.loads((results/'results.json').read_text()) if has_results else {}
    status = record.get('status', 'preflight-failed')
    if status == 'running':
        if not shutil.which('squeue'):
            raise ValueError('Campaign may still be active; check squeue before collecting')
        active = subprocess.check_output(['squeue', '-h', '-j', str(record['job_id']),
                                          '--format=%A'], text=True).strip()
        if active:
            raise ValueError('Campaign is still queued or running')
        status = 'interrupted'
    elif status not in ('complete', 'failed', 'preflight-failed'):
        raise ValueError('Unknown campaign status')
    # The build uses a private source copy under results. Preserve the build
    # log and run evidence, but do not sweep unrelated source-tree test data
    # into the review archive.
    files = {str(path.relative_to(results)): path for path in sorted(results.rglob('*'))
             if path.is_file() and path.name.endswith(('.json', '.log', '.resources'))
             and path.relative_to(results).parts[0] != 'instrumented-source'} if has_results else {}
    if slurm_log is not None and slurm_log.is_file():
        files['slurm-job.log'] = slurm_log
    if state is not None and state.is_file():
        files['submission.state'] = state
    manifest = {name: digest(path) for name, path in files.items()}
    with tarfile.open(archive, 'w:gz') as tar:
        for name, path in files.items():
            tar.add(path, arcname=str(Path(results.name)/name), recursive=False)
        data = (json.dumps(manifest, indent=2)+'\n').encode()
        info = tarfile.TarInfo(str(Path(results.name)/'REVIEW_SHA256.json'))
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    report = verify(archive)
    report['result_status'] = status
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--slurm-log', type=Path)
    parser.add_argument('--state', type=Path)
    args = parser.parse_args()
    if args.verify and len(args.paths) == 1:
        print(json.dumps(verify(args.paths[0]), indent=2))
    elif not args.verify and len(args.paths) == 2:
        collect(*args.paths, slurm_log=args.slurm_log, state=args.state)
    else:
        parser.error('Use RESULTS ARCHIVE, or --verify ARCHIVE')
