"""Build a pinned block archive; never compile in or modify the working checkout."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
from build_legacy import tracked_entries, verify
from experiment import digest


def build(repo, out, ref):
    repo, out = repo.resolve(strict=True), out.resolve()
    if out.is_relative_to(repo) or repo.is_relative_to(out):
        raise ValueError('Use a new directory outside the checkout')
    revision = subprocess.check_output(['git','-C',str(repo),'rev-parse',ref+'^{commit}'],text=True).strip()
    entries = tracked_entries(repo,revision)
    out.mkdir(parents=True,exist_ok=False)
    archive = subprocess.Popen(['git','-C',str(repo),'archive',revision],stdout=subprocess.PIPE)
    try:
        subprocess.run(['tar','-x','-C',str(out)],stdin=archive.stdout,check=True)
    finally:
        archive.stdout.close()
    if archive.wait(): raise RuntimeError('Archive failed')
    verify(out,entries)
    physics = out/'src/physics/simple_physics/phyparam'
    names = sorted(p.stem for p in (physics/'physics').glob('*.F90'))
    combined = subprocess.check_output(['bash','bash/concatenate_all_code.sh',*names],cwd=physics)
    dependencies = subprocess.check_output(['bash','bash/makedeps.sh','/dev/stdin'],cwd=physics,input=combined)
    (physics/'Makefile.inc').write_bytes(dependencies)
    command = ['make','-j1','DEBUG=false','PARAM=param_J5','TEST_CASE=climate']
    with (out/'build.log').open('x') as log:
        log.write(f'Pinned block revision: {revision}\nCommand: {command}\n');log.flush()
        result = subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT)
    verify(out,entries)
    if result.returncode: raise RuntimeError(f'Build failed: {out}/build.log')
    identity = {'role':'block','revision':revision,'tracked_contents_unchanged':True,
                'binary_sha256':digest(out/'bin/climate'),'build_log_sha256':digest(out/'build.log'),
                'command':command,'host':platform.uname()._asdict()}
    (out/'build-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    print(json.dumps(identity,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--ref',default='1b01510ae655eee0e2aab1a8d911e55c147173f2')
    a=p.parse_args();build(a.repo,a.out,a.ref)
