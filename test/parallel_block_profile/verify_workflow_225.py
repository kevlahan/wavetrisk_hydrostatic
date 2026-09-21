"""Exercise a packaged kit in an isolated home with fake Slurm/scp commands.

Never submits or transfers anything. It uses the real source extraction,
integrity checks, state scripts and archive collector.
"""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile


def verify(kit):
    kit=kit.resolve(strict=True);sys.path.insert(0,str(kit))
    from capability_223 import digest
    from candidate_capability_223 import candidate
    from build_capability_225 import KEYS
    from collect_cluster_logs_225 import verify_archive
    with tempfile.TemporaryDirectory(prefix='stage225-workflow-') as tmp:
        base=Path(tmp);sandbox=base/'test home';sandbox.mkdir();bins=base/'commands';bins.mkdir()
        (sandbox/'stage187-fixture-20260911-145307').mkdir()
        (sandbox/'stage187-fixture-20260911-145307/simple.in').write_text('fixture presence stub\n')
        for name,body in {
            'sbatch':'printf "%s\\n" "$@" >> "$STAGE223_TEST_CALLS"\nprintf "12345;testcluster\\n"\n',
            'squeue':'exit 0\n','sha256sum':'exit 0\n',
            'scp':'printf "%s\\n" "$@" >> "$STAGE223_TEST_TRANSFERS"\n'
        }.items():
            p=bins/name;p.write_text('#!/bin/bash\n'+body);p.chmod(0o755)
        env=dict(os.environ,HOME=str(sandbox),USER='workflow-test',PATH=str(bins)+os.pathsep+os.environ['PATH'],
                 STAGE223_TEST_CALLS=str(base/'submissions'),STAGE223_TEST_TRANSFERS=str(base/'transfers'))
        def script(name):return subprocess.run(['bash',str(kit/name)],env=env,capture_output=True,text=True)
        (sandbox/'stage225-state.sh').write_text('previous state\n')
        r=script('bootstrap_225.sh');assert r.returncode==0,r.stderr
        tree=sandbox/'wt-stage225'
        assert (tree/'reference/src/main.f90').is_file()
        assert next(sandbox.glob('stage225-state.previous.*')).read_text()=='previous state\n'
        state=(sandbox/'stage225-state.sh').read_bytes()
        r=script('bootstrap_225.sh');assert r.returncode!=0
        assert (sandbox/'stage225-state.sh').read_bytes()==state
        assert (base/'submissions').read_text().count('build_225.sbatch')==1
        shutil.copytree(tree/'reference',tree/'candidate',symlinks=True);candidate(tree/'candidate')
        m=json.loads((kit/'stage225-sources.json').read_text());builds={}
        for key in KEYS:
            p=tree/(key+'.fixture-binary');p.write_text('nonexecuted test binary '+key)
            builds[key]=dict(binary=p.name,sha256=digest(p),sources=m['reference_sources' if key.endswith('baseline') else 'candidate_sources'])
        record=dict(status='passed',fp_contract='off',manifest_sha256=digest(kit/'stage225-sources.json'),
                    builds=builds,build_host={'node':'bb13.example.test'})
        (tree/'capability-225-builds.json').write_text(json.dumps(record))
        r=script('launch_225.sh');assert r.returncode==0,r.stderr
        assert '--nodelist=bb13' in (base/'submissions').read_text()
        r=script('launch_225.sh');assert r.returncode!=0
        assert (base/'submissions').read_text().count('cluster_225.sbatch')==1
        out=Path(subprocess.check_output(['bash','-c','source "$HOME/stage225-state.sh"; printf "%s" "$out"'],env=env,text=True))
        out.mkdir()
        (out/'results.json').write_text(json.dumps(dict(stage=225,status='failed',slurm_environment={'SLURM_JOB_ID':'12345'})))
        (out/'capability223-rank-0.txt').write_text('partial failure evidence\n')
        r=script('collect_225.sh');assert r.returncode==0,r.stderr
        archive=sandbox/'stage225-review-12345.tar.gz';before=digest(archive);verify_archive(archive)
        r=script('collect_225.sh');assert r.returncode==0,r.stderr
        assert digest(archive)==before
        assert (base/'transfers').read_text().count('kevlahan@murier.math.mcmaster.ca:~/Desktop/')==2
    print('PASS: packaged bootstrap, state preservation, duplicate guards, same-node launch, failure archive and verified transfer reuse.')


if __name__=='__main__':verify(Path(sys.argv[1]))
