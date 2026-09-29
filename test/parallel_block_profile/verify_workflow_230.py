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
    from build_capability_230 import materialize,source_key,source_dir
    from build_capability_230 import KEYS
    from collect_cluster_logs_230 import verify_archive
    with tempfile.TemporaryDirectory(prefix='stage230-workflow-') as tmp:
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
        (sandbox/'stage230-state.sh').write_text('previous state\n')
        r=script('bootstrap_230.sh');assert r.returncode==0,r.stderr
        tree=sandbox/'wt-stage230'
        assert (tree/'reference/src/main.f90').is_file()
        assert next(sandbox.glob('stage230-state.previous.*')).read_text()=='previous state\n'
        state=(sandbox/'stage230-state.sh').read_bytes()
        r=script('bootstrap_230.sh');assert r.returncode!=0
        assert (sandbox/'stage230-state.sh').read_bytes()==state
        assert (base/'submissions').read_text().count('build_230.sbatch')==1
        for kind in ('candidate','checked'):
            materialize(kit,tree,kind)
        m=json.loads((kit/'stage230-sources.json').read_text());builds={}
        for key in KEYS:
            p=tree/(key+'.fixture-binary');p.write_text('nonexecuted test binary '+key)
            builds[key]=dict(binary=p.name,sha256=digest(p),sources=m[source_key(key.split('-',1)[1])])
        record=dict(status='passed',fp_contract='off',manifest_sha256=digest(kit/'stage230-sources.json'),
                    builds=builds,build_host={'node':'bb13.example.test'})
        (tree/'capability-230-builds.json').write_text(json.dumps(record))
        r=script('launch_230.sh');assert r.returncode==0,r.stderr
        assert '--nodelist=bb13' in (base/'submissions').read_text()
        r=script('launch_230.sh');assert r.returncode!=0
        assert (base/'submissions').read_text().count('cluster_230.sbatch')==1
        out=Path(subprocess.check_output(['bash','-c','source "$HOME/stage230-state.sh"; printf "%s" "$out"'],env=env,text=True))
        out.mkdir()
        (out/'results.json').write_text(json.dumps(dict(stage=230,status='failed',slurm_environment={'SLURM_JOB_ID':'12345'})))
        (out/'capability223-rank-0.txt').write_text('partial failure evidence\n')
        r=script('collect_230.sh');assert r.returncode==0,r.stderr
        archive=sandbox/'stage230-review-12345.tar.gz';before=digest(archive);verify_archive(archive)
        r=script('collect_230.sh');assert r.returncode==0,r.stderr
        assert digest(archive)==before
        assert (base/'transfers').read_text().count('kevlahan@murier.math.mcmaster.ca:~/Desktop/')==2
    print('PASS: packaged bootstrap, state preservation, duplicate guards, same-node launch, failure archive and verified transfer reuse.')


if __name__=='__main__':verify(Path(sys.argv[1]))
