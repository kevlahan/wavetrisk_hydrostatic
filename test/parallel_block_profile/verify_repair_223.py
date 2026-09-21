"""Exercise the placement repair in an isolated home with fake Slurm/scp commands.

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


def verify(kit, repair):
    kit=kit.resolve(strict=True);sys.path.insert(0,str(kit))
    from capability_223 import digest
    from candidate_capability_223 import candidate
    from build_capability_223 import KEYS
    from collect_cluster_logs_223 import verify_archive
    with tempfile.TemporaryDirectory(prefix='stage223-workflow-') as tmp:
        base=Path(tmp);sandbox=base/'test home';sandbox.mkdir();bins=base/'commands';bins.mkdir()
        (sandbox/'stage187-fixture-20260911-145307').mkdir()
        (sandbox/'stage187-fixture-20260911-145307/simple.in').write_text('fixture presence stub\n')
        for name,body in {
            'sbatch':'printf "%s\\n" "$@" >> "$STAGE223_TEST_CALLS"\nprintf "54724;testcluster\\n"\n',
            'squeue':'exit 0\n','sha256sum':'exit 0\n',
            'scp':'printf "%s\\n" "$@" >> "$STAGE223_TEST_TRANSFERS"\n'
        }.items():
            p=bins/name;p.write_text('#!/bin/bash\n'+body);p.chmod(0o755)
        env=dict(os.environ,HOME=str(sandbox),USER='workflow-test',PATH=str(bins)+os.pathsep+os.environ['PATH'],
                 STAGE223_TEST_CALLS=str(base/'submissions'),STAGE223_TEST_TRANSFERS=str(base/'transfers'))
        def script(name):return subprocess.run(['bash',str(kit/name)],env=env,capture_output=True,text=True)
        (sandbox/'stage223-state.sh').write_text('previous state\n')
        r=script('bootstrap_223.sh');assert r.returncode==0,r.stderr
        tree=sandbox/'wt-stage223'
        assert (tree/'reference/src/main.f90').is_file()
        assert next(sandbox.glob('stage223-state.previous.*')).read_text()=='previous state\n'
        state=(sandbox/'stage223-state.sh').read_bytes()
        r=script('bootstrap_223.sh');assert r.returncode!=0
        assert (sandbox/'stage223-state.sh').read_bytes()==state
        assert (base/'submissions').read_text().count('build_223.sbatch')==1
        shutil.copytree(tree/'reference',tree/'candidate',symlinks=True);candidate(tree/'candidate')
        m=json.loads((kit/'stage223-sources.json').read_text());builds={}
        for key in KEYS:
            p=tree/(key+'.fixture-binary');p.write_text('nonexecuted test binary '+key)
            builds[key]=dict(binary=p.name,sha256=digest(p),sources=m['reference_sources' if key.endswith('baseline') else 'candidate_sources'])
        record=dict(status='passed',fp_contract='off',manifest_sha256=digest(kit/'stage223-sources.json'),
                    builds=builds,build_host={'node':'bb13.example.test'})
        (tree/'capability-223-builds.json').write_text(json.dumps(record))
        r=script('launch_223.sh');assert r.returncode==0,r.stderr
        assert '--nodelist=bb13' in (base/'submissions').read_text()
        r=script('launch_223.sh');assert r.returncode!=0
        assert (base/'submissions').read_text().count('cluster_223.sbatch')==1
        out=Path(subprocess.check_output(['bash','-c','source "$HOME/stage223-state.sh"; printf "%s" "$out"'],env=env,text=True))
        out.mkdir()
        (out/'results.json').write_text(json.dumps(dict(stage=223,status='failed',runs={},schedule=[],
            error="CalledProcessError: --probe returned 1",slurm_environment={'SLURM_JOB_ID':'54724'})))
        old_results=(out/'results.json').read_bytes()
        old_manifest=(kit/'stage223-sources.json').read_bytes()
        state=(sandbox/'stage223-state.sh').read_bytes()
        retry=lambda: subprocess.run(['bash',str(repair/'retry_223.sh')],env=env,capture_output=True,text=True)
        # Active old job must prevent changes.
        (bins/'squeue').write_text('#!/bin/bash\nprintf "54724\\n"\n')
        r=retry();assert r.returncode!=0,r.stdout
        assert (sandbox/'stage223-state.sh').read_bytes()==state
        (bins/'squeue').write_text('#!/bin/bash\nexit 0\n')
        (bins/'sbatch').write_text('#!/bin/bash\nprintf "%s\\n" "$@" >> "$STAGE223_TEST_CALLS"\nprintf "54725;testcluster\\n"\n')
        # Ensure a new output directory even when both fake submissions are instantaneous.
        env['PATH']=str(bins)+os.pathsep+os.environ['PATH']
        (bins/'date').write_text('#!/bin/bash\nprintf "20260919-120000\\n"\n');(bins/'date').chmod(0o755)
        r=retry();assert r.returncode==0,r.stderr+r.stdout
        assert (sandbox/'stage223-state.before-placement-repair.sh').read_bytes()==state
        assert (out/'results.json').read_bytes()==old_results
        assert (kit/'stage223-sources.json').read_bytes()==old_manifest
        assert (base/'submissions').read_text().count('cluster_223.sbatch')==2
        r=retry();assert r.returncode==0,r.stderr
        assert (base/'submissions').read_text().count('cluster_223.sbatch')==2
        # A repeat installer is harmless and rechecks the complete kit.
        r=subprocess.run([sys.executable,str(repair/'repair_223.py'),str(kit)],capture_output=True,text=True)
        assert r.returncode==0,r.stderr
        assert 'already installed' in r.stdout
        newout=Path(subprocess.check_output(['bash','-c','source "$HOME/stage223-state.sh"; printf "%s" "$out"'],env=env,text=True))
        assert newout!=out
        newout.mkdir();(newout/'results.json').write_text(json.dumps(dict(stage=223,status='failed',slurm_environment={'SLURM_JOB_ID':'54725'})))
        r=script('collect_223.sh');assert r.returncode==0,r.stderr
        archive=sandbox/'stage223-review-54725.tar.gz';before=digest(archive);verify_archive(archive)
        r=script('collect_223.sh');assert r.returncode==0,r.stderr
        assert digest(archive)==before
        assert (base/'transfers').read_text().count('kevlahan@murier.math.mcmaster.ca:~/Desktop/')==2
    print('PASS: repair, active-job guard, source/build identity, state/results preservation, same-node retry, duplicate guard, idempotence and archive transfer.')


if __name__=='__main__':verify(Path(sys.argv[1]),Path(sys.argv[2]))
