#!/usr/bin/env python3
"""Checked production integration: real 41-rank restart without prototype hooks."""
import argparse,json,shutil,subprocess,time
from pathlib import Path
from capability_223 import digest,coverage,numerical
from profile_protocol import environment,parameters,checkpoint_name
from run_local_220 import SELECTORS

def run(tree,fixture,reference,out):
    out.mkdir(exist_ok=False)
    text=(fixture/'simple.in').read_text();p=parameters(text)
    for n in ['simple.in',checkpoint_name(text)]:shutil.copy2(fixture/n,out/n)
    (out/'grids').symlink_to((fixture/'grids').resolve(),target_is_directory=True)
    binary=tree/'bin-226-checked/climate';shutil.copy2(binary,out/'climate')
    env=environment(oracles=True,detail=False);env.update(SELECTORS)
    assert not any(k.startswith('WAVETRISK_STAGE223') for k in env)
    command=['mpirun','--oversubscribe','--bind-to','none','-n','41','./climate','simple.in']
    t=time.monotonic()
    with (out/'run.log').open('x') as log:
        r=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
    assert r.returncode==0, str(out/'run.log')
    a,c=coverage(out,'J4',p['run_id']);assert len(a['step_records'])==12 and c['remaps']==2
    gate=numerical(reference,out,c['new_checkpoint'],40,True);assert gate['passed'],gate
    from analyze import analyze
    ref=analyze(reference/'run.log')
    assert [s['state'] for s in a['step_records']]==[s['state'] for s in ref['step_records']]
    log=(out/'run.log').read_text()
    for marker in ['Resident block execution: MPI tasks = 41, geometric Domains = 40','Block dynamics validation oracle enabled','Block adaptation validation oracle enabled','Block vertical-remap validation oracle enabled']:
        assert marker in log,marker
    assert not list(out.glob('capability223-rank-*.txt')) and 'Stage223 capability prototype' not in log
    result=dict(status='passed',ranks=41,domains=40,steps=12,coverage=c,oracles=True,prototype_hooks=False,
                command=command,elapsed_seconds=time.monotonic()-t,binary_sha256=digest(binary),log_sha256=digest(out/'run.log'),
                numerical=gate,production_sources=json.loads((tree.parent/'sources.json').read_text()),
                environment={k:v for k,v in env.items() if k.startswith('WAVETRISK_')},
                build_log_sha256=digest(tree.parent/'build.log'),local_flags=(tree/'local-flags.mk').read_text())
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: production 41/40 checked full-oracle run; 12 steps, two reloads/remaps, exact checkpoint; no prototype hooks')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['tree','fixture','reference','out']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();run(a.tree.resolve(),a.fixture.resolve(),a.reference.resolve(),a.out.resolve())
