#!/usr/bin/env python3
"""Link the actual old/new communication modules; compare complete bit patterns."""
import argparse,json,subprocess,hashlib
from pathlib import Path

def run(a):
    a.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('test_boundary_215.f90').resolve()
    result=dict(stage=215,status='running',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),runs={})
    for name,build in [('baseline',a.baseline),('candidate',a.candidate),('checked',a.checked)]:
        build=build.resolve(strict=True);root=a.out/name;root.mkdir();exe=root/'test'
        # Only the modules needed by comm_mpi; no production program or physics.
        modules=['kind','param_J4','shared','coord_arithmetic','geom','patch','dyn_array','arch','domain',
                 'parallel_block_profile','domain_ops','comm','comm_mpi']
        cmd=['mpif90','-O0','-g','-fcheck=all','-ffp-contract=off','-I'+str(build/'mod'),str(source),
             *[str(build/(m+'.o')) for m in modules],'-o',str(exe)]
        built=subprocess.run(cmd,capture_output=True,text=True)
        (root/'build.log').write_text(built.stdout+built.stderr)
        if built.returncode:raise ValueError(built.stderr)
        p=subprocess.run(['mpirun','-np','4',str(exe)],cwd=root,capture_output=True,text=True,timeout=60)
        (root/'run.log').write_text(p.stdout+p.stderr)
        if p.returncode:raise ValueError(name+': '+p.stderr)
        hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in root.glob('boundary-*.txt')}
        if len(hashes)!=4:raise ValueError('Missing rank dumps')
        if name!='baseline' and hashes!=result['runs']['baseline']['hashes']:raise ValueError('Halo bits differ: '+name)
        result['runs'][name]=dict(hashes=hashes,binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    bad=a.out/'stale';bad.mkdir()
    p=subprocess.run(['mpirun','-np','4',str(a.out/'checked/test'),'stale'],cwd=bad,capture_output=True,text=True,timeout=60)
    (bad/'run.log').write_text(p.stdout+p.stderr)
    if p.returncode==0 or 'Stage215 stale send' not in p.stderr:raise ValueError('Stale-cache adversarial test failed')
    result.update(status='passed',stale_cache_rejected=True,scope='All overloads, mixed fields, range/single/empty levels, reuse, signed zero/edges, ordered local aliases, three active ranks and one idle; exact baseline bits.')
    (a.out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(result['scope'])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('baseline','candidate','checked','out'):p.add_argument('--'+k,type=Path,required=True)
    run(p.parse_args())
