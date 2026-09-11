#!/usr/bin/env python3
"""Sequential isolated J7 restart diagnostics, never a performance measurement."""
import argparse,json,os,re,shutil,subprocess,time
from pathlib import Path
from experiment import SWITCHES,digest
from profile_protocol import checkpoint_name
parser=argparse.ArgumentParser(description=__doc__)
for name in ('legacy','block','fixture','out'):
    parser.add_argument('--'+name,type=Path,required=True)
parser.add_argument('--time-end', default='0.3105')
parser.add_argument('--domain', type=int, help='Restrict observations to this zero-based global Domain')
parser.add_argument('--physics-all', action='store_true', help='Observe packed physics on all Domains')
parser.add_argument('--only', choices=('legacy','block'))
parser.add_argument('--oracle', action='store_true')
parser.add_argument('--remap-oracle', action='store_true',
                    help='Enable the remap oracle independently of dynamics/adaptation')
args=parser.parse_args()
args.out.mkdir(parents=True,exist_ok=False)
env={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
env.update({key:'0' for key in SWITCHES})
env['WAVETRISK_RESTART_PROBE']='1'
if args.domain is not None: env['WAVETRISK_PROBE_DOMAIN']=str(args.domain)
if args.physics_all: env['WAVETRISK_PHYSICS_PROBE_ALL']='1'
if args.oracle:
    env['WAVETRISK_VALIDATE_BLOCK_DYNAMICS']='1'
    env['WAVETRISK_VALIDATE_BLOCK_ADAPTATION']='1'
if args.oracle or args.remap_oracle:
    env['WAVETRISK_VALIDATE_BLOCK_REMAP']='1'
env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS')})
for label,binary in (('legacy',args.legacy),('block',args.block)):
    if args.only and label!=args.only: continue
    root=args.out.resolve()/label
    root.mkdir()
    shutil.copy2(binary,root/'climateJ5')
    shutil.copy2(args.fixture/checkpoint_name((args.fixture/'simple.in').read_text()),root)
    (root/'grids').symlink_to((args.fixture/'grids').resolve(),target_is_directory=True)
    text,count=re.subn(r'(?m)^time_end\s+\S+',f'time_end           {args.time_end}',
                      (args.fixture/'simple.in').read_text())
    if count!=1: raise ValueError('Expected exactly one time_end setting')
    (root/'simple.in').write_text(text)
    recorded={k:v for k,v in env.items() if k.startswith('WAVETRISK_') or k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'}
    probe_manifest=binary.resolve().parent.parent/'probe-identity.json'
    identity={'binary':str(binary),'sha256':digest(binary),'env':recorded,
              'seed_sha256':digest(root/checkpoint_name(text)),'input_sha256':digest(root/'simple.in')}
    if probe_manifest.exists():identity['probe_build']=json.loads(probe_manifest.read_text())
    (root/'identity.json').write_text(json.dumps(identity,indent=2))
    print('Starting',label,root,flush=True)
    start=time.monotonic()
    with (root/'run.log').open('x') as log:
        result=subprocess.run(['mpirun','-n','4','./climateJ5','simple.in'],cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
    identity.update(returncode=result.returncode,elapsed_seconds=time.monotonic()-start)
    (root/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    if result.returncode:raise RuntimeError(f'{label} failed: {result.returncode}, see {root}/run.log')
    if 'Total cpu time' not in (root/'run.log').read_text():raise RuntimeError('Missing normal completion marker')
    print('Completed',label,flush=True)
