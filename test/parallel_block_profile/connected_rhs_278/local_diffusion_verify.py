"""Qualify compact diffusion/halo and pressure schedules on the real J4 restart.

Validation overlay only. Domain ownership and the accepted RHS/RK path remain.
Authoritative diffusion substitution proves consumers/science, not owner cutover.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--history',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();H=a.history.resolve(strict=True);O=a.output.resolve()
R=Path(__file__).resolve().parents[3];T=Path(__file__).resolve().parent
O.mkdir(exist_ok=False);S=O/'source'
shutil.copytree(H/'.stage266-validation/local-J4/source',S,symlinks=True,
 ignore=shutil.ignore_patterns('build-*','bin-*'))
modules=['parallel_block_mass','parallel_block_rhs','parallel_block_rhs_diffusion',
 'parallel_block_primitive_walk','parallel_block_rhs_compile','parallel_block_rhs_exchange','parallel_block_rhs_boundary']
for name in modules:shutil.copyfile(R/f'src/{name}.f90',S/f'src/{name}.f90')
subprocess.run(['patch','-p1','-f','-i',str(T/'diffusion-validation.patch')],cwd=S,check=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
notes=dict(status='building',production_active=False,ownership_credit=False,
 source_sha256={f'src/{n}.f90':sha(R/f'src/{n}.f90') for n in modules},
 overlay_sha256=sha(T/'diffusion-validation.patch'),builds={},cases={})
def save(): (O/'results.json').write_text(json.dumps(notes,indent=2)+'\n')
save()
sys.path[:0]=[str(H/'test/parallel_block_profile/primitive_276'),str(H/'test/parallel_block_profile')]
from common import runtime
from build import objects
from run_local_230 import validate
fixture=H/'.stage276-validation/rk276-primitive-prototype-kit-v2/fixture-J4'
physics=objects(S)
try:
 for mode in ('optimized','checked'):
  cmd=['make','-f','Makefile','-f','stage264-flags.mk','-j1','TEST_CASE=climate','PARAM=param_J4',
   'DEBUG='+('check' if mode=='checked' else 'false'),'BUILD_DIR=build-diff-'+mode,
   'BIN_DIR=bin-diff-'+mode,'FP_CONTRACT=off','STAGE202=on','STAGE203=on']
  with (O/f'build-{mode}.log').open('x') as f:
   subprocess.run(cmd,cwd=S,stdout=f,stderr=subprocess.STDOUT,check=True)
  assert objects(S)==physics,'Physics objects changed'
  notes['builds'][mode]={'command':cmd,'binary_sha256':sha(S/f'bin-diff-{mode}/climate')};save()
  variants=[('shadow',False,False),('installed-oracle',True,True)] if mode=='optimized' else [('installed',True,False)]
  for suffix,install,oracle in variants:
   label=mode+'-'+suffix;D=O/label;D.mkdir()
   for n in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copyfile(fixture/n,D/n)
   (D/'grids').symlink_to(fixture/'grids',target_is_directory=True)
   shutil.copyfile(S/f'bin-diff-{mode}/climate',D/'climate');(D/'climate').chmod(0o755)
   env=runtime(oracle);env['OMPI_MCA_mpi_yield_when_idle']='1'
   env['WAVETRISK_RHS278_DIFFUSION_REPLACE']='1' if install else '0'
   cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
   with (D/'run.log').open('x') as f:
    subprocess.run(cmd,cwd=D,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=1800)
   row=validate(D,H/'.stage266-validation/local-J4/opt',4,oracle,cmd,env)
   row['authoritative_diffusion']=install;row['block_ownership']=False
   rows=[]
   for rank in range(4):
    pressure=D/f'pressure278-rank-{rank}.txt'
    records=[[int(x) for x in line.split()] for line in pressure.read_text().splitlines()]
    assert len(records)==48,(label,rank,len(records))
    record={'rank':rank,'pressure_transactions':len(records),'pressure_nodes':records[-1][1],
     'pressure_sha256':sha(pressure)}
    if install:
     file=D/f'diffusion278-installed-rank-{rank}.txt'
     data=[[int(x) for x in line.split()] for line in file.read_text().splitlines()]
     counts=[sum(line[0]==m for line in data) for m in range(1,5)]
     assert counts==[1440,1440,1440,4320],counts
     record['installed_sweeps']=counts;record['installed_sha256']=sha(file)
    rows.append(record)
   row['numerical_checks']=rows;notes['cases'][label]=row;save()
 notes['status']='passed';save()
 print('PASS compact diffusion/pressure qualification; block-owned RHS/RK remains unqualified')
except BaseException as e:
 notes['status']='failed';notes['error']=repr(e);save();raise
