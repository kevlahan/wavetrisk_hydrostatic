"""Full climate replays of the resident input-layout compiler using validation-only reference packets."""
from pathlib import Path
import json,shutil,subprocess,sys,os,hashlib
import argparse
parser=argparse.ArgumentParser(description='Local actual-layout validation only; no production ownership or cluster qualification')
parser.add_argument('--history',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--producer-schedule',action='store_true')
args=parser.parse_args()
H=args.history.resolve(strict=True)
R=Path(__file__).resolve().parents[3]
OUT=args.output.resolve()
if OUT.exists(): raise ValueError('Use a fresh evidence directory')
HERE=H/'test/parallel_block_profile/primitive_276'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
src=OUT/'source'
if not OUT.exists():
 OUT.mkdir(exist_ok=False)
 src=OUT/'source'
 shutil.copytree(H/'.stage276-validation/local-J4/sources/candidate',src,symlinks=True,
  ignore=shutil.ignore_patterns('build-*','bin-*'))
 shutil.copyfile(R/'src/parallel_block_rhs_layout.f90',src/'src/parallel_block_rhs_layout.f90')
 p=src/'Makefile';s=p.read_text().replace('      parallel_block.f90 \\\n','      parallel_block.f90 \\\n      parallel_block_rhs_layout.f90 \\\n')
 s+='\n$(BUILD_DIR)/parallel_block_rhs_layout.o: $(BUILD_DIR)/parallel_block.o\n'
 s+='$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/parallel_block_rhs_layout.o\n';p.write_text(s)
 p=src/'src/multi_level.f90';s=p.read_text()
 s=s.replace('module multi_level_mod\n','module multi_level_mod\n  use parallel_block_rhs_layout_mod\n',1)
 s=s.replace('  type :: Primitive276_Plan\n','  type :: Primitive276_Plan\n     type(Rhs_Input_Layout) :: resident\n',1)
 a='    real(dp) :: start,ref(9),x(16)\n';assert a in s
 s=s.replace(a,a+'    real(dp) :: layout_values(16)\n    integer :: layout_classes(16)\n    logical :: rebuilt_layout\n',1)
 a='    associate(plan=>p276_plan(b))\n';assert a in s
 s=s.replace(a,a+'    call compile_rhs_input_layout(block,p276_generation,plan%resident,rebuilt_layout)\n',1)
 a='''       nn=nodes(ii)
       do t=1,16''';assert a in s
 s=s.replace(a,'''       nn=nodes(ii)
        p=(nn-1)/PATCH_SIZE**2+1
        call rhs_input_values(block,plan%resident,p276_generation,p,mod(nn-1,PATCH_SIZE**2)+1,k, &
             layout_values,layout_classes)
       do t=1,16''',1)
 a='          values(t,ii)=p276_value(block,s,a,k,v)\n';assert a in s
 s=s.replace(a,a+'''          if (layout_classes(t)/=s.or.transfer(layout_values(t),0_int64)/=transfer(values(t,ii),0_int64)) &
                error stop 'RHS278 actual resident operand differs'
           values(t,ii)=layout_values(t)
 ''',1)
 p.write_text(s)
 if args.producer_schedule:
  from schedule_oracle import source as trace_source
  shutil.copyfile(R/'src/parallel_block_rhs_schedule.f90',src/'src/parallel_block_rhs_schedule.f90')
  shutil.copyfile(R/'test/parallel_block_profile/connected_rhs_278/schedule_check.f90',src/'src/rhs_schedule_check.f90')
  accepted_ops=subprocess.check_output(['git','show','ff740744608731a05d53771c32a3055309f3a586:src/ops.f90'],cwd=R,text=True)
  (src/'src/rhs_schedule_oracle.f90').write_text(trace_source(accepted_ops))
  p=src/'Makefile';s=p.read_text()
  s=s.replace('      parallel_block_rhs_layout.f90 '+chr(92)+chr(10),
      '      parallel_block_rhs_layout.f90 '+chr(92)+chr(10)+
      '      parallel_block_rhs_schedule.f90 '+chr(92)+chr(10)+
      '      rhs_schedule_oracle.f90 '+chr(92)+chr(10)+
      '      rhs_schedule_check.f90 '+chr(92)+chr(10))
  s+='\n$(BUILD_DIR)/parallel_block_rhs_schedule.o: $(BUILD_DIR)/parallel_block_rhs_layout.o\n'
  s+='$(BUILD_DIR)/rhs_schedule_oracle.o: $(BUILD_DIR)/patch.o\n'
  s+='$(BUILD_DIR)/rhs_schedule_check.o: $(BUILD_DIR)/parallel_block_rhs_schedule.o $(BUILD_DIR)/rhs_schedule_oracle.o\n'
  s+='$(BUILD_DIR)/multi_level.o: $(BUILD_DIR)/rhs_schedule_check.o\n';p.write_text(s)
  p=src/'src/multi_level.f90';s=p.read_text().replace('module multi_level_mod\n',
      'module multi_level_mod\n  use rhs_schedule_check_mod, only : check_schedule\n',1)
  s=s.replace('  type :: Primitive276_Plan\n',
      '  integer(int64), save :: schedule_counts(5)=0_int64\n  type :: Primitive276_Plan\n',1)
  s=s.replace('    call prepare_native_mass\n','    call prepare_native_mass\n    call check_all_rhs_schedules\n',1)
  s=s.replace('end module multi_level_mod',"""
  subroutine check_all_rhs_schedules
    integer :: context,u
    character(len=80) :: path
    schedule_counts=0_int64;context=0
    call p276_consumer(check_one_rhs_schedule,context)
    write(path,'(a,i0,a)') 'schedule278-rank-',rank,'.txt'
    open(newunit=u,file=trim(path),position='append',status='unknown')
    write(u,*) p276_transaction,native_velocity_plan_generation(),schedule_counts
    close(u)
  end subroutine
  subroutine check_one_rhs_schedule(b,block,context)
    integer, intent(in) :: b
    type(Block276), intent(in) :: block
    class(*), intent(inout) :: context
    integer(int64) :: counts(5)
    call check_schedule(block,native_velocity_plan_generation(),counts)
    schedule_counts=schedule_counts+counts
  end subroutine
end module multi_level_mod""");p.write_text(s)
 (src/'layout-flags.mk').write_text('FFLAGS += -ffp-contract=off\n')
sys.path.insert(0,str(HERE));sys.path.insert(1,str(HERE.parent))
from common import runtime
from run_local_230 import validate as validate_block
from build import objects
from contract import evidence
physics=objects(src)
notes=dict(status='building',production_active=False,source_sha256=sha(R/'src/parallel_block_rhs_layout.f90'),
 validation_overlay_sha256=sha(src/'src/multi_level.f90'),builds={},cases={},scope='actual resident layout only; accepted Stage276 independent reference packets remain validation-only')
def save(): (OUT/'results.json').write_text(json.dumps(notes,indent=2)+'\n')
save()
fixture=H/'.stage276-validation/rk276-primitive-prototype-kit-v2/fixture-J4'
if not fixture.exists():
 candidates=list((H/'.stage276-validation').glob('*/fixture-J4'))
 assert len(candidates)==1,candidates
 fixture=candidates[0]
try:
 for mode in ('optimized','checked'):
  cmd=['make','-f','Makefile','-f','layout-flags.mk','-j1','TEST_CASE=climate','PARAM=param_J4',
       'DEBUG='+('check' if mode=='checked' else 'false'),'BUILD_DIR=build-layout-'+mode,'BIN_DIR=bin-layout-'+mode,
       'FP_CONTRACT=off','STAGE202=on','STAGE203=on']
  print('Serial full climate build',mode,flush=True)
  with (OUT/f'build-{mode}.log').open('x') as log:subprocess.run(cmd,cwd=src,stdout=log,stderr=subprocess.STDOUT,check=True)
  binary=src/f'bin-layout-{mode}/climate'
  notes['builds'][mode]={'command':cmd,'binary_sha256':sha(binary),'log_sha256':sha(OUT/f'build-{mode}.log')};save()
  assert objects(src)==physics,'Physics objects changed'
  for oracle in ([False,True] if mode=='optimized' else [False]):
   label=mode+('-oracle' if oracle else '');directory=OUT/label;directory.mkdir()
   for name in ('simple.in','localJ4J6_checkpoint_0002.bin.zst'):shutil.copyfile(fixture/name,directory/name)
   (directory/'grids').symlink_to(fixture/'grids',target_is_directory=True)
   shutil.copyfile(binary,directory/'climate');(directory/'climate').chmod(0o755)
   env=runtime(oracle);env['OMPI_MCA_mpi_yield_when_idle']='1'
   cmd=['mpirun','--oversubscribe','--bind-to','none','-n','4','./climate','simple.in']
   print('Four-rank exact layout replay',label,flush=True)
   with (directory/'run.log').open('x') as log:subprocess.run(cmd,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800,check=True)
   row=validate_block(directory,H/'.stage266-validation/local-J4/opt',4,oracle,cmd,env)
   row['contract']=evidence(directory,4,True,domains=40)
   if args.producer_schedule:
    records=[]
    for rank in range(4):
     f=directory/f'schedule278-rank-{rank}.txt'
     data=[[int(x) for x in line.split()] for line in f.read_text().splitlines()]
     assert len(data)==48 and [x[0] for x in data]==list(range(1,49)),(rank,len(data))
     assert all(len(x)==7 and x[2]>0 and x[3]>0 and x[4]>0 for x in data)
     records.append({'rank':rank,'rhs_calls':len(data),'sha256':sha(f),
       'direct_events':sum(x[2] for x in data),'sw_events':sum(x[3] for x in data),
       'positive_boundary_events':sum(x[4] for x in data),
       'max_schedule_bytes':max(x[5] for x in data),
       'generations':sorted(set(x[1] for x in data))})
    row['producer_schedule']=records
   notes['cases'][label]=row;save()
 notes['status']='passed';save();print('PASS full climate builds and three actual resident layout replays',flush=True)
except BaseException as exc:
 notes['status']='failed';notes['error']=repr(exc);save();raise
