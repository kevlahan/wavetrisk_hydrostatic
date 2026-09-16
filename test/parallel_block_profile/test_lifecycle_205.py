"""Actual Fortran accounting, trace corruption, calibration and campaign guards."""
import os,subprocess,tempfile,unittest,copy
from pathlib import Path
from analyze_lifecycle_205 import collect,summarize,NR,PHASES
from cluster_205 import schedule
from instrument_lifecycle_205 import instrument
BASE=Path(__file__).parent
class Lifecycle205(unittest.TestCase):
 def test_actual_fortran_observer(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);(d/'clock.c').write_text('''double workflow_clock199(void){static double t=0;return t+=.001;}
int workflow_affinity200(void){return 0;} void workflow_hw_init201(void){} void workflow_hw_begin201(void){}
void workflow_hw_tick201(int s,int c){} void workflow_hw_finish201(int r){}
''')
   (d/'driver.f90').write_text('''program test
use workflow_profile_195_mod
use iso_fortran_env,only:int64
implicit none
integer :: i
call wf_init
call measure_enter205(1)
call measure_leave205(1)
do i=1,2
 call wf_begin
 call wf_phase(1)
 ! MPI outside selected scopes is deliberately not timed.
 call measure_enter205(26)
 call measure_leave205(26)
 call measure_enter205(18)
 call wf_enter(2)
 call lifecycle_event205(8,0,1_int64,1_int64,4_int64,0_int64)
 call measure_enter205(10)
 call lifecycle_event205(1,1,1_int64,1_int64,1_int64,0_int64)
 call lifecycle_producer_begin205(10_int64,2_int64)
 call lifecycle_producer_end205(1_int64,30_int64,3_int64)
 call measure_enter205(27)
 if(lifecycle_selected205()) call lifecycle_event205(40,1,0_int64,4_int64,4_int64,0_int64)
 call measure_leave205(27)
 call measure_leave205(10)
 call wf_leave(2)
 call measure_leave205(18)
 call wf_phase(2)
 call wf_phase(3)
 call wf_phase(4)
 call wf_flag(2)
 call wf_phase(5)
 if(i==2) call wf_flag(1)
 call wf_end
end do
call wf_finish(0)
end program
''')
   subprocess.run(['cc','-c','clock.c'],cwd=d,check=True)
   subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-132',str(BASE.resolve()/'workflow_profile_205.f90'),'driver.f90','clock.o','-o','probe'],cwd=d,check=True)
   for mode in (0,1):
    run=d/str(mode);run.mkdir();env=os.environ.copy();env.update(WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200='0',WAVETRISK_LIFECYCLE205=str(mode));env.pop('WAVETRISK_EXPECT_CPU200',None)
    subprocess.run([str(d/'probe')],cwd=run,env=env,check=True)
    q=collect(run,1,mode);self.assertEqual(q['ordinary_steps'],[1])
    if mode:
     self.assertEqual(q['lifecycle']['plans']['1']['reused'],1)
     self.assertAlmostEqual(q['lifecycle']['collective_sites']['1']['seconds_per_rank_step'],.001)
     self.assertEqual(q['lifecycle']['collective_sites']['1']['owners'],[10])
     self.assertEqual(q['lifecycle']['event_counts']['13'],1)
     trace=run/'lifecycle205-rank-0.txt';original=trace.read_text();trace.write_text(original.replace(' 40 1 0 4 4 0',' 40 1 0 -5 4 0',1))
     with self.assertRaisesRegex(ValueError,'collective counts'):collect(run,1,mode)
     trace.write_text(original)
     p=run/'operations205-rank-0.txt';p.write_text(p.read_text().replace('self 1 2 27 1','self 1 2 27 -1',1))
     with self.assertRaises(ValueError):collect(run,1,mode)
    else:self.assertEqual(q['lifecycle']['event_counts'],{})
 def fixture(self):
  r={'runs':{}}
  for side,t in [('baseline',.9),('candidate',.85),('legacy',.5)]:
   ledger=[{'context':c,'operation':p,'name':str(p),'seconds':t if (c,p)==(1,0) else 0.,'inclusive_seconds':0.} for c in range(9) for p in range(NR+1)]
   h={'quality_flags':[],'phases':{k:{'wall_seconds_rank_mean':t if k in ('whole step','dynamics WT') else 0.} for k in ('whole step',*PHASES)}}
   phases=('bare',) if side=='candidate' else ('bare','coarse','operations')
   for phase in phases:
    for rep in (1,2):
     label=f'{side}-{phase}-{rep}';r['runs'][label]={'label':label,'side':side,'phase':phase,'round':rep,'ordinary_step_median':t,'hardware':copy.deepcopy(h),'operations':{'wall_rank_mean':t,'ledger':ledger,'lifecycle':{}}}
  return r
 def test_calibration_and_signed_closure(self):
  r=self.fixture();s=summarize(r);self.assertEqual(s['decision'],'ready_for_consumer_lifetime_audit');self.assertFalse(s['optimization_authorized'])
  self.assertAlmostEqual(sum(x['difference'] for x in s['signed_locations']),.4)
  r['runs']['baseline-operations-2']['hardware']['phases']['dynamics WT']['wall_seconds_rank_mean']+=.1
  self.assertEqual(summarize(r)['decision'],'observer_effects_require_review')
  r=self.fixture();r['runs']['legacy-operations-1']['ordinary_step_median']=.6
  self.assertEqual(summarize(r)['decision'],'observer_effects_require_review')
  del r['runs']['legacy-coarse-2']
  with self.assertRaises(ValueError):summarize(r)
 def test_schedule_adjacent_reversed_pairs(self):
  jobs=schedule({});self.assertEqual(len(jobs),19)
  self.assertNotIn('candidate_profile',[x[1] for x in jobs])
  observer=[x for x in jobs if x[2] in ('coarse','operations')]
  self.assertEqual([x[2] for x in observer],['coarse','operations','coarse','operations','operations','coarse','operations','coarse'])
  self.assertEqual([x[1] for x in observer],['baseline_profile']*2+['legacy_profile']*4+['baseline_profile']*2)
 def test_conditional_instrumentation_is_inside_predicate(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);(root/'src').mkdir();p=root/'src/example.f90'
   p.write_text('module example\ncontains\nsubroutine sample(flag)\nlogical :: flag\nif (flag) call build_block_writeback_plan\nend subroutine\nend module\n')
   q=instrument(root);s=p.read_text();self.assertIn('if (flag) then\ncall measure_enter205(9)',s);self.assertTrue(q['src/example.f90']['original_bytes_recoverable'])
if __name__=='__main__':unittest.main()
