"""Check the actual Fortran observer, nested accounting and rejection paths."""
import os,subprocess,tempfile,unittest
from pathlib import Path
from analyze_operations_204 import collect,summarize,NR,PHASES
from cluster_204 import schedule
BASE=Path(__file__).parent
class Observer204(unittest.TestCase):
 def test_fortran_nested_contexts_and_corruption(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);(d/'clock.c').write_text('''double workflow_clock199(void){static double t=0;return t+=.001;}
int workflow_affinity200(void){return 0;}
void workflow_hw_init201(void){} void workflow_hw_begin201(void){}
void workflow_hw_tick201(int s,int c){} void workflow_hw_finish201(int r){}
''')
   (d/'driver.f90').write_text('''program test
use workflow_profile_195_mod
implicit none
integer :: i
call wf_init
! Initialization hooks must be inert.
call measure_enter204(1)
call measure_leave204(1)
do i=1,2
 call wf_begin
 call wf_phase(1)
 call measure_enter204(2)
 call wf_enter(1)
 call measure_enter204(3)
 call measure_enter204(47)
 call measure_leave204(47)
 call measure_leave204(3)
 call wf_leave(1)
 call measure_leave204(2)
 call wf_phase(2)
 call wf_phase(3)
 call measure_enter204(36)
 call measure_leave204(36)
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
   subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-132',str(BASE.resolve()/'workflow_profile_204.f90'),'driver.f90','clock.o','-o','probe'],cwd=d,check=True)
   for mode in (0,1):
    run=d/str(mode);run.mkdir();env=os.environ.copy();env.update(WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200='0',WAVETRISK_OPERATIONS204=str(mode))
    env.pop('WAVETRISK_EXPECT_CPU200',None)
    subprocess.run([str(d/'probe')],cwd=run,env=env,check=True)
    r=collect(run,1,mode);self.assertEqual(r['ordinary_steps'],[1])
    if mode:
     lookup={(x['context'],x['operation']):x for x in r['ledger']}
     self.assertAlmostEqual(lookup[1,47]['seconds'],.001)
     self.assertGreater(lookup[1,2]['inclusive_seconds'],lookup[1,2]['seconds'])
     self.assertEqual(r['mpi_owners'][0]['owner'],3)
     p=run/'operations204-rank-0.txt';p.write_text(p.read_text().replace('self 1 1 47 1','self 1 1 47 -1',1))
     with self.assertRaises(ValueError):collect(run,1,mode)
 def test_signed_ranking_and_calibration(self):
  result={'runs':{}}
  for side,t in [('baseline',.9),('candidate',.85),('legacy',.5)]:
   ledger=[{'context':c,'operation':p,'name':str(p),'seconds':t if (c,p)==(1,3) else 0.,
            'inclusive_seconds':t if (c,p)==(1,3) else 0.,'actual_slow_rank_path_seconds':t if (c,p)==(1,3) else 0.}
           for c in range(9) for p in range(NR+1)]
   h={'quality_flags':[],'phases':{k:{'wall_seconds_rank_mean':t if k in ('whole step','dynamics WT') else 0.} for k in ('whole step',*PHASES)}}
   for phase,repeats in [('bare',(1,2)),('coarse',(0,)),('operations',(1,2))]:
    for rep in repeats:
     label=f'{side}-{phase}-{rep}'
     result['runs'][label]={'label':label,'side':side,'phase':phase,'round':rep,'ordinary_step_median':t,
      'hardware':h,'operations':{'wall_rank_mean':t,'ledger':ledger,'mpi_owners':[]}}
  s=summarize(result)
  self.assertEqual(s['decision'],'review_operation_evidence')
  q=s['comparisons']['Stage203 minus legacy']
  self.assertAlmostEqual(sum(r['excess_seconds'] for r in q['ranking']),.35)
  result['runs']['candidate-operations-1']['ordinary_step_median']=1.2
  self.assertEqual(summarize(result)['decision'],'observer_effects_require_review')
  del result['runs']['candidate-operations-2']
  with self.assertRaises(ValueError):summarize(result)
 def test_schedule(self):
  jobs=schedule({'n83-s2':{}});self.assertEqual(len(jobs),20)
  self.assertEqual([v[1] for v in jobs if v[2]=='operations'],['baseline_profile','candidate_profile','legacy_profile','legacy_profile','candidate_profile','baseline_profile'])
  self.assertEqual(sum(v[2]=='validation' for v in jobs),2)
if __name__=='__main__':unittest.main()
