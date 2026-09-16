import copy,csv,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_operations_198 import collect,summarize,NR
from analyze_workflow_195 import collect as workflow
from cluster_198 import schedule,operation_mode,write_rankings
DRIVER='''program probe
use workflow_profile_195_mod
use iso_fortran_env, only:int64
implicit none
integer::s,i
logical::take
call wf_init
do s=1,2
 call wf_begin
 call wf_phase(1)
 call op_enter(2)
 call wf_enter(2)
 call op_enter(26)
 if(count198) leaf198(35)=leaf198(35)+11_int64
 call op_leave(26)
 call wf_leave(2)
 call wf_enter(1)
 call op_enter(65)
 call op_count(29,1_int64)
 call op_dependency(1,1,1,5)
 call op_post(2,100,100,19)
 call op_wait_begin
 call op_enter(19)
 call op_leave(19)
 call op_wait_end
 call op_leave(65)
 call wf_leave(1)
 call op_leave(2)
 call wf_phase(2)
 call op_enter(90)
 call op_leave(90)
 call wf_phase(3)
 call op_enter(4)
 call wf_enter(1)
 call op_enter(65)
 call op_count(36,7_int64)
 call op_leave(65)
 call wf_leave(1)
 call op_enter(83)
 do i=1,333
  take=sample_begin198(1)
  call sample_end198(1,take)
 end do
 call op_leave(83)
 call op_leave(4)
 call wf_phase(4)
 call wf_flag(2)
 call op_enter(89)
 call op_leave(89)
 call wf_phase(5)
 if(s==2)call wf_flag(1)
 call wf_end
end do
call wf_finish(0)
end program
'''
class Operations198(unittest.TestCase):
 @unittest.skipUnless(shutil.which('cc'),'requires C compiler')
 def test_monotonic_clock_stress(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d)
   (root/'stress.c').write_text('extern double workflow_clock198(void);\nint main(void) {double previous=workflow_clock198(); for(int i=0;i<1000000;i++){double t=workflow_clock198(); if(t<previous || t<0) return 1; previous=t;} return 0;}\n')
   subprocess.run(['cc','-O2','-Wall','-Wextra','-Werror',str(Path(__file__).with_name('workflow_clock_198.c')),str(root/'stress.c'),'-o',str(root/'stress')],check=True)
   subprocess.run([str(root/'stress')],check=True)
 def test_schedule(self):
  jobs=schedule(2,True);self.assertEqual(len(jobs),27)
  self.assertEqual(len([j for j in jobs if j[4]=='timing-bare']),6)
  self.assertEqual([j[0] for j in jobs if j[4]=='diagnostic-timers'],['candidate-timers-1','legacy-timers-1','legacy-timers-2','candidate-timers-2'])
  self.assertEqual([operation_mode(j[0],j[4]) for j in jobs[-4:]],[0,3,3,0])
 @unittest.skipUnless(shutil.which('gfortran'),'requires gfortran')
 def test_compiled_full_step_modes_and_corruption(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);(root/'driver.f90').write_text(DRIVER)
   subprocess.run(['cc','-O2','-Wall','-Wextra','-Werror','-c',str(Path(__file__).with_name('workflow_clock_198.c')),'-o',str(root/'clock.o')],check=True)
   for opt in ('-O0','-O2'):
    subprocess.run(['gfortran',opt,'-ffp-contract=off','-fcheck=all','-Wall','-Wextra','-Werror',str(Path(__file__).with_name('workflow_profile_198.f90')),str(root/'driver.f90'),str(root/'clock.o'),'-o',str(root/'probe')],cwd=root,check=True,capture_output=True)
    for mode in range(4):
     run=root/(opt+str(mode));run.mkdir()
     env={**os.environ,'WAVETRISK_WORKFLOW195':'1','WAVETRISK_OPERATIONS198':str(mode),'WAVETRISK_SAMPLE_OFFSET198':'127','WAVETRISK_PROFILE_BLOCK_DETAIL':'0','WAVETRISK_PROFILE_PARALLEL_BLOCKS':'0'}
     subprocess.run([str(root/'probe')],cwd=run,env=env,check=True)
     workflow(run,1,1);r=collect(run,1,mode)
     self.assertEqual(len(r['steps']),2);self.assertEqual(r['steps'][1]['restart'],1)
     self.assertEqual(r['work_per_rank_step'][2][34],11 if mode>=2 else 0)
     self.assertEqual(r['work_per_rank_step'][3][35],7 if mode>=2 else 0)
     self.assertEqual(r['dependency_records'],2 if mode>=2 else 0)
     self.assertEqual(sum(v['calls'] for v in r['storage']),333 if mode==3 else 0)
     if mode:self.assertAlmostEqual(sum(v['mean'] for v in r['ledger']),r['ordinary_wall_mean'],places=6)
     f=run/'operations198-rank-0.txt';raw=f.read_text();lines=raw.splitlines();idx=next(i for i,l in enumerate(lines) if l.startswith('self 1 4 83 '));parts=lines[idx].split();parts[-1]='100';lines[idx]=' '.join(parts);f.write_text('\n'.join(lines)+'\n')
     with self.assertRaises(ValueError):collect(run,1,mode)
     f.write_text(raw)
     with self.assertRaises(ValueError):collect(run,2,mode)
     f.write_text('\n'.join(l for l in raw.splitlines() if not l.startswith('self 1 4 83 '))+'\n')
     with self.assertRaises(ValueError):collect(run,1,mode)
     f.write_text(raw)
     # The complete ledger includes setup, adaptation, physics, remap and output,
     # and signed differences conserve the same statistic, never a median.
     if mode==1:
      runs={}
      for side,timing in [('baseline',2.),('candidate',2.),('legacy',1.)]:
       runs[side+'-1']={'kind':side,'phase':'timing-bare','ordinary_step_median':timing}
      for side in ('candidate','legacy'):
       for name in ('minimal','timers','work'):
        for n in (1,2):runs[f'{side}-{name}-{n}']={'kind':side+'_profile','phase':'diagnostic-'+name,'ordinary_step_median':2. if side=='candidate' else 1.,'operations':copy.deepcopy(r)}
      out=summarize(runs);self.assertAlmostEqual(sum(v['excess_seconds'] for v in out['ranking']),0.)
      self.assertFalse(out['optimization_ready']);self.assertIn('review_flags',out['calibration']['legacy']['timers'])
      # JSON round trips must preserve the ability to regenerate rankings.
      self.assertEqual(summarize(json.loads(json.dumps(runs)))['instrumented_mean_excess'],0.)
      for side,index in [('legacy',0),('candidate',101)]:
       for n in (1,2):
        op=runs[f'{side}-timers-{n}']['operations']
        op['ledger'][index]['mean']+=0.01;op['ordinary_wall_mean']+=0.01
      ranked=summarize(runs)
      self.assertAlmostEqual(ranked['ranking'][0]['excess_seconds'],0.01)
      self.assertAlmostEqual(ranked['ranking'][-1]['excess_seconds'],-0.01)
      self.assertAlmostEqual(sum(v['excess_seconds'] for v in ranked['parent_ranking']),0.)
      write_rankings(run,ranked,names={str(i):'operation '+str(i) for i in range(101)})
      with (run/'ranking.csv').open() as stream:csvrows=list(csv.DictReader(stream))
      self.assertEqual(len(csvrows),909)
      self.assertLess(float(csvrows[-1]['excess_seconds']),0.)
      self.assertIn('Full-step accounting',(run/'ranking.md').read_text())

if __name__=='__main__':unittest.main()
