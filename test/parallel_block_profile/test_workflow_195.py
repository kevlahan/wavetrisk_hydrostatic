import copy,json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_workflow_195 import collect,summarize
from cluster_195 import schedule,observer_mode
DRIVER='''program probe
use workflow_profile_195_mod
implicit none
integer :: s,p
call wf_init
do s=1,2
 call wf_begin
 do p=1,5
  call wf_phase(p)
  if(p==1) then
   call wf_enter(1)
   call wf_leave(1)
  end if
  if(p==4) call wf_flag(2)
  if(p==5.and.s==2) call wf_flag(1)
 end do
 call wf_end
end do
call wf_finish(0)
end program
'''

class Workflow195(unittest.TestCase):
 def test_schedule(self):
  jobs=schedule(6,True)
  self.assertEqual(len(jobs),58)
  labels=[j[0] for j in jobs];self.assertEqual(len(labels),len(set(labels)))
  for mode in ('bare','minimal','coarse','detail'):
   selected=[j for j in jobs if j[4]=='timing-'+mode]
   self.assertEqual(len(selected),12)
   for i in range(0,12,2):self.assertEqual({j[0].split('-')[0] for j in selected[i:i+2]},{'block','legacy'})
  self.assertEqual(observer_mode('legacy-detail-1','timing-detail'),2)
  self.assertEqual(observer_mode('block-coarse-1','timing-coarse'),1)
  self.assertEqual(observer_mode('candidate-warmup','warmup'),0)

 @unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
 def test_actual_observer_modes_and_malformed_evidence(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);(root/'driver.f90').write_text(DRIVER)
   source=Path(__file__).with_name('workflow_profile_195.f90')
   for optimization in ('-O0','-O2'):
    subprocess.run(['gfortran',optimization,'-fcheck=all','-Wall','-Wextra','-Werror','-ffp-contract=off',str(source),str(root/'driver.f90'),'-o',str(root/'probe')],cwd=root,check=True,capture_output=True)
    for mode in (0,1,2):
     directory=root/(optimization+str(mode));directory.mkdir()
     subprocess.run([str(root/'probe')],cwd=directory,env={**os.environ,'WAVETRISK_WORKFLOW195':str(mode)},check=True)
     r=collect(directory,1,mode)
     self.assertEqual(len(r['steps']),2)
     self.assertEqual(r['steps'][1]['restart'],1)
     if mode: self.assertEqual(r['steps'][0]['transaction_mean']['1']['calls_global'],1)
     path=directory/'workflow195-rank-0.txt';raw=path.read_text()
     with self.assertRaises(ValueError):collect(directory,2,mode)
     with self.assertRaises(ValueError):collect(directory,1,(mode+1)%3)
     path.write_text(raw.replace('phase','badphase',1))
     with self.assertRaises(ValueError):collect(directory,1,mode)
     path.write_text(raw+'schema 1 0 0\n')
     with self.assertRaises(ValueError):collect(directory,1,mode)
     path.write_text(raw)
     if mode==2:
      lines=raw.splitlines();lines=[l for l in lines if not l.startswith('event')]
      path.write_text('\n'.join(lines)+'\n')
      with self.assertRaises(ValueError):collect(directory,1,mode)
   for selector in ('3','x','10'):
    result=subprocess.run([str(root/'probe')],cwd=root,env={**os.environ,'WAVETRISK_WORKFLOW195':selector},capture_output=True)
    self.assertNotEqual(result.returncode,0)

 def test_paired_summary_conservation_and_overhead(self):
  runs={}
  for mode in ('bare','minimal','coarse','detail'):
   for side in ('block','legacy'):
    factor=2 if side=='block' else 1
    phases={n:factor/6 for n in ('setup','advance','physics','adapt_reconstruct','remap','output_restart')}
    runs[mode+side]={'side':side,'round':1,'phase':'timing-'+mode,'elapsed_seconds':factor*10,'ordinary_step_median':factor,
      'workflow':{'ordinary_wall_max_median':factor,'ordinary_rank_mean_wall':factor,'ordinary_phase_mean':phases,
                  'steps':[{'sequence':1,'restart':0,'remap':0,'transaction_mean':{'1':{'calls_global':5},'2':{'calls_global':4}}}]}}
  r=summarize(runs)
  self.assertEqual(r['observer_increase_percent']['legacy']['coarse'],0)
  self.assertAlmostEqual(r['paired_differences']['coarse'][0]['residual'],0)
  runs['coarseblock']['workflow']['ordinary_phase_mean']['advance']+=1
  with self.assertRaisesRegex(ValueError,'reconcile'):summarize(runs)

if __name__=='__main__':unittest.main()
