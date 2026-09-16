import json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from analyze_hardware_201 import read_rank,parse_samples,collect,EVENTS
from cluster_201 import schedule,launch
from rank_launch_201 import perf_command
BASE=Path(__file__).resolve().parent

class Protocol(unittest.TestCase):
    def test_reversed_pairs_fixed_placement(self):
        jobs=schedule({'n83-s2':{}});self.assertEqual(len(jobs),20)
        self.assertTrue(all(j[0]=='n83-s2' for j in jobs))
        for phase in ('bare','core','memory'):
            for repeat in (1,2):self.assertEqual(len([j for j in jobs if j[2]==phase and j[3]==repeat]),2)
        cmd=launch(Path('/tmp/map'),83,'--probe');self.assertIn('--ntasks=83',cmd)
        self.assertFalse(any(x.startswith(('--exact','--exclusive','--cpus-per-task')) for x in cmd))
        cmd=perf_command(Path('/tmp/perf.data'),['solver']);self.assertIn('--no-inherit',cmd)
        self.assertIn('mono',cmd);self.assertIn('dwarf,8192',cmd);self.assertNotIn('-F',cmd)
    def test_parser_preserves_intervals_counts_and_rejects_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'hardware201-rank-0.txt'
            text='schema 1 1 0 123\nevents '+' '.join(EVENTS[1])+'\ninterval 1 1 100 200 80 100 99 200 300 20 10\n'
            f.write_text(text);r=read_rank(f);self.assertEqual(r['intervals'][0]['counts'],[200,300,20,10])
            for bad in (text+'interval 1 0 190 210 20 20 20 40 60 4 2\n',text.replace('100 99','100 101'),text.replace('instructions','other')):
                f.write_text(bad)
                with self.assertRaises(ValueError):read_rank(f)
    def test_sample_header_stacks_and_lost(self):
        text=' 123 234.000000001: cycles:u: abc solver (climateJ5)\n    abc solver (climateJ5)\n    def parent (climateJ5)\n\n 123 234.100000001: instructions:u:\n    abc solver (climateJ5)\n    def parent (climateJ5)\n'
        rows,lost=parse_samples(text);self.assertFalse(lost);self.assertEqual(len(rows),2)
        self.assertEqual(rows[0]['ns'],234000000001);self.assertEqual(rows[1]['event'],'instructions:u')
        self.assertEqual(rows[1]['frames'][0],('solver','climateJ5'))
        self.assertTrue(parse_samples(text+'PERF_RECORD_LOST 10')[1])

@unittest.skipUnless(shutil.which('gfortran') and shutil.which('cc'),'Compilers needed')
class Observer(unittest.TestCase):
    def test_compiled_contexts_and_restart_filter(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'main.f90').write_text('''program main
use workflow_profile_195_mod
integer :: step,i
call wf_init
 do step=1,2
 call wf_begin
 if(step==2) call wf_flag(1)
 do i=1,5
 call wf_phase(i)
 if(i==1) then
 call wf_enter(1)
 call wf_leave(1)
 call wf_enter(2)
 call wf_leave(2)
 endif
 enddo
 call wf_end
 enddo
call wf_finish(0)
end program
''')
            subprocess.run(['cc','-std=c11','-Wall','-Wextra','-Werror','-c',str(BASE/'workflow_clock_201.c'),'-o','clock.o'],cwd=p,check=True)
            subprocess.run(['gfortran','-ffree-line-length-none','-fcheck=all',str(BASE/'workflow_profile_201.f90'),'main.f90','clock.o','-o','probe'],cwd=p,check=True)
            env={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')};env.update(WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200='0',WAVETRISK_PMU201='0')
            subprocess.run([str(p/'probe')],cwd=p,env=env,check=True)
            a=collect(p,1,0);self.assertEqual({c['step'] for c in a['cells']},{1});self.assertFalse(a['quality_flags'])
            self.assertIsNone(a['phases']['dynamics WT']['ipc'])
            # Inject known counts in the real observed intervals to test PMU coverage logic.
            f=p/'hardware201-rank-0.txt';lines=f.read_text().splitlines();lines[0]=lines[0].replace('schema 1 0','schema 1 1');lines[1]='events '+' '.join(EVENTS[1])
            for i in range(2,len(lines)):
                row=lines[i].split();row[6:]=['100','99','200','300','20','10'];lines[i]=' '.join(row)
            f.write_text('\n'.join(lines)+'\n');a=collect(p,1,1);self.assertFalse(a['quality_flags']);self.assertEqual(a['phases']['whole step']['ipc'],1.5)
            f.write_text(f.read_text().replace('100 99','100 50'));self.assertTrue(collect(p,1,1)['quality_flags'])

class Attribution(unittest.TestCase):
    def test_sample_timestamp_filter_and_scope(self):
        from analyze_hardware_201 import samples
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'perf-rank-0.txt').write_text('42 1.000000020: cycles:u:\n  abc expand+0x12 (climateJ5)\n  def caller (climateJ5)\n\n42 2.000000020: cycles:u:\n  abc restart (climateJ5)\n\n42 3.000000000: cycles:u:\n  abc outside (climateJ5)\n')
            timeline={'pid':42,'intervals':[dict(step=1,context=1,start=1000000000,end=1100000000),dict(step=2,context=1,start=2000000000,end=2100000000)]}
            with patch('analyze_hardware_201.read_rank',return_value=timeline),patch('analyze_hardware_201.common',return_value={1:{'restart':0},2:{'restart':1}}):
                a=samples(root,[0])['ranks']['0']
            self.assertEqual(a['leaf'],[dict(phase='dynamics WT',event='cycles:u',symbol='expand',dso='climateJ5',count=1)])
            self.assertEqual(a['outside_step_samples'],1)
    def test_cycle_accounting_and_calibration_remain_separate(self):
        from analyze_hardware_201 import summarize,PHASES
        runs={}
        for side,t,instructions,cycles in [('candidate',1.0,200.,800.),('legacy',.5,100.,200.)]:
            for phase in ('bare','coarse','core','memory'):
                for rep in ((1,2) if phase!='coarse' else (0,)):
                    label=f'{side}-{phase}-{rep}';r={'label':label,'side':side,'phase':phase,'round':rep,'ordinary_step_median':t}
                    if phase!='bare':
                        mode=0 if phase=='coarse' else 1 if phase=='core' else 2
                        parents={p:{'wall_seconds_rank_mean':t/7,'counts_rank_mean_per_step':dict(zip(EVENTS[mode],[cycles,instructions,20.,10.]))} for p in ('whole step',*PHASES)}
                        r['hardware']={'phases':parents,'quality_flags':[]}
                    runs[label]=r
        a=summarize({'runs':runs});r=a['counter_comparisons']['core'][0];parts=r['symmetric_cycle_accounting']
        self.assertEqual(parts['instruction_count_component']+parts['cycles_per_instruction_component'],600.)
        self.assertFalse(a['optimization_ready']);self.assertEqual(r['instruction_ratio'],2)

class Collection(unittest.TestCase):
    def test_archive_keeps_raw_and_exported_samples_excludes_solver(self):
        import tarfile
        from collect_cluster_logs_201 import collect as archive
        with tempfile.TemporaryDirectory() as d:
            home=Path(d);root=home/'results';root.mkdir();tree=home/'tree';tree.mkdir()
            (root/'results.json').write_text(json.dumps({'slurm_environment':{'SLURM_JOB_ID':'123'}}))
            (root/'perf-rank-0.data').write_bytes(b'raw sample fixture')
            (root/'perf-rank-0.txt').write_text('symbol export')
            (root/'climateJ5').write_bytes(b'executable fixture')
            with tarfile.open(archive(201,'123',home,root,tree)) as t:
                self.assertIn('results/perf-rank-0.data',t.getnames())
                self.assertIn('results/perf-rank-0.txt',t.getnames())
                self.assertNotIn('results/climateJ5',t.getnames())
                inventory=json.load(t.extractfile('collection.json'))
                self.assertEqual(inventory['files']['results/perf-rank-0.data']['bytes'],18)

if __name__=='__main__':unittest.main()
