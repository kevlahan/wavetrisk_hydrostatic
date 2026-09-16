import json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
from placement_200 import plan,verify
from cluster_200 import schedule,launch
from analyze_scaling_200 import collect,summarize,PHASES
BASE=Path(__file__).resolve().parent

class Protocol(unittest.TestCase):
    def topology(self):
        rows=[]
        for s in range(2):
            for c in range(128):
                rows.append(dict(cpu=s*128+c,core=c,socket=s,node=s,llc=list(range(s*128+(c//8)*8,s*128+(c//8+1)*8))))
        return {'host':'bb13','rows':rows}
    def test_nested_cache_spread_and_socket_pair(self):
        p=plan(self.topology())
        for n in (21,42):self.assertEqual(p[f'n{n}-s1']['mapping'],p['n83-s1']['mapping'][:n])
        self.assertEqual([r['cpu'] for r in p['n21-s1']['mapping'][:4]],[0,8,16,24])
        self.assertEqual(sum(r['socket']==0 for r in p['n83-s2']['mapping']),42)
        self.assertEqual(sum(r['socket']==1 for r in p['n83-s2']['mapping']),41)
        for case in p.values():
            self.assertEqual(len(set(r['cpu'] for r in case['mapping'])),case['ranks'])
            rows=[dict(rank=r['rank'],hostname='bb13',affinity=[r['cpu']]) for r in case['mapping']]
            self.assertTrue(verify(rows,case));rows[0]['affinity']=[999]
            with self.assertRaises(ValueError):verify(rows,case)
    def test_reject_smt_incomplete_and_unsupported_topology(self):
        t=self.topology();t['rows'][1]['core']=0
        with self.assertRaises(ValueError):plan(t)
        t=self.topology();t['rows']=t['rows'][:128]
        with self.assertRaises(ValueError):plan(t)
    def test_schedule_and_step_resource_semantics(self):
        p=plan(self.topology());jobs=schedule(p);self.assertEqual(len(jobs),42)
        for case in p:
            for side in ('candidate','legacy'):
                bare=[j for j in jobs if j[:2]==(case,side) and j[2]=='bare'];self.assertEqual([j[3] for j in bare],[1,2])
                self.assertEqual(sum(j[:3]==(case,side+'_profile','probe') for j in jobs),1)
        cmd=launch(Path('/tmp/map.json'),21,'--probe')
        self.assertIn('--ntasks=21',cmd);self.assertIn('--cpu-bind=none',cmd)
        self.assertFalse(any(x.startswith(('--exclusive','--exact','--cpus-per-task')) for x in cmd))

class Reporting(unittest.TestCase):
    def test_preserve_bare_timings_and_observer_bias(self):
        cases={name:{'ranks':n} for name,n in [('n21-s1',21),('n42-s1',42),('n83-s1',83),('n83-s2',83)]}
        result={'cases':cases,'runs':{}}
        for name,case in cases.items():
            for side,baseline in [('candidate',1.25),('legacy',.56)]:
                t=baseline*83/case['ranks']
                for phase in ('bare','coarse','probe'):
                    for repeat in (range(2) if phase=='bare' else range(1)):
                        times=t*(1.10 if phase=='probe' and side=='legacy' else 1.)
                        cells=[]
                        for ph in PHASES:
                            q=times/len(PHASES)
                            cells.append({'phase':ph,'outside_mpi_mean':q,'outside_mpi_max':q,'mpi_mean':0.,'wall_mean':q,
                              'phase_on_full_step_slowest_rank':q,'work35_global':1.,'outside_mpi_max_over_mean':1.,
                              'rank_rows':[{'counts':[1]*42,'mpi_calls':3}]})
                        result['runs'][f'{name}-{side}-{phase}-{repeat}']={'case':name,'side':side,'phase':phase,'ordinary_step_median':times,'scaling':{'cells':cells}}
        report=summarize(result)
        self.assertAlmostEqual(report['cases']['n83-s1']['gap_seconds'],.69)
        self.assertFalse(report['optimization_ready'])
        self.assertTrue(report['cases']['n83-s1']['observer']['legacy']['probe']['review_flags'])
        self.assertIsNone(report['cases']['n83-s1']['observer']['legacy']['coarse']['parents']['dynamics WT']['mpi_mean'])
        result['runs'].pop('n83-s1-legacy-bare-1')
        with self.assertRaises(ValueError):summarize(result)

@unittest.skipUnless(shutil.which('gfortran') and shutil.which('cc'),'Fortran/C compilers required')
class Observer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.root=Path(cls.tmp.name)
        (cls.root/'driver.f90').write_text('''program driver
use workflow_profile_195_mod
use iso_fortran_env, only: int64
integer :: step,p
call wf_init
 do step=1,2
 call wf_begin
 if(step==2) call wf_flag(1)
 do p=1,5
 call wf_phase(p)
 if(p==1) then
 call wf_enter(1)
 call op_enter(18)
 call op_enter(19)
 call op_count(35,6_int64)
 call op_leave(19)
 call op_leave(18)
 call op_enter(105)
 call op_leave(105)
 call wf_leave(1)
 call wf_enter(2)
 leaf199(31)=leaf199(31)+3_int64
 call wf_leave(2)
 endif
 enddo
 call wf_end
 enddo
call wf_finish(0)
end program
''')
        subprocess.run(['cc','-c',str(BASE/'workflow_clock_200.c'),'-o','clock.o'],cwd=cls.root,check=True)
        subprocess.run(['gfortran','-ffree-line-length-none','-fcheck=all',str(BASE/'workflow_profile_200.f90'),'driver.f90','clock.o','-o','probe'],cwd=cls.root,check=True)
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_modes_conservation_nested_mpi_and_counters(self):
        for mode in (0,1):
            root=self.root/str(mode);root.mkdir();env={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
            env.update(WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200=str(mode))
            subprocess.run([str(self.root/'probe')],cwd=root,env=env,check=True)
            a=collect(root,1,mode);self.assertEqual(a['ordinary_steps'],[1])
            c=next(c for c in a['cells'] if c['phase']=='dynamics WT')
            self.assertEqual(c['rank_rows'][0]['counts'][34],6*mode)
            self.assertEqual(c['rank_rows'][0]['mpi_calls'],2*mode)
            self.assertFalse(list(root.glob('batches*')))
            f=root/'scaling200-rank-0.txt';original=f.read_text();f.write_text(original.replace('self 1 0 0 0','self 1 0 0 0',1)+'\nself 1 0 0 0 0.0\n')
            with self.assertRaises(ValueError):collect(root,1,mode)
            f.write_text(original)
            with self.assertRaises(ValueError):collect(root,2,mode)
if __name__=='__main__':unittest.main()
