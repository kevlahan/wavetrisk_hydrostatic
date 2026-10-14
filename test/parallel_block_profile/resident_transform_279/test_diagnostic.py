import contextlib,io,json,os,shutil,subprocess,tarfile,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import common,diagnostic,run as runner,build as builder,control
HERE=Path(__file__).resolve().parent

def profile(exchange=.1):
    return dict(regions={str(i):dict(name=str(i),mean_on_slowest_step_rank_seconds=exchange if i==78 else .001) for i in range(1,85)},
                remap_groups={'remap':dict(step_indices=[4,9])})
def measurements(exchange=.1):
    return [dict(ranks=n,kind=k,profile=p,phase='measurement',status='passed',oracle=False,catalogue=False,
                 ordinary=dict(median=.5 if k=='main' else .75 if k=='stage266' else .8,
                  included=[dict(index=i,seconds=.8) for i in range(12) if i not in (0,7,8)]),
                 resources={},attribution=profile(exchange)) for n,k,p in diagnostic.schedule()]

class Diagnostic(unittest.TestCase):
    def test_bounded_symmetric_schedule(self):
        rows=diagnostic.schedule();self.assertEqual(len(rows),40)
        for n in diagnostic.RANKS:
            seq=[k for r,k,p in rows if r==n and k in ('previous','candidate') and p=='off']
            self.assertEqual(seq,['previous','candidate','candidate','previous','candidate','previous','previous','candidate'])
            self.assertEqual(sum(r==n and k=='stage266' for r,k,p in rows),2)
            self.assertEqual(sum(r==n and k=='main' for r,k,p in rows),2 if n<=160 else 0)
        self.assertEqual(sum(p=='coarse' for n,k,p in rows),6)
    def test_hostile_environment_cannot_enable_oracle_or_fine(self):
        with patch.dict(os.environ,{'WAVETRISK_PROFILE_RHS_FINE':'1','WAVETRISK_VALIDATE_BLOCK_DYNAMICS':'1','WAVETRISK_CONNECTED_RHS_AUDIT':'1'}):
            e=diagnostic.environment('previous');self.assertEqual(e['WAVETRISK_PROFILE_RHS_FINE'],'0');self.assertEqual(e['WAVETRISK_VALIDATE_BLOCK_DYNAMICS'],'0');self.assertEqual(e['WAVETRISK_CONNECTED_RHS_AUDIT'],'0');self.assertEqual(e['WAVETRISK_CONNECTED_RHS'],'1')
            self.assertEqual(diagnostic.environment('candidate','coarse')['WAVETRISK_PROFILE_BLOCK_DETAIL'],'1')
        for args in [('candidate','fine',True),('checked','fine'),('stage266','coarse'),('candidate','unknown')]:
            with self.assertRaises(ValueError):diagnostic.environment(*args)
    def test_no_fine_campaign_for_any_exchange_budget(self):
        for exchange in (.001,.1,10.):
            report=diagnostic.summarize_diagnostic(measurements(exchange))
            self.assertNotIn('fine_decisions',report)
            self.assertTrue(all(p!='fine' for n,k,p in diagnostic.schedule()))
    def test_coarse_previous_and_repair_have_identical_selectors(self):
        a=diagnostic.environment('previous','coarse');b=diagnostic.environment('candidate','coarse')
        self.assertEqual({k:v for k,v in a.items() if k!='WAVETRISK_RESIDENT_TRANSFORM'}, {k:v for k,v in b.items() if k!='WAVETRISK_RESIDENT_TRANSFORM'})
    def test_mean_totals_and_medians_remain_distinct(self):
        rows=measurements()
        for r in rows:
            if r['kind']=='candidate' and r['profile']=='off':r['ordinary']['included'][-1]['seconds']=8.
        s=diagnostic.summarize_diagnostic(rows)['by_ranks']['160']
        self.assertEqual(s['candidate_minus_previous_percent'],0.)
        self.assertGreater(s['ordinary_mean_candidate_minus_previous_percent'],0.)
        self.assertEqual(len(s['variants']['candidate-off']['ordinary_run_totals_seconds']),4)
    def test_missing_oracle_and_reordered_measurements_fail(self):
        for mutation in ('missing','oracle','reorder'):
            r=measurements()
            if mutation=='missing':r.pop()
            elif mutation=='oracle':r[-1]['oracle']=True
            else:r[0],r[1]=r[1],r[0]
            with self.assertRaises(ValueError):diagnostic.summarize_diagnostic(r)
    def test_reports_absolute_scaling_overhead_and_subgroups(self):
        s=diagnostic.summarize_diagnostic(measurements())
        self.assertAlmostEqual(s['by_ranks']['160']['gap_to_stage266_seconds'],.05)
        self.assertEqual(s['scaling']['candidate']['relative_efficiency_percent'],62.5)
        self.assertEqual(s['by_ranks']['83']['candidate_minus_previous_percent'],0)
        self.assertEqual(len(s['by_ranks']['83']['ordinary_subgroups']['candidate']['remap']['run_medians_seconds']),4)
        self.assertFalse(s['automatic_promotion']);self.assertEqual(s['performance_acceptance'],'requires_review')
    def test_scope_rejects_fallback_repair_and_control_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for n in ('IDENTITY.json','source-files.json','parent-source-files.json','overlay-files.json'):shutil.copyfile(HERE/n,root/n)
            manifest=json.loads((root/'source-files.json').read_text())
            common.save(root/'LOCAL_VALIDATION.json',dict(status='passed',full_source_sha256=manifest['candidate']))
            self.assertEqual(len(common.source_scope(root)),len(json.loads((HERE/'overlay-files.json').read_text())))
            for family,name in [('candidate','src/parallel_block_rhs_routes.f90'),('previous','Makefile')]:
                mutated=json.loads((HERE/'source-files.json').read_text());mutated[family][name]='wrong';common.save(root/'source-files.json',mutated)
                with self.assertRaises(ValueError):common.source_scope(root)
    def test_unsealed_preflight_rejected(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(common,'verify'),patch.object(common,'verify_packaged_grid'):
            with self.assertRaisesRegex(ValueError,'Unsealed'):common.preflight(tmp,tmp)
    def test_build_five_binaries_serial_shared_physics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);kit=root/'kit';out=root/'out';src=root/'src';kit.mkdir();out.mkdir();src.mkdir();manifests={}
            for kind in ('stage266','main','previous','candidate'):
                tree=src/kind;(tree/'src/physics/simple_physics/phyparam').mkdir(parents=True);(tree/'Makefile').write_text(kind);manifests[kind]={'Makefile':common.sha(tree/'Makefile')}
            with tarfile.open(kit/'sources.tar.gz','w:gz') as t:
                for p in src.iterdir():t.add(p,arcname=p.name)
            common.save(kit/'source-files.json',manifests);common.save(kit/'IDENTITY.json',{});common.save(kit/'files-sha256.json',{})
            def invoke(cmd,**kw):
                if cmd[0]=='make':
                    t=Path(kw['cwd']);b=next(s.split('=',1)[1] for s in cmd if s.startswith('BIN_DIR='));(t/b).mkdir();(t/b/'climate').write_text(b)
                    obj=t/'src/physics/simple_physics/phyparam/obj';obj.mkdir(exist_ok=True)
                    for i in range(23):(obj/f'{i}.o').write_text(str(i))
                return subprocess.CompletedProcess(cmd,0)
            with patch.object(builder.subprocess,'run',side_effect=invoke),patch.object(builder.subprocess,'check_output',return_value='stub compiler'):
                b=builder.build(kit,out)
            self.assertEqual(set(b['builds']),{'stage266','main','previous','candidate','checked'})
            for k,r in b['builds'].items():self.assertIn('-j1',r['command']);self.assertEqual(r['physics_objects'],b['physics_objects']);self.assertEqual('DEBUG=check' in r['command'],k=='checked')
    def test_failed_job_collection_and_duplicate_submission_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            home=Path(tmp);kit=home/'kit';kit.mkdir();calls=[]
            def invoke(cmd,**kw):calls.append(cmd);return subprocess.CompletedProcess(cmd,0,'','')
            with patch.object(Path,'home',return_value=home),patch.object(control,'preflight'),patch.object(control.shutil,'which',return_value='/stub'),patch.object(control.subprocess,'check_output',return_value='12345\n'),patch.object(control.subprocess,'run',side_effect=invoke),patch.object(control,'capture',return_value={'job':None}):
                control.submit(kit);state=json.loads(control.state_path().read_text());Path(state['log']).write_text('failed before first solver')
                with self.assertRaises(ValueError):control.submit(kit)
                control.collect(state);archive=home/'rk279-review-12345.tar.gz';h=common.sha(archive);self.assertGreater(control.verify_archive(archive),0)
                control.collect(state);self.assertEqual(common.sha(archive),h)
                self.assertEqual(sum(c[0]=='scp' for c in calls),2)
                with tarfile.open(archive) as t:self.assertEqual(t.extractfile('review/slurm.log').read(),b'failed before first solver')

class Runner(unittest.TestCase):
    def exercise(self,exchange,fail=False):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);kit=root/'kit';fixture=root/'fixture';out=root/'out';kit.mkdir();fixture.mkdir();(kit/'grids').mkdir();(kit/'files-sha256.json').write_text('{}')
            for n in ('simple.in','test_checkpoint_0003.bin.zst'):(fixture/n).write_text('fixture')
            (kit/'grids/HR95JT_WT_004').write_text('grid')
            top=dict(host='bb13',rows=[dict(cpu=c,core=c%128,socket=c//128,node=c//128,llc=list(range(c//8*8,c//8*8+8))) for c in range(256)]);cases=runner.plan(top)
            def build(k,o):
                b={}
                for kind in ('main','stage266','previous','candidate','checked'):
                    p=o/kind;p.write_text(kind);b[kind]=dict(binary=kind,binary_sha256=common.sha(p))
                return dict(host='bb13',builds=b)
            def invoke(cmd,**kw):
                n=int(next(c.split('=')[1] for c in cmd if c.startswith('--ntasks=')))
                records=[dict(rank=i,hostname='bb13',affinity=[cases[n]['mapping'][i]['cpu']],returncode=0,resource_collector='os.wait4',max_rss_kib=100,major_faults=0,minor_faults=1,user_cpu_seconds=1.,system_cpu_seconds=.1) for i in range(n)]
                if '--probe' in cmd:return subprocess.CompletedProcess(cmd,0,'\n'.join(json.dumps(r) for r in records),'')
                d=Path(kw['cwd']);env=kw['env']
                for i,r in enumerate(records):common.save(d/f'rank-{i}.json',r)
                log='Stage203 inverse workspace = on\nStage203 scalar fields = on\n'+f'Resident block execution: MPI tasks = {n}, geometric Domains = 160\n'
                if env['WAVETRISK_CONNECTED_RHS']=='1':log+='Stage278 connected resident RHS = on\n'
                if env['WAVETRISK_RESIDENT_TRANSFORM']=='1':log+='Stage279 resident provisional transforms = on\n'
                if env['WAVETRISK_VALIDATE_BLOCK_DYNAMICS']=='1':
                    for x in ('dynamics','adaptation','vertical-remap'):log+=f'Block {x} validation oracle enabled\n'
                if env['WAVETRISK_PROFILE_BLOCK_DETAIL']=='1':log+='Parallel-block profiling summary (inclusive):\n'
                kw['stdout'].write(log);kw['stdout'].flush();(d/'test_log').write_text(d.name)
                return subprocess.CompletedProcess(cmd,7 if fail else 0)
            def sha(p):
                p=Path(p)
                if p.name=='test_checkpoint_0004.bin.zst':return common.CHECKPOINTS['main' if '-main-' in p.parent.name else 'stage266']
                return common.sha(p)
            def timing(log,a,name):
                t=.5 if '-main-' in name else .75 if '-stage266-' in name else .8
                return dict(median=t,included=[dict(index=i,seconds=t) for i in range(12) if i not in (0,7,8)],excluded=[dict(index=i,seconds=t) for i in (0,7,8)])
            def attribute(d):common.save(d/'attribution.json',profile(exchange))
            with patch.dict(os.environ,dict(SLURM_JOB_ID='12345',SLURM_JOB_PARTITION='nk',SLURM_JOB_NUM_NODES='1')),patch.object(runner,'preflight',return_value={}),patch.object(runner.socket,'gethostname',return_value='bb13'),patch.object(runner.shutil,'which',return_value='/stub'),patch.object(runner,'build',side_effect=build),patch.object(runner,'topology',return_value=top),patch.object(runner.subprocess,'run',side_effect=invoke),patch.object(runner,'sha',side_effect=sha),patch.object(runner,'J5_HASHES',{n:common.sha(fixture/n) for n in ('simple.in','test_checkpoint_0003.bin.zst')}),patch.object(runner,'coverage',return_value=(dict(step_records=[dict(state=str(i)) for i in range(12)]),dict(remaps=2,new_checkpoint='test_checkpoint_0004.bin.zst'))),patch.object(runner,'vtk_members',return_value={('archive',str(i)):str(i) for i in range(62)}),patch.object(runner,'numerical',return_value=dict(passed=True)),patch.object(runner,'publication_evidence',return_value={}),patch.object(runner,'storage_evidence',return_value={}),patch.object(runner,'transform_evidence',return_value=dict(enabled=True,global_kernel_records_by_transaction=[[1,1,1]]*36)),patch.object(runner,'contract_evidence',return_value=dict(ownership_gate_passed=True)),patch.object(runner,'ordinary',side_effect=timing),patch.object(runner,'attribute',side_effect=attribute):
                if fail:
                    with self.assertRaises(RuntimeError):runner.run(kit,fixture,out)
                else:runner.run(kit,fixture,out)
            return json.loads((out/'results.json').read_text())
    def test_complete_49_run_path(self):
        r=self.exercise(.1);self.assertEqual(r['status'],'complete');self.assertEqual(len(r['schedule']),49);self.assertTrue(r['ownership_checks_passed'])
    def test_small_exchange_cannot_change_fixed_49_run_path(self):
        r=self.exercise(.001);self.assertEqual(r['status'],'complete');self.assertEqual(len(r['schedule']),49)
    def test_failed_solver_retains_active_schedule(self):
        r=self.exercise(.1,True);self.assertEqual(r['status'],'failed');self.assertEqual(len(r['schedule']),1);self.assertIn('solver/binary failure',r['error'])

if __name__=='__main__':unittest.main()
