"""Exercise submission, duplicate guards, collection and the full runner schedule.

These are orchestration tests. The separate real four-rank replay establishes
local numerical validity; these stubs do not establish cluster performance.
"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tarfile
import build as builder
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(1, str(HERE.parent))
import common
import control
import run as runner
from evidence import publication_evidence


class Workflow(unittest.TestCase):
    def test_four_builds_share_physics_and_check_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, out, sources = root/'kit', root/'results', root/'input-sources'
            kit.mkdir(); out.mkdir(); sources.mkdir()
            manifests = {}
            for kind in ('stage266', 'main', 'candidate'):
                tree = sources/kind
                (tree/'src/physics/simple_physics/phyparam').mkdir(parents=True)
                (tree/'Makefile').write_text('source-' + kind)
                manifests[kind] = {'Makefile':common.sha(tree/'Makefile')}
            with tarfile.open(kit/'sources.tar.gz', 'w:gz') as archive:
                for tree in sources.iterdir():
                    archive.add(tree,arcname=tree.name)
            common.save(kit/'source-files.json', manifests)
            common.save(kit/'IDENTITY.json', {'stage':273})
            common.save(kit/'files-sha256.json', {name:common.sha(kit/name) for name in
                        ('sources.tar.gz','source-files.json','IDENTITY.json')})
            calls = []
            def invoke(command, **kwargs):
                calls.append(command)
                if command[0] == 'make':
                    tree = Path(kwargs['cwd'])
                    binary_dir = next(x.split('=',1)[1] for x in command if x.startswith('BIN_DIR='))
                    (tree/binary_dir).mkdir()
                    (tree/binary_dir/'climate').write_text(binary_dir)
                    objects = tree/'src/physics/simple_physics/phyparam/obj'
                    objects.mkdir(exist_ok=True)
                    for i in range(23):
                        (objects/f'physics-{i}.o').write_text(str(i))
                return subprocess.CompletedProcess(command,0)
            with patch.object(builder.subprocess,'run',side_effect=invoke), \
                 patch.object(builder.subprocess,'check_output',return_value='stub compiler'):
                notes = builder.build(kit,out)
                self.assertEqual(len(notes['builds']),4)
            self.assertEqual(set(notes['builds']), {'stage266','main','candidate','checked'})
            self.assertEqual(len(notes['physics_objects']),23)
            for kind,row in notes['builds'].items():
                self.assertIn('-j1',row['command'])
                self.assertEqual(row['physics_objects'], notes['physics_objects'])
                self.assertEqual('DEBUG=check' in row['command'],kind=='checked')
            self.assertTrue(notes['builds']['checked']['binary'].startswith('sources/candidate/'))

    def test_publication_and_installation_evidence(self):
        header = 'Stage265 direct inverse SOL return = on\nStage265 tiled canonical publication = on\n'
        candidate = header + 'Stage266 provisional first installation = deferred; oracle = retained\n'
        self.assertEqual(publication_evidence(header, 'stage265t', False)['provisional_installations'], 2)
        self.assertEqual(publication_evidence(candidate, 'candidate', False)['production_provisional_installations'], 1)
        proof = ''.join(f'Stage265 exact inverse return: final={f} patches/values= 100 3000\n'
                        f'Stage266 inverse installations: final={f} count={1 if f=="T" else 2}\n'
                        for f in ['F', 'F', 'F', 'T'] * 12)
        self.assertEqual(publication_evidence(candidate + proof, 'candidate', True)['oracle_installations'], 84)
        for log, kind, oracle in ((candidate, 'candidate', True), (candidate + proof, 'candidate', False),
                                  (header, 'main', False), (candidate, 'stage265t', False),
                                  (candidate + proof.replace('3000', '0'), 'candidate', True),
                                  (candidate + proof.replace('count=2', 'count=1'), 'candidate', True),
                                  (header, 'candidate', False)):
            with self.assertRaises(ValueError):
                publication_evidence(log, kind, oracle)

    def test_seal_rejects_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'source').write_text('original')
            common.save(root / 'files-sha256.json', {'source': common.sha(root / 'source')})
            common.verify(root)
            (root / 'source').write_text('changed')
            with self.assertRaises(ValueError):
                common.verify(root)

    def test_launch_duplicate_failed_collection_and_purged_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            home = Path(tmp)
            kit = home / 'kit'
            kit.mkdir()
            calls = []
            def check_output(command, **kwargs):
                calls.append(command)
                return '12345\n'
            def invoke(command, **kwargs):
                calls.append(command)
                if command[0] == 'sacct':
                    return subprocess.CompletedProcess(command, 1, '', 'purged')
                return subprocess.CompletedProcess(command, 0, '', '')
            with patch.object(Path, 'home', return_value=home), \
                 patch.object(control, 'preflight'), \
                 patch.object(control.shutil, 'which', return_value='/stub/sbatch'), \
                 patch.object(control.subprocess, 'check_output', side_effect=check_output), \
                 patch.object(control.subprocess, 'run', side_effect=invoke):
                control.submit(kit)
                state = json.loads((home / 'rk273-latest.json').read_text())
                command = next(c for c in calls if c[0] == 'sbatch')
                self.assertEqual(command[-3:], [str(kit.resolve()), str(home / 'stage187-fixture-20260911-145307'), state['out']])
                with self.assertRaises(ValueError):
                    control.submit(kit)
                self.assertEqual(sum(c[0] == 'sbatch' for c in calls), 1)
                Path(state['log']).write_text('early build failure\n')
                control.status(state)
                result_dir=Path(state['out']);result_dir.mkdir()
                (result_dir/'cost-rank-000000.txt').write_text('partial observer failure')
                with self.assertRaises(ValueError):control.submit(kit,repeat=True)
                control.collect(state)
                archive = home / 'rk273-review-12345.tar.gz'
                self.assertGreater(control.verify_archive(archive), 0)
                with tarfile.open(archive) as stream:self.assertEqual(stream.extractfile('review/results/cost-rank-000000.txt').read(),b'partial observer failure')
                control.collect(state)
                self.assertEqual(sum(c[0] == 'scp' for c in calls), 2)

    def test_runner_reaches_all_twenty_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, fixture, out = root / 'kit', root / 'fixture', root / 'results'
            kit.mkdir()
            fixture.mkdir()
            (kit / 'files-sha256.json').write_text('{}')
            identity=dict(stage=273,revision='balanced-partition',main_commit=common.MAIN,stage266_commit=common.BASELINE,changed_source=['src/parallel_block_mpi.f90'],baseline_source_sha256={'src/parallel_block_mpi.f90':'old'},candidate_source_sha256={'src/parallel_block_mpi.f90':'new'})
            common.save(kit/'IDENTITY.json',identity)
            common.save(kit/'source-files.json',{'stage266':{'src/parallel_block_mpi.f90':'old'},'candidate':{'src/parallel_block_mpi.f90':'new'},'main':{}})
            (kit / 'grids').mkdir()
            for name in ('simple.in', 'test_checkpoint_0003.bin.zst'):
                (fixture / name).write_text('fixture')
            (kit / 'grids/HR95JT_WT_004').write_text('grid')
            top = dict(host='bb13', rows=[])
            for cpu in range(256):
                top['rows'].append(dict(cpu=cpu, core=cpu % 128, socket=cpu // 128,
                                       node=cpu // 128, llc=list(range(cpu // 8 * 8, cpu // 8 * 8 + 8))))
            cases = runner.plan(top)
            sources = {}
            def fake_build(kit_arg, out_arg):
                builds = {}
                for kind in ('stage266', 'main', 'candidate', 'checked'):
                    binary = out_arg / ('binary-' + kind)
                    binary.write_text(kind)
                    sources[kind] = common.sha(binary)
                    log=out_arg/('build-'+kind+'.log');log.write_text('serial build')
                    builds[kind] = dict(binary=binary.name, binary_sha256=common.sha(binary),physics_objects={str(i):str(i) for i in range(23)},log_sha256=common.sha(log),command=['make','-j1','DEBUG='+('check' if kind=='checked' else 'false')])
                notes=dict(host='bb13', builds=builds,status='passed',level='J5',physics_objects={str(i):str(i) for i in range(23)},source_manifest_sha256=common.sha(kit/'source-files.json'))
                common.save(out_arg/'builds.json',notes)
                return notes
            def invoke(command, **kw):
                self.assertEqual(command[0], 'srun')
                n = int(next(c.split('=')[1] for c in command if c.startswith('--ntasks=')))
                if '--probe' in command:
                    rows = [dict(rank=i, hostname='bb13', affinity=[cases[n]['mapping'][i]['cpu']]) for i in range(n)]
                    return subprocess.CompletedProcess(command, 0, '\n'.join(json.dumps(r) for r in rows), '')
                directory = Path(kw['cwd'])
                self.assertEqual(kw['env']['WAVETRISK_COST_AUDIT'],'0')
                self.assertEqual(kw['env']['OMP_NUM_THREADS'],'1')
                expected='1' if '-validation-' in directory.name and directory.name.endswith(('-candidate','-checked')) else '0'
                self.assertEqual(kw['env']['WAVETRISK_BLOCK_DETAILED_DIAGNOSTICS'],expected)
                for rank in range(n):
                    common.save(directory / f'rank-{rank}.json', dict(rank=rank, hostname='bb13', returncode=0,
                                affinity=[cases[n]['mapping'][rank]['cpu']]))
                log = f'number of processors = {n}\nRestarting from checkpoint 3\nRestarting from checkpoint 4\nRemapping vertical coordinates\nRemapping vertical coordinates\nStage203 inverse workspace = on\nStage203 scalar fields = on\n'
                log += f'Resident block execution: MPI tasks = {n}, geometric Domains = 160\n'
                for name in ('dynamics', 'adaptation', 'vertical-remap'):
                    log += f'Block {name} validation oracle enabled\n'
                kw['stdout'].write(log)
                kw['stdout'].flush()
                (directory / 'test_log').write_text('precise timing')
                return subprocess.CompletedProcess(command, 0)
            def fake_sha(path):
                path = Path(path)
                if path.name == 'test_checkpoint_0004.bin.zst':
                    return common.CHECKPOINTS['main' if path.parent.name.endswith('-main') else 'stage266']
                return common.sha(path)
            def timing(log, analysis, precise):
                return dict(median=.75, excluded=[dict(index=i) for i in (0, 7, 8)])
            env = dict(SLURM_JOB_ID='12345', SLURM_JOB_PARTITION='nk', SLURM_JOB_NUM_NODES='1')
            with patch.dict(os.environ, env), \
                 patch.object(runner, 'preflight', return_value=identity), \
                 patch.object(runner.socket, 'gethostname', return_value='bb13'), \
                 patch.object(runner.shutil, 'which', return_value='/stub'), \
                 patch.object(runner, 'build', side_effect=fake_build) as builds, \
                 patch.object(runner, 'topology', return_value=top), \
                 patch.object(runner.subprocess, 'run', side_effect=invoke), \
                 patch.object(runner, 'sha', side_effect=fake_sha), \
                 patch.object(runner, 'J5_HASHES', {'simple.in':common.sha(fixture/'simple.in'),
                     'test_checkpoint_0003.bin.zst':common.sha(fixture/'test_checkpoint_0003.bin.zst'),
                     'grids/HR95JT_WT_004':common.sha(kit/'grids/HR95JT_WT_004')}), \
                 patch.object(runner, 'coverage', return_value=(dict(step_records=[dict(state=str(i)) for i in range(12)]),
                      dict(remaps=2, new_checkpoint='test_checkpoint_0004.bin.zst'))), \
                 patch.object(runner, 'vtk_members', return_value={('archive', str(i)): str(i) for i in range(62)}), \
                 patch.object(runner, 'numerical', return_value=dict(passed=True,exact_required=True,report={})), \
                 patch.object(runner, 'storage_evidence', return_value={}), \
                 patch.object(runner, 'publication_evidence', return_value={}), \
                 patch.object(runner, 'catalogue_evidence', return_value={}), \
                 patch.object(runner, 'ordinary', side_effect=timing):
                runner.run(kit, fixture, out)
                initial = json.loads((out / 'results.json').read_text())
                self.assertEqual(len(initial['schedule']),20)
                self.assertFalse(initial['comparison']['automatic_promotion'])
                with self.assertRaises(ValueError):runner.run(kit, fixture, out, repeat=True)
                self.assertEqual(builds.call_count,1)
            result = json.loads((out / 'results.json').read_text())
            self.assertEqual(result['status'],'complete')
            self.assertEqual([(r['ranks'],r['kind']) for r in result['runs'].values() if r['phase']=='timing'],list(common.ORDER))
            self.assertEqual([(r['ranks'],r['kind'],r['oracle']) for r in result['runs'].values() if r['phase']=='validation'],list(common.VALIDATIONS))
            self.assertGreater(result['elapsed_seconds'],0)
            import review as reviewer
            state=dict(job=12345,kit=str(kit),out=str(out),log=str(root/'slurm.log'))
            with patch.object(Path,'home',return_value=root),patch.object(control,'queued',return_value=False),patch.object(control.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'','')):
                control.collect(state)
            archive=root/'rk273-review-12345.tar.gz'
            with patch.object(reviewer,'J5_HASHES',result['input_hashes']), \
                 patch.object(reviewer,'analyze',return_value=dict(completed=True,step_records=[dict(state=str(i)) for i in range(12)])), \
                 patch.object(reviewer,'accept',return_value=[]), \
                 patch.object(reviewer,'publication_evidence',return_value={}), \
                 patch.object(reviewer,'storage_evidence',return_value={}), \
                 patch.object(reviewer,'catalogue_evidence',return_value={}), \
                 patch.object(reviewer,'ordinary',side_effect=timing):
                checked=reviewer.review(archive)
                self.assertEqual(checked['status'],'verified')
                self.assertEqual(checked['comparison'],result['comparison'])


    def test_scaling_pairing_and_no_legacy_above_domains(self):
        rows=[dict(ranks=n,kind=k,ordinary=dict(median=(.70 if k=='candidate' else .75)*83/n)) for n,k in common.ORDER]
        result=common.summarize(rows)
        self.assertEqual(len(rows),14)
        self.assertIsNone(result['by_ranks']['256']['means_seconds']['main'])
        for n in common.TIMING_RANKS:
            record=result['by_ranks'][str(n)]
            self.assertEqual(record['assessment'],'both_pairs_faster')
            self.assertAlmostEqual(result['scaling']['candidate']['83'][str(n)]['relative_efficiency'],1.)
        rows[3]['ordinary']['median']=1.
        self.assertEqual(common.summarize(rows)['by_ranks']['83']['assessment'],'mixed_pair_signs')
        with self.assertRaises(ValueError):common.summarize(rows[:-1])

    def test_accounting_preserves_exact_runtime_and_unavailability(self):
        from accounting import parse_accounting,capture
        line='55467|COMPLETED|2026-10-07T17:40:43|2026-10-07T17:56:21|938|00:15:38\n'
        self.assertEqual(parse_accounting(line,'55467')['elapsed_seconds'],938)
        self.assertIsNone(parse_accounting(line,'other'))
        with self.assertRaises(ValueError):parse_accounting(line+line,'55467')
        with patch('accounting.subprocess.run',side_effect=FileNotFoundError('no sacct')):
            result=capture(55467)
            self.assertIsNone(result['job']);self.assertEqual(result['returncode'],-1)

    def test_source_scope_rejects_additional_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit=Path(tmp);name='src/parallel_block_mpi.f90'
            manifests={'stage266':{name:'old','Makefile':'same'},'candidate':{name:'new','Makefile':'same'},'main':{}}
            common.save(kit/'source-files.json',manifests)
            common.save(kit/'IDENTITY.json',dict(changed_source=[name],baseline_source_sha256={name:'old'},candidate_source_sha256={name:'new'}))
            self.assertEqual(common.source_scope(kit),[name])
            manifests['candidate']['Makefile']='different'
            common.save(kit/'source-files.json',manifests)
            with self.assertRaises(ValueError):common.source_scope(kit)

    def test_placement_covers_256_unique_physical_cores(self):
        from placement_273 import plan
        top=dict(host='bb13',rows=[dict(cpu=c,core=c%128,socket=c//128,node=c//128,llc=[c//8]) for c in range(256)])
        cases=plan(top)
        self.assertEqual(set(cases),{16,83,160,256})
        for n,case in cases.items():
            self.assertEqual(len({(r['socket'],r['core']) for r in case['mapping']}),n)
        self.assertEqual(cases[256]['socket_rank_counts'],[128,128])

    def test_catalogue_parser_and_timing_exclusion(self):
        from partition_evidence import catalogue_evidence,balanced
        text='Stage273 partition catalogue: 4 2\nStage273 partition block: 0 0 1 64 0 0\nStage273 partition block: 1 1 1 64 1 1\n'
        result=catalogue_evidence(text,4,True,domains=2)
        self.assertEqual(result['catalogues'][0]['candidate']['zero_block_ranks'],2)
        self.assertIsNone(result['catalogues'][0]['baseline_counterfactual'])
        for log,p,enabled in ((text,4,False),('',4,True),(text,3,True),
                              (text.replace('64 1 1','64 0 1'),4,True),
                              (text.replace('block: 1','block: 0'),4,True)):
            with self.assertRaises(ValueError):catalogue_evidence(log,p,enabled,domains=2)
        self.assertEqual(catalogue_evidence('',4,False)['catalogues'],[])
        owners=balanced([1]*640,256)
        counts=[owners.count(r) for r in range(256)]
        self.assertEqual(counts.count(2),128);self.assertEqual(counts.count(3),128)

if __name__=='__main__':unittest.main()
