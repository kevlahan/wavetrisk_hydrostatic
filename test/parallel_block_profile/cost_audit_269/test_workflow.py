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


from cost_analysis import NAMES, calibrate, compare_costs, read_costs


def fake_costs(value=.75):
    metrics=dict(self_wall=0.,self_cpu=0.,inclusive_wall=0.,inclusive_cpu=0.,calls=0)
    phases={name:dict(metrics) for name in NAMES}
    phases['step'].update(self_wall=value,inclusive_wall=value)
    return dict(phases=phases)

class Workflow(unittest.TestCase):
    def test_five_builds_share_physics_and_check_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, out, sources = root/'kit', root/'results', root/'input-sources'
            kit.mkdir(); out.mkdir(); sources.mkdir()
            manifests = {}
            for kind in ('stage266', 'main', 'observer', 'main_observer'):
                tree = sources/kind
                (tree/'src/physics/simple_physics/phyparam').mkdir(parents=True)
                (tree/'Makefile').write_text('source-' + kind)
                manifests[kind] = {'Makefile':common.sha(tree/'Makefile')}
            with tarfile.open(kit/'sources.tar.gz', 'w:gz') as archive:
                for tree in sources.iterdir():
                    archive.add(tree,arcname=tree.name)
            common.save(kit/'source-files.json', manifests)
            common.save(kit/'IDENTITY.json', {'stage':269})
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
                self.assertEqual(len(notes['builds']),5)
            self.assertEqual(set(notes['builds']), {'stage266','main','observer','main_observer','checked'})
            self.assertEqual(len(notes['physics_objects']),23)
            for kind,row in notes['builds'].items():
                self.assertIn('-j1',row['command'])
                self.assertEqual(row['physics_objects'], notes['physics_objects'])
                self.assertEqual('DEBUG=check' in row['command'],kind=='checked')
            self.assertTrue(notes['builds']['checked']['binary'].startswith('sources/observer/'))

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
                state = json.loads((home / 'rk269-latest.json').read_text())
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
                archive = home / 'rk269-review-12345.tar.gz'
                self.assertGreater(control.verify_archive(archive), 0)
                with tarfile.open(archive) as stream:self.assertEqual(stream.extractfile('review/results/cost-rank-000000.txt').read(),b'partial observer failure')
                control.collect(state)
                self.assertEqual(sum(c[0] == 'scp' for c in calls), 2)

    def test_runner_reaches_all_fifteen_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, fixture, out = root / 'kit', root / 'fixture', root / 'results'
            kit.mkdir()
            fixture.mkdir()
            (kit / 'files-sha256.json').write_text('{}')
            identity=dict(stage=269,revision='bounded-cost-audit',main_commit=common.MAIN,stage266_commit=common.BASELINE,changed_source=dict(observer=[],main_observer=[]))
            common.save(kit/'IDENTITY.json',identity)
            common.save(kit/'source-files.json',{kind:{} for kind in ('stage266','main','observer','main_observer')})
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
                for kind in ('stage266', 'main', 'observer', 'main_observer', 'checked'):
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
                observed=directory.name.endswith(('-observer','-main_observer','-checked'))
                self.assertEqual(kw['env']['WAVETRISK_COST_AUDIT'],'1' if observed else '0')
                self.assertEqual(kw['env']['OMP_NUM_THREADS'],'1')
                for rank in range(n):
                    common.save(directory / f'rank-{rank}.json', dict(rank=rank, hostname='bb13', returncode=0,
                                affinity=[cases[n]['mapping'][rank]['cpu']]))
                log = f'number of processors = {n}\nRestarting from checkpoint 3\nRestarting from checkpoint 4\nRemapping vertical coordinates\nRemapping vertical coordinates\nStage203 inverse workspace = on\nStage203 scalar fields = on\n'
                log += 'Resident block execution: MPI tasks = 167, geometric Domains = 160\n'
                for name in ('dynamics', 'adaptation', 'vertical-remap'):
                    log += f'Block {name} validation oracle enabled\n'
                kw['stdout'].write(log)
                kw['stdout'].flush()
                (directory / 'test_log').write_text('precise timing')
                return subprocess.CompletedProcess(command, 0)
            def fake_sha(path):
                path = Path(path)
                if path.name == 'test_checkpoint_0004.bin.zst':
                    return common.CHECKPOINTS['main' if path.parent.name.endswith(('-main','-main_observer')) else 'stage266']
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
                 patch.object(runner, 'ordinary', side_effect=timing), \
                 patch.object(runner, 'read_costs', return_value=fake_costs()):
                runner.run(kit, fixture, out)
                initial = json.loads((out / 'results.json').read_text())
                self.assertEqual(len(initial['schedule']),15)
                self.assertEqual(initial['comparison']['attribution']['status'],'calibrated_diagnostic_only')
                with self.assertRaises(ValueError):runner.run(kit, fixture, out, repeat=True)
                self.assertEqual(builds.call_count,1)
            result = json.loads((out / 'results.json').read_text())
            self.assertEqual(result['status'],'complete')
            self.assertEqual([r['kind'] for r in result['runs'].values() if r['phase']=='timing'],list(common.ORDER))
            self.assertEqual([(r['ranks'],r['kind'],r['oracle']) for r in result['runs'].values() if r['phase']=='validation'],list(common.VALIDATIONS))
            self.assertEqual(sum(r['observer'] for r in result['runs'].values()),9)
            import review as reviewer
            state=dict(job=12345,kit=str(kit),out=str(out),log=str(root/'slurm.log'))
            with patch.object(Path,'home',return_value=root),patch.object(control,'queued',return_value=False),patch.object(control.subprocess,'run',return_value=subprocess.CompletedProcess([],0)):
                control.collect(state)
            archive=root/'rk269-review-12345.tar.gz'
            with patch.object(reviewer,'J5_HASHES',result['input_hashes']), \
                 patch.object(reviewer,'analyze',return_value=dict(completed=True,step_records=[dict(state=str(i)) for i in range(12)])), \
                 patch.object(reviewer,'accept',return_value=[]), \
                 patch.object(reviewer,'publication_evidence',return_value={}), \
                 patch.object(reviewer,'storage_evidence',return_value={}), \
                 patch.object(reviewer,'ordinary',side_effect=timing), \
                 patch.object(reviewer,'read_costs',return_value=fake_costs()):
                checked=reviewer.review(archive)
                self.assertEqual(checked['status'],'verified')
                self.assertEqual(checked['comparison'],result['comparison'])


    def test_calibration_failure_retained_and_signed_offsets(self):
        rows=[]
        for kind in common.ORDER:
            family='main' if kind.startswith('main') else 'stage266'
            observed=kind.endswith('observer')
            costs=fake_costs(.54 if family=='main' else .75)
            if family=='stage266':
                costs['phases']['step']['self_wall']=.70
                costs['phases']['wavelet transaction']['self_wall']=.03
                costs['phases']['tendency context']['self_wall']=.02
            if family=='main':
                costs['phases']['step']['self_wall']=.50
                costs['phases']['tendency context']['self_wall']=.04
            rows.append(dict(kind=kind,family=family,observer=observed,costs=costs,
                             ordinary=dict(median=.54 if family=='main' else .75)))
        result=common.summarize(rows)
        self.assertAlmostEqual(result['candidate_minus_main_seconds'],.21)
        self.assertAlmostEqual(result['attribution']['additive_mean_rank_difference_seconds'],.21)
        self.assertEqual(result['attribution']['signed_self_wall_seconds']['tendency context'],-.02)
        self.assertIsNone(result['attribution']['recoverable_budget_seconds'])
        rows[1]['ordinary']['median']+=.012
        self.assertEqual(common.summarize(rows)['attribution']['status'],'observer_costs_unqualified')
        with self.assertRaises(ValueError):common.summarize(rows[:-1])

    def test_overlay_scope_rejects_unsealed_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit=Path(tmp);base={'src/main.f90':'a','Makefile':'m'}
            observer=dict(base,**{'src/main.f90':'new','src/cost_observer.f90':'cost'})
            manifests=dict(stage266=base,main=base,observer=observer,main_observer=observer)
            common.save(kit/'source-files.json',manifests)
            common.save(kit/'IDENTITY.json',dict(changed_source={kind:['src/cost_observer.f90','src/main.f90'] for kind in ('observer','main_observer')}))
            common.source_scope(kit)
            manifests['observer']=dict(observer,**{'src/other.f90':'unsealed'})
            common.save(kit/'source-files.json',manifests)
            with self.assertRaises(ValueError):common.source_scope(kit)

    def test_raw_records_conservation_and_one_rank_composition(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            for rank in range(2):
                lines=[f'STAGE269_COST_V1 {rank} 12 27']
                for step in range(1,13):
                    for phase in range(1,28):
                        calls={1:1,2:1,3:4,4:4}.get(phase,0)
                        total=1.+rank
                        self_time=total if phase==(3 if rank else 4) else 0.
                        inclusive=total if phase in (1,2,3,4) else 0.
                        lines.append(f'{step} {phase} {calls} {self_time} {self_time} {inclusive} {inclusive} 0 0 0')
                (directory/f'cost-rank-{rank:06d}.txt').write_text('\n'.join(lines)+'\n')
            result=read_costs(directory,2,False)
            self.assertAlmostEqual(result['mean_rank_step_seconds'],1.5)
            self.assertAlmostEqual(result['mean_longest_local_step_seconds'],2.)
            self.assertEqual(result['longest_local_step_rank_phases']['wavelet transaction'],0.)
            path=directory/'cost-rank-000001.txt'
            path.write_text(path.read_text().replace('1 3 4 2.0 2.0','1 3 4 1.9 1.9',1))
            with self.assertRaisesRegex(ValueError,'conservation'):read_costs(directory,2,False)

    def test_wrapper_balances_early_return_and_nested_routine(self):
        from instrument import wrap, inject
        source='subroutine outer\n integer :: i\n i=1\n if(i==1) return\n contains\n subroutine child\n return\n end subroutine\nend subroutine outer\n'
        result=wrap(source,'outer','INSTALL')
        self.assertEqual(result.count('call cost_enter'),1)
        self.assertEqual(result.count('call cost_leave'),2)
        self.assertIn('subroutine child\n return\n',result)
        with self.assertRaises(ValueError):inject(source,'absent','new')

if __name__=='__main__':unittest.main()
