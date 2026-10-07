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
            for kind in ('stage230', 'main', 'candidate'):
                tree = sources/kind
                (tree/'src/physics/simple_physics/phyparam').mkdir(parents=True)
                (tree/'Makefile').write_text('source-' + kind)
                manifests[kind] = {'Makefile':common.sha(tree/'Makefile')}
            with tarfile.open(kit/'sources.tar.gz', 'w:gz') as archive:
                for tree in sources.iterdir():
                    archive.add(tree,arcname=tree.name)
            common.save(kit/'source-files.json', manifests)
            common.save(kit/'IDENTITY.json', {'stage':265})
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
            self.assertEqual(set(notes['builds']), {'stage230','main','candidate','checked'})
            self.assertEqual(len(notes['physics_objects']),23)
            for kind,row in notes['builds'].items():
                self.assertIn('-j1',row['command'])
                self.assertEqual(row['physics_objects'], notes['physics_objects'])
                self.assertEqual('DEBUG=check' in row['command'],kind=='checked')
            self.assertTrue(notes['builds']['checked']['binary'].startswith('sources/candidate/'))

    def test_pairing(self):
        rows = [dict(kind=kind, ordinary=dict(median=.54 if kind == 'main' else .755))
                for kind in common.ORDER]
        result = common.summarize(rows)
        self.assertEqual(result['independent_runs_per_version'], 4)
        self.assertAlmostEqual(result['stage230_minus_main_seconds'], .215)
        self.assertAlmostEqual(result['candidate_minus_stage230_seconds'], 0.0)
        with self.assertRaises(ValueError):
            common.summarize(rows[:-1])

    def test_source_scope_and_order_balanced_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            kit = Path(tmp)
            baseline = {'src/parallel_block_inverse.f90':'a', 'src/parallel_block_mpi.f90':'b', 'src/other.f90':'c'}
            candidate = dict(baseline, **{'src/parallel_block_inverse.f90':'new-a', 'src/parallel_block_mpi.f90':'new-b'})
            identity = dict(changed_source=['src/parallel_block_inverse.f90','src/parallel_block_mpi.f90'],
                            candidate_source_sha256={name:candidate[name] for name in candidate if candidate[name] != baseline[name]})
            common.save(kit / 'IDENTITY.json', identity)
            manifest = dict(stage230=baseline, candidate=candidate)
            common.save(kit / 'source-files.json', manifest)
            self.assertEqual(len(common.source_scope(kit)), 2)
            candidate['src/other.f90'] = 'unexpected'
            common.save(kit / 'source-files.json', manifest)
            with self.assertRaises(ValueError):
                common.source_scope(kit)
        medians = [.54, .76, .74, .745, .765, .545] * 2
        rows = [dict(kind=k, ordinary=dict(median=v)) for k,v in zip(common.ORDER,medians)]
        result = common.summarize(rows)
        for gap in result['paired_candidate_gaps_seconds']:
            self.assertAlmostEqual(gap, -.020)

    def test_publication_evidence(self):
        header = 'Stage265 direct inverse SOL return = on\n'
        self.assertEqual(publication_evidence(header, 'candidate', False)['exact_transactions'], 0)
        proof = ''.join(f'Stage265 exact inverse return: final={f} patches/values= 100 3000\n'
                        for f in ['F', 'F', 'F', 'T'] * 12)
        self.assertEqual(publication_evidence(header + proof, 'candidate', True)['exact_transactions'], 48)
        for log, kind, oracle in ((header, 'candidate', True), (header + proof, 'candidate', False),
                                  (header, 'main', False), (header + proof.replace('3000', '0'), 'candidate', True)):
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
                state = json.loads((home / 'rk265-latest.json').read_text())
                command = next(c for c in calls if c[0] == 'sbatch')
                self.assertEqual(command[-3:], [str(kit.resolve()), str(home / 'stage187-fixture-20260911-145307'), state['out']])
                with self.assertRaises(ValueError):
                    control.submit(kit)
                self.assertEqual(sum(c[0] == 'sbatch' for c in calls), 1)
                Path(state['log']).write_text('early build failure\n')
                control.status(state)
                control.collect(state)
                archive = home / 'rk265-review-12345.tar.gz'
                self.assertGreater(control.verify_archive(archive), 0)
                control.collect(state)
                self.assertEqual(sum(c[0] == 'scp' for c in calls), 2)

    def test_runner_reaches_all_nineteen_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            kit, fixture, out = root / 'kit', root / 'fixture', root / 'results'
            kit.mkdir()
            fixture.mkdir()
            (kit / 'files-sha256.json').write_text('{}')
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
                for kind in ('stage230', 'main', 'candidate', 'checked'):
                    binary = out_arg / ('binary-' + kind)
                    binary.write_text(kind)
                    sources[kind] = common.sha(binary)
                    builds[kind] = dict(binary=binary.name, binary_sha256=common.sha(binary))
                return dict(host='bb13', builds=builds)
            def invoke(command, **kw):
                self.assertEqual(command[0], 'srun')
                n = int(next(c.split('=')[1] for c in command if c.startswith('--ntasks=')))
                if '--probe' in command:
                    rows = [dict(rank=i, hostname='bb13', affinity=[cases[n]['mapping'][i]['cpu']]) for i in range(n)]
                    return subprocess.CompletedProcess(command, 0, '\n'.join(json.dumps(r) for r in rows), '')
                directory = Path(kw['cwd'])
                for rank in range(n):
                    common.save(directory / f'rank-{rank}.json', dict(rank=rank, hostname='bb13', returncode=0,
                                affinity=[cases[n]['mapping'][rank]['cpu']]))
                log = 'Stage203 inverse workspace = on\nStage203 scalar fields = on\n'
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
                    return common.CHECKPOINTS['main' if path.parent.name.endswith('-main') else 'stage230']
                return common.sha(path)
            def timing(log, analysis, precise):
                return dict(median=.75, excluded=[dict(index=i) for i in (0, 7, 8)])
            env = dict(SLURM_JOB_ID='12345', SLURM_JOB_PARTITION='nk', SLURM_JOB_NUM_NODES='1')
            with patch.dict(os.environ, env), \
                 patch.object(runner, 'preflight', return_value={}), \
                 patch.object(runner.socket, 'gethostname', return_value='bb13'), \
                 patch.object(runner.shutil, 'which', return_value='/stub'), \
                 patch.object(runner, 'build', side_effect=fake_build), \
                 patch.object(runner, 'topology', return_value=top), \
                 patch.object(runner.subprocess, 'run', side_effect=invoke), \
                 patch.object(runner, 'runtime', return_value={}), \
                 patch.object(runner, 'sha', side_effect=fake_sha), \
                 patch.object(runner, 'J5_HASHES', {'simple.in':common.sha(fixture/'simple.in'),
                     'test_checkpoint_0003.bin.zst':common.sha(fixture/'test_checkpoint_0003.bin.zst'),
                     'grids/HR95JT_WT_004':common.sha(kit/'grids/HR95JT_WT_004')}), \
                 patch.object(runner, 'coverage', return_value=(dict(step_records=[dict(state=str(i)) for i in range(12)]),
                      dict(remaps=2, new_checkpoint='test_checkpoint_0004.bin.zst'))), \
                 patch.object(runner, 'vtk_members', return_value={('archive', str(i)): str(i) for i in range(62)}), \
                 patch.object(runner, 'numerical', return_value=dict(passed=True)), \
                 patch.object(runner, 'storage_evidence', return_value={}), \
                 patch.object(runner, 'publication_evidence', return_value={}), \
                 patch.object(runner, 'ordinary', side_effect=timing):
                runner.run(kit, fixture, out)
            result = json.loads((out / 'results.json').read_text())
            self.assertEqual(result['status'], 'complete')
            self.assertEqual(len(result['schedule']), 19)
            self.assertEqual([r['kind'] for r in result['runs'].values() if r['phase'] == 'timing'], list(common.ORDER))


if __name__ == '__main__':
    unittest.main()
