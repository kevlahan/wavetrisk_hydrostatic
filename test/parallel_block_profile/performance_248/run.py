"""One-node bbserv campaign: build, verify, then compare separate full-run timings."""
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import tarfile
import hashlib

from build import run as build, save, sha, verify_seal

STEP_RESTRICTIONS = ('SLURM_CPUS_PER_TASK', 'SRUN_CPUS_PER_TASK', 'SLURM_EXACT',
                     'SLURM_EXCLUSIVE', 'SLURM_TRES_PER_TASK', 'SLURM_HINT',
                     'SLURM_THREADS_PER_CORE')
J5_HASHES = dict(simple='ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
                 seed='8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d',
                 grid='d4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d')


def verify_performance_kit(directory):
    for name, expected in json.loads((directory/'files-sha256.json').read_text()).items():
        if sha(directory/name) != expected:
            raise ValueError('Performance script differs: '+name)


def input_snapshot(kit, fixture, out, level, parameters, checkpoint_name):
    fixture = fixture.resolve(strict=True)
    p = parameters((fixture/'simple.in').read_text())
    expected = dict(test_case='climate', physics_type='Simple', max_level='6' if level=='J4' else '7',
                    zlevels='30', Nsoil='10', NCAR_topo='F', sso='F', CP_EVERY='1',
                    resume='2' if level=='J4' else '3')
    if any(p.get(key) != value for key, value in expected.items()):
        raise ValueError('Fixture parameters differ from the qualified climate restart')
    if (float(p['time_end']), float(p['dt_write'])) != ((0.33, 0.02) if level=='J4' else (0.315, 0.01)):
        raise ValueError('Fixture timing/restart window differs')
    seed = checkpoint_name((fixture/'simple.in').read_text())
    grid = 'HR95JT_WT_003' if level=='J4' else 'HR95JT_WT_004'
    source_grid = kit/'fixture-J4/grids'/grid
    if level=='J5' and (sha(fixture/'simple.in'), sha(fixture/seed), sha(source_grid)) != (
            J5_HASHES['simple'], J5_HASHES['seed'], J5_HASHES['grid']):
        raise ValueError('J5 fixture differs from the established Stage230 input')
    out.mkdir()
    for name in ('simple.in', seed):
        shutil.copy2(fixture/name, out/name)
    (out/'grids').symlink_to((kit/'fixture-J4/grids').resolve(), target_is_directory=True)
    return dict(parameters=p, seed=seed, grid=grid,
                simple_sha256=sha(out/'simple.in'), seed_sha256=sha(out/seed), grid_sha256=sha(source_grid))


def runtime_environment(kit):
    env = {k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
    for line in (kit/'runtime.env').read_text().splitlines():
        match = re.fullmatch(r'export ([A-Za-z_][A-Za-z_0-9]*)=([A-Za-z0-9_+-]+)', line)
        if not match:
            raise ValueError('Invalid sealed runtime selector: '+line)
        env[match[1]] = match[2]
    env.update(MV2_ENABLE_AFFINITY='0', OMPI_MCA_hwloc_base_binding_policy='none')
    return env


def vtk_members(root):
    archives = sorted(root.glob('*.vtk.tgz'))
    if len(archives) != 2:
        raise ValueError('Require both VTK output events')
    members = {}
    for archive in archives:
        with tarfile.open(archive) as stream:
            for entry in stream:
                if entry.isfile():
                    key = (archive.name, entry.name)
                    if key in members:
                        raise ValueError('Duplicate VTK member')
                    members[key] = hashlib.sha256(stream.extractfile(entry).read()).hexdigest()
    if len(members) != 62:
        raise ValueError('Expected exactly 62 VTK members')
    return members


def storage_evidence(log, kind, ranks, steps):
    absent = [int(x) for x in re.findall(r'Domain prognostic storage absent rank=(\d+)', log)]
    output = [int(x) for x in re.findall(r'Native output alias adapter rank=(\d+)', log)]
    scratch = [int(x) for x in re.findall(r'Legacy scalar staging and RK snapshots absent rank=(\d+)', log)]
    if kind=='candidate':
        for rank in range(ranks):
            if absent.count(rank) != steps+output.count(rank) or scratch.count(rank) != output.count(rank) or output.count(rank)<1:
                raise ValueError(f'Native storage absence incomplete on rank {rank}')
    elif absent or scratch:
        raise ValueError('Mirrored baseline unexpectedly dropped Domain storage')
    return dict(absence_checks=len(absent), bounded_output_checks=len(scratch))


def run(kit, out, level, comparison, fixture=None):
    kit, out = kit.resolve(strict=True), out.resolve()
    perf = Path(__file__).resolve().parent
    verify_performance_kit(perf)
    verify_seal(kit)
    if level not in ('J4','J5') or comparison not in ('mirrored','feature-off') or (level=='J5') != (fixture is not None):
        raise ValueError('Choose J4 or J5 and mirrored or feature-off; J5 requires its established fixture')
    fixture = (kit/'fixture-J4') if level=='J4' else fixture
    if out.exists() or out.is_relative_to(kit) or out.is_relative_to(perf):
        raise ValueError('Use a fresh results directory outside both kits')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:
        raise ValueError('Run inside one allocated Slurm node')
    if socket.gethostname().split('.')[0]!='bb13' or os.environ.get('SLURM_JOB_PARTITION')!='nk':
        raise ValueError('Require the 256-core bb13 node through partition nk')
    for name in ('srun','numactl','lscpu','zstd','mpif90','gfortran','nf-config'):
        if not shutil.which(name):
            raise ValueError('Missing cluster dependency: '+name)
    removed = {key:os.environ.pop(key) for key in STEP_RESTRICTIONS if key in os.environ}
    sys.path.insert(0, str(kit/'source/test/parallel_block_profile'))
    from analyze import analyze
    from capability_223 import coverage, numerical
    from placement_200 import topology, verify
    from profile_protocol import parameters, checkpoint_name
    from placement import plan
    from timing import ordinary, screen, confirm
    low, high = ((16,41) if level=='J4' else (16,167))
    timed_high = 41 if level=='J4' else 83
    ranks_needed = sorted({low,high,timed_high})
    out.mkdir()
    result = dict(status='running', level=level, comparison=comparison, host=socket.gethostname(), slurm_job_id=os.environ['SLURM_JOB_ID'], slurm_submit_dir=os.environ.get('SLURM_SUBMIT_DIR'),
                  source_identity=dict(frozen_candidate=json.loads((kit/'LOCAL_VALIDATION.json').read_text())['frozen_candidate'],
                                       local_validation_sha256=sha(kit/'LOCAL_VALIDATION.json'),
                                       correction_sha256=sha(kit/'CORRECTION.json'),
                                       visitor_correction_sha256=sha(kit/'VISITOR_CORRECTION.json'),
                                       pentagon_correction_sha256=sha(kit/'PENTAGON_CORRECTION.json')),
                  removed_step_restrictions=removed,
                  policy=dict(regression_guard=0.03, near_threshold=0.005, maximum_abba_batches=2,
                              timing='median ordinary steps per independent full run'), runs={}, gates={}, decisions={}, schedule=[])
    save(out/'results.json',result)
    try:
        identity = input_snapshot(kit, fixture, out/'input', level, parameters, checkpoint_name)
        result['input']=identity
        result['build']=build(kit,out,level,comparison)
        env=runtime_environment(kit)
        top=topology()
        if top['host']!=socket.gethostname():
            raise ValueError('Build and placement host differ')
        save(out/'topology.json', top)
        cases={n:plan(top,n) for n in ranks_needed}
        result['placements']={str(n):dict(socket_rank_counts=case['socket_rank_counts'],host=case['host']) for n,case in cases.items()}
        rank_launcher=kit/'source/test/parallel_block_profile/rank_launch_223.py'
        def launch(n, *args):
            return ['srun','--nodes=1',f'--ntasks={n}','--cpu-bind=none','--kill-on-bad-exit=1',
                    sys.executable,str(rank_launcher),str(out/f'placement-{n}.json'),*map(str,args)]
        for n,case in cases.items():
            save(out/f'placement-{n}.json',case)
            probe=subprocess.run(launch(n,'--probe'),env=env,capture_output=True,text=True)
            (out/f'placement-{n}-probe.stdout.log').write_text(probe.stdout)
            (out/f'placement-{n}-probe.stderr.log').write_text(probe.stderr)
            if probe.returncode:
                raise RuntimeError(f'Placement probe {n} exited {probe.returncode}')
            rows=[json.loads(line) for line in probe.stdout.splitlines() if line.startswith('{')]
            verify(rows,case)
            save(out/f'placement-{n}-verified.json',rows)
        save(out/'results.json',result)
        anchors={}
        def execute(n,kind,phase,suffix):
            label=f'{level}-n{n}-{phase}-{suffix}-{kind}'
            directory=out/label
            directory.mkdir()
            for name in ('simple.in',identity['seed']):
                shutil.copy2(out/'input'/name,directory/name)
            (directory/'grids').symlink_to((kit/'fixture-J4/grids').resolve(),target_is_directory=True)
            binary=Path(result['build']['builds'][kind]['binary'])
            if sha(binary)!=result['build']['builds'][kind]['sha256']:
                raise ValueError('Built binary changed')
            shutil.copy2(binary,directory/'climate')
            result['active']=label
            result['schedule'].append(label)
            save(out/'results.json',result)
            print('Starting '+label,flush=True)
            command=launch(n,directory,'./climate','simple.in')
            with (directory/'run.log').open('x') as log:
                code=subprocess.run(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
            row=dict(kind=kind,phase=phase,ranks=n,returncode=code,command=command,
                     binary_sha256=sha(directory/'climate'),log_sha256=sha(directory/'run.log'))
            result['runs'][label]=row
            save(out/'results.json',result)
            if code or row['binary_sha256']!=result['build']['builds'][kind]['sha256']:
                raise RuntimeError(f'{label} solver/binary failure {code}')
            rank_rows=[]
            for rank in range(n):
                record=json.loads((directory/f'rank-{rank}.json').read_text())
                if record['rank']!=rank or record['returncode'] or record['hostname']!=top['host'] or record['affinity']!=[cases[n]['mapping'][rank]['cpu']]:
                    raise ValueError(f'{label}: rank {rank} placement/exit mismatch')
                rank_rows.append(record)
            row['max_rank_peak_rss_kib']=max(x['max_rss_kib'] for x in rank_rows)
            row['rank_peak_rss_sum_kib']=sum(x['max_rss_kib'] for x in rank_rows)
            analysis,cov=coverage(directory,level,identity['parameters']['run_id'])
            if len(analysis['step_records'])!=12 or cov['remaps']!=2:
                raise ValueError(label+': incomplete adaptive/restart trajectory')
            row.update(steps=12,remaps=2,storage=storage_evidence((directory/'run.log').read_text(),kind,n,12))
            states=[r['state'] for r in analysis['step_records']]
            vtk=vtk_members(directory)
            if n not in anchors:
                if kind!='baseline' or phase!='validation':
                    raise ValueError('Baseline validation must anchor each rank count')
                anchors[n]=(directory,states,vtk)
            else:
                reference,reference_states,reference_vtk=anchors[n]
                if states!=reference_states or vtk!=reference_vtk:
                    raise ValueError(label+': state/VTK differs from same-rank baseline')
                gate=numerical(reference,directory,cov['new_checkpoint'],cov['domains'],False)
                result['gates'][label]=gate
                if not gate['passed']:
                    raise ValueError(label+': checkpoint numerical gate failed')
            if n==high and low in anchors and states!=anchors[low][1]:
                raise ValueError('Adaptive states differ across rank counts')
            precise=(directory/(identity['parameters']['run_id']+'_log')).read_text()
            row['ordinary']=ordinary((directory/'run.log').read_text(),analysis,precise)
            row.update(state_sha256=hashlib.sha256('\n'.join(states).encode()).hexdigest(),vtk_members=len(vtk),
                       checkpoint_sha256=sha(directory/cov['new_checkpoint']))
            save(directory/'analysis.json',dict(step_records=analysis['step_records'],coverage=cov,
                                                 ordinary=row['ordinary'],storage=row['storage']))
            save(out/'results.json',result)
            print(f"PASS {label}: ordinary median {row['ordinary']['median']:.9g}s",flush=True)
            return row['ordinary']['median']
        for n in (low,high):
            for kind in ('baseline','candidate'):
                execute(n,kind,'validation',kind)
        for n,name in ((low,'screen'),(timed_high,'transfer')):
            if n not in anchors:
                execute(n,'baseline','validation','baseline')
                execute(n,'candidate','validation','candidate')
            a=execute(n,'baseline',name,'A')
            b=execute(n,'candidate',name,'B')
            decision=screen(a,b)
            result['decisions'][name]=decision
            save(out/'results.json',result)
            if decision['decision']=='review_regression':
                result.update(status='complete',decision='review_regression')
                result.pop('active',None)
                save(out/'results.json',result)
                print('COMPLETE: review_regression; confirmation stopped.',flush=True)
                return
        batches=[]
        while True:
            batch=[]
            for i,kind in enumerate(('baseline','candidate','candidate','baseline'),1):
                batch.append((kind,execute(timed_high,kind,'confirm',f'{len(batches)+1}-{i}')))
            batches.append(batch)
            decision=confirm(batches)
            result['decisions']['confirmation']=decision
            save(out/'results.json',result)
            if decision['decision']!='repeat_abba':
                break
        result.update(status='complete',decision=decision['decision'])
        result.pop('active',None)
        save(out/'results.json',result)
        print('COMPLETE: '+decision['decision']+'; no automatic promotion.',flush=True)
    except BaseException as error:
        result.update(status='failed',error=repr(error))
        save(out/'results.json',result)
        raise


if __name__=='__main__':
    if len(sys.argv) not in (5,6):
        raise SystemExit('Usage: run.py CANDIDATE_KIT FRESH_RESULTS J4|J5 mirrored|feature-off [J5_FIXTURE]')
    run(Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],sys.argv[4],Path(sys.argv[5]) if len(sys.argv)==6 else None)
