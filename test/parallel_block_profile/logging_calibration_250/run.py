"""Matched bb13 calibration of candidate logging and the residual main gap."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import socket
import statistics
import subprocess
import sys
import tarfile

from build import run as build, save, sha, verify_seal

J5_HASHES = dict(simple='ed173774b8400630657243fa67b4065faff1faf409ba93c1712bb80174d6ae8d',
                 seed='8c7b651c2d6bc6f481838163591565771e1f9a95fd54ec9adf0ce9b74a117e2d',
                 grid='d4d8393bf53d6275432f6f954ef3178a3ae2e3447a17160f5a3aa029ac39098d')
PRIOR_CANDIDATE_STATE = '20dc4b302215e0fbac234f29bc0e6e4a18dd108cf3665d500d6249b6398f6b88'
PRIOR_CANDIDATE_CHECKPOINT = '1c944c62270f6ea5e97eca5af0666b0e6cd1d612884d73192143d07dc9f23b53'
STEP_RESTRICTIONS = ('SLURM_CPUS_PER_TASK','SRUN_CPUS_PER_TASK','SLURM_EXACT',
                     'SLURM_EXCLUSIVE','SLURM_TRES_PER_TASK','SLURM_HINT','SLURM_THREADS_PER_CORE')


def verify_scripts(perf):
    files = json.loads((perf/'files-sha256.json').read_text())
    for name, expected in files.items():
        if sha(perf/name) != expected:
            raise ValueError('Performance kit file differs: '+name)


def snapshot_input(candidate_kit, fixture, out, parameters, checkpoint_name):
    fixture = fixture.resolve(strict=True)
    text = (fixture/'simple.in').read_text()
    p = parameters(text)
    needed = dict(test_case='climate',physics_type='Simple',run_id='test',max_level='7',
                  zlevels='30',Nsoil='10',NCAR_topo='F',sso='F',CP_EVERY='1',resume='3')
    if any(p.get(key) != value for key,value in needed.items()):
        raise ValueError('J5 fixture parameters differ from established campaign')
    if (float(p['time_end']),float(p['dt_write'])) != (.315,.01):
        raise ValueError('J5 timing or restart window differs')
    seed = checkpoint_name(text)
    grid = candidate_kit/'fixture-J4/grids/HR95JT_WT_004'
    if (sha(fixture/'simple.in'),sha(fixture/seed),sha(grid)) != tuple(J5_HASHES.values()):
        raise ValueError('J5 input, seed or grid differs from established campaign')
    out.mkdir()
    for name in ('simple.in',seed):
        shutil.copy2(fixture/name,out/name)
    (out/'grids').symlink_to((candidate_kit/'fixture-J4/grids').resolve(),target_is_directory=True)
    return dict(parameters=p,seed=seed,grid='HR95JT_WT_004',
                simple_sha256=sha(out/'simple.in'),seed_sha256=sha(out/seed),grid_sha256=sha(grid))


def runtime_environment(candidate_kit):
    env = {key:value for key,value in os.environ.items() if not key.startswith('WAVETRISK_')}
    for line in (candidate_kit/'runtime.env').read_text().splitlines():
        match = re.fullmatch(r'export ([A-Za-z_][A-Za-z_0-9]*)=([A-Za-z0-9_+-]+)',line)
        if not match:
            raise ValueError('Invalid sealed runtime selector: '+line)
        env[match[1]]=match[2]
    env.update(MV2_ENABLE_AFFINITY='0',OMPI_MCA_hwloc_base_binding_policy='none')
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
                    key = (archive.name,entry.name)
                    if key in members:
                        raise ValueError('Duplicate VTK member')
                    members[key]=hashlib.sha256(stream.extractfile(entry).read()).hexdigest()
    if len(members) != 62:
        raise ValueError('Expected exactly 62 VTK members')
    return members


def candidate_storage(log,ranks,quiet):
    absent = [int(x) for x in re.findall(r'Domain prognostic storage absent rank=(\d+)',log)]
    output = [int(x) for x in re.findall(r'Native output alias adapter rank=(\d+)',log)]
    scratch = [int(x) for x in re.findall(r'Legacy scalar staging and RK snapshots absent rank=(\d+)',log)]
    for rank in range(ranks):
        if (absent.count(rank) != (0 if quiet else 12+output.count(rank)) or
                scratch.count(rank) != output.count(rank) or output.count(rank)!=2):
            raise ValueError(f'Candidate storage absence incomplete on rank {rank}')
    return dict(absence_checks=len(absent),bounded_output_checks=len(scratch),quiet=quiet)


def run(perf, candidate_kit, main_tree, out, fixture):
    perf, candidate_kit, main_tree = (Path(p).resolve(strict=True) for p in (perf,candidate_kit,main_tree))
    out, fixture = Path(out).resolve(), Path(fixture).resolve(strict=True)
    verify_scripts(perf)
    verify_seal(candidate_kit)
    verify_seal(main_tree)
    if out.exists() or any(out.is_relative_to(p) for p in (perf,candidate_kit,main_tree,fixture)):
        raise ValueError('Use a fresh external result directory')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:
        raise ValueError('Run inside one allocated Slurm node')
    if socket.gethostname().split('.')[0]!='bb13' or os.environ.get('SLURM_JOB_PARTITION')!='nk':
        raise ValueError('Require bb13 through partition nk')
    for name in ('srun','numactl','lscpu','zstd','mpif90','gfortran','nf-config','git'):
        if not shutil.which(name):
            raise ValueError('Missing cluster dependency: '+name)
    removed = {key:os.environ.pop(key) for key in STEP_RESTRICTIONS if key in os.environ}
    sys.path.insert(0,str(candidate_kit/'source/test/parallel_block_profile'))
    from capability_223 import coverage,numerical
    from placement_200 import topology,verify
    from profile_protocol import parameters,checkpoint_name
    from placement import plan
    from timing import ordinary
    out.mkdir()
    result = dict(status='running',comparison='logging-on versus logging-off candidate; quiet candidate versus main',
                  level='J5',slurm_job_id=os.environ['SLURM_JOB_ID'],
                  slurm_submit_dir=os.environ.get('SLURM_SUBMIT_DIR'),host=socket.gethostname(),
                  main_identity=json.loads((main_tree/'MAIN_IDENTITY.json').read_text()),
                  candidate_identity=dict(local_validation_sha256=sha(candidate_kit/'LOCAL_VALIDATION.json'),
                                          correction_sha256=sha(candidate_kit/'CORRECTION.json'),
                                          visitor_correction_sha256=sha(candidate_kit/'VISITOR_CORRECTION.json'),
                                          pentagon_correction_sha256=sha(candidate_kit/'PENTAGON_CORRECTION.json')),
                  removed_step_restrictions=removed,runs={},gates={},schedule=[],decisions={},
                  policy=dict(timing='median of nine ordinary steps per independent full run',
                              comparison='J5 83-rank two A-B-B-A logging batches and one A-B-B-A quiet-main batch; same patched candidate binary in both logging modes'))
    save(out/'results.json',result)
    try:
        identity=snapshot_input(candidate_kit,fixture,out/'input',parameters,checkpoint_name)
        result['input']=identity
        result['build']=build(candidate_kit,main_tree,out,perf)
        env=runtime_environment(candidate_kit)
        top=topology()
        if top['host']!=socket.gethostname():
            raise ValueError('Build and placement host differ')
        save(out/'topology.json',top)
        cases={83:plan(top,83)}
        result['placements']={str(n):dict(socket_rank_counts=case['socket_rank_counts'],host=case['host']) for n,case in cases.items()}
        rank_launcher=candidate_kit/'source/test/parallel_block_profile/rank_launch_223.py'
        def launch(n,*args):
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
        anchors={}
        def execute(mode,phase,suffix):
            n=83
            kind='main' if mode=='main' else 'candidate'
            label=f'J5-n{n}-{phase}-{suffix}-{mode}'
            directory=out/label
            directory.mkdir()
            for name in ('simple.in',identity['seed']):
                shutil.copy2(out/'input'/name,directory/name)
            (directory/'grids').symlink_to((candidate_kit/'fixture-J4/grids').resolve(),target_is_directory=True)
            binary=Path(result['build']['builds'][kind]['binary'])
            if sha(binary)!=result['build']['builds'][kind]['sha256']:
                raise ValueError('Built binary changed')
            shutil.copy2(binary,directory/'climate')
            result['active']=label
            result['schedule'].append(label)
            save(out/'results.json',result)
            print('Starting '+label,flush=True)
            command=launch(n,directory,'./climate','simple.in')
            run_env=env.copy()
            run_env['WAVETRISK_QUIET_STORAGE_ABSENCE']='1' if mode=='quiet' else '0'
            with (directory/'run.log').open('x') as log:
                code=subprocess.run(command,cwd=directory,env=run_env,stdout=log,stderr=subprocess.STDOUT).returncode
            row=dict(kind=kind,mode=mode,phase=phase,ranks=n,returncode=code,command=command,
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
            analysis,cov=coverage(directory,'J5',identity['parameters']['run_id'])
            if len(analysis['step_records'])!=12 or cov['remaps']!=2:
                raise ValueError(label+': incomplete adaptive/restart trajectory')
            states=[r['state'] for r in analysis['step_records']]
            vtk=vtk_members(directory)
            row.update(steps=12,remaps=2,vtk_members=len(vtk),
                       state_sha256=hashlib.sha256('\n'.join(states).encode()).hexdigest(),
                       checkpoint_sha256=sha(directory/cov['new_checkpoint']))
            if kind=='candidate':
                row['storage']=candidate_storage((directory/'run.log').read_text(),n,mode=='quiet')
                if (row['state_sha256'],row['checkpoint_sha256']) != (PRIOR_CANDIDATE_STATE,PRIOR_CANDIDATE_CHECKPOINT):
                    raise ValueError(label+': candidate differs from qualified v7 J5 result')
                if states != anchors['main'][1]:
                    raise ValueError(label+': main/candidate adaptive trajectories differ')
                gate=numerical(anchors['main'][0],directory,cov['new_checkpoint'],160,False)
                result['gates'][label]=gate
                if not gate['passed']:
                    raise ValueError(label+': main/candidate numerical gate failed')
            if mode not in anchors:
                anchors[mode]=(directory,states,vtk,row['checkpoint_sha256'])
            else:
                reference,reference_states,reference_vtk,reference_checkpoint=anchors[mode]
                if states!=reference_states or vtk!=reference_vtk or row['checkpoint_sha256']!=reference_checkpoint:
                    raise ValueError(label+': repeated run differs from same-version reference')
            if mode=='quiet' and 'loud' in anchors:
                _,loud_states,loud_vtk,loud_checkpoint=anchors['loud']
                if states!=loud_states or vtk!=loud_vtk or row['checkpoint_sha256']!=loud_checkpoint:
                    raise ValueError(label+': disabling logging changed candidate output')
            row['ordinary']=ordinary((directory/'run.log').read_text(),analysis,
                                     (directory/(identity['parameters']['run_id']+'_log')).read_text())
            save(directory/'analysis.json',dict(step_records=analysis['step_records'],coverage=cov,
                                                 ordinary=row['ordinary'],storage=row.get('storage')))
            save(out/'results.json',result)
            print(f"PASS {label}: ordinary median {row['ordinary']['median']:.9g}s",flush=True)
            return row['ordinary']['median']
        for mode in ('main','loud','quiet'):
            execute(mode,'validation',mode)
        logging_batches=[]
        for batch_id in (1,2):
            batch=[]
            for i,mode in enumerate(('loud','quiet','quiet','loud'),1):
                batch.append((mode,execute(mode,'logging',f'{batch_id}-{i}')))
            logging_batches.append(batch)
            save(out/'results.json',result)
        loud_values=[value for batch in logging_batches for mode,value in batch if mode=='loud']
        quiet_values=[value for batch in logging_batches for mode,value in batch if mode=='quiet']
        logging_pairs=[(batch[0][1]-batch[1][1],batch[3][1]-batch[2][1]) for batch in logging_batches]
        loud_mean,quiet_mean=statistics.mean(loud_values),statistics.mean(quiet_values)
        result['decisions']['logging']=dict(loud_mean_seconds=loud_mean,quiet_mean_seconds=quiet_mean,
                                           saving_seconds=loud_mean-quiet_mean,
                                           saving_fraction=(loud_mean-quiet_mean)/loud_mean,
                                           paired_savings_seconds=[value for pair in logging_pairs for value in pair],
                                           independent_runs_per_mode=4,
                                           reference_stage249_loud_seconds=1.6776144205)
        save(out/'results.json',result)
        residual=[]
        for i,mode in enumerate(('main','quiet','quiet','main'),1):
            residual.append((mode,execute(mode,'residual',str(i))))
        main_mean=statistics.mean(value for mode,value in residual if mode=='main')
        quiet_residual_mean=statistics.mean(value for mode,value in residual if mode=='quiet')
        result['decisions']['residual']=dict(main_mean_seconds=main_mean,
                                            quiet_candidate_mean_seconds=quiet_residual_mean,
                                            candidate_minus_main_seconds=quiet_residual_mean-main_mean,
                                            candidate_to_main_ratio=quiet_residual_mean/main_mean,
                                            paired_regression_fractions=[residual[1][1]/residual[0][1]-1,
                                                                         residual[2][1]/residual[3][1]-1],
                                            independent_runs_per_version=2)
        result.update(status='complete',decision='review_logging_and_residual')
        result.pop('active',None)
        save(out/'results.json',result)
        print('COMPLETE: logging and residual comparisons; no automatic promotion.',flush=True)
    except BaseException as error:
        result.update(status='failed',error=repr(error))
        save(out/'results.json',result)
        raise


if __name__=='__main__':
    if len(sys.argv)!=6:
        raise SystemExit('Usage: run.py PERF_KIT CANDIDATE_KIT MAIN_TREE FRESH_RESULTS J5_FIXTURE')
    run(*(Path(p) for p in sys.argv[1:]))
