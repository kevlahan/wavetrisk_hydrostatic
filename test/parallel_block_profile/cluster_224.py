#!/usr/bin/env python3
"""Stage224: bounded 41/40 and 167/160 capability gates, no performance promotion."""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from capability_223 import digest, save, verify_kit, verify_sources, schedule, coverage, events, numerical
from build_capability_224 import KEYS
from profile_protocol import parameters, checkpoint_name, environment
from placement_223 import topology, plan, verify
from rank_metrics_190 import rank_resources
from run_local_220 import SELECTORS


# Explicit -c implies --exact and narrows the CPU set, defeating whole-node placement.
STEP_RESTRICTIONS = ('SLURM_CPUS_PER_TASK', 'SRUN_CPUS_PER_TASK', 'SLURM_EXACT',
                     'SLURM_EXCLUSIVE', 'SLURM_TRES_PER_TASK', 'SLURM_HINT',
                     'SLURM_THREADS_PER_CORE')


def clear_step_restrictions(env):
    return {key: env.pop(key) for key in STEP_RESTRICTIONS if key in env}


def probe_placement(kit, out, n, case):
    path = out / f'placement-{n}.json'
    save(path, case)
    command = launch(kit, path, n, '--probe')
    probe = subprocess.run(command, capture_output=True, text=True)
    (out / f'placement-{n}-probe.stdout.log').write_text(probe.stdout)
    (out / f'placement-{n}-probe.stderr.log').write_text(probe.stderr)
    save(out / f'placement-{n}-probe.json', dict(command=command, returncode=probe.returncode))
    if probe.returncode:
        print(probe.stderr, file=sys.stderr, flush=True)
        raise RuntimeError(f'Placement probe for {n} ranks exited {probe.returncode}; '
                           f'see {out}/placement-{n}-probe.stderr.log')
    rows = [json.loads(line) for line in probe.stdout.splitlines() if line.startswith('{')]
    verify(rows, case)
    save(out / f'placement-{n}-verified.json', rows)


def launch(kit, placement, ranks, *args):
    return ['srun','--nodes=1',f'--ntasks={ranks}','--cpu-bind=none',
            '--kill-on-bad-exit=1','--time=00:15:00',sys.executable,str(kit/'rank_launch_223.py'),
            str(placement),*map(str,args)]


def run(kit, tree, fixture, out):
    kit=kit.resolve(strict=True);tree=tree.resolve(strict=True);fixture=fixture.resolve(strict=True);out=out.resolve()
    verify_kit(kit)
    if out.exists() or any(out.is_relative_to(p) or p.is_relative_to(out) for p in (kit,tree,fixture)):
        raise ValueError('Fresh separate result directory required')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:
        raise ValueError('Use one exclusive nk node')
    for command in ('srun','numactl','lscpu','zstd'):
        if not shutil.which(command):raise ValueError('Missing '+command)
    removed_restrictions=clear_step_restrictions(os.environ)
    out.mkdir();m=json.loads((kit/'stage224-sources.json').read_text())
    notes=json.loads((tree/'capability-224-builds.json').read_text())
    result=dict(stage=224,status='running',runs={},gates={},schedule=[],build_notes=notes,
                removed_step_restrictions=removed_restrictions,
                slurm_environment={k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
                scope='Capability and correctness only. Lifecycle evidence is not performance instrumentation; no speedup claim or automatic promotion.')
    def persist():save(out/'results.json',result)
    persist()
    try:
        if notes['status']!='passed' or set(notes['builds'])!=set(KEYS) or notes['fp_contract']!='off' or notes['manifest_sha256']!=digest(kit/'stage224-sources.json'):
            raise ValueError('Build identity incomplete/different')
        for kind in ('reference','candidate'):verify_sources(tree/kind,m[kind+'_sources'])
        top=topology();cases=plan(top);save(out/'topology.json',top)
        if notes['build_host']['node']!=top['host']:raise ValueError('Build/test native-compilation host differs')
        result['cases']=cases;persist()
        for n,case in cases.items():
            probe_placement(kit,out,n,case)
        inputs=out/'inputs';inputs.mkdir();result['inputs']={}
        for level,source in (('J4',kit/'fixture-J4'),('J5',fixture)):
            dst=inputs/level;dst.mkdir();text=(source/'simple.in').read_text();p=parameters(text)
            needed={'max_level':'6' if level=='J4' else '7','zlevels':'30','Nsoil':'10',
                    'physics_type':'Simple','NCAR_topo':'F','sso':'F','CP_EVERY':'1','resume':'2' if level=='J4' else '3'}
            if any(p.get(k)!=v for k,v in needed.items()):raise ValueError('Unexpected '+level+' fixture parameters')
            expected=(.33,.02) if level=='J4' else (.315,.01)
            if (float(p['time_end']),float(p['dt_write']))!=expected:raise ValueError('Changed fixture trajectory')
            seed=checkpoint_name(text)
            for name in ('simple.in',seed):shutil.copy2(source/name,dst/name)
            (dst/'grids').mkdir();grid='HR95JT_WT_003' if level=='J4' else 'HR95JT_WT_004'
            if level=='J5' and digest(source/'grids'/grid)!=digest(kit/'grids'/grid):
                raise ValueError('Cluster J5 grid differs from frozen grid')
            shutil.copy2(kit/'grids'/grid,dst/'grids'/grid)
            result['inputs'][level]=dict(parameters=p,seed=seed,input_sha256=digest(dst/'simple.in'),
                                        seed_sha256=digest(dst/seed),grid=grid,grid_sha256=digest(dst/'grids'/grid))
        (inputs/'binaries').mkdir()
        for key,record in notes['builds'].items():
            binary=tree/record['binary']
            expected_sources=m['reference_sources' if key.endswith('-baseline') else 'candidate_sources']
            if record['sources']!=expected_sources or digest(binary)!=record['sha256']:raise ValueError('Build binary/source mismatch')
            shutil.copy2(binary,inputs/'binaries'/key)
        persist();anchors={};same_count={};states={}
        for level,n,kind,oracle,layout in schedule():
            label=f'{level}-n{n}-{kind}-'+('oracles' if oracle else 'production')+'-'+layout
            directory=out/label;directory.mkdir();identity=result['inputs'][level];source=inputs/level
            result['active']=label;result['schedule'].append(label);persist();print('Starting '+label,flush=True)
            for name in ('simple.in',identity['seed']):shutil.copy2(source/name,directory/name)
            shutil.copytree(source/'grids',directory/'grids')
            binary=inputs/'binaries'/f'{level}-{kind}';shutil.copy2(binary,directory/'climate')
            env=environment(oracles=oracle,detail=False);env.update(SELECTORS)
            env.update(MV2_ENABLE_AFFINITY='0',OMPI_MCA_hwloc_base_binding_policy='none',WAVETRISK_SAMPLE_RANKS201='[]')
            if kind!='baseline':env.update(WAVETRISK_STAGE223_ALLOW='1',WAVETRISK_STAGE223_LAYOUT=layout)
            command=launch(kit,out/f'placement-{n}.json',n,directory,'./climate','simple.in')
            start=time.monotonic()
            with (directory/'run.log').open('x') as f:code=subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
            row=dict(ranks=n,level=level,kind=kind,oracle=oracle,layout=layout,returncode=code,
                     elapsed_seconds=time.monotonic()-start,command=command,binary_sha256=digest(directory/'climate'),
                     environment={k:v for k,v in env.items() if k.startswith(('WAVETRISK_','MV2_','OMPI_')) or k.endswith('NUM_THREADS')})
            result['runs'][label]=row;persist()
            if code:raise RuntimeError(label+f' exited {code}; later gates not run')
            if row['binary_sha256']!=notes['builds'][f'{level}-{kind}']['sha256']:raise ValueError('Executable snapshot changed')
            resources=rank_resources(directory);verify(resources,cases[n]);row['rank_metrics']=resources
            a,cov=coverage(directory,level,identity['parameters']['run_id']);save(directory/'analysis.json',a);row.update(cov)
            row['steps']=len(a['step_records']);row['log_sha256']=digest(directory/'run.log')
            log=(directory/'run.log').read_text()
            for pattern,expected in ((r'Stage 191 requested install/mass-overlap =\s*([01])\s+([01])',[('1','0')]),
                                     (r'Stage 197 requested resident wavelets = ([01])',['1'])):
                if re.findall(pattern,log)!=expected:raise ValueError('Runtime selector evidence differs')
            if kind!='baseline':row['capability']=events(directory,n,cov['domains'],row['steps'],layout)
            sequence=[r['state'] for r in a['step_records']]
            if level not in anchors:
                if kind!='baseline':raise ValueError('Baseline must anchor each fixture')
                anchors[level]=directory;states[level]=sequence
            if sequence!=states[level]:raise ValueError('Adaptive trajectory differs')
            exact=n==cov['domains']
            gate=numerical(anchors[level],directory,cov['new_checkpoint'],cov['domains'],exact)
            result['gates'][label]=gate
            if not gate['passed']:raise ValueError('Cross-count/normal-layout numerical gate failed: '+label)
            key=(level,n,layout)
            if key in same_count:
                gate=numerical(same_count[key],directory,cov['new_checkpoint'],cov['domains'],True)
                result['gates']['same-count-'+label]=gate
                if not gate['passed']:raise ValueError('Same-count exact gate failed')
            else:same_count[key]=directory
            row['checkpoint_sha256']=digest(directory/cov['new_checkpoint'])
            persist();print('PASS '+label,flush=True)
        for level,identity in result['inputs'].items():
            for name,expected in [('simple.in',identity['input_sha256']),(identity['seed'],identity['seed_sha256']),
                                  ('grids/'+identity['grid'],identity['grid_sha256'])]:
                if digest(inputs/level/name)!=expected:raise ValueError('Input snapshot changed')
        for key,record in notes['builds'].items():
            if digest(inputs/'binaries'/key)!=record['sha256']:raise ValueError('Binary snapshot changed')
        result.update(status='passed',decision='capability_gates_passed_pending_review');result.pop('active',None);persist()
        print('PASS: 41/40 and 167/160 capability gates; review before lifting the production guard.',flush=True)
    except BaseException as e:
        result.update(status='failed',error=repr(e));persist();raise


if __name__=='__main__':run(*(Path(p) for p in sys.argv[1:]))
