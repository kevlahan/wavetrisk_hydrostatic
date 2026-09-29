#!/usr/bin/env python3
"""Qualify removal of Domain tendency storage; bounded performance regression check."""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from capability_223 import digest, save, verify_kit, verify_sources, coverage, numerical
from protocol_230 import POLICY, ordinary, screen, confirm
from run_local_230 import vtk_members
from validation_230 import storage_evidence
from build_capability_230 import KEYS, source_key, source_dir
from profile_protocol import parameters, checkpoint_name, environment
from placement_230 import topology, plan, verify
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


def campaign(execute,checkpoint):
    for n in (16,167):
        versions=[('baseline',False),('checked',False),('candidate',True)]
        if n==167:versions.append(('candidate',False))
        for kind,oracle in versions:
            execute('J5',n,kind,oracle,'validation',kind+('-oracle' if oracle else ''))
    decisions={}
    for n,phase in ((16,'screen'),(83,'transfer')):
        a=execute('J5',n,'baseline',False,phase,'A')
        b=execute('J5',n,'candidate',False,phase,'B')
        decisions[phase]=screen(a,b);checkpoint(decisions)
        if decisions[phase]['decision']=='review_regression':return decisions
    batches=[]
    while True:
        batch=[]
        for i,kind in enumerate(('baseline','candidate','candidate','baseline'),1):
            batch.append((kind,execute('J5',83,kind,False,'confirm',f'{len(batches)+1}-{i}')))
        batches.append(batch);decisions['confirmation']=confirm(batches);checkpoint(decisions)
        if decisions['confirmation']['decision']!='repeat_abba':return decisions


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
    out.mkdir();m=json.loads((kit/'stage230-sources.json').read_text())
    notes=json.loads((tree/'capability-230-builds.json').read_text())
    result=dict(stage=230,status='running',runs={},gates={},schedule=[],build_notes=notes,
                removed_step_restrictions=removed_restrictions,
                slurm_environment={k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
                policy=POLICY, scope='Stage229 versus deletion of Domain tendency storage. Correctness at16/167; bare16/83 screens and83 ABBA. Regression guard3%; no minimum gain and no automatic promotion.')
    def persist():save(out/'results.json',result)
    persist()
    try:
        if notes['status']!='passed' or set(notes['builds'])!=set(KEYS) or notes['fp_contract']!='off' or notes['manifest_sha256']!=digest(kit/'stage230-sources.json'):
            raise ValueError('Build identity incomplete/different')
        for kind in ('baseline','candidate','checked'):
            verify_sources(tree/source_dir(kind),m[source_key(kind)])
        top=topology();cases=plan(top);save(out/'topology.json',top)
        if notes['build_host']['node']!=top['host']:raise ValueError('Build/test native-compilation host differs')
        result['cases']=cases;persist()
        for n,case in cases.items():
            probe_placement(kit,out,n,case)
        inputs=out/'inputs';inputs.mkdir();result['inputs']={}
        for level,source in (('J5',fixture),):
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
            expected_sources=m[source_key(key.split('-',1)[1])]
            if record['sources']!=expected_sources or digest(binary)!=record['sha256']:raise ValueError('Build binary/source mismatch')
            shutil.copy2(binary,inputs/'binaries'/key)
        persist();anchors={};same_count={};states={}
        def execute(level,n,kind,oracle,phase,suffix):
            label=f'{level}-n{n}-{phase}-{suffix}-{kind}'
            directory=out/label;directory.mkdir();identity=result['inputs'][level];source=inputs/level
            result['active']=label;result['schedule'].append(label);persist();print('Starting '+label,flush=True)
            for name in ('simple.in',identity['seed']):shutil.copy2(source/name,directory/name)
            shutil.copytree(source/'grids',directory/'grids')
            binary=inputs/'binaries'/f'{level}-{kind}';shutil.copy2(binary,directory/'climate')
            env=environment(oracles=oracle,detail=False);env.update(SELECTORS)
            env.update(MV2_ENABLE_AFFINITY='0',OMPI_MCA_hwloc_base_binding_policy='none',WAVETRISK_SAMPLE_RANKS201='[]')
            command=launch(kit,out/f'placement-{n}.json',n,directory,'./climate','simple.in')
            start=time.monotonic()
            with (directory/'run.log').open('x') as f:code=subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
            row=dict(ranks=n,level=level,kind=kind,oracle=oracle,phase=phase,returncode=code,
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
                if not kind.startswith('legacy') and re.findall(pattern,log)!=expected:raise ValueError('Runtime selector evidence differs')
            if row['steps']!=12 or cov['remaps']!=2:raise ValueError('Incomplete adaptive/restart trajectory')
            expected203=True
            row['storage']=storage_evidence(log,n,kind,oracle)
            for marker in ('Stage203 inverse workspace = on','Stage203 scalar fields = on'):
                if (marker in log)!=expected203:raise ValueError('Stage203 build/runtime selector differs')
            if kind=='candidate':
                if list(directory.glob('capability223-rank-*.txt')) or 'Stage223 capability prototype' in log:
                    raise ValueError('Production candidate contains prototype hooks')
                if any(k.startswith('WAVETRISK_STAGE223') for k in env):raise ValueError('Prototype environment leaked')
                if n>cov['domains'] and f'Resident block execution: MPI tasks = {n}, geometric Domains = {cov["domains"]}' not in log:
                    raise ValueError('Missing production rank-limit qualification')
            if oracle:
                for marker in ('Block dynamics validation oracle enabled','Block adaptation validation oracle enabled','Block vertical-remap validation oracle enabled'):
                    if marker not in log:raise ValueError('Missing oracle evidence: '+marker)
                if log.count('remap validation oracle max time')!=cov['remaps']:raise ValueError('Missing remap oracle completion')

            sequence=[r['state'] for r in a['step_records']]
            if level not in anchors:
                if kind!='baseline':raise ValueError('Baseline must anchor each fixture')
                anchors[level]=directory;states[level]=sequence
            if sequence!=states[level]:raise ValueError('Adaptive trajectory differs')
            exact=(anchors[level]==directory)
            gate=numerical(anchors[level],directory,cov['new_checkpoint'],cov['domains'],exact)
            result['gates'][label]=gate
            if not gate['passed']:raise ValueError('Cross-count/normal-layout numerical gate failed: '+label)
            key=(level,n)
            if key not in same_count:
                if kind!='baseline':raise ValueError('Same-rank baseline must run first')
                same_count[key]=directory
            reference=same_count[key]
            gate=numerical(reference,directory,cov['new_checkpoint'],cov['domains'],True)
            result['gates']['same-count-'+label]=gate
            if not gate['passed']:raise ValueError('Same-count exact checkpoint gate failed')
            actual,expected=vtk_members(directory),vtk_members(reference)
            save(directory/'vtk-members.json',[dict(archive=a,member=m,sha256=h) for (a,m),h in sorted(actual.items())])
            row['vtk']=dict(reference=reference.name,members=len(actual),exact=actual==expected)
            if len(actual)!=62 or actual!=expected:raise ValueError('Same-count exact VTK gate failed')
            row['checkpoint_sha256']=digest(directory/cov['new_checkpoint'])
            precise=directory/(identity['parameters']['run_id']+'_log')
            timing=ordinary(log,a,precise.read_text());row['ordinary']=timing
            if [v['index'] for v in timing['excluded']]!=[0,7,8]:
                raise ValueError('Timing window differs from the validated 12-step fixture')
            row['ordinary_step_median']=timing['median']
            if phase!='validation' and (oracle or kind=='checked'):raise ValueError('Instrumented timing run')
            persist();print(f"PASS {label}: ordinary median {timing['median']:.9g}s",flush=True)
            return timing['median']
        def checkpoint(decisions):result['decisions']=decisions;persist()
        decisions=campaign(execute,checkpoint)
        for level,identity in result['inputs'].items():
            for name,expected in [('simple.in',identity['input_sha256']),(identity['seed'],identity['seed_sha256']),
                                  ('grids/'+identity['grid'],identity['grid_sha256'])]:
                if digest(inputs/level/name)!=expected:raise ValueError('Input snapshot changed')
        for key,record in notes['builds'].items():
            if digest(inputs/'binaries'/key)!=record['sha256']:raise ValueError('Binary snapshot changed')
        decision=decisions.get('confirmation',decisions.get('transfer',decisions['screen']))['decision']
        result.update(status='passed',decision=decision,decisions=decisions);result.pop('active',None);persist()
        print('COMPLETE: '+decision+'; no automatic promotion.',flush=True)
    except BaseException as e:
        result.update(status='failed',error=repr(e));persist();raise


if __name__=='__main__':run(*(Path(p) for p in sys.argv[1:]))
