#!/usr/bin/env python3
"""Stage215: correctness, 16-task AB, early 83-task AB, then fresh ABBA."""
import argparse,json,os,re,shutil,struct,subprocess,sys,time
from pathlib import Path
from cluster_199 import save,digest,copy_grids,validate_log,workload_gate,exact_gate,legacy_gate,stage191_evidence,stage194_evidence
from profile_protocol import parameters,checkpoint_name,grid_identity,environment
from rank_metrics_190 import memory,rank_resources
from placement_200 import topology,plan,verify
from protocol_215 import POLICY,ordinary,screen,confirm

KINDS=('baseline','candidate','checked')

def launch(case_file,ranks,*command):
    return ['srun','--nodes=1',f'--ntasks={ranks}','--cpu-bind=none','--kill-on-bad-exit=1',
            sys.executable,str(Path(__file__).with_name('rank_launch_201.py').resolve()),str(case_file),*map(str,command)]

def cases_for(top):
    old=plan(top);small=dict(old['n21-s1']);small.update(ranks=16,mapping=old['n21-s1']['mapping'][:16])
    return {'n16-s1':small,'n83-s2':old['n83-s2']}

def campaign(execute,checkpoint):
    """Injectable sequential scheduler; no later runs after an early stop."""
    for kind,oracle,label in [('baseline',False,'validate-A'),('checked',False,'validate-B-audit'),
                              ('candidate',True,'validate-B-oracles'),('candidate',False,'validate-B')]:
        execute('n16-s1',kind,'validation',label,oracle)
    decisions={}
    for case,phase in [('n16-s1','screen'),('n83-s2','transfer')]:
        a=execute(case,'baseline',phase,'A',False)
        b=execute(case,'candidate',phase,'B',False)
        decisions[phase]=screen(a,b)
        if phase=='transfer' and decisions[phase]['decision']=='advance' and a-b<POLICY['acceptance_seconds']-POLICY['near_threshold_seconds']:
            decisions[phase]['decision']='stop_below_target'
        checkpoint(decisions)
        if decisions[phase]['decision']!='advance':return decisions
    batches=[]
    while True:
        batch=[]
        for i,kind in enumerate(('baseline','candidate','candidate','baseline'),1):
            value=execute('n83-s2',kind,'confirm',f'{len(batches)+1}-{i}',False)
            batch.append((kind,value))
        batches.append(batch);decisions['confirmation']=confirm(batches);checkpoint(decisions)
        if decisions['confirmation']['decision']!='repeat_abba':return decisions


def run(args):
    fixture=args.fixture.resolve(strict=True);out=args.out.resolve()
    binaries={k:getattr(args,k).resolve(strict=True) for k in KINDS}
    if out.exists() or out.is_relative_to(fixture) or fixture.is_relative_to(out):raise ValueError('Require fresh output outside fixture')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:raise ValueError('Require one exclusive Slurm node')
    for command in ('srun','numactl','zstd','lscpu'):
        if not shutil.which(command):raise ValueError('Missing '+command)
    for key in ('SLURM_CPUS_PER_TASK','SRUN_CPUS_PER_TASK','SLURM_EXACT','SLURM_EXCLUSIVE'):os.environ.pop(key,None)
    text=(fixture/'simple.in').read_text();p=parameters(text)
    for key,value in {'max_level':'7','zlevels':'30','Nsoil':'10','physics_type':'Simple','NCAR_topo':'F','sso':'F','resume':'3','CP_EVERY':'1'}.items():
        if p.get(key)!=value:raise ValueError(f'Require {key}={value}')
    if float(p['time_end'])!=.315 or float(p['dt_write'])!=.01:raise ValueError('Require validated .315/.01 fixture')
    seed=checkpoint_name(text)
    with subprocess.Popen(['zstd','-dc',str(fixture/seed)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL) as process:
        header=process.stdout.read(24);process.stdout.close();process.wait()
    if len(header)!=24 or struct.unpack('<3q',header)!=(0x5741564554524953,1,160):raise ValueError('Require J5 160-Domain checkpoint')
    notes=json.loads(args.build_notes_file.read_text())
    if notes.get('status')!='passed':raise ValueError('Builds have not passed')
    for kind,binary in binaries.items():
        if digest(binary)!=notes['builds'][kind]['sha256']:raise ValueError('Binary differs from build record')
    out.mkdir(parents=True)
    result=dict(stage=215,status='running',policy=POLICY,build_notes=notes,parameters=p,runs={},gates={},schedule=[],decisions={},
                slurm_environment={k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
                scope='Compile generic Domain boundary integer routes per topology and level interval; direct packing; no profiled timing runs.')
    def persist():save(out/'results.json',result)
    persist()
    try:
        top=topology();cases=cases_for(top);result['cases']=cases;save(out/'topology.json',top)
        if notes['build_host']['node']!=top['host']:raise ValueError('Build/test host differs for native compilation')
        result['memory_policy']='CPU pinned before MPI init; numactl --localalloc; MV2_ENABLE_AFFINITY=0'
        for name,case in cases.items():
            path=out/(name+'-placement.json');save(path,case)
            probe=subprocess.run(launch(path,case['ranks'],'--probe'),capture_output=True,text=True,check=True)
            rows=[json.loads(line) for line in probe.stdout.splitlines() if line.startswith('{')];verify(rows,case)
            save(out/(name+'-verified-placement.json'),rows)
        inputs=out/'inputs';inputs.mkdir();grids=copy_grids(fixture/'grids',inputs/'grids')
        for name in ('simple.in',seed):shutil.copy2(fixture/name,inputs/name)
        if (inputs/'simple.in').read_text()!=text:raise ValueError('Input changed while snapshotting')
        for kind,binary in binaries.items():
            shutil.copy2(binary,inputs/kind)
            if digest(inputs/kind)!=notes['builds'][kind]['sha256']:raise ValueError('Executable snapshot differs')
        result.update(input_sha256=digest(inputs/'simple.in'),seed_sha256=digest(inputs/seed),grids=grids,
                      executables={k:dict(source=str(b),sha256=digest(inputs/k)) for k,b in binaries.items()})
        anchors={};analyses={};ordinary_states={}
        def execute(case_name,kind,phase,suffix,oracle):
            case=cases[case_name];label=f'{case_name}-{phase}-{suffix}'
            directory=out/label;directory.mkdir();result['active']=label
            result['schedule'].append(dict(label=label,kind=kind,oracle=oracle));persist();print('Starting '+label,flush=True)
            for name in ('simple.in',seed):shutil.copy2(inputs/name,directory/name)
            shutil.copy2(inputs/kind,directory/'climateJ5');copy_grids(inputs/'grids',directory/'grids',grids)
            env=environment(oracles=oracle,detail=False)
            env.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',WAVETRISK_STAGE192_SPARSE_OUTER='0',
                       WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1',
                       MV2_ENABLE_AFFINITY='0',WAVETRISK_SAMPLE_RANKS201='[]')
            command=launch(out/(case_name+'-placement.json'),case['ranks'],directory,'./climateJ5','simple.in')
            before=memory();start=time.monotonic()
            with (directory/'run.log').open('x') as f:code=subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
            record=dict(label=label,case=case_name,ranks=case['ranks'],kind=kind,phase=phase,oracle=oracle,returncode=code,
                        elapsed_seconds=time.monotonic()-start,command=command,memory_before=before,memory_after=memory(),
                        environment={k:v for k,v in env.items() if k.startswith(('WAVETRISK_','MV2_')) or k.endswith('NUM_THREADS')},
                        binary_sha256=digest(directory/'climateJ5'))
            result['runs'][label]=record;persist()
            if record['binary_sha256']!=result['executables'][kind]['sha256']:raise ValueError('Run binary differs from frozen snapshot')
            if code:raise RuntimeError(f'{label} exited {code}')
            record['rank_metrics']=rank_resources(directory);verify(record['rank_metrics'],case)
            analysis,coverage=validate_log(directory,p['run_id'],3,[4],True);save(directory/'analysis.json',analysis)
            log=(directory/'run.log').read_text();timing=ordinary(log,analysis)
            record.update(coverage,step_count=len(analysis['step_records']),ordinary=timing,ordinary_step_median=timing['median'])
            record['stage191']=stage191_evidence(directory,analysis,'1','0',False);record['stage194']=stage194_evidence(directory,'1')
            if re.findall(r'Stage 197 requested resident wavelets = ([01])',log)!=['1']:raise ValueError('Resident selector differs')
            if case_name not in anchors:
                if kind!='baseline':raise ValueError('Baseline must anchor each task count')
                anchors[case_name]=directory;analyses[case_name]=analysis
                ordinary_states[case_name]=[r['state'] for r in timing['included']]
            gate=exact_gate(anchors[case_name],directory,p['run_id'],[4]);result['gates'][label]=gate
            if not gate['passed']:raise ValueError('Exact same-count checkpoint gate failed: '+label)
            work=workload_gate(analyses['n16-s1'],analysis);result['gates']['work-'+label]=work
            if not work['passed']:raise ValueError('Adaptive workload differs')
            if [r['state'] for r in timing['included']]!=ordinary_states[case_name]:raise ValueError('Ordinary timestep selection differs')
            if case_name=='n83-s2' and directory==anchors[case_name]:
                cross=legacy_gate(anchors['n16-s1'],directory,p['run_id'],[4]);cross['scope']='Fixed cross-rank-count tolerance; not bitwise.'
                result['gates']['cross-count']=cross
                if not cross['passed']:raise ValueError('Cross-count numerical screen failed')
            persist();print(f"PASS {label}: ordinary median {timing['median']:.6g}s",flush=True)
            return timing['median']
        def checkpoint(decisions):result['decisions']=decisions;persist()
        decisions=campaign(execute,checkpoint)
        if digest(inputs/'simple.in')!=result['input_sha256'] or digest(inputs/seed)!=result['seed_sha256'] or grid_identity(inputs/'grids')!=grids:
            raise ValueError('Input snapshot changed')
        for kind in KINDS:
            if digest(inputs/kind)!=result['executables'][kind]['sha256']:raise ValueError('Executable snapshot changed')
        decision=decisions.get('confirmation',decisions.get('transfer',decisions['screen']))['decision']
        result.update(status='passed',measurement_status=decision,decisions=decisions);result.pop('active',None);persist()
        save(out/'decision.json',dict(status=decision,policy=POLICY,decisions=decisions))
        print('COMPLETE: '+decision+'; correctness PASS alone is not optimization acceptance.',flush=True)
    except BaseException as e:result.update(status='failed',error=repr(e));persist();raise

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in (*KINDS,'fixture','out','build-notes-file'):p.add_argument('--'+key,type=Path,required=True)
    run(p.parse_args())
