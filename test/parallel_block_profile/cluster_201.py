#!/usr/bin/env python3
"""Stage 201: distinguish representation work from memory stalls at fixed placement."""
import argparse,json,os,re,shutil,socket,struct,subprocess,sys,time
from pathlib import Path
from cluster_199 import (save,digest,copy_grids,validate_log,workload_gate,exact_gate,legacy_gate,stage191_evidence,stage194_evidence)
from profile_protocol import parameters,checkpoint_name,grid_identity,environment
from rank_metrics_190 import memory,rank_resources
from placement_200 import topology,plan,verify
from analyze_scaling_200 import collect
from analyze_hardware_201 import collect as hardware, samples, summarize, write_reports
from preflight_201 import preflight

KINDS=('candidate','checked','candidate_profile','legacy','legacy_profile')

def schedule(cases):
    case='n83-s2'
    jobs=[(case,'checked','validation',0,False),(case,'checked','validation',0,True)]
    jobs.extend((case,s,'warmup',0,False) for s in ('candidate','legacy'))
    for repeat in (1,2):
        jobs.extend((case,s,'bare',repeat,False) for s in (('candidate','legacy') if repeat==1 else ('legacy','candidate')))
    jobs.extend((case,s+'_profile','coarse',0,False) for s in ('candidate','legacy'))
    for repeat in (1,2):
        for phase in ('core','memory'):
            jobs.extend((case,s+'_profile',phase,repeat,False) for s in (('legacy','candidate') if repeat==1 else ('candidate','legacy')))
    jobs.extend((case,s+'_profile','sample',0,False) for s in ('legacy','candidate'))
    return jobs


def launch(case_file,ranks,*command):
    # Do not use -c/--exact/step --exclusive: they restrict the step CPU subset.
    # The job owns the node exclusively and all steps execute sequentially.
    return ['srun','--nodes=1',f'--ntasks={ranks}','--cpu-bind=none','--kill-on-bad-exit=1',
            sys.executable,str(Path(__file__).with_name('rank_launch_201.py').resolve()),str(case_file),*map(str,command)]

def system_snapshot(case):
    value=memory()
    value['loadavg']=Path('/proc/loadavg').read_text().strip()
    value['cpu_stat']=Path('/proc/stat').read_text()
    value['cpu_frequency_snapshots']={}
    for row in case['mapping']:
        policy=Path(f"/sys/devices/system/cpu/cpu{row['cpu']}/cpufreq")
        fields={}
        for name in ('scaling_driver','scaling_governor','scaling_min_freq','scaling_max_freq','scaling_cur_freq'):
            try:fields[name]=(policy/name).read_text().strip()
            except OSError:fields[name]=None
        value['cpu_frequency_snapshots'][str(row['cpu'])]=fields
    value['scope']='Before/after subprocess, node-wide counters and instantaneous frequencies; not timed-step hardware counters.'
    return value


def run(args):
    fixture=args.fixture.resolve(strict=True);out=args.out.resolve();kit=Path(__file__).resolve().parent
    binaries={k:getattr(args,k).resolve(strict=True) for k in KINDS}
    if out.exists() or out.is_relative_to(fixture) or fixture.is_relative_to(out):raise ValueError('Require a new output directory outside fixture')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=1:raise ValueError('Require one exclusive Slurm node')
    for command in ('srun','numactl','zstd','lscpu'):
        if not shutil.which(command):raise ValueError('Missing required command: '+command)
    for binary in binaries.values():
        if not os.access(binary,os.X_OK) or binary.is_relative_to(out):raise ValueError('Invalid executable')
    for key in ('SLURM_CPUS_PER_TASK','SRUN_CPUS_PER_TASK','SLURM_EXACT','SLURM_EXCLUSIVE'):
        os.environ.pop(key,None)
    text=(fixture/'simple.in').read_text();p=parameters(text)
    for key,value in {'max_level':'7','zlevels':'30','Nsoil':'10','physics_type':'Simple','NCAR_topo':'F','sso':'F','resume':'3','CP_EVERY':'1'}.items():
        if p.get(key)!=value:raise ValueError(f'Require {key}={value}')
    if float(p['time_end'])!=.315 or float(p['dt_write'])!=.01:raise ValueError('Use the validated time_end=0.315, dt_write=0.01 fixture')
    seed=checkpoint_name(text)
    with subprocess.Popen(['zstd','-dc',str(fixture/seed)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL) as process:
        header=process.stdout.read(24);process.stdout.close();process.wait()
    if len(header)!=24 or struct.unpack('<3q',header)!=(0x5741564554524953,1,160):raise ValueError('Require the J5 160-Domain checkpoint')
    notes=json.loads(args.build_notes_file.read_text())
    out.mkdir(parents=True)
    result={'stage':201,'status':'running','slurm_environment':{k:v for k,v in os.environ.items() if k.startswith('SLURM_')},
            'parameters':p,'build_notes':notes,'runs':{},'gates':{},'cases':{},'purpose':'Diagnose substantial legacy/block gap; no numerical optimization.'}
    save(out/'results.json',result)
    try:
        top=topology();cases={'n83-s2':plan(top)['n83-s2']};result['cases']=cases;save(out/'topology.json',top)
        if notes['build_host']['node']!=top['host']:raise ValueError('Build host and test host differ; native compilation is frozen to the tested node')
        result['memory_policy']='numactl --localalloc; CPU pinned before MPI initialization; MV2_ENABLE_AFFINITY=0'
        result['numactl_hardware']=subprocess.check_output(['numactl','--hardware'],text=True)
        for name,case in cases.items():
            path=out/(name+'-placement.json');save(path,case)
            probe=subprocess.run(launch(path,case['ranks'],'--probe'),capture_output=True,text=True,check=True)
            rows=[json.loads(line) for line in probe.stdout.splitlines() if line.startswith('{')];verify(rows,case)
            save(out/(name+'-verified-placement.json'),rows)
        result['capabilities']=preflight(kit,out,out/'n83-s2-placement.json',launch)
        result['sample_cohort']=json.loads((kit/'sample-cohort-201.json').read_text())
        inputs=out/'inputs';inputs.mkdir();grids=copy_grids(fixture/'grids',inputs/'grids')
        for name in ('simple.in',seed):shutil.copy2(fixture/name,inputs/name)
        for kind,binary in binaries.items():shutil.copy2(binary,inputs/kind)
        result.update(input_sha256=digest(inputs/'simple.in'),seed_sha256=digest(inputs/seed),grids=grids,
                      executables={k:{'source':str(b),'sha256':digest(inputs/k)} for k,b in binaries.items()})
        shutil.copy2(kit/'operation-instrumentation.json',out/'operation-instrumentation.json')
        jobs=schedule(cases);result['schedule']=[list(j) for j in jobs]
        for case_name,kind,phase,repeat,oracle in jobs:
            case=cases[case_name];side='legacy' if kind.startswith('legacy') else 'candidate'
            suffix=f'{phase}-{repeat}' if phase in ('bare','core','memory') else f'checked-{int(oracle)}' if phase=='validation' else phase
            label=f'{case_name}-{side}-{suffix}';directory=out/label;directory.mkdir()
            result['active']=label;save(out/'results.json',result);print('Starting '+label,flush=True)
            for name in ('simple.in',seed):shutil.copy2(inputs/name,directory/name)
            shutil.copy2(inputs/kind,directory/'climateJ5');copy_grids(inputs/'grids',directory/'grids',grids)
            env=environment(oracles=oracle,detail=False)
            env.update(WAVETRISK_STAGE191_INSTALL='1',WAVETRISK_STAGE191_MASS_OVERLAP='0',WAVETRISK_STAGE192_SPARSE_OUTER='0',
                       WAVETRISK_STAGE193_SCALAR_COLUMNS='1',WAVETRISK_STAGE194_INNER_GEOMETRY='1',WAVETRISK_STAGE197_RESIDENT_WAVELETS='1',
                       WAVETRISK_WORKFLOW195='1',WAVETRISK_SCALING200='0',MV2_ENABLE_AFFINITY='0',
                       WAVETRISK_PMU201=str({'core':1,'memory':2}.get(phase,0)),
                       WAVETRISK_SAMPLE_RANKS201=json.dumps(result['sample_cohort']['ranks'] if phase=='sample' else []))
            command=launch(out/(case_name+'-placement.json'),case['ranks'],directory,'./climateJ5','simple.in')
            before=system_snapshot(case);start=time.monotonic()
            with (directory/'run.log').open('x') as stream:code=subprocess.run(command,cwd=directory,env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
            record={'label':label,'case':case_name,'ranks':case['ranks'],'side':side,'kind':kind,'phase':phase,'round':repeat,'oracle':oracle,'returncode':code,
                    'elapsed_seconds':time.monotonic()-start,'command':command,'memory_before':before,'memory_after':system_snapshot(case),
                    'environment':{k:v for k,v in env.items() if k.startswith(('WAVETRISK_','MV2_')) or k.endswith('NUM_THREADS')},'binary_sha256':digest(directory/'climateJ5')}
            result['runs'][label]=record;save(out/'results.json',result)
            if code:raise RuntimeError(f'{label} exited {code}; see {directory}/run.log')
            record['rank_metrics']=rank_resources(directory);verify(record['rank_metrics'],case)
            analysis,coverage=validate_log(directory,p['run_id'],3,[4],True);save(directory/'analysis.json',analysis)
            record.update(coverage,step_count=len(analysis['step_records']),ordinary_step_median=analysis['non_checkpoint_median'])
            if record['ordinary_step_median'] is None:raise ValueError('Missing ordinary timing')
            if side=='candidate':
                record['stage191']=stage191_evidence(directory,analysis,'1','0',False);record['stage194']=stage194_evidence(directory,'1')
                log=(directory/'run.log').read_text()
                if re.findall(r'Stage 197 requested resident wavelets = ([01])',log)!=['1']:raise ValueError('Resident selector differs')
            reference=out/(f'{case_name}-{side}-warmup' if phase!='validation' else 'n83-s2-candidate-checked-0')
            if phase=='warmup':reference=directory
            gate=exact_gate(reference,directory,p['run_id'],[4]);gate.update(reference=str(reference),self_check=reference==directory)
            result['gates'][label]=gate
            if not gate['passed']:raise ValueError('Same-configuration exact gate failed: '+label)
            # One common work/time/grid sequence across every implementation and placement.
            anchor=out/'n83-s2-candidate-checked-0'
            work=workload_gate(json.loads((anchor/'analysis.json').read_text()),analysis);result['gates']['work-'+label]=work
            if not work['passed']:raise ValueError('Different adaptive workload: '+label)
            if phase=='warmup':
                if side=='candidate' and case_name=='n83-s2':cross=exact_gate(anchor,directory,p['run_id'],[4])
                elif side=='candidate' and case_name=='n83-s1':cross=exact_gate(out/'n83-s2-candidate-warmup',directory,p['run_id'],[4])
                elif side=='legacy':cross=legacy_gate(directory,out/(case_name+'-candidate-warmup'),p['run_id'],[4])
                else:
                    cross=legacy_gate(directory,out/'n83-s2-candidate-warmup',p['run_id'],[4])
                    cross['scope']='Cross-rank-count fixed numerical screen; identical printed workload still required. Not bitwise equality.'
                result['gates']['cross-'+label]=cross
                if not cross['passed']:raise ValueError('Cross-configuration/mode numerical gate failed: '+label)
                if side=='legacy' and case_name!='n83-s2':
                    ref=out/'n83-s2-legacy-warmup'
                    other=exact_gate(ref,directory,p['run_id'],[4]) if case_name=='n83-s1' else legacy_gate(ref,directory,p['run_id'],[4])
                    if case_name!='n83-s1':other['scope']='Same-implementation cross-rank-count fixed numerical screen; not bitwise equality.'
                    result['gates']['cross-count-'+label]=other
                    if not other['passed']:raise ValueError('Legacy cross-count/placement gate failed: '+label)
            if phase in ('coarse','core','memory','sample'):
                record['hardware']=hardware(directory,case['ranks'],{'core':1,'memory':2}.get(phase,0))
                if record['hardware']['quality_flags']:raise ValueError('Invalid PMU scheduling/count coverage; inspect hardware quality_flags')
                if phase=='sample':record['samples']=samples(directory,result['sample_cohort']['ranks'])
                record['scaling']=collect(directory,case['ranks'],int(phase=='probe'))
                if len(record['scaling']['workflow']['steps'])!=record['step_count']:raise ValueError('Observer/log step count differs')
                steps=record['scaling']['workflow']['steps']
                if sum(s['restart'] for s in steps)!=1:raise ValueError('Observer restart flag differs')
                if sum(s['remap'] for s in steps)!=(directory/'run.log').read_text().count('Remapping vertical coordinates'):raise ValueError('Observer remap flag differs')
            save(out/'results.json',result)
            print(f"PASS {label}: ordinary-step median {record['ordinary_step_median']:.4g}s",flush=True)
        if digest(inputs/'simple.in')!=result['input_sha256'] or digest(inputs/seed)!=result['seed_sha256'] or grid_identity(inputs/'grids')!=grids:raise ValueError('Input snapshot changed')
        result['summary']=summarize(result);result['status']='passed';result.pop('active',None);save(out/'results.json',result)
        save(out/'hardware-summary.json',result['summary'])
        write_reports(out,result)
        print('PASS: Stage 201 hardware campaign; results.json and hardware-summary.json ready for review.',flush=True)
    except BaseException as error:
        result.update(status='failed',error=repr(error));save(out/'results.json',result);raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in (*KINDS,'fixture','out','build-notes-file'):p.add_argument('--'+key,type=Path,required=True)
    run(p.parse_args())
if __name__=='__main__':main()
