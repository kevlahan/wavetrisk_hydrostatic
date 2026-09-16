"""Fail early if the frozen AMD events or DWARF sampling cannot be measured."""
import json,os,shutil,subprocess,sys
from pathlib import Path
from rank_launch_201 import perf_command,export_perf
from analyze_hardware_201 import read_rank,parse_samples

def preflight(kit,out,case_file,launch):
    dest=out/'preflight';dest.mkdir();record={'status':'running','required_cpu':'AuthenticAMD family25 model160','groups':{},'sampling':{}}
    def save(): (dest/'capabilities.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        for cmd in ('mpicc','perf'):
            if not shutil.which(cmd):raise ValueError('Missing '+cmd+'; ask cluster support for matching compiler/perf tools. No solver runs started.')
        for cmd,name in ((['perf','version'],'perf-version.txt'),(['perf','list','--details'],'perf-events.txt'),(['uname','-a'],'kernel.txt')):
            with (dest/name).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        for p in ('/proc/sys/kernel/perf_event_paranoid','/proc/sys/kernel/nmi_watchdog','/proc/sys/kernel/perf_event_max_sample_rate','/sys/bus/event_source/devices/cpu/type'):
            try:record[p]=Path(p).read_text().strip()
            except OSError:record[p]=None
        exe=dest/'pmu-probe'
        cmd=['mpicc','-O2','-g','-DSTANDALONE201',str(kit/'workflow_clock_201.c'),'-o',str(exe)]
        record['compile']=cmd
        with (dest/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        # Test simultaneous per-task groups on the actual map before solver execution.
        for mode in (1,2):
            directory=dest/f'group-{mode}';directory.mkdir()
            # Each process needs a separate cwd because the standalone report is rank0.
            command=launch(case_file,83,directory,sys.executable,str(kit/'preflight_201.py'),'--child',str(exe),str(directory),str(mode))
            env=os.environ.copy();env.pop('WAVETRISK_SAMPLE_RANKS201',None)
            with (directory/'run.log').open('w') as f:subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
            reports=[read_rank(directory/f'rank-{rank}'/'hardware201-rank-0.txt') for rank in range(83)]
            record['groups'][str(mode)]={'ranks':83,'status':'passed','reports':reports};save()
        # Sample a CPU-bound test on a prescribed CPU; prove both events and unwind frames exist.
        case=json.loads(case_file.read_text());cpu=case['mapping'][0]['cpu']
        directory=dest/'sampling';directory.mkdir();path=directory/'perf-probe.data'
        env=os.environ.copy();env['WAVETRISK_PMU201']='0'
        with (directory/'record.log').open('w') as f:
            subprocess.run(perf_command(path,['taskset','-c',str(cpu),str(exe)]),cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        export_perf(path);samples,lost=parse_samples(path.with_suffix('.txt').read_text())
        if lost or any(not any(x['event']==e and len(x['frames'])>=2 for x in samples) for e in ('cycles:u','instructions:u')):
            raise ValueError('Preflight needs both cycles/instructions and usable DWARF stacks; inspect sampling/*.txt')
        record['sampling']={'status':'passed','samples':len(samples),'period':20000000,'stack_bytes':8192}
        record['status']='passed';save();return record
    except BaseException as error:
        record.update(status='failed',error=repr(error));save();raise

if __name__=='__main__':
    if len(sys.argv)!=5 or sys.argv[1]!='--child':raise SystemExit('Internal preflight child only')
    exe,dest,mode=sys.argv[2:];directory=Path(dest)/('rank-'+os.environ['SLURM_PROCID']);directory.mkdir()
    env=os.environ.copy();env['WAVETRISK_PMU201']=mode
    raise SystemExit(subprocess.run([exe],cwd=directory,env=env).returncode)
