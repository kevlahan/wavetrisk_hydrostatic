"""Fail early if the frozen AMD core PMU group cannot run on all 83 ranks."""
import json,os,shutil,subprocess,sys
from pathlib import Path
from analyze_hardware_201 import read_rank

def preflight(kit,out,case_file,launch):
    dest=out/'preflight';dest.mkdir();record={'status':'running','required_cpu':'AuthenticAMD family25 model160','groups':{}}
    def save(): (dest/'capabilities.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        for cmd in ('mpicc',):
            if not shutil.which(cmd):raise ValueError('Missing '+cmd+'; ask cluster support for matching compiler tools. No solver runs started.')
        for p in ('/proc/sys/kernel/perf_event_paranoid','/proc/sys/kernel/nmi_watchdog','/proc/sys/kernel/perf_event_max_sample_rate','/sys/bus/event_source/devices/cpu/type'):
            try:record[p]=Path(p).read_text().strip()
            except OSError:record[p]=None
        exe=dest/'pmu-probe'
        cmd=['mpicc','-O2','-g','-DSTANDALONE201',str(kit/'workflow_clock_201.c'),'-o',str(exe)]
        record['compile']=cmd
        with (dest/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        # Test simultaneous per-task groups on the actual map before solver execution.
        for mode in (1,):
            directory=dest/f'group-{mode}';directory.mkdir()
            # Each process needs a separate cwd because the standalone report is rank0.
            command=launch(case_file,83,directory,sys.executable,str(kit/'preflight_202.py'),'--child',str(exe),str(directory),str(mode))
            env=os.environ.copy();env.pop('WAVETRISK_SAMPLE_RANKS201',None)
            with (directory/'run.log').open('w') as f:subprocess.run(command,cwd=directory,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
            reports=[read_rank(directory/f'rank-{rank}'/'hardware201-rank-0.txt') for rank in range(83)]
            record['groups'][str(mode)]={'ranks':83,'status':'passed','reports':reports};save()
        record['status']='passed';save();return record
    except BaseException as error:
        record.update(status='failed',error=repr(error));save();raise

if __name__=='__main__':
    if len(sys.argv)!=5 or sys.argv[1]!='--child':raise SystemExit('Internal preflight child only')
    exe,dest,mode=sys.argv[2:];directory=Path(dest)/('rank-'+os.environ['SLURM_PROCID']);directory.mkdir()
    env=os.environ.copy();env['WAVETRISK_PMU201']=mode
    raise SystemExit(subprocess.run([exe],cwd=directory,env=env).returncode)
