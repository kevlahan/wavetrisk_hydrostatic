"""One native verification attempt, supervised outside its owned child group.

Only saved-capture Python verification is executable here. No compiler, solver,
resubmission, resume or performance entry point exists. Aggregate memory/disk
guards are sampled; each process also inherits a hard16GiB address-space limit.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import selectors
import signal
import socket
import stat
import subprocess
import sys
import time

sys.dont_write_bytecode = True
GiB = 1024**3
MiB = 1024**2
CAP_BYTES325 = 4*GiB
CAP_MEMORY325 = 16*GiB
WORKFLOW_BYTES325 = CAP_BYTES325-16*MiB
WALL_SECONDS325 = 1800
WORKER_SECONDS325 = 1740
SUPERVISOR_SECONDS325 = 1750
GRACE_SECONDS325 = 10
LOG_BYTES325 = 2*MiB
PRIOR_SECONDS325 = 2656


def require325(ok, reason):
    if not ok: raise ValueError('native.runner325.'+reason)


def sha325(path):
    with Path(path).open('rb') as source:
        return hashlib.file_digest(source,'sha256').hexdigest()


def save325(path, value):
    path=Path(path);data=(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    require325(len(data)<=MiB,'report_bound')
    temp=path.with_suffix(path.suffix+'.tmp')
    require325(not path.is_symlink() and not temp.is_symlink(),'report_symlink')
    with temp.open('wb') as output:
        output.write(data);output.flush();os.fsync(output.fileno())
    os.replace(temp,path)
    descriptor=os.open(path.parent,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(descriptor)
    finally:os.close(descriptor)


def package_check325(kit):
    kit=Path(kit).resolve(strict=True)
    manifest_path=kit/'native-files-sha256.json'
    manifest=json.loads(manifest_path.read_bytes())
    require325(type(manifest) is dict and len(manifest)>50,'package_manifest')
    actual={}
    for base,dirs,files in os.walk(kit,followlinks=False):
        for name in dirs+files:
            path=Path(base)/name;info=path.lstat()
            require325(not stat.S_ISLNK(info.st_mode),'package_symlink')
            require325(stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode),'package_type')
            if stat.S_ISREG(info.st_mode):actual[str(path.relative_to(kit))]=path
    require325(set(actual)==set(manifest)|{'native-files-sha256.json'},'package_inventory')
    for name,expected in manifest.items():
        rel=Path(name)
        require325(not rel.is_absolute() and '..' not in rel.parts and str(rel)==name,'package_path')
        require325(type(expected) is str and len(expected)==64 and sha325(actual[name])==expected,
                   'package_bytes:'+name)
    return sha325(manifest_path)


def allocation325(env=None, hostname=None, cpus=None, limits=None):
    env=os.environ if env is None else env
    hostname=socket.gethostname() if hostname is None else hostname
    cpus=sorted(os.sched_getaffinity(0)) if cpus is None else cpus
    limits=resource.getrlimit(resource.RLIMIT_AS) if limits is None else limits
    require325(env.get('SLURM_JOB_ID','').isdigit(),'job_identity')
    for key,want in [('SLURM_JOB_PARTITION','nk'),('SLURM_JOB_NUM_NODES','1'),
                     ('SLURM_NTASKS','1'),('SLURM_CPUS_PER_TASK','16')]:
        require325(env.get(key)==want,'allocation_'+key)
    require325(hostname.split('.')[0]=='bb13' and len(cpus)==len(set(cpus))==16,'host_or_cpu_mask')
    require325(tuple(limits)==(CAP_MEMORY325,CAP_MEMORY325),'hard_and_soft_address_limit')
    require325(platform.system()=='Linux','linux_required')
    return dict(job=int(env['SLURM_JOB_ID']),hostname=hostname,cpus=cpus,
                address_space_limits=list(limits),python=platform.python_version(),
                executable=sys.executable,sqlite=__import__('sqlite3').sqlite_version)


def install_address_limit325():
    """Set both bounds inside the actual step, rather than assuming inheritance."""
    before=resource.getrlimit(resource.RLIMIT_AS)
    require325(before[1]==resource.RLIM_INFINITY or before[1]>=CAP_MEMORY325,
               'inherited_hard_limit_below_required_cap:'+repr(before))
    resource.setrlimit(resource.RLIMIT_AS,(CAP_MEMORY325,CAP_MEMORY325))
    after=resource.getrlimit(resource.RLIMIT_AS)
    require325(after==(CAP_MEMORY325,CAP_MEMORY325),'installed_address_limit:'+repr(after))
    return dict(before=list(before),after=list(after),
                installation='inside Python verification step before allocation/capture checks')


def tree_bytes325(root):
    total=0;seen=set()
    for base,dirs,files in os.walk(root,followlinks=False):
        for name in dirs+files:
            path=Path(base)/name
            try:info=path.lstat()
            except FileNotFoundError:continue  # A committed rename can race this sample.
            require325(not stat.S_ISLNK(info.st_mode),'derived_symlink')
            require325(stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode),'derived_type')
            if stat.S_ISREG(info.st_mode):
                require325(info.st_nlink==1,'derived_hardlink')
                identity=(info.st_dev,info.st_ino)
                # A phase rename may be seen under both old and new paths by
                # a concurrent walk. Count the same regular file only once.
                if identity not in seen:total+=info.st_size;seen.add(identity)
    return total


def processes325(root_pid, *, group=None, proc=Path('/proc')):
    rows={}
    for entry in proc.iterdir():
        if not entry.name.isdigit():continue
        try:
            fields={line.split(':',1)[0]:line.split(':',1)[1].strip()
                    for line in (entry/'status').read_text().splitlines() if ':' in line}
            # /proc/stat fields after the last ')' begin at process state.
            bits=(entry/'stat').read_text().rsplit(')',1)[1].split()
            rows[int(entry.name)]=dict(ppid=int(fields['PPid']),group=int(bits[2]),
                virtual_bytes=int(fields.get('VmSize','0 kB').split()[0])*1024,
                rss_bytes=int(fields.get('VmRSS','0 kB').split()[0])*1024,state=bits[0])
        except (FileNotFoundError,ProcessLookupError):continue
        except PermissionError:continue  # Unrelated users need not expose /proc.
    owned={root_pid}
    if group is not None:owned.update(pid for pid,row in rows.items() if row['group']==group)
    while True:
        added={pid for pid,row in rows.items() if row['ppid'] in owned}
        if added<=owned:break
        owned.update(added)
    return {pid:rows[pid] for pid in sorted(owned) if pid in rows}


def kill_group325(pid, sig):
    try:os.killpg(pid,sig)
    except ProcessLookupError:pass


def stop_group325(child, *, grace=None):
    """Bound cleanup to this exact child session, including surviving grandchildren."""
    if grace is None:grace=GRACE_SECONDS325
    kill_group325(child.pid,signal.SIGTERM)
    deadline=time.monotonic()+grace
    while time.monotonic()<deadline:
        child.poll()
        rows=processes325(os.getpid(),group=child.pid)
        live=[pid for pid,row in rows.items() if row['group']==child.pid and row['state']!='Z']
        if not live:break
        time.sleep(.05)
    kill_group325(child.pid,signal.SIGKILL)
    try:code=child.wait(timeout=2)
    except subprocess.TimeoutExpired:
        # An uninterruptible I/O task may not reap even after SIGKILL. Preserve
        # that fact; never lose the supervisor report to another wait timeout.
        return None
    return code


def worker325(kit,out,*,address_limit_proof=None):
    limits=address_limit_proof or dict(after=list(resource.getrlimit(resource.RLIMIT_AS)),installation='caller supplied; verified below')
    package_check325(kit);environment=allocation325()
    save325(Path(out)/'worker-environment.json',dict(address_limit=limits,environment=environment))
    from resume_native325 import NativeBackend325
    from resume_workflow325 import run325
    backend=NativeBackend325(Path(kit)/'recovery_v4p1')
    # No explicit_resume/discard_partial path is offered by this native runner.
    report=run325(backend,Path(out)/'state',max_seconds=WORKER_SECONDS325,
                  cap_bytes=WORKFLOW_BYTES325)
    require325(report['status']=='passed' and report['native_gate_passed'] is True,
               'complete_workflow_required')
    save325(Path(out)/'worker-complete.json',dict(passed=True,phase_certificates=report['certificate_sha256'],
            native_qualification=False,independent_review_required=True,source_job=55569,
            solver_runs=0,compiler_runs=0,performance_runs=0))


def protected_paths325(kit):
    from resume_native325 import NativeBackend325
    return (Path(kit),*NativeBackend325(Path(kit)/'recovery_v4p1').protected_roots)


def protect_output325(work,out,protected):
    for root in protected:
        root=Path(root).resolve()
        require325(all(not candidate.is_relative_to(root) and not root.is_relative_to(candidate)
                       for candidate in (work,out)),'protected_output_overlap')


def completion_check325(out):
    from resume_state325 import PHASES325
    out=Path(out)
    completion=json.loads((out/'worker-complete.json').read_text())
    attempt=json.loads((out/'state.attempt.json').read_text())
    receipts=completion.get('phase_certificates',{})
    require325(completion.get('passed') is True and attempt.get('status')=='passed' and
               attempt.get('native_gate_passed') is True and set(receipts)==set(PHASES325) and
               receipts==attempt.get('certificate_sha256'),'complete_gate_certificate')
    for phase,digest in receipts.items():
        path=out/'state/phases'/phase/'certificate.json'
        require325(path.is_file() and not path.is_symlink() and path.stat().st_size<=MiB and
                   sha325(path)==digest,'phase_certificate_bytes:'+phase)
    return completion


def supervise325(kit,out,*,address_limit_proof=None):
    started=time.monotonic();kit=Path(kit).resolve(strict=True)
    out=Path(out).resolve();work=out.parent
    require325(not out.exists() and work.is_dir(),'fresh_output_required')
    first_signal=[];old_handlers={}
    def notice(signum,frame):
        if not first_signal:first_signal.append(dict(number=signum,name=signal.Signals(signum).name,
            monotonic=time.monotonic(),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    for sig in (signal.SIGTERM,signal.SIGINT,signal.SIGALRM):old_handlers[sig]=signal.signal(sig,notice)
    child=None;selector=selectors.DefaultSelector();log=None
    report=dict(stage=325,revision='native-resumable-v5p1',status='preflight',source_job=55569,
        started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),phase='allocation_and_package',
        prior_elapsed_seconds=PRIOR_SECONDS325,new_allocation_cap_seconds=WALL_SECONDS325,
        cumulative_maximum_seconds=PRIOR_SECONDS325+WALL_SECONDS325,worker_alarm_seconds=WORKER_SECONDS325,
        supervisor_deadline_seconds=SUPERVISOR_SECONDS325,derived_logical_byte_cap=CAP_BYTES325,
        per_process_address_space_cap=CAP_MEMORY325,sampled_aggregate_virtual_cap=CAP_MEMORY325,
        sampling_policy='0.25second periodic /proc and tree samples; peaks are sampled, not hard aggregate quotas',
        sampled_peak_derived_bytes=0,sampled_peak_virtual_bytes=0,sampled_peak_rss_bytes=0,
        stdout_bytes=0,retained_log_limit=LOG_BYTES325,builds=0,solver_runs=0,performance_runs=0,
        native_qualification=False,performance_qualified=False,automatic_retry=False,signals=first_signal)
    try:
        report['address_limit_before_installation']=(address_limit_proof or {}).get('before',list(resource.getrlimit(resource.RLIMIT_AS)))
        report['address_limit_installation']=address_limit_proof or dict(
            after=list(resource.getrlimit(resource.RLIMIT_AS)),installation='caller supplied; verified below')
        environment=allocation325();manifest=package_check325(kit)
        protect_output325(work,out,protected_paths325(kit))
        require325(work.stat().st_uid==os.getuid() and stat.S_IMODE(work.stat().st_mode)==0o700,'private_work_root')
        out.mkdir(mode=0o700)
        report['package_manifest_sha256']=manifest;report['environment']=environment
        save325(out/'environment.json',environment)
        save325(out/'supervisor-report.json',report)
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                 MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',TMPDIR=str(out))
        child=subprocess.Popen([sys.executable,'-B','-u',str(kit/'native_runner325.py'),'worker',str(kit),str(out)],
            cwd=kit,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
        report.update(status='running',phase='native_workflow',child_pid=child.pid,child_group=child.pid)
        os.set_blocking(child.stdout.fileno(),False);selector.register(child.stdout,selectors.EVENT_READ)
        log=(out/'native-worker.log').open('xb');retained=0;last_save=0
        while True:
            for key,events in selector.select(timeout=.25):
                data=os.read(key.fileobj.fileno(),65536)
                if data:
                    report['stdout_bytes']+=len(data);keep=data[:max(0,LOG_BYTES325-retained)]
                    log.write(keep);retained+=len(keep);log.flush()
                else:selector.unregister(key.fileobj)
            rows=processes325(os.getpid(),group=child.pid)
            vm=sum(row['virtual_bytes'] for row in rows.values());rss=sum(row['rss_bytes'] for row in rows.values())
            used=tree_bytes325(work)
            report['sampled_peak_derived_bytes']=max(report['sampled_peak_derived_bytes'],used)
            report['sampled_peak_virtual_bytes']=max(report['sampled_peak_virtual_bytes'],vm)
            report['sampled_peak_rss_bytes']=max(report['sampled_peak_rss_bytes'],rss)
            report['observed_pids']=sorted(rows);elapsed=time.monotonic()-started
            require325(vm<=CAP_MEMORY325,'aggregate_virtual_memory_cap')
            require325(used<=CAP_BYTES325,'derived_byte_cap')
            if first_signal:raise InterruptedError('supervisor_received_'+first_signal[0]['name'])
            if elapsed>=SUPERVISOR_SECONDS325:raise TimeoutError('supervisor_monotonic_deadline')
            if child.poll() is not None and not selector.get_map():break
            if elapsed-last_save>=5:
                report['elapsed_seconds']=elapsed;save325(out/'supervisor-report.json',report);last_save=elapsed
        report['child_returncode']=child.wait()
        require325(child.returncode==0,'verification_child_failed')
        completion=completion_check325(out)
        require325(package_check325(kit)==manifest,'package_changed_after_verification')
        report.update(status='passed',phase='complete',gate_passed=True,independent_review_required=True,
                      phase_certificates=completion['phase_certificates'],retention='all derived certificates/databases retained')
    except BaseException as error:
        report.update(status='failed',error=dict(type=type(error).__name__,message=str(error)[:4096]),
                      termination_initiator='supervisor' if child and child.poll() is None else 'child_or_preflight')
        raise
    finally:
        if child is not None:
            try:
                report['child_returncode']=stop_group325(child)
                report['cleanup_wait_timed_out']=report['child_returncode'] is None
                remaining=processes325(os.getpid(),group=child.pid)
                report['remaining_live_child_group_pids']=[pid for pid,row in remaining.items()
                    if row['group']==child.pid and row['state']!='Z']
                if report['remaining_live_child_group_pids'] or report['cleanup_wait_timed_out']:
                    report['status']='failed'
            except BaseException as cleanup_error:
                report.update(status='failed',cleanup_error=repr(cleanup_error)[:2048])
        attempt=out/'state.attempt.json'
        if attempt.is_file() and not attempt.is_symlink() and attempt.stat().st_size<=MiB:
            try:
                value=json.loads(attempt.read_bytes())
                report['worker_attempt_summary']={key:value.get(key) for key in ('status','phase','signal','error')}
            except (ValueError,OSError) as read_error:report['attempt_summary_error']=repr(read_error)[:512]
        if log is not None:log.flush();os.fsync(log.fileno());log.close()
        selector.close();report['elapsed_seconds']=time.monotonic()-started
        report['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
        if out.exists():save325(out/'supervisor-report.json',report)
        else:sys.stderr.write(json.dumps(report)+'\n')
        for sig,handler in old_handlers.items():signal.signal(sig,handler)
    require325(report['status']=='passed','supervisor_final_verdict')
    return report


if __name__=='__main__':
    require325(len(sys.argv)==4 and sys.argv[1] in ('supervise','worker'),'usage_supervise_or_worker_KIT_OUT')
    before=list(resource.getrlimit(resource.RLIMIT_AS))
    sys.stderr.write(json.dumps(dict(address_limit_bootstrap=dict(before=before,target=[CAP_MEMORY325]*2)))+'\n')
    sys.stderr.flush()
    try:
        proof=install_address_limit325()
    except BaseException as error:
        sys.stderr.write(json.dumps(dict(address_limit_bootstrap_failed=dict(before=before,
            actual=list(resource.getrlimit(resource.RLIMIT_AS)),error=repr(error))))+'\n')
        sys.stderr.flush();raise
    sys.stderr.write(json.dumps(dict(address_limit_bootstrap_complete=proof))+'\n');sys.stderr.flush()
    (supervise325 if sys.argv[1]=='supervise' else worker325)(*sys.argv[2:],address_limit_proof=proof)
