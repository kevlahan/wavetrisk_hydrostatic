"""Shared, fail-closed local experiment inputs and process/memory observations."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import threading
import time
from experiment import SWITCHES, digest
from analyze import analyze


def parameters(text):
    result = {}
    for line in text.splitlines():
        words=line.split('!',1)[0].split()
        if not words: continue
        if len(words)!=2 or words[0] in result: raise ValueError('Ambiguous input parameter: '+line)
        result[words[0]]=words[1]
    return result


def replace_parameter(text,name,value):
    text,n=re.subn(r'(?m)^(\s*'+re.escape(name)+r'\s+)\S+',lambda m:m[1]+str(value),text)
    if n!=1: raise ValueError('Expected one '+name)
    return text


def checkpoint_name(text):
    p=parameters(text)
    if int(p['resume'])<0: raise ValueError('An evolved checkpoint is required')
    run_id=p['run_id']
    if not re.fullmatch(r'[A-Za-z0-9_-]+',run_id): raise ValueError('Unsafe run_id')
    return f"{run_id}_checkpoint_{int(p['resume']):04d}.bin.zst"


def grid_identity(root):
    root=root.resolve(strict=True)
    h=hashlib.sha256();count=0
    def visit(directory, logical, ancestors):
        resolved=directory.resolve(strict=True)
        if resolved in ancestors: raise ValueError('Grid directory symlink cycle: '+str(directory))
        for p in sorted(directory.iterdir()):
            relative=logical/p.name
            p.resolve(strict=True)  # Reject dangling links, including directory assets.
            if p.is_dir():
                yield from visit(p,relative,ancestors|{resolved})
            elif p.is_file():
                yield relative,p
            else: raise ValueError('Unsupported grid asset: '+str(p))
    for relative,p in visit(root,Path(),set()):
        h.update(str(relative).encode()+b'\0'+digest(p).encode()+b'\n');count+=1
    if not count: raise ValueError('Empty grid assets')
    return {'files':count,'sha256':h.hexdigest()}


def fixture_identity(root):
    text=(root/'simple.in').read_text();p=parameters(text)
    p.pop('time_end',None)
    return {'parameters':p,'checkpoint':checkpoint_name(text),
            'checkpoint_sha256':digest(root/checkpoint_name(text)), 'grids':grid_identity(root/'grids')}


def environment(oracles=False,detail=False):
    env={k:v for k,v in os.environ.items() if not k.startswith('WAVETRISK_')}
    env.update({key:'0' for key in SWITCHES})
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key]='1'
    if oracles:
        for name in ('DYNAMICS','ADAPTATION','REMAP'): env['WAVETRISK_VALIDATE_BLOCK_'+name]='1'
    if detail:
        env['WAVETRISK_PROFILE_BLOCK_DETAIL']='1';env['WAVETRISK_PROFILE_PARALLEL_BLOCKS']='1'
    return env


def vm_snapshot():
    result=subprocess.run(['vm_stat'],capture_output=True,text=True)
    if result.returncode: raise RuntimeError(result.stderr)
    page=re.search(r'page size of (\d+) bytes',result.stdout)
    if not page: raise ValueError('Unknown vm_stat format')
    fields={m[1]:int(m[2]) for m in re.finditer(r'^([^:\n]+):\s+(\d+)\.',result.stdout,re.M)}
    return {'page_bytes':int(page[1]),'counters':fields}


def child_rss(text,launcher):
    rows=[]
    for line in text.splitlines():
        f=line.split(None,3)
        if len(f)==4 and all(x.isdigit() for x in f[:3]): rows.append((int(f[0]),int(f[1]),int(f[2]),f[3]))
    descendants={launcher}
    while True:
        new=descendants|{pid for pid,parent,_,_ in rows if parent in descendants}
        if new==descendants: break
        descendants=new
    return {str(pid):rss for pid,_,rss,cmd in rows if pid in descendants and Path(cmd).name=='climateJ5'}


class MemoryMonitor:
    """System-wide VM deltas are a confounder screen, not application DRAM traffic."""
    def __init__(self,root):
        self.root=root;self.stop=threading.Event();self.samples=[];self.errors=[];self.thread=None
        try: self.before=vm_snapshot()
        except (OSError,ValueError,RuntimeError) as e: self.before=None;self.errors.append(str(e))

    def start(self,pid):
        def observe():
            while not self.stop.is_set():
                try:
                    p=subprocess.run(['ps','-axo','pid,ppid,rss,comm'],capture_output=True,text=True)
                    if p.returncode: raise RuntimeError(p.stderr)
                    rss=child_rss(p.stdout,pid)
                    self.samples.append({'time':time.monotonic(),'rss_kib':rss})
                except (OSError,RuntimeError) as e:
                    self.errors.append(str(e));return
                self.stop.wait(2)
        self.thread=threading.Thread(target=observe,daemon=True);self.thread.start()

    def finish(self):
        self.stop.set()
        if self.thread: self.thread.join()
        after=None
        try: after=vm_snapshot()
        except (OSError,ValueError,RuntimeError) as e: self.errors.append(str(e))
        deltas={}
        if self.before and after:
            for key in ('Decompressions','Compressions','Swapins','Swapouts'):
                if key not in self.before['counters'] or key not in after['counters']:
                    self.errors.append('Missing VM counter '+key);continue
                delta=after['counters'][key]-self.before['counters'][key]
                if delta<0: self.errors.append('VM counter reset '+key)
                deltas[key]=delta*after['page_bytes']
        reasons=list(self.errors)
        if deltas.get('Decompressions',0)>1024**3: reasons.append('More than 1 GiB system-wide decompression during run')
        if deltas.get('Swapouts',0)>64*1024**2: reasons.append('More than 64 MiB system-wide swapout during run')
        peak=max((sum(s['rss_kib'].values()) for s in self.samples),default=0)
        if not peak: reasons.append('No child RSS samples captured')
        result={'scope':'whole-run system-wide VM deltas; sampled simulation-child RSS, excludes compressed footprint',
                'vm_delta_bytes':deltas,'sampled_simultaneous_rss_peak_kib':peak,
                'usable_for_local_timing':not reasons,'reasons':reasons,'samples':self.samples}
        (self.root/'memory.json').write_text(json.dumps(result,indent=2)+'\n')
        return result


def execute(binary,fixture,out,ranks,time_end,oracles=False,detail=False,monitor=False,overrides=None):
    binary=binary.resolve(strict=True);fixture=fixture.resolve(strict=True);out=out.resolve()
    if out.is_relative_to(fixture) or fixture.is_relative_to(out):
        raise ValueError('Execution output must be outside the immutable fixture')
    out.mkdir(parents=True,exist_ok=False)
    text=(fixture/'simple.in').read_text()
    name=checkpoint_name(text)
    for k,v in (overrides or {}).items(): text=replace_parameter(text,k,v)
    text=replace_parameter(text,'time_end',time_end)
    if checkpoint_name(text)!=name: raise ValueError('Do not relabel the seed checkpoint')
    (out/'simple.in').write_text(text)
    shutil.copy2(fixture/name,out/name)
    (out/'grids').symlink_to((fixture/'grids').resolve(strict=True),target_is_directory=True)
    shutil.copy2(binary,out/'climateJ5')
    env=environment(oracles,detail)
    identity={'binary':str(binary),'sha256':digest(binary),'ranks':ranks,'time_end':str(time_end),
              'parameters':parameters(text),'seed_sha256':digest(out/name),
              'environment':{k:v for k,v in env.items() if k.startswith('WAVETRISK_') or k.endswith('NUM_THREADS') or k=='VECLIB_MAXIMUM_THREADS'}}
    print('Starting',out.name,flush=True)
    memory=MemoryMonitor(out) if monitor else None
    start=time.monotonic()
    try:
        with (out/'run.log').open('x') as log:
            process=subprocess.Popen(['mpirun','-n',str(ranks),'./climateJ5','simple.in'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
            if memory: memory.start(process.pid)
            code=process.wait()
    finally:
        observation=memory.finish() if memory else None
    identity.update(returncode=code,elapsed_seconds=time.monotonic()-start)
    (out/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    report=analyze(out/'run.log')
    if code or not report['completed']: raise RuntimeError(f'Run failed: {out}/run.log')
    print('Completed',out.name,flush=True)
    return {'directory':str(out),'identity':identity,'analysis':report,'memory':observation}
