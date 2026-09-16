"""Discover OS CPU topology; use nested, LLC-spread physical-core subsets."""
import json, os, socket, subprocess
from pathlib import Path

CASES=(('n21-s1',21,1),('n42-s1',42,1),('n83-s1',83,1),('n83-s2',83,2))

def expand(text):
    out=[]
    for item in text.strip().split(','):
        a,_,b=item.partition('-');out.extend(range(int(a),int(b or a)+1))
    return out

def topology():
    text=subprocess.check_output(['lscpu','--parse=CPU,CORE,SOCKET,NODE'],text=True)
    rows=[];allowed=os.sched_getaffinity(0)
    for line in text.splitlines():
        if not line or line.startswith('#'):continue
        cpu,core,sock,node=map(int,line.split(','))
        if cpu not in allowed:continue
        caches=[]
        for cache in Path(f'/sys/devices/system/cpu/cpu{cpu}/cache').glob('index*'):
            if (cache/'level').read_text().strip()=='3':caches.append(tuple(expand((cache/'shared_cpu_list').read_text())))
        if len(caches)!=1 or node<0:raise ValueError('Require an identified L3 cache and NUMA node for every CPU')
        rows.append(dict(cpu=cpu,core=core,socket=sock,node=node,llc=list(caches[0])))
    return {'host':socket.gethostname(),'rows':rows,'lscpu':subprocess.check_output(['lscpu'],text=True),'allowed_cpus':sorted(allowed)}

def plan(top):
    rows=top['rows'];sockets=sorted({r['socket'] for r in rows})
    if len(sockets)!=2:raise ValueError('This experiment requires exactly two allocated sockets')
    if len({(r['socket'],r['core']) for r in rows})!=len(rows):raise ValueError('Require one hardware thread per physical core')
    if len({r['cpu'] for r in rows})!=len(rows):raise ValueError('Duplicate topology CPU')
    orders={}
    for sock in sockets:
        groups={}
        for r in rows:
            if r['socket']==sock:groups.setdefault((r['node'],tuple(r['llc'])),[]).append(r)
        groups=[sorted(v,key=lambda r:r['cpu']) for _,v in sorted(groups.items())]
        orders[sock]=[g[i] for i in range(max(map(len,groups))) for g in groups if i<len(g)]
        if len(orders[sock])<83:raise ValueError('Need at least 83 allocated physical cores on each socket')
    result={}
    for name,n,ns in CASES:
        mapping=(orders[sockets[0]][:n] if ns==1 else [orders[sockets[i%2]][i//2] for i in range(n)])
        result[name]={'ranks':n,'sockets':ns,'host':top['host'],'mapping':[dict(r,rank=i) for i,r in enumerate(mapping)],
                      'policy':'Nested prefixes, round-robin across discovered LLC groups within each socket; two-socket ranks alternate sockets.'}
    return result

def verify(rows,case):
    if len(rows)!=case['ranks'] or {r['rank'] for r in rows}!=set(range(case['ranks'])):raise ValueError('Incomplete/duplicate rank placement')
    for row in rows:
        expected=case['mapping'][row['rank']]
        if row['hostname']!=case['host'] or row['affinity']!=[expected['cpu']]:raise ValueError('Rank affinity differs from prescribed OS CPU map')
        if row.get('returncode',0)!=0:raise ValueError('Nonzero child exit')
    return True
