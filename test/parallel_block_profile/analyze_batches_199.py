"""Ordered, inclusive batch records. No cross-rank clock subtraction or causal speedup claim."""
from collections import defaultdict
from pathlib import Path
import math
import heapq
import statistics as st

EXCHANGES={108,109,113,114,115}
PRODUCERS={101,102,103,104,107,110,111,112,130}


def collect(root,ranks,mode,side=None):
    root=Path(root)
    files=list(root.glob('batches199-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing Stage199 batch rank files')
    records=[];seen=set();coverage=defaultdict(int);ordinary_steps=None
    for path in files:
        rank=int(path.stem.rsplit('-',1)[1])
        if rank not in range(ranks) or rank in seen:raise ValueError('Duplicate/invalid batch rank')
        seen.add(rank);lines=path.read_text().splitlines()
        if lines[0].split()!=['schema','1',str(mode),str(rank)]:raise ValueError('Invalid batch schema')
        from analyze_operations_199 import common
        steps=common(root/f'workflow195-rank-{rank}.txt')
        ordinary=[n for n,v in steps.items() if not v['restart']]
        if ordinary_steps is not None and ordinary!=ordinary_steps:raise ValueError('Batch rank ordinary steps differ')
        ordinary_steps=ordinary
        byseq={};active=[];previous_start=(-1,-1.)
        for seq,line in enumerate(lines[1:],1):
            f=line.split()
            if len(f)!=20:raise ValueError('Invalid batch field count')
            ints=list(map(int,f[:17]));start,end,mpi=map(float,f[17:])
            n,c,index,parent,op,k,l,pos,closed=ints[:9]
            if index!=seq or n not in steps or c not in range(9) or op not in range(101,136) or closed!=1 or parent>=seq or parent<0:
                raise ValueError('Invalid batch identity')
            if any(v<0 for v in ints[9:]) or not all(math.isfinite(v) and v>=0 for v in (start,end,mpi)) or end<start or end>steps[n]['wall']+1e-6 or mpi>end-start+1e-6:
                raise ValueError('Invalid batch duration/work')
            if (n,start)<previous_start:raise ValueError('Batch start order differs')
            previous_start=(n,start)
            while active and active[-1]['sequence']!=parent:
                if active[-1]['step']==n and active[-1]['end']>start+1e-9:raise ValueError('Overlapping sibling batches')
                active.pop()
            if parent!=(active[-1]['sequence'] if active else 0):raise ValueError('Broken batch parent tree')
            if parent:
                p=byseq[parent]
                if (p['step'],p['context'])!=(n,c) or end>p['end']+1e-9:raise ValueError('Batch escapes parent')
            r=dict(rank=rank,step=n,restart=steps[n]['restart'],context=c,sequence=seq,parent=parent,operation=op,
                   layer=k,level=l,position=pos,work=ints[9],storage_counts=ints[10:17],start=start,end=end,
                   seconds=end-start,mpi_seconds=mpi,cpu_seconds=max(0.,end-start-mpi))
            active.append(r);byseq[seq]=r;records.append(r);coverage[op]+=1
        if mode<2 and byseq:raise ValueError('Disabled batch recording produced records')
    if mode>=2 and side:
        required={101,105,106}
        if side=='block':required|={104,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,129,130,131,132,133,134,135}
        missing=required-set(coverage)
        if missing:raise ValueError('Missing Stage199 batch coverage: '+repr(sorted(missing)))
    # Aggregate each operation independently. These rows overlap by design.
    cells=defaultdict(lambda:defaultdict(lambda:[0.,0.,0,0]))
    ordered=defaultdict(list)
    for r in records:
        if r['restart']:continue
        key=(r['step'],r['context'],r['operation'],r['layer'],r['level'],r['position'])
        a=cells[key][r['rank']];a[0]+=r['cpu_seconds'];a[1]+=r['mpi_seconds'];a[2]+=r['work'];a[3]+=1
        ordered[r['rank'],r['step']].append(r)
    rows=[]
    for key,rankrows in sorted(cells.items()):
        vals=[rankrows.get(i,[0.,0.,0,0]) for i in range(ranks)]
        cpu=[v[0] for v in vals];wait=[v[1] for v in vals];work=[v[2] for v in vals]
        rows.append(dict(zip(('step','context','operation','layer','level','position'),key),
            cpu_mean=st.mean(cpu),cpu_max=max(cpu),cpu_max_minus_mean=max(cpu)-st.mean(cpu),
            mpi_mean=st.mean(wait),cpu_mpi_correlation=correlation(cpu,wait),cpu_work_correlation=correlation(cpu,work),
            rank_rows=[{'rank':i,'cpu':v[0],'mpi':v[1],'work':v[2],'calls':v[3]} for i,v in enumerate(vals)]))
    # Chronological association is explicit, not proof of data dependency.
    links=[]
    for key,stream in ordered.items():
        pending=[];completed={}
        for r in sorted(stream,key=lambda v:(v['start'],v['sequence'])):
            while pending and pending[0][0]<=r['start']:
                _,_,p=heapq.heappop(pending)
                completed[p['operation'],p['layer'],p['level']]=p
            if r['operation'] in EXCHANGES:
                eligible=[p for p in completed.values() if p['operation'] in PRODUCERS and p['end']<=r['start'] and
                          (r['layer']==0 or p['layer'] in (0,r['layer'])) and
                          (r['level']<=0 or p['level'] in (0,-1,r['level'],r['level']+1))]
                p=max(eligible,key=lambda p:p['end']) if eligible else None
                links.append({'rank':r['rank'],'step':r['step'],'exchange_sequence':r['sequence'],
                    'exchange_operation':r['operation'],'producer_sequence':p['sequence'] if p else None,
                    'producer_operation':p['operation'] if p else None,
                    'producer_cpu_seconds':p['cpu_seconds'] if p else None,
                    'intervening_seconds':r['start']-p['end'] if p else None,'exchange_mpi_seconds':r['mpi_seconds']})
            if r['operation'] in PRODUCERS:heapq.heappush(pending,(r['end'],r['sequence'],r))
    return {'records':len(records),'ordinary_steps':ordinary_steps,'coverage':dict(coverage),'rows':rows,'preceding_producers':links,
            'scope':'Batch times are inclusive and must not be added. MPI subset uses existing disjoint operation clocks. CPU includes unclassified transport/preparation. Keys 0/-1 mean all or unspecified layers/levels. Links identify nearest compatible preceding batches, not proven data dependencies. Rank clocks have independent origins.'}


def correlation(a,b):
    ma,mb=st.mean(a),st.mean(b)
    aa=sum((x-ma)**2 for x in a);bb=sum((x-mb)**2 for x in b)
    if aa==0 or bb==0:return None
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(aa*bb)


def campaign_summary(out,runs):
    """Compare independent inclusive batches, preserving calibration separation."""
    import json
    groups=defaultdict(lambda:defaultdict(list))
    for label,r in runs.items():
        if r.get('phase')!='diagnostic-work':continue
        report=json.loads((Path(out)/label/'batches-summary.json').read_text())
        side=r['side']
        bykey=defaultdict(list)
        for row in report['rows']:
            key=tuple(row[k] for k in ('context','operation','layer','level','position'))
            bykey[key].append(row)
        for key,rows in bykey.items():
            groups[key][side].append({name:sum(v[name] for v in rows)/len(report['ordinary_steps']) for name in ('cpu_mean','mpi_mean','cpu_max_minus_mean')})
    result=[]
    for key,sides in sorted(groups.items()):
        r=dict(zip(('context','operation','layer','level','position'),key))
        for side in ('block','legacy'):
            vals=sides.get(side,[])
            r[side]={name:st.mean(v[name] for v in vals) if vals else 0. for name in ('cpu_mean','mpi_mean','cpu_max_minus_mean')}
        r['cpu_excess']=r['block']['cpu_mean']-r['legacy']['cpu_mean']
        r['mpi_excess']=r['block']['mpi_mean']-r['legacy']['mpi_mean']
        result.append(r)
    return {'rows':sorted(result,key=lambda r:r['cpu_excess'],reverse=True),
            'scope':'Inclusive independent batch comparisons averaged over ordinary steps and two work passes. Do not sum rows. Refer to operation-only calibration before interpretation. No inferred optimization saving.'}
