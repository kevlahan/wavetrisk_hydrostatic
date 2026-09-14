"""Read all experimental legacy rank profiles; fail closed on missing records."""
from pathlib import Path
import statistics
import math


def collect(root,ranks):
    files=list(Path(root).glob('legacy-detail-rank-*.txt'))
    if len(files)!=ranks: raise ValueError('Missing legacy rank profile files')
    data={}; names=None
    for path in files:
        rank=int(path.stem.rsplit('-',1)[1]); regions={}; steps={}; ordinal={}; labels={}
        if rank in data or not 0<=rank<ranks: raise ValueError('Invalid legacy profile rank')
        for line in path.read_text().splitlines():
            f=line.split()
            if not f or f[0]=='#': continue
            if f[0]=='name': labels[int(f[1])]=' '.join(f[2:])
            elif f[0]=='region': regions[int(f[1])]=list(map(float,f[2:6]))
            elif f[0]=='step':
                seq=int(f[1])
                if seq in steps: raise ValueError('Duplicate legacy step')
                steps[seq]={'restarts':int(f[2]),'remaps':int(f[3]),'wall':float(f[4]),'cpu':float(f[5]),'self':{}}
                ordinal[len(ordinal)+1]=seq
            elif f[0]=='self':
                sequence=ordinal[int(f[1])]; i=int(f[2])
                if i in steps[sequence]['self']: raise ValueError('Duplicate legacy region/step')
                steps[sequence]['self'][i]=[float(f[3]),float(f[4])]
        expected=set(range(1,66))
        if set(regions)!=expected or set(labels)!=expected or not steps: raise ValueError('Incomplete legacy profile schema')
        if names is not None and labels!=names: raise ValueError('Legacy labels differ by rank')
        names=labels
        for s in steps.values():
            values=[s['wall'],s['cpu'],*(v for row in s['self'].values() for v in row)]
            if not all(math.isfinite(v) and v>=0 for v in values): raise ValueError('Nonfinite or negative legacy time')
            if set(s['self'])!=expected: raise ValueError('Incomplete legacy step')
            if abs(sum(v[0] for v in s['self'].values())-s['wall'])>1e-4: raise ValueError('Legacy exclusive times do not conserve step')
        data[rank]=steps
    sequences=set(data[0])
    if any(set(s)!=sequences for s in data.values()): raise ValueError('Legacy step coverage differs by rank')
    result=[]
    for seq in sorted(sequences):
        current=[data[r][seq] for r in range(ranks)]
        slow=max(range(ranks),key=lambda r:current[r]['wall'])
        result.append({'sequence':seq,'restart':any(s['restarts'] for s in current),'remap':any(s['remaps'] for s in current),
                       'wall_mean':statistics.mean(s['wall'] for s in current),'wall_max':current[slow]['wall'],'max_rank':slow,
                       'region_mean_self_wall':[statistics.mean(s['self'][i][0] for s in current) for i in range(1,66)],
                       'region_self_wall_on_max_rank':[current[slow]['self'][i][0] for i in range(1,66)]})
    return {'ranks':ranks,'regions':names,'steps':result,
            'scope':'Experimental legacy phase profile; labels share IDs with block code but some boundaries differ. Excluded from timing medians.'}
