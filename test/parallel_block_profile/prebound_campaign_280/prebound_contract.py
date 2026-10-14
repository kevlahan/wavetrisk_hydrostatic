"""Validate compact publication bindings from audited Stage280 preparation."""
from pathlib import Path

def evidence(directory,ranks,enabled):
    root=Path(directory);records={rank:[] for rank in range(ranks)}
    for line in (root/'run.log').read_text().splitlines():
        if not line.startswith('Transform280 publication:'):continue
        fields=[int(x) for x in line.split(':',1)[1].split()]
        if len(fields)!=4:raise ValueError('Compact publication width')
        rank,outputs,imports,local=fields
        if rank not in records or min(fields)<0 or local>outputs or imports>outputs-local:
            raise ValueError('Compact publication extent')
        records[rank].append([outputs,imports,local])
    if not enabled:
        if any(records.values()):raise ValueError('Compact publication audit in a timing/control run')
        return {'enabled':False}
    totals=[0,0,0]
    for rank,rows in records.items():
        cost=[[int(x) for x in line.split()[:13]] for line in (root/f'transform279-cost-rank-{rank}.txt').read_text().splitlines()]
        if len(rows)!=12 or len(cost)!=12:raise ValueError('Compact publication requires all12 preparations')
        for row,plan in zip(rows,cost):
            if row[0]!=plan[6]:raise ValueError('Publication output census differs')
            totals=[a+b for a,b in zip(totals,row)]
    return dict(enabled=True,epoch_rank_samples=12*ranks,logical_outputs=totals[0],unique_remote_rows=totals[1],
                local_canonical_outputs=totals[2],avoided_numerical_rows=totals[0]-totals[1],
                scope='Summed across12 epochs; not simultaneous memory. Actual RSS and named plan/numeric buffers remain separate evidence.')
