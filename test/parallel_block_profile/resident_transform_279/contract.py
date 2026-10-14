"""Check authoritative RHS/RK work, bounded adapter volumes and all rank participation."""
from pathlib import Path
import math

def evidence(directory,ranks,enabled,domains=160):
    directory=Path(directory);paths=sorted(directory.glob('rhs278-production-rank-*.txt'))
    if not enabled:
        if paths or list(directory.glob('rhs278-cost-rank-*.txt')):raise ValueError('RHS audit output in bare control')
        return dict(enabled=False)
    if len(paths)!=ranks:raise ValueError('Missing RHS rank evidence')
    if len(list(directory.glob('rhs278-cost-rank-*.txt')))!=ranks:raise ValueError('Missing replacement cost evidence')
    ids=set();zero_ranks=zero_work=0;ndom=0;records=[]
    for p in paths:
        rank=int(p.stem.rsplit('-',1)[1]);ids.add(rank)
        rows=[]
        for line in p.read_text().splitlines():
            fields=line.split()
            if len(fields)!=18:raise ValueError('RHS audit width differs')
            row=[int(x) for x in fields[:16]]+[float(x) for x in fields[16:]]
            if any(x<0 for x in row) or any(not math.isfinite(x) for x in row):raise ValueError('Invalid RHS counter')
            rows.append(row)
        if len(rows)!=48 or [r[0] for r in rows]!=list(range(1,49)) or any(r[0]!=r[1] for r in rows):
            raise ValueError('Incomplete authoritative RHS/RK transaction sequence')
        cost=cost_evidence(directory,rank,rows)
        previous=None;work=[0]*8
        for row in rows:
            same=previous is not None and row[2:4]==previous[2:4]
            for j,column in enumerate((6,7,8,9,10,11,12,13)):
                delta=row[column]-(previous[column] if same else 0)
                if delta<0:raise ValueError('Counter regressed within RHS plan epoch')
                work[j]+=delta
            if row[5]>0 and (row[6]==0 or row[7]==0 or row[8]==0):raise ValueError('Block owner performed no connected numerical work')
            previous=row
        if work[3]!=48:raise ValueError('RHS not retained for every real RK stage')
        ndom+=rows[-1][4]
        if rows[-1][4]==0:
            zero_ranks+=1;zero_work+=work[0]
            if rows[-1][5]>0 and work[0]==0:raise ValueError('Zero-Domain block owner did no RHS work')
        records.append(dict(rank=rank,replacement_cost=cost,domains=rows[-1][4],blocks=rows[-1][5],epochs=len(set(tuple(r[2:4]) for r in rows)),
                            primitive_events=work[0],scalar_targets=work[1],velocity_targets=work[2],
                            input_adapter_values=work[4],input_resident_values=work[5],derived_sent_values=work[6],
                            derived_messages=work[7],rk_adapter_values=rows[-1][14],rk_resident_values=rows[-1][15],evaluate_seconds=rows[-1][16],rk_adapter_seconds=rows[-1][17]))
    for kind in ('input','rk'):
        if sum(r['replacement_cost'][kind+'_sent_values'] for r in records)!=sum(r['replacement_cost'][kind+'_received_values'] for r in records):
            raise ValueError('Global replacement traffic does not balance: '+kind)
    if ids!=set(range(ranks)) or ndom!=domains:raise ValueError('Rank/Domain census differs')
    if ranks>domains and (not zero_ranks or not zero_work):raise ValueError('No above-Domain numerical participation')
    return dict(enabled=True,ownership_gate_passed=True,rhs_per_rank=48,layers=30,zero_domain_ranks=zero_ranks,
                zero_domain_rank_primitive_events=zero_work,ranks=sorted(records,key=lambda r:r['rank']),
                scope='Authoritative RHS retained by real resident RK; old Domain RHS and capture entry guards active. Local/cluster science gates are separate. RK adapter counts include raw five-field bounded reader values and five-field resident deliveries (only temperature installed). Input/derived counters exclude topology preparation and metadata; whole-step timing includes all costs. No speedup claim.')


def cost_evidence(directory,rank,stages):
    path=directory/f'rhs278-cost-rank-{rank}.txt';plans={}
    for line in path.read_text().splitlines():
        fields=line.split()
        if len(fields)!=18:raise ValueError('Replacement cost width differs')
        row=[int(x) for x in fields[:17]]+[float(fields[17])]
        if any(x<0 for x in row) or any(not math.isfinite(x) for x in row):raise ValueError('Invalid replacement cost')
        key=tuple(row[:2])
        if key in plans:raise ValueError('Duplicate replacement cost epoch')
        plans[key]=row
    if set(plans)!={tuple(r[2:4]) for r in stages}:raise ValueError('Replacement cost epoch coverage differs')
    selected=[plans[tuple(r[2:4])] for r in stages]
    return dict(input_sent_values=sum(p[6] for p in selected),input_received_values=sum(p[7] for p in selected),
                input_send_messages=sum(p[8] for p in selected),rk_sent_values=sum(p[9] for p in selected),
                rk_received_values=sum(p[10] for p in selected),rk_send_messages=sum(p[11] for p in selected),
                max_route_payload_bytes=max(p[12] for p in plans.values()),
                max_numerical_buffer_bytes=max(p[13] for p in plans.values()),
                max_rounding_column_bytes=max(p[14] for p in plans.values()),
                placement_sent_bytes=sum(p[15] for p in plans.values()),
                placement_received_bytes=sum(p[16] for p in plans.values()),
                backend_preparation_seconds=sum(p[17] for p in plans.values()),
                scope='Numerical buffers, rounding columns and resident/derived/RK route payload allocations only; excludes descriptors, compilation peak, geometry-only transform/metric adapters, retained accepted workspaces and MPI storage. Placement metadata is only its named topology stream. Stage traffic excludes 8-byte per-peer epoch headers. Full-step time and rank RSS include all costs.')
