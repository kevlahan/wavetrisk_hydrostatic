"""Reproduce bounded fused-producer evidence; never qualify shadow timings."""
from pathlib import Path
import math

def evidence(directory,ranks,enabled,domains=160):
    paths=sorted(Path(directory).glob('primitive276-rank-*.txt'))
    if not enabled:
        if paths:raise ValueError('Prototype output in control run')
        return dict(enabled=False)
    if len(paths)!=ranks:raise ValueError('Missing prototype rank evidence')
    total=dict(nodes=0,outputs=0,output_differences=0,reference_differences=0,skipped_domain_nodes=0,
               inputs=[0,0,0],input_differences=[0,0,0],zero_domain_rank_nodes=0)
    rank_ids=set();phases=set();byrank=[];balance={};phase_rows={}
    for path in paths:
        lines=path.read_text().splitlines();header=lines[0].split()
        if header[0]!='RANK' or list(map(int,header[2:4]))!=[ranks,domains]:raise ValueError('Rank header differs')
        rank,nd=int(header[1]),int(header[4])
        if rank in rank_ids or path.name!=f'primitive276-rank-{rank}.txt':raise ValueError('Duplicate rank')
        if not 0<=nd<=domains:raise ValueError('Invalid Domain count')
        rank_ids.add(rank);seen=set();timings=[0.0]*7;nodes=0
        for line in lines[1:]:
            row=line.split()
            if row[0]!='PHASE' or len(row)!=24:raise ValueError('Incomplete prototype phase')
            ints=list(map(int,row[1:17]));times=list(map(float,row[17:]))
            transaction,k,l,generation,sent,n,skipped,*counts=ints
            key=(transaction,k,l)
            if key in seen or transaction not in [1,2,3,9,17,33,41,48] or k not in [1,30]:raise ValueError('Duplicate/unexpected phase')
            seen.add(key)
            balance[key]=balance.get(key,0)+sent-n
            if min(ints)<0 or any(not math.isfinite(t) or t<0 for t in times):raise ValueError('Invalid count/time')
            inputs,diffs=counts[:3],counts[3:6];outputs,odiff,rdiff=counts[6:]
            if sum(inputs)!=16*n or outputs!=9*n:raise ValueError('Prototype coverage inconsistent')
            if any(d>n for d,n in zip(diffs,inputs)) or max(odiff,rdiff)>outputs:raise ValueError('Difference counts exceed coverage')
            total['nodes']+=n;nodes+=n;total['outputs']+=outputs;total['output_differences']+=odiff
            total['reference_differences']+=rdiff;total['skipped_domain_nodes']+=skipped
            for j in range(3):total['inputs'][j]+=inputs[j];total['input_differences'][j]+=diffs[j]
            timings=[a+b for a,b in zip(timings,times)]
            phase_rows.setdefault(key,[]).append(dict(rank=rank,domains=nd,nodes=n,seconds=times))
            if nd==0:total['zero_domain_rank_nodes']+=n
        if {(t,k) for t,k,l in seen}!={(t,k) for t in (1,2,3,9,17,33,41,48) for k in (1,30)}:
            raise ValueError('Missing sampled RHS/layer')
        if not phases:phases=seen
        if seen!=phases:raise ValueError('Ranks did not observe common phases')
        byrank.append(dict(rank=rank,domains=nd,nodes=nodes,seconds=timings))
    if rank_ids!=set(range(ranks)) or not total['nodes']:raise ValueError('Incomplete rank coverage')
    if any(balance.values()):raise ValueError('Reference send/receive node counts differ')
    if sum(r['domains'] for r in byrank)!=domains:raise ValueError('Domain ownership census differs')
    if ranks>domains and not total['zero_domain_rank_nodes']:raise ValueError('No prototype work on zero-Domain ranks')
    total['zero_domain_ranks']=sum(r['domains']==0 for r in byrank)
    total['phase_costs']=[dict(transaction=key[0],layer=key[1],level=key[2],nodes=sum(r['nodes'] for r in rows),
                             max_rank_seconds=[max(r['seconds'][i] for r in rows) for i in range(7)],
                             max_rank_gather_compute_seconds=max(sum(r['seconds'][i] for i in (4,5)) for r in rows),
                             zero_domain_nodes=sum(r['nodes'] for r in rows if r['domains']==0))
                          for key,rows in sorted(phase_rows.items())]
    total.update(enabled=True,phases_per_rank=len(phases),ranks=sorted(byrank,key=lambda r:r['rank']),
                 resident_gate_passed=not any(total['input_differences']) and not total['output_differences'],
                 reference_gate_passed=not total['reference_differences'],performance_qualified=False,
                 timing_fields=['map_prepare','reference_capture','reference_exchange','address_prepare','resident_gather','fused_compute_and_store','validation'],
                 scope='Positive block patch interiors only, selected RHS/layers; Domain physics and geometry reference packets. No restriction, alias publication, B/Exner or source consumer replacement.')
    if total['reference_differences']:raise ValueError('Fused formula differs from original step1')
    return total
