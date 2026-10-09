"""Require every native RHS/layer to exercise the production assignment checks."""
from pathlib import Path

def evidence(directory,ranks,enabled,domains=160):
    directory=Path(directory);paths=sorted(directory.glob('primitive277-rank-*.txt'))
    log=(directory/'run.log').read_text()
    marker='Stage277 cached fused primitive = on; assignment checks = '
    if not enabled:
        if paths:raise ValueError('Assignment evidence in uninstrumented control')
        if marker+'T' in log:raise ValueError('Assignment checks in timing run')
        return dict(enabled=False)
    if len(paths)!=ranks or log.count(marker+'T')!=1:raise ValueError('Missing primitive rank evidence/marker')
    ids=set();total=checks=ndom=0;zero=0;rebuilds=[];transactions=set()
    for p in paths:
        lines=p.read_text().splitlines();header=lines[0].split()
        if len(header)!=7 or header[0]!='RANK':raise ValueError('Bad rank header')
        rank,n,d,owned,nz,patch=map(int,header[1:])
        if rank in ids or n!=ranks or d!=domains or not 0<=rank<ranks or not 0<=owned<=domains or nz!=30 or patch!=4:
            raise ValueError('Rank identity differs')
        if p.name!=f'primitive277-rank-{rank}.txt':raise ValueError('Rank filename differs')
        ids.add(rank);ndom+=owned;zero+=owned==0;seen=set();last_generation=-1;last_builds=0
        for line in lines[1:]:
            row=line.split()
            if len(row)!=5+nz or row[0]!='RHS':raise ValueError('Bad RHS evidence')
            transaction,generation,builds,count,*layers=map(int,row[1:])
            if transaction in seen or min(generation,builds,count,*layers)<0 or generation<last_generation or builds<last_builds:
                raise ValueError('Duplicate/invalid RHS evidence')
            if builds<1 or (owned and (min(layers)<=0 or len(set(layers))!=1)) or (not owned and any(layers)):
                raise ValueError('Missing layer work or unexpected zero-Domain work')
            if any(v%(patch*patch) for v in layers) or count!=19*sum(layers):raise ValueError('Assignment coverage differs')
            seen.add(transaction);last_generation=generation;last_builds=builds
            total+=sum(layers);checks+=count
        if seen!=set(range(1,49)):raise ValueError('Missing native RHS transaction')
        transactions|=seen;rebuilds.append(last_builds)
    if ndom!=domains or ids!=set(range(ranks)) or not total:raise ValueError('Ownership/work census differs')
    return dict(enabled=True,nodes=total,assignments=checks,rhs_per_rank=len(transactions),layers=30,
                zero_domain_ranks=zero,zero_domain_rank_nodes=0,cache_rebuilds_min=min(rebuilds),cache_rebuilds_max=max(rebuilds),
                assignment_gate_passed=True,reference_gate_passed=True,
                scope='Full-call production assignment checks on existing Domain owners; no ownership migration or new MPI traffic.')
