"""Validate the selected provisional transform, its per-epoch work and bounded costs."""
from pathlib import Path
import math

def evidence(directory,ranks,enabled,oracle=False,domains=160):
    root=Path(directory)
    files=list(root.glob('transform279-*-rank-*.txt'))
    if not enabled:
        if files:raise ValueError('Transform audit in an unaudited run')
        return {'enabled':False}
    records=[];sends=[0]*36;receives=[0]*36;domain_count=0;global_records=[[0,0,0] for _ in range(36)]
    if len(files)!=3*ranks:raise ValueError('Incomplete transform rank evidence')
    for rank in range(ranks):
        rows=[[int(x) for x in line.split()] for line in (root/f'transform279-production-rank-{rank}.txt').read_text().splitlines()]
        if len(rows)!=36 or any(len(row)!=19 for row in rows) or [r[0] for r in rows]!=list(range(1,37)):
            raise ValueError('Require all 36 provisional transform transactions')
        plans={}
        for line in (root/f'transform279-cost-rank-{rank}.txt').read_text().splitlines():
            f=line.split()
            if len(f)!=14:raise ValueError('Transform cost width')
            r=[int(x) for x in f[:13]]+[float(f[13])]
            if any(not math.isfinite(x) or x<0 for x in r):raise ValueError('Invalid transform cost')
            key=tuple(r[:2])
            if key in plans:raise ValueError('Transform epoch rebuilt without invalidation')
            plans[key]=r
        work={key:[0,0,0,0] for key in plans};blocks={key:set() for key in plans}
        for line in (root/f'transform279-work-rank-{rank}.txt').read_text().splitlines():
            f=[int(x) for x in line.split()]
            if len(f)!=7 or any(x<0 for x in f):raise ValueError('Transform operation vector width')
            key=tuple(f[:2]);b=f[2]
            if key not in work or b in blocks[key]:raise ValueError('Transform work epoch or block duplicate')
            blocks[key].add(b);work[key]=[a+b for a,b in zip(work[key],f[3:])]
        if set(plans)!={tuple(r[1:3]) for r in rows}:raise ValueError('Transform cost epoch coverage')
        for j,r in enumerate(rows):
            p=plans[tuple(r[1:3])]
            if any(x<0 for x in r) or r[8]!=int(oracle) or r[15]!=1 or p[11]!=1:
                raise ValueError('Transform reference-path or context lifetime contract')
            if r[13:15]!=p[7:9] or r[16]!=p[12] or r[4]>p[3]:raise ValueError('Transform layout costs disagree')
            if work[tuple(r[1:3])]!=[p[3],*r[5:8]]:raise ValueError('Per-block operation census disagrees')
            # Placement follows each operation's output/dependency component.
            # A numerical owner may legitimately receive just one kernel class.
            if r[4]>0 and sum(r[5:8])<=0:raise ValueError('Canonical owner has no transform records')
            global_records[j]=[a+b for a,b in zip(global_records[j],r[5:8])]
            sends[j]+=r[9];receives[j]+=r[10]
        domain_count+=rows[-1][3]
        records.append(dict(rank=rank,domains=rows[-1][3],transactions=36,plan_builds=len(plans),
            kernel_records=sum(sum(r[5:8]) for r in rows),sent_values=sum(r[9] for r in rows),
            received_values=sum(r[10] for r in rows),send_messages=sum(r[11] for r in rows),
            route_invocations=sum(r[12] for r in rows),max_route_payload_bytes=max(r[13] for r in rows),
            max_numerical_buffer_bytes=max(r[14] for r in rows),max_live_contexts=1,
            compatibility_publication_values=sum(r[17] for r in rows),block_boundary_values=sum(r[18] for r in rows),
            preparation_seconds=sum(p[13] for p in plans.values())))
    if any(min(row)<=0 for row in global_records):raise ValueError('Missing global transform kernel class')
    if sends!=receives:raise ValueError('Transform global traffic does not balance by transaction')
    if domain_count!=domains:raise ValueError('Transform Domain census differs')
    zero=[r for r in records if r['domains']==0]
    if ranks>domains and (len(zero)!=ranks-domains or any(r['kernel_records']<=0 for r in zero)):
        raise ValueError('No useful above-Domain transform participation')
    return dict(enabled=True,global_kernel_records_by_transaction=global_records,transactions_per_rank=36,reference_path=oracle,production_retirement=not oracle,
        zero_domain_ranks=len(zero),zero_domain_kernel_records=sum(r['kernel_records'] for r in zero),ranks=records,
        scope='Three provisional transitions per step. Fixed kernel records are work descriptors, not FLOPs or time weights. '
        'Production bypasses the old complete WT path; oracle runs compare old outputs and restore candidate SOL for real consumers. '
        'Traffic excludes metadata and epoch headers; route invocations are not global synchronization rounds. '
        'Memory is named numerical/route payload, excludes kernel/geometry metadata, descriptors, compilation peak, retained final-transform workspace and MPI. '
        'Rank RSS and whole-step time remain separate gates; final lifecycle and granularity are unchanged.')


def check_global_census(records):
    """Same fixed trajectory must retain every compiled record across placements."""
    active=[r['global_kernel_records_by_transaction'] for r in records if r.get('enabled')]
    if active and any(r!=active[0] for r in active[1:]):
        raise ValueError('Global transform census differs across validation runs')
