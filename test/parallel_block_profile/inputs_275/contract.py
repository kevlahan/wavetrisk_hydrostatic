"""Exact component budget and live resident-input epoch evidence; no timing claims."""
from collections import Counter,defaultdict
from pathlib import Path
import ownership
from ownership import read_plan,compile_snapshot

READS=((1,1),(2,1),(3,1),(1,2),(2,2),(3,2),(1,3),(2,3),(3,3),
       (1,4),(2,4),(3,4),(1,5),(2,5),(2,6),(2,7),(3,7),(3,8),(1,9))
TRANSACTIONS=(1,2,3,9,17,33,41,48)
FIELDS=('mass_state','mass_mean','rho','temperature_state','temperature_mean','rho_theta')


def footprint(snapshot,rank):
    nodes=snapshot['nodes'];batches=defaultdict(set)
    def owner(n):return nodes[n][0] if nodes[n][0]>=0 else rank
    for a in snapshot['actions']:
        if a['kind']!='DIRECT' or not a['active']:continue
        target=a['nodes'][0];key=('source',a['level'],owner(target))
        for field in ('flux','pv'):
            for edge,position in READS:batches[key].add((field,a['nodes'][position-1],edge))
        for edge in (1,2,3):batches[key].add(('physics',target,edge))
    for a in snapshot['gradients']:
        if not a['active']:continue
        key=('gradient',0,owner(a['nodes'][0]))
        for n in a['nodes']:
            for field in ('bernoulli','exner','rho','rho_theta'):batches[key].add((field,n,0))
    counts=defaultdict(Counter)
    for key,values in batches.items():
        for field,n,edge in values:
            c=counts[field];c['logical_values']+=1
            c['Domain_remote_bytes']+=8*(rank!=key[-1])
            c['hypothetical_upstream_remote_bytes']+=8*(owner(n)!=key[-1])
            c['unmapped_values']+=nodes[n][0]<0
    return {k:dict(v) for k,v in sorted(counts.items())}


def residency(directory,ranks,expected_domains=None):
    files=list(Path(directory).glob('inputs275-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing input probe rank files')
    resident={};wanted={};checks={};layers=None;domains=None;seen=set()
    for path in sorted(files):
        with path.open() as stream:
            h=stream.readline().split()
            if len(h)!=5 or h[0]!='RANK':raise ValueError('Invalid input header')
            rank,n,d,z=map(int,h[1:])
            if rank in seen or not 0<=rank<ranks or n!=ranks or path.name!=f'inputs275-rank-{rank}.txt':
                raise ValueError('Invalid input rank identity')
            if layers is not None and (z!=layers or d!=domains):raise ValueError('Input header extent differs')
            if z<1 or d<1:raise ValueError('Invalid input extent')
            if expected_domains is not None and d!=expected_domains:raise ValueError('Input Domain extent differs')
            seen.add(rank);layers=z;domains=d;rows=[]
            for line in stream:
                parts=line.split();tag=parts[0];v=tuple(map(int,parts[1:]))
                if tag=='CHECK':
                    if len(v)!=4 or min(v)<0:raise ValueError('Invalid check record')
                    rows.append(v);continue
                if tag not in ('D','B') or len(v)!=13:raise ValueError('Invalid input record')
                if v[0] not in TRANSACTIONS or v[1] not in (1,z) or min(v[2:4])<0:
                    raise ValueError('Invalid input epoch/key')
                key=v[:7];table=resident if tag=='B' else wanted
                if key in table:raise ValueError('Duplicate input identity')
                table[key]=(rank,v[7:])
        if not rows:raise ValueError('Missing shadow execution evidence')
        expected=[(tx,k) for tx in range(1,rows[-1][0]+1) for k in range(1,z+1)]
        if [v[:2] for v in rows]!=expected:raise ValueError('Incomplete gradient phase sequence')
        if any(a[2]>b[2] or a[3]>b[3] for a,b in zip(rows,rows[1:])):raise ValueError('Decreasing shadow checks')
        checks[rank]=rows
    sequences={tuple(v[:2] for v in rows) for rows in checks.values()}
    if len(sequences)!=1:raise ValueError('Rank transaction sequences differ')
    source_checks=sum(v[-1][2] for v in checks.values());gradient_checks=sum(v[-1][3] for v in checks.values())
    if min(source_checks,gradient_checks)<=0 or not wanted or not resident:raise ValueError('Empty input diagnostic')
    if gradient_checks!=domains*layers*checks[0][-1][0]:raise ValueError('Incomplete global gradient shadow coverage')
    epochs=sorted({(tx,k) for tx,k in next(iter(sequences)) if tx in TRANSACTIONS and k in (1,layers)})
    if sorted({k[:2] for k in wanted})!=epochs or sorted({k[:2] for k in resident})!=epochs:
        raise ValueError('Missing selected input epoch')
    if any(k[:2] not in epochs for k in wanted|resident):raise ValueError('Input epoch not executed')
    reports=[];examples=[]
    for epoch in epochs:
        keys=[k for k in wanted if k[:2]==epoch];different=Counter({f:0 for f in FIELDS})
        missing=matched=remote=0
        for key in keys:
            dr,dv=wanted[key]
            if key not in resident:missing+=1;continue
            br,bv=resident[key];matched+=1;remote+=br!=dr
            for field,a,b in zip(FIELDS,dv,bv):different[field]+=a!=b
            if dv!=bv and len(examples)<8:examples.append(dict(key=key,Domain_rank=dr,block_rank=br,Domain_bits=dv,block_bits=bv))
        reports.append(dict(transaction=epoch[0],layer=epoch[1],required_mapped_nodes=len(keys),matched_nodes=matched,
            missing_nodes=missing,cross_rank_matches=remote,different_values=dict(different)))
    return dict(source_shadow_calls=source_checks,gradient_shadow_calls=gradient_checks,transactions=checks[0][-1][0],
        layers=layers,epochs=reports,mismatch_examples=examples,
        sampled_owner_interior_reuse_matches=all(v['missing_nodes']==0 and not any(v['different_values'].values()) for v in reports),
        scope='Exact bit comparisons of mapped gradient inputs against resident owner INTERIORS at selected global call ordinals, first/last layer. No halo freshness, intervening calls/layers, PV/flux/physics or Bernoulli/Exner residency proof. Missing or stale values are diagnostic findings, not silently accepted as reusable.')


def evidence(directory,ranks,enabled,domains=160):
    result=ownership.evidence(directory,ranks,enabled,domains)
    if not enabled:
        if list(Path(directory).glob('inputs275-rank-*.txt')):raise ValueError('Unexpected input probe')
        return result
    sampled=defaultdict(lambda:defaultdict(Counter))
    for path in sorted(Path(directory).glob('ownership275-rank-*.txt')):
        plan=read_plan(path)
        for snapshot in plan['snapshots']:
            for field,c in footprint(snapshot,plan['rank']).items():sampled[str(snapshot['sample'])][field].update(c)
    result['minimal_input_budget']={sample:{field:dict(c) for field,c in sorted(fields.items())} for sample,fields in sorted(sampled.items())}
    result['input_epochs']=residency(directory,ranks,domains)
    if result['input_epochs']['transactions']!=48:raise ValueError('Require all 48 RHS transactions')
    result['schema_version']=3
    result['scope']+=' Stage275 minimal component budget still assumes Domain-fed inputs; owner-interior epoch comparisons are separate samples and cannot be subtracted as proven halo traffic savings.'
    return result
