"""Compile sampled real direct-kernel plans into an ownership/communication budget.

This offline design compiler never evaluates numerical expressions or introduces
an execution-depth interpreter. Counts are per layer and per sampled plan, not
measured timings or sums over an RK trajectory. Unmapped nodes remain explicit
Domain adapters; they are never silently assigned to an invented block owner.
"""
from collections import Counter,defaultdict
from pathlib import Path
import json,re


def read_plan(path):
    lines=iter(Path(path).read_text().splitlines())
    def row(tag,width=None):
        try:parts=next(lines).split()
        except StopIteration:raise ValueError('Truncated '+tag) from None
        if not parts or parts[0]!=tag:raise ValueError('Expected '+tag)
        try:values=list(map(int,parts[1:]))
        except ValueError:raise ValueError('Noninteger '+tag) from None
        if width is not None and len(values)!=width:raise ValueError('Invalid width '+tag)
        return values
    rank,ranks,domains,layers=row('RANK',4)
    if not (0<=rank<ranks and domains>0 and layers>0):raise ValueError('Invalid rank header')
    snapshots=[];seen=Counter()
    while True:
        try:line=next(lines)
        except StopIteration:break
        parts=line.split()
        if not parts or parts[0]!='SNAP' or len(parts)!=7:raise ValueError('Invalid snapshot header')
        dom,sample,generation,n,lo,hi=map(int,parts[1:])
        seen[dom]+=1
        if not (0<=dom<domains and sample==seen[dom] and sample<=3 and n>0 and 0<=lo<=hi):
            raise ValueError('Invalid snapshot identity/sequence')
        nodes={}
        for i in range(1,n+1):
            node,owner,block,level=row('NODE',4)
            if node!=i or not (-1<=owner<ranks):raise ValueError('Invalid node owner/identity')
            if (owner==-1)!=(block==-1) or (owner==-1)!=(level==-1):raise ValueError('Partial node mapping')
            if owner>=0 and (block<0 or level<0):raise ValueError('Invalid mapped node')
            nodes[node]=(owner,block,level)
        def address(node):
            if node not in nodes:raise ValueError('Invalid stencil address')
        actions=[]
        for expected in range(hi,lo-1,-1):
            level,count=row('LEVEL',2)
            if level!=expected or count<0:raise ValueError('Invalid level schedule')
            for index in range(1,count+1):
                parts=next(lines).split();kind=parts[0];v=list(map(int,parts[1:]))
                if not v or v[0]!=index:raise ValueError('Invalid ordered action')
                if kind=='DIRECT':
                    if len(v)!=11 or v[1] not in (0,1):raise ValueError('Invalid direct action')
                    address(v[2])
                    if v[1]:
                        for node in v[2:]:address(node)
                    actions.append(dict(kind=kind,level=level,index=index,active=v[1],nodes=v[2:]))
                elif kind=='RESTRICT':
                    if len(v)!=9 or any(x not in (0,1) for x in v[6:]):raise ValueError('Invalid restriction')
                    address(v[1])
                    for e,active in enumerate(v[6:]):
                        if active:address(v[2]);address(v[3+e])
                    actions.append(dict(kind=kind,level=level,index=index,target=v[1],child=v[2],neighbor=v[3:6],edges=v[6:]))
                else:raise ValueError('Unknown action')
        count,=row('GRADIENTS',1);gradients=[]
        if count<0:raise ValueError('Invalid gradient extent')
        for index in range(1,count+1):
            v=row('GRADIENT',6)
            if v[0]!=index or v[1] not in (0,1):raise ValueError('Invalid gradient')
            address(v[2])
            if v[1]:
                for node in v[2:]:address(node)
            gradients.append(dict(active=v[1],nodes=v[2:]))
        row('END',0)
        snapshots.append(dict(domain=dom,sample=sample,generation=generation,nodes=nodes,actions=actions,gradients=gradients))
    return dict(rank=rank,ranks=ranks,domains=domains,layers=layers,snapshots=snapshots)


def compile_snapshot(snapshot,domain_owner,ranks):
    nodes=snapshot['nodes'];actions=snapshot['actions'];counts=Counter()
    target_work=Counter();adapters=Counter();inputs=defaultdict(set);writers={};restriction_routes=set()
    same_level=[];unresolved=[];source_reads=0
    def owner(node):
        destination=nodes[node][0]
        return domain_owner if destination<0 else destination
    def target(node,kind):
        destination=owner(node);target_work[(destination,kind)]+=1
        if nodes[node][0]<0:adapters[kind]+=1
        return destination
    for action in actions:
        level,index=action['level'],action['index']
        if action['kind']=='DIRECT':
            node=action['nodes'][0];destination=target(node,'direct')
            counts['direct']+=1;counts['active_direct']+=action['active']
            if action['active']:
                # Proposed full 3x9 argument packing: deliberately includes the
                # whole array arguments. This is a concrete conservative design
                # budget, not an assertion that the kernel loads every element.
                for dep in action['nodes']:
                    for field in ('flux','pv'):
                        for e in range(3):inputs[(level,destination)].add((field,dep,e))
                for e in range(3):inputs[(level,destination)].add(('physics',node,e))
            for e in range(3):writers[(node,e)]=(destination,level,index)
        else:
            node=action['target'];destination=target(node,'restriction')
            counts['restriction_actions']+=1
            # Read both operands before updating each target component, as in
            # execute_velocity_sources. Never collapse repeated writes by node.
            for e,active in enumerate(action['edges']):
                if not active:continue
                counts['restricted_edges']+=1
                for dep in (action['child'],action['neighbor'][e]):
                    source_reads+=1;producer=writers.get((dep,e))
                    if producer is None:
                        unresolved.append([level,index,dep,e]);continue
                    src,source_level,source_index=producer
                    if src!=destination:
                        restriction_routes.add((source_level,source_index,level,destination,src,dep,e))
                        if source_level==level:same_level.append([level,source_index,index,src,destination,dep,e])
                writers[(node,e)]=(destination,level,index)
    gradient_inputs=defaultdict(set);gradient_source_routes=set()
    for index,g in enumerate(snapshot['gradients'],1):
        node=g['nodes'][0];destination=target(node,'gradient')
        counts['gradient']+=1;counts['active_gradient']+=g['active']
        if not g['active']:continue
        for dep in g['nodes']:
            for field in ('bernoulli','exner','rho','rho_theta'):gradient_inputs[destination].add((field,dep,0))
        for e in range(3):
            producer=writers.get((node,e))
            if producer is None:unresolved.append(['gradient',index,node,e]);continue
            if producer[0]!=destination:gradient_source_routes.add((producer[0],destination,node,e))
    def volume(batches):
        logical=remote_current=remote_fused=adapter_values=0;peers=set()
        for key,values in batches.items():
            destination=key[-1] if isinstance(key,tuple) else key
            logical+=len(values)
            for field,node,e in values:
                if domain_owner!=destination:remote_current+=1
                assumed=owner(node)
                if assumed!=destination:remote_fused+=1;peers.add((assumed,destination))
                if nodes[node][0]<0:adapter_values+=1
        return dict(logical_packed_values=logical,retained_Domain_producer_remote_bytes=8*remote_current,
                    hypothetical_upstream_block_producer_remote_bytes=8*remote_fused,
                    unresolved_adapter_input_values=adapter_values,upstream_block_peer_pairs=len(peers))
    return dict(domain=snapshot['domain'],sample=snapshot['sample'],generation=snapshot['generation'],
        observing_rank=domain_owner,
        counts=dict(counts),target_work=[dict(rank=r,kind=k,count=n) for (r,k),n in sorted(target_work.items())],
        unmapped_target_adapter_work=dict(adapters),source_input_budget=volume(inputs),gradient_input_budget=volume(gradient_inputs),
        restriction_source_reads=source_reads,restriction_remote_versioned_values=len(restriction_routes),
        same_level_cross_owner_dependencies=len(same_level),same_level_examples=same_level[:12],
        gradient_remote_source_values=len(gradient_source_routes),unresolved_source_reads=len(unresolved),unresolved_examples=unresolved[:12],
        no_detected_same_level_or_uninitialized_source_blocker=not same_level and not unresolved,
        current_native_tendency_export_remote_bytes_per_layer=24*sum(v[0]>=0 and v[0]!=domain_owner for v in nodes.values()),
        production_ready=False)


def evidence(directory,ranks,enabled,domains=160):
    files=sorted(Path(directory).glob('ownership274-rank-*.txt'))
    if not enabled:
        if files:raise ValueError('Unexpected ownership diagnostic')
        return dict(enabled=False)
    if len(files)!=ranks:raise ValueError('Missing rank plan evidence')
    plans=[read_plan(p) for p in files]
    if any(p.name!=f'ownership274-rank-{plan["rank"]}.txt' for p,plan in zip(files,plans)):
        raise ValueError('Rank filename/header differs')
    if len({p['layers'] for p in plans})!=1:raise ValueError('Layer count differs across ranks')
    if {p['rank'] for p in plans}!=set(range(ranks)) or any(p['ranks']!=ranks or p['domains']!=domains for p in plans):
        raise ValueError('Rank/domain identity differs')
    seen=set();compiled=[];samples=Counter();observed_domains=set();domain_ranks=defaultdict(set)
    for plan in plans:
        for s in plan['snapshots']:
            key=(plan['rank'],s['domain'],s['sample'])
            if key in seen:raise ValueError('Duplicate rank-local Domain snapshot')
            seen.add(key);samples[s['sample']]+=1
            observed_domains.add(s['domain']);domain_ranks[s['domain']].add(plan['rank'])
            compiled.append(compile_snapshot(s,plan['rank'],ranks))
    if observed_domains!=set(range(domains)):raise ValueError('Incomplete Domain coverage')
    compiled.sort(key=lambda x:(x['sample'],x['domain'],x['observing_rank']))
    aggregate={}
    for sample in sorted(samples):
        rows=[c for c in compiled if c['sample']==sample];counts=Counter();work=Counter()
        for c in rows:
            counts.update(c['counts'])
            for x in c['target_work']:work[(x['rank'],x['kind'])]+=x['count']
        aggregate[str(sample)]=dict(domains=len({c['domain'] for c in rows}),observations=len(rows),counts=dict(counts),
            work_by_rank=[dict(rank=r,direct=work[(r,'direct')],restriction=work[(r,'restriction')],gradient=work[(r,'gradient')]) for r in range(ranks)],
            same_level_cross_owner_dependencies=sum(c['same_level_cross_owner_dependencies'] for c in rows),
            unresolved_source_reads=sum(c['unresolved_source_reads'] for c in rows),
            current_native_tendency_export_remote_bytes_per_layer=sum(c['current_native_tendency_export_remote_bytes_per_layer'] for c in rows),
            source_input_budget={k:sum(c['source_input_budget'][k] for c in rows) for k in rows[0]['source_input_budget']},
            gradient_input_budget={k:sum(c['gradient_input_budget'][k] for c in rows) for k in rows[0]['gradient_input_budget']},
            unmapped_target_adapter_work=dict(sum((Counter(c['unmapped_target_adapter_work']) for c in rows),Counter())))
    return dict(enabled=True,rank_files=ranks,domains=domains,sampled_plans=len(compiled),by_sample=aggregate,plans=compiled,
        schema_version=2,domains_observed_on_multiple_ranks={str(d):sorted(rs) for d,rs in sorted(domain_ranks.items()) if len(rs)>1},
        aggregation='Sums of rank-local observations, including repeated Domain observations after ownership changes; not a global mesh census.',
        performance_qualified=False,production_ready=False,
        scope='First three actual plan builds per Domain per observing rank; per-layer structural counts, not RK call-weighted totals. Rank-local sample ordinals and generations do not establish simultaneous global epochs. Repeated Domain observations on different ranks are retained, not deduplicated. Full kernel argument packing is conservative. Boundary aliases remain Domain-local identities; hypothetical upstream ownership is conditional, not an implemented producer. No timing or speedup inference.')
