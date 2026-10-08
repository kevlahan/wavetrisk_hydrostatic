"""Catalogue validation and exact-integer reference cuts; no performance inference."""
from fractions import Fraction
import bisect
import statistics


def balanced(weights,ranks):
    if ranks<1 or any(w<0 for w in weights):raise ValueError('Invalid partition inputs')
    if not weights:return []
    effective=weights if sum(weights) else [1]*len(weights)
    prefix=[0]
    for w in effective:prefix.append(prefix[-1]+w)
    active=min(ranks,len(weights));previous=0;owners=[]
    n=len(weights)
    def suffix(cap):
        # Dynamic programming reference; independent of the Fortran two-pointer scan.
        need=[n+1]*(n+1);need[n]=0
        for i in range(n-1,-1,-1):
            for j in range(i+1,n+1):
                if prefix[j]-prefix[i]>cap:break
                need[i]=min(need[i],1+need[j])
        return need
    low=max(effective);high=prefix[-1]
    while low<high:
        mid=(low+high)//2
        if suffix(mid)[0]<=active:high=mid
        else:low=mid+1
    need=suffix(low)
    for r in range(1,active):
        target=Fraction(prefix[-1]*r,active)
        feasible=[k for k in range(previous+1,n-(active-r)+1)
                  if prefix[k]-prefix[previous]<=low and need[k]<=active-r]
        cut=min(feasible,key=lambda k:(abs(prefix[k]-target),k))
        owners.extend([r-1]*(cut-previous));previous=cut
    owners.extend([active-1]*(len(weights)-previous))
    return owners


def legacy(weights,ranks):
    if ranks<1 or len(weights)<ranks or any(w<=0 for w in weights):
        return None
    goal=.05;average=sum(weights)/ranks
    for attempt in range(1000):
        loads=[0]*ranks;owners=[-1]*len(weights);assigned=0
        for r in range(ranks):
            while loads[r]<average and len(weights)-assigned>ranks-r-1:
                owners[assigned]=r;loads[r]+=weights[assigned];assigned+=1
            if assigned>0 and loads[r]>average*(1+goal):
                loads[r]-=weights[assigned-1];owners[assigned-1]=-1;assigned-=1
        if assigned==len(weights):return owners
        goal*=1.2
    raise ValueError('Legacy assignment did not converge')


def metrics(weights,owners,ranks,domain_owners):
    loads=[0]*ranks;counts=[0]*ranks
    for w,r in zip(weights,owners):loads[r]+=w;counts[r]+=1
    average=sum(loads)/ranks
    return dict(weight_by_rank=loads,blocks_by_rank=counts,weight_min=min(loads),
        weight_max=max(loads),weight_mean=average,
        weight_max_over_mean=None if average==0 else max(loads)/average,
        weight_variance=statistics.pvariance(loads),
        Domain_owner_weight_share=None if sum(loads)==0 else sum(loads[r] for r in domain_owners)/sum(loads),
        Domain_free_weight_share=None if sum(loads)==0 else sum(loads[r] for r in range(ranks) if r not in domain_owners)/sum(loads),
        zero_block_ranks=counts.count(0))


def catalogue_evidence(log,ranks,enabled,domains=160):
    rows=[];current=None
    for line in log.splitlines():
        if line.startswith('Stage273 partition catalogue:'):
            if current is not None:rows.append(current)
            p,n=map(int,line.split(':')[1].split())
            if p!=ranks or n<domains:raise ValueError('Catalogue extent differs')
            current=dict(ranks=p,blocks=n,records=[])
        elif line.startswith('Stage273 partition block:'):
            if current is None:raise ValueError('Catalogue row without header')
            values=list(map(int,line.split(':')[1].split()))
            if len(values)!=6:raise ValueError('Catalogue record width differs')
            current['records'].append(values)
    if current is not None:rows.append(current)
    if bool(rows)!=enabled:raise ValueError('Validation catalogue selector/output differs')
    for row in rows:
        records=row['records'];n=row['blocks']
        if len(records)!=n or [v[0] for v in records]!=list(range(n)):
            raise ValueError('Incomplete/duplicate catalogue records')
        if any(not (0<=v[1]<domains and v[2]>=0 and v[3]>0 and 0<=v[4]<ranks and 0<=v[5]<ranks) for v in records):
            raise ValueError('Invalid catalogue identity/weight/owner')
        roots={v[1] for v in records}
        if roots!=set(range(domains)):raise ValueError('Missing geometric Domain')
        if any(sum(v[1]==d for v in records) not in (1,4) for d in roots):raise ValueError('Invalid subtree split')
        if any(len({v[5] for v in records if v[1]==d})!=1 for d in roots):raise ValueError('Geometric owner differs within root')
        weights=[v[3] for v in records];owners=[v[4] for v in records]
        if owners!=balanced(weights,ranks):raise ValueError('Actual candidate owners differ from exact reference')
        old=legacy(weights,ranks);geom={v[5] for v in records}
        if old is not None:
            assert max(metrics(weights,owners,ranks,geom)['weight_by_rank'])<=max(metrics(weights,old,ranks,geom)['weight_by_rank'])
        row.update(candidate=metrics(weights,owners,ranks,geom),
            baseline_counterfactual=None if old is None else metrics(weights,old,ranks,geom),
            reassigned_blocks=None if old is None else sum(a!=b for a,b in zip(old,owners)))
    return dict(catalogues=rows,scope='Actual candidate catalogue records. Baseline ownership is an offline counterfactual on those exact weights; not measured baseline timing. Leaf weights are a proxy, not executed work or communication cost.')
