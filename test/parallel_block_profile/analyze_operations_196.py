"""Conserved dynamics accounting by WT caller context, with explicit residuals."""
from pathlib import Path
import math,statistics as st
from analyze_workflow_195 import summarize as workflow_summary

GROUPS={
 'Domain boundary CPU preparation/copies':[28,58,59,61,62,63],
 'Domain boundary MPI posting/wait':[38,60],
 'Inverse transfer planning/packing/install':[15,16,17,20],
 'Inverse transfer MPI post/self-copy/wait':[18,19],
 'Block ghost preparation/install':[33], 'Block ghost MPI':[34],
 'Scalar transport/install':[24], 'Scalar transport MPI':[37],
 'Native mass exchange preparation/post/install':[68,69], 'Native mass MPI wait':[39],
 'Domain writeback preparation/install':[32], 'Domain writeback MPI':[40],
 'Domain full import preflight/install':[70],
 'Domain full import pack/zero/local traversal':[71],
 'Domain full import MPI':[72],
 'Domain boundary import preflight':[73],
 'Domain boundary import zero/pack':[74],
 'Domain boundary import MPI':[75],
 'Domain boundary import install/local copies':[76],
 'Domain publication commit/residual':[77],
 'Scalar input storage/geometry/capture':[22,23,41,51],
 'Inverse scaffold staging':[42],
 'Inverse numerical kernels and traversal':[21,36,66,67],
 'Scalar restriction/divergence and traversal':[25],
 'Wavelet/compression self work':[43,44],
 'Shared producers/basic operators/native mass/velocity':[26,27,29,30,31,52],
 'RK assembly/retention/tendency/compatibility':[45,46,47,53],
 'State publication/refresh':[48,49,50],
 'Other named or unclassified self work':[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,35,54,55,56,57,64,65],
}
BOUNDARY_GROUPS=list(GROUPS)[:20]
assert sorted(i for ids in GROUPS.values() for i in ids)==list(range(78))


def collect(root,ranks,mode,side=None):
    root=Path(root);files=list(root.glob('operations196-rank-*.txt'));data={}
    if len(files)!=ranks:raise ValueError('Missing operation rank files')
    for file in files:
        rank=int(file.stem.rsplit('-',1)[1]);rows={};seen=set();schema=None
        if rank in data or not 0<=rank<ranks:raise ValueError('Bad operation rank')
        for line in file.read_text().splitlines():
            f=line.split()
            if f[0]=='schema':
                if schema is not None:raise ValueError('Duplicate schema')
                schema=list(map(int,f[1:]));continue
            if f[0]=='step':
                seq,restart,remap=map(int,f[1:4])
                if seq in rows or restart not in (0,1) or remap not in (0,1):raise ValueError('Invalid step')
                rows[seq]={'restart':restart,'remap':remap,'wall':float(f[4]),'self':{},'calls':{},'count':{},'volume':{}}
                continue
            if f[0] not in ('self','count','volume'):raise ValueError('Unknown record')
            seq,c,p=map(int,f[1:4]);key=f[0],seq,c,p
            if key in seen or seq not in rows or not 0<=c<=2:raise ValueError('Duplicate/invalid record')
            seen.add(key)
            if f[0]=='self':
                if not 0<=p<=77:raise ValueError('Bad operation id')
                calls=int(f[4]);seconds=float(f[5])
                if calls<0 or not math.isfinite(seconds) or seconds<0:raise ValueError('Invalid operation values')
                rows[seq]['self'][c,p]=seconds;rows[seq]['calls'][c,p]=calls
            elif f[0]=='count':
                if not 1<=p<=27 or int(f[4])<0:raise ValueError('Invalid counter')
                rows[seq]['count'][c,p]=int(f[4])
            else:
                if not 1<=p<=56 or len(f)!=6 or min(map(int,f[4:]))<0:raise ValueError('Invalid volume')
                rows[seq]['volume'][c,p]=list(map(int,f[4:]))
        if schema!=[1,mode,rank] or not rows or set(rows)!=set(range(1,len(rows)+1)):raise ValueError('Incomplete schema/steps')
        # Read independent common-clock envelopes for this exact rank.
        envelope={};transactions={};flags={};transaction_calls={}
        for line in (root/f'workflow195-rank-{rank}.txt').read_text().splitlines():
            f=line.split()
            if f[0]=='step':flags[int(f[1])]=(int(f[2]),int(f[3]))
            if f[0]=='phase' and f[2]=='1':envelope[int(f[1])]=float(f[3])
            if f[0]=='transaction' and f[2]=='1':
                transactions[int(f[1]),int(f[3])]=float(f[5]);transaction_calls[int(f[1]),int(f[3])]=int(f[4])
        if set(envelope)!=set(rows):raise ValueError('Common workflow step coverage differs')
        for seq,row in rows.items():
            if not math.isfinite(row['wall']) or row['wall']<0:raise ValueError('Invalid dynamics wall')
            if (row['restart'],row['remap'])!=flags[seq] or abs(row['wall']-envelope[seq])>1e-6:raise ValueError('Envelope/flags disagree')
            if len(row['self'])!=234 or len(row['count'])!=81:raise ValueError('Missing dense records')
            expected={1:transactions[seq,1],2:transactions[seq,2]}
            expected[0]=row['wall']-expected[1]-expected[2]
            row['context_envelopes']={'WT':expected[1],'outside_WT':expected[0]+expected[2]}
            if min(expected.values()) < -1e-6:raise ValueError('Overlapping workflow transactions')
            if mode==2 and row['calls'][1,65]!=transaction_calls[seq,1]:raise ValueError('WT operation hooks do not cover common WT calls')
            if mode==0:
                if any(row['self'].values()) or any(row['calls'].values()) or any(row['count'].values()) or row['volume']:raise ValueError('Disabled observer recorded operations')
            else:
                for c in range(3):
                    if abs(sum(row['self'][c,p] for p in range(78))-expected[c])>1e-6:raise ValueError('Context self times fail conservation')
                    if mode==1:
                        allowed={0,*[i for g in BOUNDARY_GROUPS for i in GROUPS[g]]}
                        if any(row['self'][c,p] or row['calls'][c,p] for p in range(78) if p not in allowed):raise ValueError('Unexpected boundary-mode scope')
        data[rank]=rows
    if any(set(rows)!=set(data[0]) for rows in data.values()):raise ValueError('Different rank steps')
    steps=[]
    for seq in data[0]:
        values=[data[rank][seq] for rank in range(ranks)]
        if len({(v['restart'],v['remap']) for v in values})!=1:raise ValueError('Rank workload differs')
        slow=max(range(ranks),key=lambda i:values[i]['wall'])
        row={'sequence':seq,'restart':values[0]['restart'],'remap':values[0]['remap'],
             'dynamics_wall_mean':st.mean(v['wall'] for v in values),'dynamics_wall_max':values[slow]['wall'],'slow_rank':slow,'context_wall_means':{c:st.mean(v['context_envelopes'][c] for v in values) for c in ('WT','outside_WT')},'contexts':{}}
        for label,contexts in [('WT',[1]),('outside_WT',[0,2]),('tendency_tag_diagnostic_only',[2])]:
            ids=[st.mean(sum(v['self'][c,p] for c in contexts) for v in values) for p in range(78)]
            row['contexts'][label]={'self_seconds':ids,'groups':{g:sum(ids[p] for p in ps) for g,ps in GROUPS.items()},
                'calls_global':[sum(v['calls'][c,p] for c in contexts for v in values) for p in range(78)],
                'counters_global':[sum(v['count'][c,p] for c in contexts for v in values) for p in range(1,28)],
                'block_volume_by_profile_id':{str(p):[sum(v['volume'].get((c,p),[0,0])[i] for c in contexts for v in values) for i in range(2)] for p in range(1,57) if any((c,p) in v['volume'] for c in contexts for v in values)}}
        steps.append(row)
    if side not in (None,'block','legacy'):raise ValueError('Invalid implementation side')
    totals=[sum(row['calls'][c,p] for rows in data.values() for row in rows.values() if not row['restart'] for c in range(3)) for p in range(78)]
    required=[]
    if mode>0:required=[28,58,38]
    if side=='block' and mode>0:required += [16,19,33,37,68,69,70,71,72,73,74,75,76,77]
    if mode==2:required += [65]
    if mode==2 and side=='block':required += [25,43,44,66,67]
    if mode==2 and side=='legacy':required += [14,27]
    if any(totals[p]==0 for p in required):raise ValueError('Required operation scope has no ordinary coverage: '+str([p for p in required if not totals[p]]))
    ordinary=[s for s in steps if not s['restart']]
    if not ordinary:raise ValueError('No ordinary steps')
    return {'schema':1,'mode':mode,'ranks':ranks,'side':side,'required_scopes':required,'steps':steps,
            'ordinary_context_wall_means':{c:st.mean(s['context_wall_means'][c] for s in ordinary) for c in ('WT','outside_WT')},
            'ordinary_context_group_means':{c:{g:st.mean(s['contexts'][c]['groups'][g] for s in ordinary) for g in GROUPS} for c in ('WT','outside_WT')},
            'scope':'Only dynamics. WT and outside_WT are disjoint; tendency tag is nested diagnostic only. Self time conserves each rank/context. MPI may include progress/arrival waits. Group names do not establish equivalent numerical work or avoidable cost.'}


def summarize(runs):
    out=workflow_summary(runs);out['operation_excess']={}
    for mode in ('coarse','detail'):
        pairs=[]
        for b in [v for v in runs.values() if v['phase']=='timing-'+mode and v['side']=='block']:
            l=next(v for v in runs.values() if v['phase']==b['phase'] and v['round']==b['round'] and v['side']=='legacy')
            pair={'round':b['round'],'contexts':{}}
            for context in ('WT','outside_WT'):
                values={side:run['operations']['ordinary_context_group_means'][context] for side,run in [('block',b),('legacy',l)]}
                diff={g:values['block'][g]-values['legacy'][g] for g in GROUPS}
                pair['contexts'][context]={'sides':values,'group_excess':diff,'total_excess':sum(diff.values()),
                    'boundary_and_transfer_excess':sum(diff[g] for g in BOUNDARY_GROUPS)}
            dynamics_excess=b['workflow']['ordinary_phase_mean']['advance']-l['workflow']['ordinary_phase_mean']['advance']
            residual=dynamics_excess-sum(v['total_excess'] for v in pair['contexts'].values())
            if abs(residual)>2e-6:raise ValueError('Operation differences fail to reconcile common dynamics excess')
            pair['reconciliation_residual']=residual
            pairs.append(pair)
        out['operation_excess'][mode]=pairs
    out['context_observer_effect_percent']={}
    for side in ('block','legacy'):
        selected={m:[v for v in runs.values() if v['phase']=='timing-'+m and v['side']==side] for m in ('minimal','coarse','detail')}
        out['context_observer_effect_percent'][side]={m:{c:100*(st.median(v['operations']['ordinary_context_wall_means'][c] for v in selected[m])/st.median(v['operations']['ordinary_context_wall_means'][c] for v in selected['minimal'])-1) for c in ('WT','outside_WT')} for m in ('coarse','detail')}
    out['scope']='Minimal=common outer observer only; coarse=boundary/transfer operations; detail=all operation scopes. Old detailed profilers disabled. Compare WT and outside-WT separately. Residual is not automatically computation; named MPI time is not network service time.'
    return out
