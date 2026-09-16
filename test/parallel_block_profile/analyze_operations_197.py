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
GROUPS.update({'Resident initialization':[78],'Direct scalar flux/traversal':[79],
 'Scalar/flux restriction traversal':[80],'Divergence/tendency traversal':[81],'Temperature boundary closure':[82]})
BOUNDARY_GROUPS=list(GROUPS)[:20]
assert sorted(i for ids in GROUPS.values() for i in ids)==list(range(83))


def collect(root,ranks,mode,side=None):
    root=Path(root);files=list(root.glob('operations197-rank-*.txt'));data={}
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
                if not 0<=p<=82:raise ValueError('Bad operation id')
                calls=int(f[4]);seconds=float(f[5])
                if calls<0 or not math.isfinite(seconds) or seconds<0:raise ValueError('Invalid operation values')
                rows[seq]['self'][c,p]=seconds;rows[seq]['calls'][c,p]=calls
            elif f[0]=='count':
                if not 1<=p<=35 or int(f[4])<0:raise ValueError('Invalid counter')
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
            if len(row['self'])!=249 or len(row['count'])!=105:raise ValueError('Missing dense records')
            expected={1:transactions[seq,1],2:transactions[seq,2]}
            expected[0]=row['wall']-expected[1]-expected[2]
            row['context_envelopes']={'WT':expected[1],'outside_WT':expected[0]+expected[2]}
            if min(expected.values()) < -1e-6:raise ValueError('Overlapping workflow transactions')
            if mode==2 and row['calls'][1,65]!=transaction_calls[seq,1]:raise ValueError('WT operation hooks do not cover common WT calls')
            if mode==0:
                if any(row['self'].values()) or any(row['calls'].values()) or any(row['count'].values()) or row['volume']:raise ValueError('Disabled observer recorded operations')
            else:
                for c in range(3):
                    if abs(sum(row['self'][c,p] for p in range(83))-expected[c])>1e-6:raise ValueError('Context self times fail conservation')
                    if mode==1:
                        allowed={0,*[i for g in BOUNDARY_GROUPS for i in GROUPS[g]]}
                        if any(row['self'][c,p] or row['calls'][c,p] for p in range(83) if p not in allowed):raise ValueError('Unexpected boundary-mode scope')
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
            ids=[st.mean(sum(v['self'][c,p] for c in contexts) for v in values) for p in range(83)]
            row['contexts'][label]={'self_seconds':ids,'groups':{g:sum(ids[p] for p in ps) for g,ps in GROUPS.items()},
                'calls_global':[sum(v['calls'][c,p] for c in contexts for v in values) for p in range(83)],
                'counters_global':[sum(v['count'][c,p] for c in contexts for v in values) for p in range(1,36)],
                'block_volume_by_profile_id':{str(p):[sum(v['volume'].get((c,p),[0,0])[i] for c in contexts for v in values) for i in range(2)] for p in range(1,57) if any((c,p) in v['volume'] for c in contexts for v in values)}}
        steps.append(row)
    if side not in (None,'block','legacy'):raise ValueError('Invalid implementation side')
    totals=[sum(row['calls'][c,p] for rows in data.values() for row in rows.values() if not row['restart'] for c in range(3)) for p in range(83)]
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
    result={'timings':{},'profiles':{},'scope':'Development screen: two rounds are not a confidence interval. Profiles excluded from timing medians. MPI wait is not network service time.'}
    for variant in ('baseline','candidate','legacy'):
        selected=[v for label,v in runs.items() if label.startswith(variant+'-') and v['phase']=='timing-bare']
        result['timings'][variant]={'ordinary_steps':[r['ordinary_step_median'] for r in selected],
          'elapsed_seconds':[r['elapsed_seconds'] for r in selected],
          'median_step':st.median(r['ordinary_step_median'] for r in selected),
          'median_elapsed':st.median(r['elapsed_seconds'] for r in selected)}
        probe=runs[variant+'-profile']
        result['profiles'][variant]={'context_groups':probe['operations']['ordinary_context_group_means'],
           'context_seconds':probe['operations']['ordinary_context_wall_means'],
           'printed_step_change_percent':100*(probe['ordinary_step_median']/result['timings'][variant]['median_step']-1),
           'dependencies':probe['operations'].get('dependencies',{}),
           'work_per_rank_step':probe['operations'].get('ordinary_counters_per_rank_step',{})}
    result['candidate_step_reduction_percent']=100*(1-result['timings']['candidate']['median_step']/result['timings']['baseline']['median_step'])
    result['candidate_to_legacy_ratio']=result['timings']['candidate']['median_step']/result['timings']['legacy']['median_step']
    return result


_base_collect=collect

def collect(root,ranks,mode,side=None):
    result=_base_collect(root,ranks,mode,side)
    files=list(Path(root).glob('dependencies197-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing dependency rank files')
    records=[];seen=set()
    for file in files:
        rank=int(file.stem.rsplit('-',1)[1]);lines=file.read_text().splitlines()
        if rank in seen or rank not in range(ranks) or lines[0].split()!=['schema','1',str(rank)]:raise ValueError('Invalid dependency rank/schema')
        seen.add(rank)
        for seq,line in enumerate(lines[1:],1):
            f=line.split()
            if len(f)!=14:raise ValueError('Invalid dependency row')
            v=list(map(int,f[:12]));times=list(map(float,f[12:]))
            if v[11]!=seq or v[0] not in range(1,len(result['steps'])+1) or v[1] not in range(3) or v[2] not in (1,2,3):raise ValueError('Invalid dependency key')
            if min(v[7:10])<0 or any(not math.isfinite(t) or t<0 for t in times):raise ValueError('Invalid dependency values')
            records.append({'rank':rank,'step':v[0],'context':v[1],'kind':v[2],'component':v[3],'family':v[4],
                'first':v[5],'last':v[6],'requests':v[7],'sent_values':v[8],'received_values':v[9],'tag':v[10],
                'post_to_wait_seconds':times[0],'completion_seconds':times[1]})
    if mode==0 and records:raise ValueError('Disabled dependency observer recorded work')
    if mode>0 and side=='block' and not records:raise ValueError('Missing block dependencies')
    ordinary=[v for v in records if not result['steps'][v['step']-1]['restart']]
    groups={}
    for v in ordinary:
        key='/'.join(str(v[k]) for k in ('context','kind','component','family','first','last','tag'))
        a=groups.setdefault(key,{'calls':0,'request_bearing':0,'requests':0,'sent_values':0,'received_values':0,'post_to_wait_seconds':0.,'completion_seconds':0.})
        a['calls']+=1;a['request_bearing']+=int(v['requests']>0)
        for k in ('requests','sent_values','received_values','post_to_wait_seconds','completion_seconds'):a[k]+=v[k]
    if side=='block' and mode>0:
        waits=sum(step['contexts'][c]['calls_global'][19] for step in result['steps'] for c in ('WT','outside_WT'))
        if len(records)!=waits:raise ValueError('Dependency records do not cover inverse completion scopes')
    selected=[step for step in result['steps'] if not step['restart']]
    result['ordinary_counters_per_rank_step']={c:[st.mean(step['contexts'][c]['counters_global'][i]/ranks for step in selected) for i in range(35)] for c in ('WT','outside_WT')}
    if mode>0 and side in ('block','legacy'):
        work=[sum(step['contexts'][c]['counters_global'][i] for step in result['steps'] for c in ('WT','outside_WT')) for i in range(35)]
        if work[33]>work[32]:raise ValueError('Divergence evaluations exceed samples')
        divergence_calls=sum(step['contexts'][c]['calls_global'][81] for step in result['steps'] for c in ('WT','outside_WT'))
        if side=='block' and divergence_calls>0 and work[34]>0 and work[32]==0:
            raise ValueError('Production divergence sample counter missing')
        if work[8]>0 and work[11]==0:raise ValueError('Boundary pack counter still disabled')
        if work[9]>0 and any(work[i]==0 for i in (12,13,14)):raise ValueError('Boundary receive/zero counters still disabled')
    result['dependencies']={'records':len(records),'ordinary_groups_global_sums':groups,
      'scope':'Rank-local elapsed intervals, not synchronized arrival or useful-computation time. Zero-request completions retained. Keys=context/kind/component/family/first/last/tag; -999 means all/unspecified levels.'}
    return result
