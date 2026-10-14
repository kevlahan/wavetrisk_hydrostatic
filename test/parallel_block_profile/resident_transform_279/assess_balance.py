"""Read-only ordinary-step balance assessment of independently verified job 55504.
No solver timing forecast. Counters and profiles come from separate exact-science replays.
Usage: python3 assess278_balance.py RETURN_RESULTS OUTPUT_JSON
"""
import hashlib,json,math,re,statistics as st,sys
from pathlib import Path

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dist(x):
    return dict(min=min(x),mean=st.mean(x),max=max(x),max_over_mean=max(x)/st.mean(x))
def corr(x,y):
    mx,my=st.mean(x),st.mean(y)
    den=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
    return None if den==0 else sum((a-mx)*(b-my) for a,b in zip(x,y))/den

def analyze(root):
    data=json.loads((root/'results.json').read_text())
    assert data['status']=='complete' and str(data['job_id'])=='55504'
    out=dict(job=55504,scope='Validation work counters joined by step/rank to a separate coarse exact-science replay; descriptive association, not causal attribution or speedup forecast.',inputs={'results.json':digest(root/'results.json')},ranks={})
    totals=[]
    for n in (83,160,256):
        matches=lambda phase,profile:[(name,r) for name,r in data['runs'].items() if r['kind']=='candidate' and r['ranks']==n and r['phase']==phase and r['profile']==profile]
        (vn,v),=matches('validation','off');(pn,p),=matches('measurement','coarse')
        assert v['status']==p['status']=='passed'
        keep=[x['index']+1 for x in p['ordinary']['included']]
        assert keep==[2,3,4,5,6,7,10,11,12]
        assert keep==[x['index']+1 for x in v['ordinary']['included']]
        work={};blocks={};domains={};whole=[];costs={}
        for r in range(n):
            path=root/vn/f'rhs278-production-rank-{r}.txt'
            out['inputs'][str(path.relative_to(root))]=digest(path)
            rows=[[int(x) for x in line.split()[:16]] for line in path.read_text().splitlines()]
            assert len(rows)==48 and [x[0] for x in rows]==list(range(1,49))
            assert all(x[0]==x[1] for x in rows)
            path=root/vn/f'rhs278-cost-rank-{r}.txt'
            out['inputs'][str(path.relative_to(root))]=digest(path)
            plans={}
            for line in path.read_text().splitlines():
                fields=[int(x) for x in line.split()[:17]]
                key=tuple(fields[:2]);assert key not in plans;plans[key]=fields
            assert set(plans)=={tuple(row[2:4]) for row in rows}
            for k,row in enumerate(rows):costs[k+1,r]=plans[tuple(row[2:4])]
            prev=None;delta=[]
            for row in rows:
                same=prev is not None and row[2:4]==prev[2:4]
                d=[row[i]-(prev[i] if same else 0) for i in (6,7,8)]
                assert min(d)>=0
                delta.append(d);prev=row
            total=[sum(d[i] for d in delta) for i in range(3)]
            cr=v['contract']['ranks'][r]
            assert cr['rank']==r and total==[cr[k] for k in ('primitive_events','scalar_targets','velocity_targets')]
            whole.append(total)
            for s in range(1,13):
                lo=4*(s-1);rs=rows[lo:lo+4]
                assert len({tuple(x[4:6]) for x in rs})==1
                work[s,r]=[sum(d[i] for d in delta[lo:lo+4]) for i in range(3)]
                domains[s,r],blocks[s,r]=rs[0][4:6]
        path=root/pn/'run.log';out['inputs'][str(path.relative_to(root))]=digest(path)
        log=path.read_text();timing={}
        for line in log.splitlines():
            m=re.match(r'\s*detail step (\d+) rank (\d+) restarts (\d+) remaps (\d+) wall\s+(\S+) CPU\s+(\S+)',line)
            if m:
                key=int(m[1]),int(m[2]);assert key not in timing
                timing[key]=dict(wall=float(m[5]),remaps=int(m[4]))
            m=re.match(r'\s*detail step (\d+) rank (\d+) self-wall = (.*)',line)
            if m:timing[int(m[1]),int(m[2])]['self']=[float(x) for x in m[3].split()]
        assert set(timing)=={(s,r) for s in range(1,13) for r in range(n)}
        assert all(len(t['self'])==84 and abs(sum(t['self'])-t['wall'])<.002*t['wall'] for t in timing.values())
        assert all(sum(blocks[s,r] for r in range(n))==640 for s in range(1,13))
        selected=[]
        for s in keep:
            w=[work[s,r][0] for r in range(n)];a=[1000*timing[s,r]['self'][74] for r in range(n)]
            heavy=max(range(n),key=lambda r:w[r]);slow=max(range(n),key=lambda r:timing[s,r]['wall'])
            singleton=[r for r in range(n) if blocks[s,r]==1]
            entry=dict(step=s,remaps=timing[s,slow]['remaps'],global_counters=[sum(work[s,r][i] for r in range(n)) for i in range(3)],primitive=dist(w),arithmetic_ms=dist(a),primitive_arithmetic_correlation=corr(w,a),max_work_rank=heavy,max_work_rank_blocks=blocks[s,heavy],max_work_rank_arithmetic_ms=a[heavy],singleton_ranks=len(singleton),max_singleton_primitive=max((w[r] for r in singleton),default=None),slowest_window_rank=slow,slowest_window_ms=timing[s,slow]['wall']*1000,groups={})
            for label,rr in [('domain',[r for r in range(n) if domains[s,r]>0]),('zero_domain',[r for r in range(n) if domains[s,r]==0])]:
                if rr:entry['groups'][label]=dict(ranks=len(rr),primitive_mean=st.mean(w[r] for r in rr),primitive_share=sum(w[r] for r in rr)/sum(w),arithmetic_mean_ms=st.mean(a[r] for r in rr),exchange_mean_ms=st.mean(timing[s,r]['self'][77]*1000 for r in rr),inverse_wait_mean_ms=st.mean(timing[s,r]['self'][18]*1000 for r in rr))
            selected.append(entry)
        phases={}
        for key,ids in [('preparation',[68,69,70,71,72,73,84]),('arithmetic',[75]),('rhs_exchange',[78]),('inverse_wait',[19]),('resident_input',[74]),('rk_adapter',[77])]:
            f=lambda s,r:sum(timing[s,r]['self'][i-1] for i in ids)*1000
            phases[key]=dict(mean_all_rank_ms=st.mean(f(s,r) for s in keep for r in range(n)),mean_step_max_ms=st.mean(max(f(s,r) for r in range(n)) for s in keep),mean_on_slowest_window_rank_ms=st.mean(f(e['step'],e['slowest_window_rank']) for e in selected))
        plan_summary={}
        stages=[4*(s-1)+j for s in keep for j in (1,2,3,4)]
        for key,col in [('program_nodes',3),('route_payload_bytes',12),('numerical_buffer_bytes',13),('rounding_column_bytes',14)]:
            plan_summary[key]=dict(mean_stage_global_total=st.mean(sum(costs[k,r][col] for r in range(n)) for k in stages),mean_rank_stage=st.mean(costs[k,r][col] for k in stages for r in range(n)))
        avgwork=[st.mean(work[s,r][0] for s in keep) for r in range(n)]
        avgarith=[st.mean(timing[s,r]['self'][74]*1000 for s in keep) for r in range(n)]
        out['ranks'][str(n)]=dict(validation=vn,profile=pn,ordinary_plan_extents=plan_summary,ordinary_steps=selected,ordinary_primitive_by_rank_summary=dist(avgwork),rank_mean_work_arithmetic_correlation=corr(avgwork,avgarith),mean_step_max_over_mean=st.mean(e['primitive']['max_over_mean'] for e in selected),mean_step_max_primitive=st.mean(e['primitive']['max'] for e in selected),phases=phases,weighted_load_summary_lines=[x.strip() for x in log.splitlines() if 'weighted load min/avg/max' in x],mean_logged_step_ms=st.mean(x['seconds']*1000 for x in p['ordinary']['included']),mean_profile_window_ms=st.mean(e['slowest_window_ms'] for e in selected),trajectory_totals=[sum(w[i] for w in whole) for i in range(3)])
        totals.append([e['global_counters'] for e in selected])
        for i in (19,75,78):
            assert abs(phases[{19:'inverse_wait',75:'arithmetic',78:'rhs_exchange'}[i]]['mean_on_slowest_window_rank_ms']/1000-p['attribution']['regions'][str(i)]['mean_on_slowest_step_rank_seconds'])<1e-10
    assert totals[0]==totals[1]==totals[2]
    out['checks']=['48 transactions/rank and epoch-aware counter deltas','Per-rank trajectory counters agree with verified contract','All nine ordinary-step global primitive/scalar/velocity totals identical at 83/160/256','Complete 12-step, 84-region rank profiles partition each window','Recomputed critical-rank phase means agree with trusted attribution']
    return out

if __name__=='__main__':
    out=analyze(Path(sys.argv[1]));Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
    for n,d in out['ranks'].items():
        print(n,'mean step max/mean',d['mean_step_max_over_mean'],'max work mean',d['mean_step_max_primitive'],'work/arith correlation',d['rank_mean_work_arithmetic_correlation'])
        for s in d['ordinary_steps']: print(s['step'], 'max/mean',round(s['primitive']['max_over_mean'],3),'maxwork',s['primitive']['max'],'rank/blocks',s['max_work_rank'],s['max_work_rank_blocks'],'singleton max',s['max_singleton_primitive'],'arith max',round(s['arithmetic_ms']['max'],3),'corr',round(s['primitive_arithmetic_correlation'],3))
        print('PHASES',json.dumps(d['phases']))
