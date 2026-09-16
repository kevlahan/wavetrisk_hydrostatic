"""Whole-phase counters and timestamp-matched samples; never infer removable time from CPI alone."""
import bisect,collections,csv,json,re,statistics as st
from pathlib import Path
from analyze_operations_199 import common,PHASES

EVENTS={0:['none']*4,1:['cycles','instructions','ex_no_retire.not_complete','ex_no_retire.load_not_complete'],2:['cycles','instructions','l2_cache_req_stat.ls_rd_blk_c','ls_dmnd_fills_from_sys.dram_io_near_or_far']}

def read_rank(path):
    lines=path.read_text().splitlines();h=lines[0].split()
    if len(h)!=5 or h[:2]!=['schema','1']:raise ValueError('Invalid PMU schema')
    mode,rank,pid=map(int,h[2:]);names=lines[1].split()[1:]
    if names!=EVENTS[mode]:raise ValueError('Event identity mismatch')
    rows=[]
    for line in lines[2:]:
        v=line.split()
        if len(v)!=12 or v[0]!='interval':raise ValueError('Invalid interval')
        step,ctx,start,end,cpu,enabled,running,*counts=map(int,v[1:])
        if min(start,end,cpu,enabled,running,*counts)<0 or step<1 or ctx not in range(9) or end<start or running>enabled:raise ValueError('Invalid counter/timestamp')
        if rows and (start<rows[-1]['end'] or step<rows[-1]['step'] or (step==rows[-1]['step'] and start!=rows[-1]['end'])):raise ValueError('Noncontiguous timeline')
        if mode==0 and any([enabled,running,*counts]):raise ValueError('Disabled PMU reported counts')
        rows.append(dict(step=step,context=ctx,start=start,end=end,cpu=cpu,enabled=enabled,running=running,counts=counts))
    if not rows:raise ValueError('Empty timeline')
    return dict(mode=mode,rank=rank,pid=pid,events=names,intervals=rows)

def collect(root,ranks,mode):
    root=Path(root);files=list(root.glob('hardware201-rank-*.txt'))
    if len(files)!=ranks:raise ValueError('Missing PMU rank reports')
    data={};cells=[];flags=[]
    for f in files:
        q=read_rank(f);rank=q['rank']
        if rank in data or rank not in range(ranks) or q['mode']!=mode or f.name!=f'hardware201-rank-{rank}.txt':raise ValueError('Rank/mode mismatch')
        wf=common(root/f'workflow195-rank-{rank}.txt')
        if {x['step'] for x in q['intervals']}!=set(wf):raise ValueError('PMU step coverage mismatch')
        data[rank]=q
        for step,w in wf.items():
            rows=[x for x in q['intervals'] if x['step']==step]
            # Boundary clocks include a few clock/read instructions beyond common timers.
            if abs(sum(x['end']-x['start'] for x in rows)*1e-9-w['wall'])>max(.002,w['wall']*.01):raise ValueError('PMU/coarse whole-step coverage differs')
            if w['restart']:continue
            for phase,contexts in {'whole step':list(range(9)),**PHASES}.items():
                selected=[x for x in rows if x['context'] in contexts]
                counts=[sum(x['counts'][i] for x in selected) for i in range(4)]
                enabled=sum(x['enabled'] for x in selected);running=sum(x['running'] for x in selected)
                wall=sum(x['end']-x['start'] for x in selected)*1e-9
                expected=sum(w['context'][c] for c in contexts)
                if abs(wall-expected)>max(.002,.01*expected):raise ValueError('PMU/coarse phase coverage differs')
                valid=mode==0 or enabled>0 and running/enabled>=.95
                if mode and phase in ('dynamics WT','dynamics outside WT') and (not counts[0] or not counts[1]):valid=False
                if not valid:flags.append({'rank':rank,'step':step,'phase':phase,'reason':'Zero count/less than 95% scheduled; do not scale or interpret.'})
                cells.append(dict(rank=rank,step=step,phase=phase,wall_seconds=wall,thread_cpu_seconds=sum(x['cpu'] for x in selected)*1e-9,counts=counts,time_enabled_ns=enabled,time_running_ns=running,valid=valid))
    phases={}
    for phase in ('whole step',*PHASES):
        cs=[x for x in cells if x['phase']==phase];counts=[st.mean(x['counts'][i] for x in cs) for i in range(4)]
        phases[phase]={'wall_seconds_rank_mean':st.mean(x['wall_seconds'] for x in cs),'thread_cpu_seconds_rank_mean':st.mean(x['thread_cpu_seconds'] for x in cs),'counts_rank_mean_per_step':dict(zip(EVENTS[mode],counts)) if mode else {},'ipc':counts[1]/counts[0] if counts[0] else None,'valid':all(x['valid'] for x in cs)}
    return dict(mode=mode,events=EVENTS[mode],ranks=ranks,phases=phases,cells=cells,quality_flags=flags,
                scope='User-space main-thread counters include MPI execution/spinning. Thread CPU time includes kernel time. Phases are disjoint; whole step overlaps them. No multiplex correction. Counts/stalls are not recoverable wall seconds.')

HEADER=re.compile(r'^\s*(\d+)\s+(\d+)\.(\d+):\s+(\S+):\s*(.*)$')
FRAME=re.compile(r'^\s*[0-9a-fA-F]+\s+(.+?)\s+\((.*)\)\s*$')
def parse_samples(text):
    samples=[];cur=None;lost=False
    for line in text.splitlines():
        if 'LOST' in line.upper():lost=True
        m=HEADER.match(line)
        if m:
            pid,sec,frac,event,rest=m.groups();cur={'pid':int(pid),'ns':int(sec)*10**9+int(frac.ljust(9,'0')[:9]),'event':event,'frames':[]};samples.append(cur)
            f=FRAME.match(rest)
            if f:cur['frames'].append(f.groups())
        elif cur:
            f=FRAME.match(line)
            if f:cur['frames'].append(f.groups())
    return samples,lost

def samples(root,ranks):
    result={}
    for rank in ranks:
        hardware=read_rank(root/f'hardware201-rank-{rank}.txt');timeline=hardware['intervals'];starts=[x['start'] for x in timeline]
        wf=common(root/f'workflow195-rank-{rank}.txt');raw,lost=parse_samples((root/f'perf-rank-{rank}.txt').read_text())
        if not raw:raise ValueError('No parsed samples for rank '+str(rank))
        counts=collections.Counter();leaf=collections.Counter();stack=collections.Counter();outside=0;unknown=0;short=0
        for s in raw:
            if s['pid']!=hardware['pid']:continue
            i=bisect.bisect_right(starts,s['ns'])-1
            if i<0 or s['ns']>=timeline[i]['end']:outside+=1;continue
            row=timeline[i]
            if wf[row['step']]['restart']:continue
            phase=next(name for name,ctx in PHASES.items() if row['context'] in ctx)
            key=phase,s['event'];counts[key]+=1
            frames=[(re.sub(r'\+0x[0-9a-fA-F]+$','',symbol),dso) for symbol,dso in s['frames']]
            if not frames:unknown+=1;continue
            if frames[0][0] in ('[unknown]','[unknown] ([unknown])') or frames[0][1]=='[unknown]':unknown+=1
            if len(frames)<2:short+=1
            leaf[key+tuple(frames[0])]+=1
            # Inclusive stack ancestry; not additive across symbols.
            for symbol,dso in set(frames):stack[key+(symbol,dso)]+=1
        result[str(rank)]={'ordinary_samples':[dict(phase=k[0],event=k[1],count=v) for k,v in sorted(counts.items())],
            'leaf':[dict(phase=k[0],event=k[1],symbol=k[2],dso=k[3],count=v) for k,v in leaf.most_common()],
            'inclusive_stack':[dict(phase=k[0],event=k[1],symbol=k[2],dso=k[3],count=v) for k,v in stack.most_common()],
            'lost_events':lost,'no_frames':unknown,'short_stacks':short,'outside_step_samples':outside,
            'quality_flags':(['Lost samples'] if lost else [])+(['More than 5% unresolved leaves'] if unknown>0.05*sum(counts.values()) else [])+(['More than 20% shallow stacks'] if short>0.2*sum(counts.values()) else [])+(['Sparse dynamics samples; inspect raw stacks'] if any(counts[p,e]<100 for p in ('dynamics WT','dynamics outside WT') for e in ('cycles:u','instructions:u')) else [])}
    return {'ranks':result,'scope':'Selected ranks only, not an unbiased node average. Fixed 20-million-event periods; compare within-event shares and retain event identity, not unweighted mixtures of cycles and instructions. Inclusive ancestors overlap. MPI stacks must be distinguished from solver maintenance. No per-symbol time saving claimed.'}

def summarize(result):
    runs=result['runs'];out={'bare':{},'calibration':{},'counter_comparisons':{},'optimization_ready':False}
    for side in ('candidate','legacy'):
        rs=[r for r in runs.values() if r['side']==side];bare=[r['ordinary_step_median'] for r in rs if r['phase']=='bare'];out['bare'][side]={'samples':bare,'median':st.median(bare)}
        coarse=next(r for r in rs if r['phase']=='coarse');base=out['bare'][side]['median']
        for r in rs:
            if 'hardware' not in r:continue
            row={'printed_change_from_bare_percent':100*(r['ordinary_step_median']/base-1),'phase_changes':{},'review_required':False}
            for p in PHASES:
                a=r['hardware']['phases'][p]['wall_seconds_rank_mean'];b=coarse['hardware']['phases'][p]['wall_seconds_rank_mean'];delta=a-b
                flagged=abs(delta)>.001 and abs(delta)>.02*b
                row['phase_changes'][p]={'seconds':delta,'percent':100*delta/b if b else None,'flag':flagged};row['review_required']|=flagged
            row['review_required']|=abs(row['printed_change_from_bare_percent'])>3 or bool(r['hardware']['quality_flags'])
            out['calibration'][r['label']]=row
    for group in ('core','memory'):
        comparisons=[]
        for repeat in (1,2):
            pair={side:next(r for r in runs.values() if r['side']==side and r['phase']==group and r['round']==repeat) for side in ('candidate','legacy')}
            for phase in ('whole step',*PHASES):
                b,l=(pair[side]['hardware']['phases'][phase] for side in ('candidate','legacy'))
                bc,lc=(x['counts_rank_mean_per_step'] for x in (b,l))
                ib,il=bc['instructions'],lc['instructions'];cb,cl=bc['cycles'],lc['cycles']
                row={'repeat':repeat,'phase':phase,'candidate':b,'legacy':l,
                     'instruction_ratio':ib/il if il else None,'cycles_ratio':cb/cl if cl else None,
                     'instruction_excess_per_rank_step':ib-il,'cycles_excess_per_rank_step':cb-cl,
                     'observer_review_required':any(out['calibration'][pair[side]['label']]['review_required'] for side in pair)}
                if ib and il:
                    cpib,cpil=cb/ib,cl/il
                    row['cpi_ratio']=cpib/cpil if cpil else None
                    row['event_rates']={side:{event:{'per_kilo_instruction':1000*counts[event]/counts['instructions'] if counts['instructions'] else None,'per_cycle':counts[event]/counts['cycles'] if counts['cycles'] else None} for event in EVENTS[1 if group=='core' else 2][2:]} for side,counts in (('candidate',bc),('legacy',lc))}
                    row['symmetric_cycle_accounting']={
                        'instruction_count_component':(ib-il)*(cpib+cpil)/2,
                        'cycles_per_instruction_component':(cpib-cpil)*(ib+il)/2,
                        'scope':'Exact algebra on rank-mean cycles/instructions, not causal wall-time savings; MPI spinning is included.'}
                comparisons.append(row)
        out['counter_comparisons'][group]=comparisons
    out['sample_quality']={r['label']:{rank:q['quality_flags'] for rank,q in r['samples']['ranks'].items()} for r in runs.values() if 'samples' in r}
    out['scope']='Use bare timing for the gap, coarse clocks for locations, paired counters for execution evidence, and selected-rank stacks for source attribution. IPC and load-wait cycles alone do not prove DRAM bandwidth saturation or removable work.'
    return out


def write_reports(root,result):
    summary=result['summary']
    with (root/'hardware-phases.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['run','side','mode','phase','rank','step','wall_seconds','thread_cpu_seconds','cycles','instructions','event2','event3','enabled_ns','running_ns','valid'])
        for label,r in result['runs'].items():
            if 'hardware' not in r:continue
            for c in r['hardware']['cells']:
                w.writerow([label,r['side'],r['hardware']['mode'],c['phase'],c['rank'],c['step'],c['wall_seconds'],c['thread_cpu_seconds'],*c['counts'],c['time_enabled_ns'],c['time_running_ns'],c['valid']])
    with (root/'sample-leaves.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['run','side','rank','phase','event','symbol','dso','samples'])
        for label,r in result['runs'].items():
            if 'samples' not in r:continue
            for rank,q in r['samples']['ranks'].items():
                for row in q['leaf']:w.writerow([label,r['side'],rank,row['phase'],row['event'],row['symbol'],row['dso'],row['count']])
    lines=['# Stage201 hardware review','', 'Numerical/placement gates passed. Review observer and sampling quality before interpretation. No production optimization is approved automatically.','', '| Version | Bare ordinary-step medians | Median seconds |','|---|---|---:|']
    for side,q in summary['bare'].items():lines.append(f"| {side} | {q['samples']} | {q['median']:.6f} |")
    lines+=['','## Complete dynamics: paired counter evidence','','Ratios below are blocks/legacy. Includes MPI execution/spinning; sampled stacks must identify its contribution. Counts are user-space main-thread events, not recoverable wall seconds.','','| Group / repeat | Phase | Instructions ratio | CPI ratio | Observer review flag |','|---|---|---:|---:|---|']
    for group,rows in summary['counter_comparisons'].items():
        for r in rows:
            if r['phase'] not in ('dynamics WT','dynamics outside WT'):continue
            lines.append(f"| {group}/{r['repeat']} | {r['phase']} | {r['instruction_ratio']} | {r.get('cpi_ratio')} | {r['observer_review_required']} |")
    lines+=['','See hardware-summary.json for event rates, paired cycle accounting, phase calibration and sampling flags; hardware-phases.csv retains per-rank/per-step counts. sample-leaves.csv contains exclusive sampled leaves; inclusive stacks remain in results.json and raw exports.','', 'A high CPI or high load-wait count does not establish bandwidth saturation. No exact function saving is inferred from sample shares.']
    (root/'review.md').write_text('\n'.join(lines)+'\n')
