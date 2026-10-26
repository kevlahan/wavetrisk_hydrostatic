#!/usr/bin/env python3
"""Fail-closed compact work ledger validation. Counts never imply elapsed cause."""
from pathlib import Path
import argparse, hashlib, json, re

MAX_FILE=8*1024*1024
TAG=re.compile(r'^[a-z][a-z0-9.]*$')
LABELS=('checked-L256','qual-C0','qual-L256')

def insist(condition,message):
    if not condition: raise ValueError('account325.'+message)

def parse(path,rank,schema):
    path=Path(path);insist(0<path.stat().st_size<=MAX_FILE,'size')
    with path.open('rb') as stream: data=stream.read(MAX_FILE+1)
    insist(0<len(data)<=MAX_FILE,'size')
    try: lines=data.decode('ascii').splitlines()
    except UnicodeDecodeError: raise ValueError('account325.ascii')
    insist(lines and lines[0]==f'RK325WORK1 {rank}','header')
    entries={}; peers={}; footer=None;rows=0
    for line in lines[1:]:
        p=line.split()
        insist(p and footer is None,'trailing_or_blank')
        if p[0]=='COMPLETE':
            insist(len(p)==3,'footer_extent')
            footer=tuple(map(int,p[1:])); continue
        insist(p[0] in ('C','P') and len(p)>3 and TAG.fullmatch(p[1]),'record')
        tag=p[1];insist(tag in schema,'unknown_tag')
        try: v=list(map(int,p[2:]))
        except ValueError: raise ValueError('account325.integer')
        insist(all(0<=x<2**63 for x in v),'range')
        phase=v[0];key=(tag,phase);rows+=1
        if p[0]=='C':
            insist(key not in entries and len(v)>=5,'duplicate_or_short_counter')
            calls,width=v[1:3]
            insist(calls>0 and width==len(schema[tag]['fields']) and len(v)==3+2*width,'counter_extent')
            total=v[3:3+width];peak=v[3+width:]
            insist(all(m<=t<=calls*m for t,m in zip(total,peak)),'peak_sum')
            entries[key]={'calls':calls,'sum':total,'max':peak,'peers':{}}
        else:
            insist(schema[tag]['peer_unit'] is not None and len(v)==6,'peer_extent')
            peer=v[1];pid=(tag,phase,peer)
            insist(pid not in peers,'duplicate_peer')
            insist(0<max(v[2],v[3]) and v[4]<=v[2] and v[5]<=v[3],'peer_peak_sum')
            peers[pid]=v[2:]
    insist(footer==(len(entries),rows),'footer')
    for (tag,phase,peer),v in peers.items():
        insist((tag,phase) in entries,'orphan_peer')
        e=entries[tag,phase]
        insist(v[0]<=e['calls']*v[2] and v[1]<=e['calls']*v[3],'peer_call_bound')
        e['peers'][peer]=v
    return {'rank':rank,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'entries':entries}

def serialized_rank(obj):
    return {'rank':obj['rank'],'bytes':obj['bytes'],'sha256':obj['sha256'],
            'records':[{'tag':tag,'phase':phase,**entry} for (tag,phase),entry in sorted(obj['entries'].items())]}

def validate_run(path,nrank=256,schema=None,expected_rebuilds=5,arm=None):
    path=Path(path)
    if schema is None:
        choices=[Path(__file__).resolve().parent/'work325-schema.json',path/'work325-schema.json',path.parent/'work325-schema.json',path.parent.parent/'work325-schema.json']
        found=next((p for p in choices if p.exists()),None)
        insist(found is not None,'schema_missing')
        schema=json.loads(found.read_text())['records']
    elif 'records' in schema: schema=schema['records']
    files=list(path.glob('work325-rank-*.tsv'))
    insist(len(files)==nrank,'rank_file_count')
    ranks=[]
    for rank in range(nrank):
        p=path/f'work325-rank-{rank}.tsv';insist(p.exists(),'rank_file_missing')
        ranks.append(parse(p,rank,schema))
    # Every actual sender/receiver pair reconciles, including self and headers.
    for r in ranks:
        for key,e in r['entries'].items():
            for peer,(sent,received,sp,rp) in e['peers'].items():
                insist(peer<nrank,'peer_range')
                reverse=ranks[peer]['entries'].get(key,{}).get('peers',{}).get(r['rank'],[0,0,0,0])
                insist(sent==reverse[1] and received==reverse[0],'peer_conservation')
    required=('migration','rhs.geometry','transform.geometry','route.exchange','transform.workspace','selector.call')
    summaries={}
    for tag in schema:
        present=[(r,key,e) for r in ranks for key,e in r['entries'].items() if key[0]==tag]
        if tag in required: insist(len({r['rank'] for r,key,e in present})==nrank,'required_'+tag)
        for phase in sorted({key[1] for r,key,e in present}):
            found=[(r,e) for r,key,e in present if key[1]==phase]
            n=len(schema[tag]['fields'])
            summaries[f'{tag}:{phase}']={'fields':schema[tag]['fields'],'active_ranks':len(found),
                'call_range':[min(e['calls'] for r,e in found),max(e['calls'] for r,e in found)],
                'sum':[sum(e['sum'][i] for r,e in found) for i in range(n)],
                'max_rank_total':[max(e['sum'][i] for r,e in found) for i in range(n)],
                'max_single_call':[max(e['max'][i] for r,e in found) for i in range(n)]}
    geo=summaries['rhs.geometry:0']
    insist(geo['call_range'][0]==geo['call_range'][1] and geo['sum'][0]==352*geo['call_range'][0],
           'geometry_source_selection_census')
    migration=[r['entries']['migration',0] for r in ranks]
    insist(all(e['calls']==expected_rebuilds for e in migration),'migration_lifecycle')
    source=[(r,key[1],e) for r in ranks for key,e in r['entries'].items() if key[0]=='migration.source']
    insist(len(source)==640 and {b for r,b,e in source}==set(range(1,641)),'source_census')
    for r,b,e in source:
        insist(r['rank']==(b-1)//4 and e['calls']==expected_rebuilds,'source_home_lifecycle')
        insist(e['max'][0]>0 and e['max'][1] in (0,1),'source_extent_remote')
        insist(e['sum'][0]==expected_rebuilds*e['max'][0] and
                e['sum'][1]==expected_rebuilds*e['max'][1],'source_fixed_case')
    remote=634 if arm=='C0' else 240 if arm=='L256' else None
    if remote is not None:
        insist(sum(e['sum'][1] for r,b,e in source)==remote*expected_rebuilds,'remote_map_census')
        for field,want in ((0,640),(1,640-remote),(2,remote),(3,remote),(7,17*nrank)):
            insist(sum(e['sum'][field] for e in migration)==want*expected_rebuilds,'migration_count_'+str(field))
        peak=4 if arm=='C0' else 3
        insist(max(e['max'][2] for e in migration)==peak and max(e['max'][3] for e in migration)==3,'migration_endpoint_peak')
    full_sent=sum(e['sum'][0] for r,b,e in source if e['max'][1])
    insist(full_sent==sum(e['sum'][4] for e in migration)==sum(e['sum'][5] for e in migration),'full_payload_bytes')
    for r in ranks:
        for key,e in r['entries'].items():
            tag=key[0]; sent=sum(v[0] for v in e['peers'].values());recv=sum(v[1] for v in e['peers'].values())
            if tag=='migration': si,ri=4,5
            elif tag=='rhs.geometry': si,ri=1,2
            elif tag=='transform.geometry': si,ri=0,1
            else: continue
            insist(sent==e['sum'][si] and recv==e['sum'][ri],'peer_counter_reconciliation')
    vectors=path/'work325-rank-vectors.json'
    vectors.write_text(json.dumps({'schema':'stage325-work-vectors-v1','per_rank':[serialized_rank(r) for r in ranks]},separators=(',',':'))+'\n')
    return {'passed':True,'status':'PASS','schema':'stage325-work-run-v1','arm':arm,'nrank':nrank,
            'expected_rebuilds':expected_rebuilds,'bytes':sum(r['bytes'] for r in ranks),
            'summary':summaries,'per_rank_artifact':{'path':str(vectors),'bytes':vectors.stat().st_size,
                'sha256':hashlib.sha256(vectors.read_bytes()).hexdigest()},
            'limitations':['Per-phase sums and maxima are counts, not elapsed attribution.',
                'Phase maxima are never added across phases or counter families.',
                'State, provider and ordered-program equivalence require the separate semantic gate.']}

def qualify(results,kit=None):
    results=Path(results);kit=Path(kit) if kit else results
    paths=[kit/'work325-schema.json',results/'work325-schema.json']
    schema_path=next((p for p in paths if p.exists()),None)
    insist(schema_path is not None,'schema_missing')
    schema=json.loads(schema_path.read_text());runs={}
    for label in LABELS:
        folder=results/label
        # Runs may put solver cwd below a label; exactly one complete ledger directory.
        candidates={p.parent for p in folder.rglob('work325-rank-0.tsv')}
        insist(len(candidates)==1,'ledger_directory_'+label)
        run=validate_run(next(iter(candidates)),schema=schema,arm='C0' if label=='qual-C0' else 'L256')
        out=results/(label+'-work325.json');out.write_text(json.dumps(run,separators=(',',':'))+'\n')
        runs[label]={'report':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
                     'bytes':run['bytes'],'summary':run['summary']}
    # Numerical route peer permutation and physics equivalence are not silently inferred here.
    extents={}
    for label in LABELS:
        native=json.loads((results/(label+'-work325.json')).read_text())
        vector=json.loads(Path(native['per_rank_artifact']['path']).read_text())
        extents[label]={e['phase']:e['max'][0] for row in vector['per_rank']
                        for e in row['records'] if e['tag']=='migration.source'}
    insist(extents['checked-L256']==extents['qual-C0']==extents['qual-L256'],
           'paired_source_serialized_extents')
    a=runs['qual-C0']['summary'];b=runs['qual-L256']['summary']
    geometry=a['rhs.geometry:0']['sum'][0],b['rhs.geometry:0']['sum'][0]
    insist(geometry[0]==geometry[1],'geometry_selection_work')
    report={'passed':True,'status':'PASS_ACCOUNTING_NOT_PERFORMANCE','runs':runs,
            'decision':'AWAIT_SEMANTIC_AND_NATIVE_QUALIFICATION; NO_AUTOMATIC_PERFORMANCE_PROMOTION',
            'schema_sha256':hashlib.sha256(schema_path.read_bytes()).hexdigest(),
            'performance_qualified':False}
    (results/'account325.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('results');p.add_argument('--kit');a=p.parse_args()
    print(json.dumps(qualify(a.results,a.kit),indent=2))
