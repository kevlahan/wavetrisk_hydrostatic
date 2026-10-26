"""Small review export of counters and actual hardware endpoint burdens.

Full peer vectors stay native and are hash-bound. Each exported counter remains
separate by phase and rank; bytes and counts are never converted into elapsed time.
"""
from pathlib import Path
import hashlib
import json


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def project(directory,placement):
    directory=Path(directory);source=directory/'work325-rank-vectors.json'
    rows=json.loads(source.read_text())['per_rank'];mapping=placement['mapping']
    n=len(mapping)
    if [r['rank'] for r in rows]!=list(range(n)) or [r['rank'] for r in mapping]!=list(range(n)):
        raise ValueError('Incomplete hardware/rank correspondence')
    projected=[];files=[]
    for rank,row in enumerate(rows):
        records=[];core=mapping[rank]
        for record in row['records']:
            # Directed send/receive volumes. They are not added together.
            totals=[0]*8
            for peer,counts in record['peers'].items():
                peer=int(peer)
                if peer<0 or peer>=n:raise ValueError('Peer outside actual placement')
                sent,received=counts[:2];totals[0]+=sent;totals[1]+=received
                if peer==rank:totals[2]+=sent;totals[3]+=received
                if mapping[peer]['socket']!=core['socket']:totals[4]+=sent;totals[5]+=received
                if mapping[peer]['llc']!=core['llc']:totals[6]+=sent;totals[7]+=received
            records.append([record['tag'],record['phase'],record['calls'],record['sum'],record['max'],totals])
        result={'schema':'stage325-review-rank-v1','rank':rank,'core':core,'native_ledger_sha256':row['sha256'],
                'columns':['tag','phase','calls','counter_sums','counter_maxima','peer_byte_totals'],
                'peer_byte_columns':['sent','received','self_sent','self_received','cross_socket_sent','cross_socket_received','cross_llc_sent','cross_llc_received'],
                'records':records}
        p=directory/f'work325-review-rank-{rank}.json';p.write_text(json.dumps(result,separators=(',',':'))+'\n')
        if p.stat().st_size>256*1024:raise ValueError('Compact rank review exceeds256KiB; native evidence preserved')
        files.append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
        projected.append(records)
    families=sorted({(r[0],r[1]) for records in projected for r in records})
    endpoint=[]
    for tag,phase in families:
        matched=[next((r for r in records if r[:2]==[tag,phase]),None) for records in projected]
        vectors=[r[5] if r else [0]*8 for r in matched]
        summed=[sum(v[i] for v in vectors) for i in range(8)]
        if summed[0]!=summed[1] or summed[2]!=summed[3] or summed[4]!=summed[5] or summed[6]!=summed[7]:
            raise ValueError('Projected directed endpoint conservation failed')
        endpoint.append({'tag':tag,'phase':phase,'global_byte_totals':summed,
                         'rank_byte_maxima':[max(v[i] for v in vectors) for i in range(8)]})
    result={'passed':True,'schema':'stage325-hardware-review-v1','full_native_vectors_sha256':sha(source),
            'rank_files':files,'endpoint_burdens':endpoint,
            'limits':['Counters retain distinct phases; no cross-phase critical-path sum.',
                      'Endpoint bytes include payloads recorded by each named hook, not MPI internal traffic.',
                      'No time, cache-miss or memory-bandwidth inference. Core/L3 grouping uses actual verified placement.']}
    p=directory/'work325-hardware-review.json';p.write_text(json.dumps(result,separators=(',',':'))+'\n')
    return {'passed':True,'file':p.name,'sha256':sha(p),'rank_review_bytes':sum(f['bytes'] for f in files)}
