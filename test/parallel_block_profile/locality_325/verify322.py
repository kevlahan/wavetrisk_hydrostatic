"""Fail-closed Stage322 numerical-identity and route verification.
No elapsed/performance claim is made here. All source streams retain exact bits.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
import gzip, hashlib, json, math, re, struct, sqlite3, tempfile
import numpy as np
from coverage322 import validate_coverage322

MAGIC = b'RK316V01'
LIMIT = 512 * 1024 * 1024
WIDTH = {1:4, 2:8, 3:8, 4:1, 10:8, 11:8, 12:8}
HEADER = struct.Struct('<8s5i')
RECORD = struct.Struct('<iiq')

class AuditError(ValueError):
    pass

@dataclass
class Entry:
    index: int
    offset: int
    name: str
    kind: int
    count: int
    data: bytes

    def descriptor(self):
        return [self.name, self.kind, self.count]

def open_stream(path):
    path = Path(path)
    return gzip.open(path, 'rb') if path.suffix == '.gz' else path.open('rb')

def file_for(directory, rank):
    stem = Path(directory)/f'execution316-rank-{rank}.bin'
    matches = [p for p in (stem, Path(str(stem)+'.gz')) if p.is_file()]
    if len(matches) != 1:
        raise AuditError(f'Expected exactly one raw or gzip stream for rank {rank}, found {len(matches)}')
    return matches[0]

def exact_read(f, n, label):
    data = f.read(n)
    if len(data) != n:
        raise AuditError(f'Truncated {label}: expected {n}, received {len(data)} bytes')
    return data

class Reader:
    def __init__(self, path, rank=None):
        self.path = Path(path)
        self.expected_rank = rank
        self.summary = None

    def __iter__(self):
        digest = hashlib.sha256()
        with open_stream(self.path) as f:
            raw = exact_read(f, HEADER.size, 'file header')
            magic, version, rank, integer_size, real_size, endian = HEADER.unpack(raw)
            if (magic, version, integer_size, real_size, endian) != (MAGIC,1,4,8,0x01020304):
                raise AuditError('Invalid header, unsupported byte order or representation')
            if rank < 0 or (self.expected_rank is not None and rank != self.expected_rank):
                raise AuditError('Audit rank differs from filename/expected rank')
            digest.update(raw)
            offset, index, done = len(raw), 0, False
            while True:
                head = f.read(RECORD.size)
                if not head:
                    break
                if done:
                    raise AuditError('Records after complete footer')
                if len(head) != RECORD.size:
                    raise AuditError('Truncated record header')
                length, kind, count = RECORD.unpack(head)
                if not 1 <= length <= 256 or kind not in WIDTH or count < 0:
                    raise AuditError('Invalid record descriptor')
                size = count * WIDTH[kind]
                if size > LIMIT or offset + RECORD.size + length + size > LIMIT:
                    raise AuditError('Audit file exceeds sealed size bound')
                tag = exact_read(f, length, 'tag')
                try:
                    name = tag.decode('ascii')
                except UnicodeError as exc:
                    raise AuditError('Non-ASCII record tag') from exc
                if any(ord(c) < 32 or ord(c) > 126 for c in name):
                    raise AuditError('Non-printable record tag')
                data = exact_read(f, size, 'record payload')
                if kind==3 and not np.isfinite(np.frombuffer(data,dtype='<f8')).all():
                    raise AuditError('record.nonfinite_real: '+name)
                if kind in (10,11) and count != 1:
                    raise AuditError('Invalid scope epoch')
                if kind == 4 and any(b not in (0,1) for b in data):
                    raise AuditError('Invalid logical encoding')
                if kind == 12:
                    if name != 'complete' or count != 2:
                        raise AuditError('Invalid complete footer')
                    if struct.unpack('<2q', data) != (index, offset):
                        raise AuditError('Footer record count/byte count mismatch')
                    done = True
                entry = Entry(index, offset, name, kind, count, data)
                for part in (head,tag,data):
                    digest.update(part)
                offset += RECORD.size + length + size
                index += 1
                yield entry
            if not done:
                raise AuditError('Missing complete footer')
            self.summary = dict(rank=rank, records=index, uncompressed_bytes=offset,
                                stream_sha256=digest.hexdigest(), file=self.path.name)

def integers(entry):
    if entry.kind not in (1,2,10,11,12):
        raise AuditError('Expected integer record: '+entry.name)
    code = 'i' if entry.kind == 1 else 'q'
    return list(struct.unpack('<'+str(entry.count)+code,entry.data))


def require(ok,code,detail=''):
    if not ok:raise AuditError(code+(': '+str(detail) if detail else ''))

def value(e):
    require(e.kind in (1,2,4,10,11,12),'record.integer_kind',e.name)
    return list(e.data) if e.kind==4 else integers(e)

def prefix(count):
    result=[];total=0
    for n in count:
        require(isinstance(n,int) and n>=0,'route.negative_count')
        result.append(total);total+=n
    return result

def digest_entry(e):
    return (e.name,e.kind,e.count,hashlib.sha256(e.data).hexdigest())

ROUTE_FIELDS=tuple('route.allocated route.epoch route.header route.flags route.sc route.sd route.rc route.rd route.source route.target route.address route.orientation route.send_peer route.receive_peer route.target_peer route.allocations'.split())

def check_route(q,n,rank,empty=False):
    require(q.get('route.allocated') in ([0],[1]),'route.allocated')
    if q['route.allocated']==[0]:
        require(set(q)=={'route.allocated'},'route.unallocated_extra')
        return
    require(set(q)==set(ROUTE_FIELDS),'route.fields')
    h=q['route.header'];require(len(h)==7 and h[3:5]==[rank,n],'route.header')
    require(h[1]>=0 and h[2]>0,'route.extent')
    for c,d in [('sc','sd'),('rc','rd')]:
        require(len(q['route.'+c])==n,'route.world_extent',c)
        require(q['route.'+d]==prefix(q['route.'+c]),'route.prefix_sum',d)
    sc,rc=q['route.sc'],q['route.rc']
    require(len(q['route.source'])==sum(sc),'route.source_extent')
    require(all(1<=x<=h[1] for x in q['route.source']),'route.source_slot')
    count=len(q['route.target'])
    require(all(len(q['route.'+k])==count for k in ('address','orientation','target_peer')),'route.target_extent')
    require(all(1<=x<=h[1] for x in q['route.target']),'route.target_slot')
    require(all(x in (-1,1) for x in q['route.orientation']),'route.orientation')
    require(q['route.send_peer']==[i+1 for i,x in enumerate(sc) if x],'route.send_order')
    require(q['route.receive_peer']==[i+1 for i,x in enumerate(rc) if x],'route.receive_order')
    for addr,peer in zip(q['route.address'],q['route.target_peer']):
        require(1<=peer<=n,'route.provider_rank')
        require(q['route.rd'][peer-1]<=addr<q['route.rd'][peer-1]+rc[peer-1],'route.provider_address')
    allocation=q['route.allocations'];request=len(q['route.send_peer'])+len(q['route.receive_peer'])
    require(allocation in ([0,0,request],[h[2]*sum(sc)+n,h[2]*sum(rc)+n,request]),'route.buffer_formula')
    if empty:
        require(sum(sc)==sum(rc)==count==0,'route.empty_rank_payload')

def reconcile_routes(routes):
    """Resolve every receiver slot from the actual producer wire slot.
    Returns ordered producer identities, never infers identity from equal values.
    """
    n=len(routes);resolved=[]
    require(all(q.get('route.allocated')==routes[0].get('route.allocated') for q in routes),'route.collective_allocation')
    if routes[0]['route.allocated']==[0]:return [[] for _ in routes]
    for receiver,q in enumerate(routes):
        check_route(q,n,receiver)
        out=[]
        for producer,p in enumerate(routes):
            require(q['route.rc'][producer]==p['route.sc'][receiver],'route.wire_count',[receiver,producer])
        for target,address,peer,orientation in zip(q['route.target'],q['route.address'],q['route.target_peer'],q['route.orientation']):
            producer=peer-1;p=routes[producer];ordinal=address-q['route.rd'][producer]
            source=p['route.source'][p['route.sd'][receiver]+ordinal]
            out.append((target,producer,source,orientation,ordinal))
        resolved.append(out)
    return resolved

def check_adapter(q,n,rank,empty=False):
    for c,d in [('sc','sd'),('rc','rd')]:
        require(len(q['plan%'+c])==n,'adapter.world_extent')
        require(q['plan%'+d]==prefix(q['plan%'+c]),'adapter.prefix_sum')
    columns=q['extent'][2];sc,rc=q['plan%sc'],q['plan%rc']
    require(q['plan%source.shape']==[3,sum(rc)] and len(q['plan%source'])==3*sum(rc),'adapter.source_extent')
    require(q['plan%sent.shape']==[columns,sum(rc)] and q['plan%received.shape']==[columns,sum(sc)],'adapter.buffer_formula')
    require(q['plan%send_peer']==[i+1 for i,x in enumerate(rc) if x],'adapter.send_order')
    require(q['plan%receive_peer']==[i+1 for i,x in enumerate(sc) if x],'adapter.receive_order')
    require(len(q['plan%target'])==len(q['plan%sign'])==len(q['plan%index']),'adapter.target_extent')
    require(all(x in (-1,1) for x in q['plan%sign']),'adapter.orientation')
    if empty:require(not sum(sc) and not sum(rc) and not q['plan%target'],'adapter.empty_rank_payload')

def reconcile_adapters(adapters,domains):
    """Transform adapters use rc for sends, sc for receives (opposite RHS)."""
    n=len(adapters);result=[]
    for receiver,q in enumerate(adapters):
        check_adapter(q,n,receiver);rows=[]
        for producer,p in enumerate(adapters):
            require(q['plan%sc'][producer]==p['plan%rc'][receiver],'adapter.wire_count')
        for target,index,sign in zip(q['plan%target'],q['plan%index'],q['plan%sign']):
            # Fortran index refers to the packed received column, one based.
            require(1<=index<=sum(q['plan%sc']),'adapter.receive_index')
            producer=next((p for p in range(n) if q['plan%sd'][p]<=index-1<q['plan%sd'][p]+q['plan%sc'][p]),None)
            require(producer is not None,'adapter.provider_rank')
            ordinal=index-1-q['plan%sd'][producer];p=adapters[producer]
            slot=p['plan%rd'][receiver]+ordinal;local,node,component=p['plan%source'][3*slot:3*slot+3]
            require(1<=local<=len(domains[producer]),'adapter.provider_domain')
            rows.append((target,producer,domains[producer][local-1]+1,node,component,sign,ordinal))
        result.append(rows)
    return result

@__import__('functools').lru_cache(maxsize=1)
def _world_fields():
    schema=json.loads((Path(__file__).with_name('schema322.json')).read_text())
    return set(schema.get('world_array_fields',[]))|{'route.sc','route.sd','route.rc','route.rd','plan%sc','plan%sd','plan%rc','plan%rd'}

GHOST_LEVEL_WORLD_FIELDS=frozenset('ghost_exchange_plan%'+direction+'_level_'+suffix
    for direction in ('send','recv') for suffix in ('count','displ'))

def normalize_entry(e,n,target=160):
    """Closed numerical shape normalization. No metadata field is ignored."""
    if e.kind in (10,11):return (e.name,e.kind,e.count,value(e))
    if e.kind==12:return ('complete',12) # Reader independently verifies original footer totals.
    # Source-specific corrections to the sealed world-extent inventory.
    # Full parsing independently proves each level matrix equals its existing
    # record counts/prefixes, and requires the exact inverse request capacity.
    # Actual metadata/storage sizes remain in the unmodified raw streams.
    if e.name in GHOST_LEVEL_WORLD_FIELDS or (e.name.endswith('.shape') and
            e.name[:-6] in GHOST_LEVEL_WORLD_FIELDS) or e.name=='requests.shape':
        require(e.kind==1,'normalization.fixed_world_kind',e.name)
        require(n==target or (n,target)==(256,160),'normalization.world_pair')
        x=value(e)
        if e.name=='requests.shape':
            require(x==[2*n],'normalization.inverse_request_capacity')
            return (e.name,e.kind,1,[2*target])
        if e.name.endswith('.shape'):
            require(x==[n,1],'normalization.ghost_level_shape',e.name)
            return (e.name,e.kind,2,[target,1])
        require(len(x)==n and all(v>=0 for v in x),'normalization.ghost_level_extent',e.name)
        if n!=target:
            if e.name.endswith('_count'):
                require(not any(x[target:]),'normalization.nonzero_extra_peer',e.name)
            else:
                require(len(set(x[target:]))==1 and x[target]>=x[target-1],
                        'normalization.displacement_tail',e.name)
        return (e.name,e.kind,target,x[:target])
    if e.kind==3:return digest_entry(e)
    x=value(e)
    if n==target:return (e.name,e.kind,len(x),x)
    require((n,target)==(256,160),'normalization.world_pair')
    name=e.name;world=_world_fields()
    if name=='route.header':
        require(len(x)==7 and x[4]==n,'normalization.route_header');x[4]=target
    elif name in world:
        require(len(x)==n,'normalization.world_extent',name)
        if name.endswith(('displ','.sd','.rd','%sd','%rd')):
            # Per-scope validation verifies full prefix sums; tail must repeat total.
            require(len(set(x[target:]))<=1 and (not x[target:] or x[target]>=x[target-1]),'normalization.displacement_tail',name)
        else:require(all(v==0 for v in x[target:]),'normalization.nonzero_extra_peer',name)
        x=x[:target]
    elif name.endswith('.shape') and name[:-6] in world:
        if x!=[-1]:require(x==[n],'normalization.world_shape',name);x=[target]
    elif name=='main.identity':
        require(len(x)>=4 and x[3]==n,'normalization.main_identity');x[3]=target
    elif name=='route.allocations':
        require(len(x)==3,'normalization.route_allocations')
        if x[:2]!=[0,0]:x[0]-=n-target;x[1]-=n-target
    elif name in ('transport.scratch.shape','exchange_scratch.shape'):
        require(len(x)==1 and x[0]>=2*n,'normalization.scratch_extent',name)
        x[0]-=2*(n-target)
    return (name,e.kind,len(x),x)

def parse_stream(path,rank,n,arm,coverage='candidate'):
    reader=Reader(path,rank);stack=[];scopes=[];driver=None
    for e in reader:
        if e.kind==10:
            parent=tuple((s['name'],s['epoch']) for s in stack)
            step=next((s['epoch'] for s in stack if s['name']=='step'),None)
            s=dict(name=e.name,epoch=value(e)[0],step=step,parent=parent,record=e.index,
                   tags=set(),values={},records=[],record_counts=Counter(),record_descriptors={})
            stack.append(s);scopes.append(s)
        elif e.kind==11:
            require(stack and (e.name,value(e)[0])==(stack[-1]['name'],stack[-1]['epoch']),'coverage.scope_end')
            stack.pop()['end_record']=e.index
        elif e.kind==12:require(not stack,'coverage.scope_footer')
        elif stack:
            # Reals are retained as exact SHA256 bit digests, not Python floats.
            stored=Entry(e.index,e.offset,e.name,e.kind,e.count,e.data)
            s=stack[-1];s['records'].append(stored);s['tags'].add(e.name);s['values'][e.name]=stored
            s['record_counts'][e.name]+=1;s['record_descriptors'][e.name]=(e.kind,e.count)
        else:
            require(e.name=='driver.steps' and driver is None,'coverage.root_record',e.name);driver=value(e)
    require(driver==[4],'coverage.driver_steps');require(reader.summary is not None,'coverage.footer')
    summary=validate_coverage322(scopes,rank,n,4,coverage,arm)
    routes={};adapters={};static={};seen=Counter();domains=None;route_scopes=[]
    for s in scopes:
        if s['name']=='catalog.executed':domains=value(s['values']['domain.local_glo_id'])
        path=tuple(s['parent'])+((s['name'],s['epoch']),)
        # Same path can legitimately repeat adapter invocations; retain order.
        seen[path]+=1;path=repr((path,seen[path]))
        q=None;serial=0
        for e in s['records']:
            if e.name=='route.allocated':
                if q is not None:
                    check_route(q,n,rank,arm in ('B','checked-B') and rank>=160)
                    routes[path+f'/route/{serial}']=q;route_scopes.append((s,q));serial+=1
                q={}
            if e.name.startswith('route.'):
                require(q is not None and e.name not in q,'route.record_order',e.name);q[e.name]=value(e)
        if q is not None:
            check_route(q,n,rank,arm in ('B','checked-B') and rank>=160)
            routes[path+f'/route/{serial}']=q;route_scopes.append((s,q))
        if s['name'] in ('transform.scalar_adapter','transform.vector_adapter'):
            q={e.name:value(e) for e in s['records']};check_adapter(q,n,rank,arm in ('B','checked-B') and rank>=160)
            adapters[path]=q
        if s['name'] in ('ghost.plan','writeback.plan','inverse.interior','inverse.level_interior','inverse.level_boundary','inverse.outer_install','inverse.alias','inverse.range_alias') or (s['name']=='inverse.boundary' and any(x[0]=='inverse.plan' for x in s['parent'])):
            static[path]={'__scope':s['name'],**{e.name:value(e) for e in s['records']}}
        validate_world_arrays(s,n,arm,rank)
        if s['name']=='rhs.prepare':
            qs=[v for k,v in routes.items() if k.startswith(path+'/route/') and v['route.allocated']==[1]]
            validate_shared_scratch(value(s['values']['transport.scratch.shape']),qs,n,'route.shared_scratch_formula')
    for plan in (s for s in scopes if s['name']=='transform.plan'):
        qs=[q for s,q in route_scopes if plan['record']<s['record']<s['end_record']<plan['end_record'] and q['route.allocated']==[1]]
        validate_shared_scratch(value(plan['values']['exchange_scratch.shape']),qs,n,'adapter.shared_scratch_formula')
    return dict(summary=summary,reader=reader.summary,routes=routes,adapters=adapters,static=static,domains=domains,scopes=scopes)

def validate_shared_scratch(shape,routes,n,code):
    expected=max((q['route.header'][2]*(sum(q['route.sc'])+sum(q['route.rc']))+2*n
                  for q in routes if q['route.allocated']==[1]),default=0)
    require(shape==[expected],code)

def validate_world_arrays(s,n,arm,rank):
    fields=_world_fields();records={e.name:e for e in s['records']}
    for name in fields:
        if name not in records:continue
        # Repeated generic route arrays are validated individually above.
        if name.startswith('route.'):continue
        x=value(records[name]);require(len(x)==n,'route.world_extent',name)
        if name.endswith('displ'):other=name[:-5]+'count'
        elif name.endswith(('sd','rd')):other=name[:-1]+'c'
        else:other=None
        if other:
            require(other in records and x==prefix(value(records[other])),'route.prefix_sum',name)
        else:require(all(v>=0 for v in x),'route.negative_count',name)
        if arm in ('B','checked-B') and not other:
            require(not any(x[160:]),'route.extra_rank_count',name)
            if rank>=160:require(not any(x),'route.empty_rank_count',name)

def compare_entries(a,b):
    for i,(x,y) in enumerate(zip_longest(a,b)):
        require(x==y,'identity.record_difference',{'record':i,'left':str(x)[:240],'right':str(y)[:240]})


def _stream_set(directory,n):
    found=list(Path(directory).glob('execution316-rank-*.bin'))+list(Path(directory).glob('execution316-rank-*.bin.gz'))
    require(len(found)==n,'coverage.rank_count',{'expected':n,'found':len(found)})
    return [file_for(directory,rank) for rank in range(n)]

def _cache_open(path):
    c=sqlite3.connect(path)
    c.execute('create table route(arm text, rank integer, kind text, key text, payload text, primary key(arm,rank,kind,key))')
    c.execute('create table domain(arm text, rank integer, payload text, primary key(arm,rank))')
    return c

def _cache_rank(c,arm,rank,obj):
    for kind in ('routes','adapters','static'):
        c.executemany('insert into route values (?,?,?,?,?)',[(arm,rank,kind,k,json.dumps(v,separators=(',',':'))) for k,v in obj[kind].items()])
    c.execute('insert into domain values (?,?,?)',(arm,rank,json.dumps(obj['domains'])))
    c.commit()

def _reconcile_cache(c,arm,n,retain):
    results={};counts=Counter()
    domains=[json.loads(row[0]) for row in c.execute('select payload from domain where arm=? order by rank',(arm,))]
    require(len(domains)==n,'route.domain_rank_count')
    keys=c.execute('select kind,key,count(*) from route where arm=? group by kind,key order by kind,key',(arm,)).fetchall()
    for kind,key,num in keys:
        require(num==n,'route.missing_rank_or_plan',{'arm':arm,'family':key,'ranks':num})
        rows=c.execute('select rank,payload from route where arm=? and kind=? and key=? order by rank',(arm,kind,key)).fetchall()
        require([r for r,_ in rows]==list(range(n)),'route.rank_order')
        plans=[json.loads(p) for _,p in rows]
        output=reconcile_routes(plans) if kind=='routes' else (reconcile_adapters(plans,domains) if kind=='adapters' else reconcile_static(plans))
        if arm in ('B','checked-B'):
            require(not any(output[160:]),'route.extra_rank_consumer')
            if kind!='static':require(all(all(x[1]<160 for x in rank_rows) for rank_rows in output[:160]),'route.extra_rank_provider')
        # Actual send-slot/provider identities and receiver order are canonical.
        results[(kind,key)]=hashlib.sha256(json.dumps(output[:retain],separators=(',',':')).encode()).hexdigest()
        counts[kind]+=1
    return results,dict(counts)

def inspect_one_run(directory,n,arm,c,coverage='candidate'):
    paths=_stream_set(directory,n);rows=[]
    for rank,path in enumerate(paths):
        obj=parse_stream(path,rank,n,arm,coverage)
        _cache_rank(c,arm,rank,obj);rows.append(obj['reader'])
    return {'arm':arm,'ranks':n,'passed':True,'streams':rows}

def compare_runs(a,b,na,nb,arm_a,arm_b,c):
    aa=_stream_set(a,na);bb=_stream_set(b,nb);rows_a=[];rows_b=[]
    require(nb>=na,'identity.rank_relation')
    for rank in range(nb):
        right=parse_stream(bb[rank],rank,nb,arm_b)
        _cache_rank(c,arm_b,rank,right);rows_b.append(right['reader']);del right
        if rank<na:
            left=parse_stream(aa[rank],rank,na,arm_a)
            _cache_rank(c,arm_a,rank,left);rows_a.append(left['reader']);del left
            # Second streaming pass compares one raw record per side. Never retain
            # a rank's full numerical values or normalize them into Python floats.
            for index,(x,y) in enumerate(zip_longest(Reader(aa[rank],rank),Reader(bb[rank],rank))):
                require(x is not None and y is not None,'identity.record_count',rank)
                require(normalize_entry(x,na,na)==normalize_entry(y,nb,na),
                        'identity.record_difference',dict(rank=rank,index=index,field=x.name))
    return {'passed':True,'left':arm_a,'right':arm_b,'compared_ranks':na,
            'left_streams':rows_a,'right_streams':rows_b}

def verify_campaign_qualification(results,kit=None,placements=None):
    """Public campaign gate. Raises AuditError at first unclassified failure.

    Streams are read only. SQLite is a disposable semantic-route cache created
    outside input evidence, never a replacement for original evidence.
    """
    results=Path(results)
    source_binding=bind_schema_sources(Path(__file__).resolve().parent if kit is None else kit)
    report={'schema':'stage322-semantic-gate-v1','passed':False,'performance_qualified':False,
            'pair_claim':'exact configured numerical state/program/route identity with closed world-size normalization',
            'checks':[],'source_binding':source_binding}
    cache_dir=tempfile.TemporaryDirectory(prefix='stage322-semantic-cache-')
    cache=Path(cache_dir.name)/'routes.sqlite'
    report['route_cache_policy']='temporary outside input evidence; removed after verification'
    c=_cache_open(cache)
    try:
        report['checks'].append(inspect_one_run(results/'checked-B',256,'checked-B',c,'checked'))
        report['checks'].append(compare_runs(results/'qual-A0',results/'qual-B',160,256,'A0','B',c))
        report['checks'].append(compare_runs(results/'qual-C0',results/'qual-Cs',256,256,'C0','Cs',c))
        proofs={};counts={}
        for arm,n,retain in [('checked-B',256,160),('A0',160,160),('B',256,160),('C0',256,256),('Cs',256,256)]:
            proofs[arm],counts[arm]=_reconcile_cache(c,arm,n,retain)
        require(proofs['A0']==proofs['B'],'identity.crossrank_provider_AB')
        require(proofs['C0']==proofs['Cs'],'identity.crossrank_provider_CCs')
        report['crossrank_reconstructed_families']=counts
        report['route_identity_digests']={a:hashlib.sha256(json.dumps(sorted(p.items())).encode()).hexdigest() for a,p in proofs.items()}
        report['passed']=True
        return report
    finally:
        c.close();cache_dir.cleanup()


def _array(q,name):
    require(name in q,'static.missing_field',name)
    return q[name]


def _static_counts(q,n,prefix_name='plan%',names=('sc','sd','rc','rd')):
    sc,sd,rc,rd=[_array(q,prefix_name+k) for k in names]
    require(len(sc)==len(rc)==n and sd==prefix(sc) and rd==prefix(rc),'static.prefix_sum')
    return sc,sd,rc,rd


def reconcile_static(plans):
    """Actual wire correspondence from final committed tables, not value equality.

    Ghost request manifests reverse direction for field payloads. Inverse
    transfer sc describes workspace requests, rc the requested block providers.
    Alias and writeback sc are ordinary producer-to-consumer traffic.
    """
    n=len(plans);scope=plans[0]['__scope'];out=[[] for _ in plans]
    require(all(q['__scope']==scope for q in plans),'static.scope_family')
    if scope=='ghost.plan':
        pre='ghost_exchange_plan%'
        counts=[_static_counts(q,n,pre,('send_record_count','send_record_displ','recv_record_count','recv_record_displ')) for q in plans]
        for consumer,q in enumerate(plans):
            sc,sd,rc,rd=counts[consumer]
            arrays=[q[pre+field] for field in ('source_block','source_local_patch','destination_block','destination_ghost','source_patch_level')]
            source_owner=q[pre+'source_owner']
            require(all(len(a)==len(source_owner) for a in arrays),'ghost.request_extent')
            requested=list(zip(*arrays))
            require(len(q[pre+'request_index'])==sum(sc) and len(q[pre+'send_data'])==5*sum(sc) and
                    len(q[pre+'recv_data'])==5*sum(rc),'ghost.wire_extent')
            require(sorted(q[pre+'request_index'])==[i+1 for i,p in enumerate(source_owner) if p!=consumer],
                    'ghost.request_index_coverage')
            for i,key in enumerate(requested):
                provider=source_owner[i];require(0<=provider<n,'ghost.provider_rank')
                if provider==consumer:out[consumer].append((i,provider,tuple(key),'local'))
            for provider,p in enumerate(plans):
                psc,psd,prc,prd=counts[provider]
                require(sc[provider]==prc[consumer],'ghost.wire_count')
                for i in range(sc[provider]):
                    sent=q[pre+'send_data'][5*(sd[provider]+i):5*(sd[provider]+i+1)]
                    received=p[pre+'recv_data'][5*(prd[consumer]+i):5*(prd[consumer]+i+1)]
                    require(len(sent)==5 and sent==received,'ghost.actual_provider_manifest')
                    request=q[pre+'request_index'][sd[provider]+i]-1
                    require(0<=request<len(requested) and tuple(sent)==tuple(requested[request]) and source_owner[request]==provider,
                            'ghost.receiver_request_binding')
                    out[consumer].append((request,provider,tuple(received),'wire',i))
                for field in ('scalar','vector'):
                    require(q[pre+field+'_recv_count'][provider]==p[pre+field+'_send_count'][consumer],
                            'ghost.field_wire_count')
            # These are wire/request rows. Local execution order remains independently
            # protected by exact local_order297/local_start297 table comparison.
        return out
    if scope=='writeback.plan':
        pre='block_writeback_plan%';counts=[_static_counts(q,n,pre,('send_count','send_displ','recv_count','recv_displ')) for q in plans]
        for receiver,q in enumerate(plans):
            sc,sd,rc,rd=counts[receiver]
            for provider,p in enumerate(plans):
                psc,psd,prc,prd=counts[provider]
                require(rc[provider]==psc[receiver],'writeback.wire_count')
                for i in range(rc[provider]):
                    send_index=psd[receiver]+i;recv_index=rd[provider]+i
                    record=[]
                    for field in ('block','patch_count','scalar_nvalue','vector_nvalue'):
                        src=p[pre+'send_'+field][send_index];dst=q[pre+'recv_'+field][recv_index]
                        require(src==dst,'writeback.actual_provider_manifest',field);record.append(src)
                    out[receiver].append((provider,i,tuple(record)))
                for field in ('scalar','vector'):
                    require(q[pre+field+'_recv_count'][provider]==p[pre+field+'_send_count'][receiver],
                            'writeback.field_wire_count')
        return out
    if scope.startswith('inverse.'):
        # Unallocated level routes are an explicit all-rank state.
        if plans[0].get('plan%sc.shape')==[-1]:
            require(all(q.get('plan%sc.shape')==[-1] for q in plans),'inverse.collective_allocation');return out
        counts=[_static_counts(q,n) for q in plans]
        alias=scope in ('inverse.alias','inverse.range_alias')
        for receiver,q in enumerate(plans):
            sc,sd,rc,rd=counts[receiver]
            for provider,p in enumerate(plans):
                psc,psd,prc,prd=counts[provider]
                if alias:
                    require(rc[provider]==psc[receiver],'inverse.alias_wire_count')
                    require(len(q['plan%dest'])==sum(rc) and len(p['plan%source'])==sum(psc),'inverse.alias_extent')
                    for i in range(rc[provider]):
                        out[receiver].append((q['plan%dest'][rd[provider]+i],provider,p['plan%source'][psd[receiver]+i],i))
                else:
                    # Receiver here is the requesting workspace. The producer is
                    # the block owner whose address is read in gather, written in scatter.
                    require(sc[provider]==prc[receiver],'inverse.transfer_wire_count')
                    require(len(q['plan%slot'])==sum(sc),'inverse.transfer_slot_extent')
                    shape=p['plan%key.shape'];require(shape==[4,sum(prc)],'inverse.transfer_key_extent')
                    address_shape=p['plan%address.shape'];require(len(address_shape)==2 and address_shape[1]==sum(prc),
                                                                'inverse.transfer_address_extent')
                    width=address_shape[0]
                    for i in range(sc[provider]):
                        request_index=sd[provider]+i;producer_index=prd[receiver]+i
                        key=p['plan%key'][4*producer_index:4*(producer_index+1)]
                        address=p['plan%address'][width*producer_index:width*(producer_index+1)]
                        require(len(key)==4 and len(address)==width,'inverse.transfer_provider_extent')
                        if q.get('plan%source_info.shape')!=[-1] and p.get('plan%target_info.shape')!=[-1]:
                            require(q['plan%source_info'][2*request_index:2*request_index+2]==p['plan%target_info'][2*producer_index:2*producer_index+2],
                                    'inverse.transfer_source_metadata')
                        out[receiver].append((q['plan%slot'][request_index],provider,tuple(key),tuple(address),i))
            if alias and q.get('plan%install_order.shape',[-1])!=[-1]:
                order=q['plan%install_order']
                require(sorted(order)==list(range(1,sum(rc)+1)),'inverse.alias_install_permutation')
                out[receiver]=[out[receiver][i-1] for i in order]
        return out
    raise AuditError('static.unknown_family: '+scope)

def bind_schema_sources(kit):
    kit=Path(kit)
    schema=json.loads((kit/'schema322.json').read_text())
    manifest=json.loads((kit/'source-files.json').read_text())
    candidate=manifest['candidate']
    for relative,expected in schema['sources'].items():
        require(candidate.get(relative)==expected,'schema.source_identity',relative)
    require(schema['configured']==dict(level_start=5,level_end=5,layers=30,steps=4,rhs_per_step=4,halo_events_per_rhs=630),
            'schema.configured_contract')
    return {'passed':True,'source_files':len(schema['sources']),
            'schema_sha256':hashlib.sha256((kit/'schema322.json').read_bytes()).hexdigest()}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('results');p.add_argument('--output',required=True)
    args=p.parse_args()
    try:out=verify_campaign_qualification(args.results)
    except (ValueError,OSError,EOFError,struct.error,sqlite3.Error) as exc:
        out={'schema':'stage322-semantic-gate-v1','passed':False,'reason':str(exc),'performance_qualified':False}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print('PASS semantic qualification' if out['passed'] else 'FAIL semantic qualification: '+out['reason'])
    raise SystemExit(0 if out['passed'] else 1)
