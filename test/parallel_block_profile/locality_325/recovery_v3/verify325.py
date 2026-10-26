"""Stage325 bounded, physical-identity qualification. No performance verdict.

The inherited v4 schema checks every original record and lifecycle. Values are
finite-checked and SHA256 hashed at their source; integers remain exact. Source
bound additional records separate moving numerical rows from fixed Domain rows.
"""
from collections import Counter, defaultdict
from pathlib import Path
import ast, copy, gzip, hashlib, json, math, sqlite3, struct, tempfile, time
import verify322 as v
import coverage322 as cov
from reconcile_boundary325 import reconcile_boundary325

PERM=cov.PERMUTATION325
INV=[PERM.index(i) for i in range(256)]
N=256
EXTRA={
 'backend.input':{'input.provider_key325.shape','input.provider_key325','input.physical_key325.shape','input.physical_key325'},
 'backend.rk_input':{'input.provider_key325.shape','input.provider_key325','input.physical_key325.shape','input.physical_key325'},
 'transform.plan':{'publication.key325.shape','publication.key325'},
 'rhs.rk_values322':{'rk.adapter.numeric.shape','rk.adapter.numeric','rk.adapter.domain.shape','rk.adapter.domain'},
 'transform.values322':{'transform.numeric.scalar.shape','transform.numeric.scalar','transform.numeric.vector.shape','transform.numeric.vector','transform.numeric.saved_scalar.shape','transform.numeric.saved_scalar','transform.numeric.saved_vector.shape','transform.numeric.saved_vector'},
}
LIMIT=128*1024*1024
COMPRESSED_LIMIT=64*1024*1024
SCRATCH_LIMIT=4*1024*1024*1024
# Source-derived physical routes for the sealed level_start=level_end=5 case.
# parallel_block_rhs_routes.f90 selects spatial components; independently,
# parallel_block_rhs_distributed.f90 selects the channel capacity. The logical
# injected-Bernoulli phase 18 reuses physical route 16 and has no own route.
RHS_PHASE325={
    1:(3,2),2:(1,2),3:(3,1),4:(2,1),5:(3,2),6:(2,1),7:(1,1),
    11:(2,1),12:(3,1),13:(3,2),14:(1,2),15:(3,1),16:(1,2),20:(3,2),
}

def capture_environment(kind,owner):
    return {'WAVETRISK_EXECUTION_AUDIT316':'1','WAVETRISK_COMPACT_AUDIT325':'1'}

def is_locality(arm):return arm in ('L256','checked-L256','qual-L256')
def group(rank,arm):
    require(0<=rank<N,'physical.numerical_rank_range',{'rank':rank,'arm':arm})
    return INV[rank] if is_locality(arm) else rank

def require(ok,code,detail=''):
    v.require(ok,code,detail)

def integers(s,name):return v.value(s['values'][name])
def cols(values,width):
    require(width>0 and len(values)%width==0,'physical.column_extent')
    return [tuple(values[i:i+width]) for i in range(0,len(values),width)]
def matrix(s,name,width=None):
    shape=integers(s,name+'.shape')
    require(len(shape)==2 and all(x>=0 for x in shape),'physical.matrix_shape',name)
    if width is not None:require(shape[0]==width,'physical.matrix_width',name)
    data=integers(s,name);require(len(data)==math.prod(shape),'physical.matrix_payload',name)
    return cols(data,shape[0]) if shape[0] else []
def descriptor(e):
    return [e.name,e.kind,e.count,e.data.hex() if e.kind==3 else v.value(e)]
def sha(obj):return hashlib.sha256(json.dumps(obj,separators=(',',':'),sort_keys=True).encode()).hexdigest()
def physical(key,component=1):return ('physical',key[0],key[1],component)
def physical_rows(keys,m=1):return [physical(key,c) for key in keys for c in range(1,m+1)]

class Reader325:
    def __init__(self,path,rank=None):self.path=Path(path);self.expected_rank=rank;self.summary=None
    def __iter__(self):
        require(self.path.stat().st_size<=COMPRESSED_LIMIT,'compact.compressed_size')
        digest=hashlib.sha256();logical=0
        with v.open_stream(self.path) as f:
            raw=v.exact_read(f,v.HEADER.size,'header325')
            magic,version,rank,isz,rsz,endian=v.HEADER.unpack(raw)
            require((magic,version,isz,rsz,endian)==(b'RK325V01',1,4,8,0x01020304),'compact.header')
            require(rank>=0 and (self.expected_rank is None or rank==self.expected_rank),'compact.rank')
            digest.update(raw);offset=len(raw);index=0;done=False;real_records=0
            while True:
                head=f.read(v.RECORD.size)
                if not head:break
                require(not done and len(head)==v.RECORD.size,'compact.truncated_or_after_footer')
                length,kind,count=v.RECORD.unpack(head)
                require(1<=length<=256 and kind in (1,2,4,10,11,12,13) and count>=0,'compact.descriptor')
                size=32 if kind==13 else count*v.WIDTH[kind]
                require(offset+16+length+size<=LIMIT and count<=2**30,'compact.size_bound')
                tag=v.exact_read(f,length,'tag325');name=tag.decode('ascii')
                require(all(32<=ord(c)<=126 for c in name),'compact.tag')
                data=v.exact_read(f,size,'payload325')
                if kind in (10,11):require(count==1,'compact.scope_epoch')
                if kind==4:require(all(x in (0,1) for x in data),'compact.logical')
                if kind==13:
                    real_records+=1
                    if count==0:require(data==hashlib.sha256(b'').digest(),'compact.empty_digest')
                if kind==12:
                    require(name=='complete' and count==2 and struct.unpack('<2q',data)==(index,offset),'compact.footer')
                    done=True
                digest.update(head);digest.update(tag);digest.update(data)
                yield v.Entry(index,offset,name,3 if kind==13 else kind,count,data)
                logical+=8*count if kind==13 else size
                index+=1;offset+=16+length+size
            require(done and real_records>0,'compact.complete_real_capture')
            self.summary={'rank':rank,'records':index,'physical_bytes':offset,'logical_payload_bytes':logical,
                          'real_digest_records':real_records,'stream_sha256':digest.hexdigest(),'file':self.path.name}

ORIGINAL_COVERAGE=v.validate_coverage322

def validate_coverage(scopes,rank,n,steps=4,coverage='candidate',arm=None):
    stripped=[]
    for s in scopes:
        allowed=EXTRA.get(s['name'],set())
        extras=[e for e in s['records'] if e.name in allowed]
        require({e.name for e in extras}==allowed and all(sum(e.name==name for e in extras)==1 for name in allowed),
                'compact.extra_coverage',s['name'])
        t=dict(s)
        for field in ('records','tags','values','record_counts','record_descriptors'):
            t[field]=([e for e in s[field] if e.name not in allowed] if field=='records' else
                      {x for x in s[field] if x not in allowed} if field=='tags' else
                      {k:x for k,x in s[field].items() if k not in allowed})
        stripped.append(t)
        if s['name'] in ('backend.input','backend.rk_input'):
            np,nd=integers(s,'input.extent')
            require(len(matrix(s,'input.provider_key325',4))==np and len(matrix(s,'input.physical_key325',4))==nd,
                    'compact.input_key_extents')
        elif s['name']=='transform.plan':
            require(len(matrix(s,'publication.key325',4))==integers(s,'extent')[3],'compact.publication_key_extent')
        elif s['name']=='rhs.rk_values322':
            whole=integers(s,'rk.adapter.shape');a=integers(s,'rk.adapter.numeric.shape');b=integers(s,'rk.adapter.domain.shape')
            require(len(a)==len(b)==2 and a[0]+b[0]==whole[0] and a[1]==b[1]==whole[1] and b==integers(s,'rk.value.shape'),
                    'compact.rk_split')
        elif s['name']=='transform.values322':
            nd=integers(s,'transform.phase')[1]
            require(integers(s,'transform.numeric.scalar.shape')==[nd,60] and
                    integers(s,'transform.numeric.vector.shape')==[3*nd,30],'compact.transform_split')
            plan=prior(scopes,s,'transform.plan');saved=len(integers(plan,'save_node'))+len(integers(plan,'fixed_node'))
            require(integers(s,'transform.numeric.saved_scalar.shape')==[saved,60] and integers(s,'transform.numeric.saved_vector.shape')==[3*saved,30],'compact.transform_saved_split')
        for e in extras:
            if e.name.endswith('.shape'):require(e.kind==1,'compact.extra_shape_kind')
            elif e.name.startswith(('input.','publication.')):
                require(e.kind==1,'compact.extra_key_kind')
            else:
                require(e.kind==3,'compact.extra_real_kind')
                shape=integers(s,e.name+'.shape');require(all(x>=0 for x in shape) and e.count==math.prod(shape),'compact.extra_real_extent')
    return ORIGINAL_COVERAGE(stripped,rank,n,steps,coverage,arm)

v.Reader=Reader325
v.validate_coverage322=validate_coverage

def scope_key(s):
    return (tuple(tuple(x) for x in s['parent']),s['name'],s['epoch'])
def ancestor(scopes,s,name):
    found=[t for t in scopes if t['name']==name and t['record']<=s['record'] and t['end_record']>=s['end_record']]
    require(len(found)==1,'physical.ancestor',name)
    return found[0]
def prior(scopes,s,name):
    found=[t for t in scopes if t['name']==name and t['record']<=s['record']]
    require(bool(found),'physical.prior_plan',name)
    return found[-1]
def rhs_route_layout(q,nd,rank,coordinate=None):
    """Validate one compiled route against the closed live source contract."""
    require(q.get('route.allocated')==[1],'physical.rhs_allocated',{'rank':rank,'coordinate':coordinate})
    h=q.get('route.header',[])
    require(len(h)==7,'physical.rhs_header',{'rank':rank,'header':h})
    encoded,nslot,ncolumn,me,nrank,tag,last_columns=h
    level,phase=divmod(encoded,100)
    detail={'rank':rank,'phase':phase,'level':level,'coordinate':coordinate,
            'nodes':nd,'slots':nslot,'components':nslot/nd if nd else None,
            'columns':ncolumn,'header':h}
    require(level==5 and phase in RHS_PHASE325,'physical.rhs_phase',detail)
    require(coordinate is None or coordinate==[phase,level],'physical.rhs_coordinate_header',detail)
    m,nc=RHS_PHASE325[phase]
    detail['expected_components']=m;detail['expected_columns']=nc
    require(nd>0 and nslot==m*nd,'physical.rhs_components',detail)
    require(ncolumn==last_columns==nc,'physical.rhs_columns',detail)
    require([me,nrank]==[rank,N] and tag==3*encoded+nc,'physical.rhs_header_contract',detail)
    require(q.get('route.flags')==[0,1],'physical.rhs_shared_variable_route',detail)
    return m

def validate_rhs_routes(s,rank):
    """Close the complete physical census and logical-to-physical bindings.

    Do not infer components from the observed slot count: two-component triangle
    routes are required, and a wrong count common to both arms must still fail.
    The coordinate preceding each actual route is checked against its header.
    """
    detail={'rank':rank,'epoch':s['epoch'],'step':s['step']}
    require(integers(s,'transport.route.bounds')==[5,5],'physical.rhs_level_bounds',detail)
    require(integers(s,'transport.bound')==[int(p in RHS_PHASE325) for p in range(1,21)],
            'physical.rhs_bound_census',detail)
    logical={**{p:p for p in RHS_PHASE325},18:16}
    require(integers(s,'transport.route_phase')==[logical.get(p,0) for p in range(1,21)],
            'physical.rhs_logical_phase_binding',detail)
    require(integers(s,'transport.route_level')==[5 if p in logical else 0 for p in range(1,21)],
            'physical.rhs_logical_level_binding',detail)
    nd=len(matrix(s,'placement.key',4));routes=[];coordinate=None;q=None
    for e in s['records']:
        if e.name=='transport.route.coordinate':
            if q is not None:routes.append((coordinate,q));q=None;coordinate=None
            require(coordinate is None,'physical.rhs_coordinate_without_route',detail)
            require(e.kind==1 and e.count==2,'physical.rhs_coordinate_format',detail)
            coordinate=v.value(e)
        elif e.name=='route.allocated':
            require(coordinate is not None and q is None,'physical.rhs_route_without_coordinate',detail)
            q={e.name:v.value(e)}
        elif e.name.startswith('route.'):
            require(q is not None and e.name not in q,'physical.rhs_route_record_order',detail)
            q[e.name]=v.value(e)
    if q is not None:routes.append((coordinate,q));coordinate=None
    require(coordinate is None,'physical.rhs_dangling_coordinate',detail)
    require([a for a,b in routes]==[[p,5] for p in RHS_PHASE325],
            'physical.rhs_physical_census',detail)
    for coordinate,q in routes:rhs_route_layout(q,nd,rank,coordinate)

def record_scope_routes(scopes,rank):
    seen=Counter();out={}
    for s in scopes:
        if s['name']=='rhs.prepare':validate_rhs_routes(s,rank)
        path=tuple(s['parent'])+((s['name'],s['epoch']),);seen[path]+=1
        path=repr((path,seen[path]));serial=0
        for e in s['records']:
            if e.name=='route.allocated':out[path+f'/route/{serial}']=s;serial+=1
    return out

def block_storage(scopes,s):
    installed=prior(scopes,s,'storage.installed')
    blocks=[q for q in scopes if q['name']=='block' and installed['record']<q['record']<installed['end_record']]
    return blocks

def address_key_fields(address,blocks):
    require(len(address)==6,'physical.block_address_width')
    local,storage,q,n,offset,layers=address
    require(1<=local<=len(blocks),'physical.block_address_local')
    block=blocks[local-1];catalog=block['catalog'];fields=block['fields']
    ident=fields['identity']
    require((offset,layers)==(1-ident[7],ident[8]),'physical.block_address_layers')
    if storage==1:
        patches=cols(fields['patch'],fields['patch.shape'][0])
        matches=[(i+1,q-row[1]) for i,row in enumerate(patches) if row[0]==0 and 0<=q-row[1]<16]
        require(n==fields['node.shape'][1] and len(matches)==1,'physical.patch_address')
    elif storage==2:
        entries=cols(fields['bdry_storage'],7)
        matches=[(i+1,q-row[6]) for i,row in enumerate(entries) if 0<=q-row[6]<row[5]]
        require(n==fields['bdry_node.shape'][1] and len(matches)==1,'physical.boundary_address')
    else:raise v.AuditError('physical.unsupported_block_storage')
    return (catalog,storage,*matches[0])

BLOCK_ADDRESS_FIELDS={"identity","patch.shape","patch","node.shape","bdry_storage","bdry_node.shape"}
def compact_blocks(blocks):
    return [{"catalog":b["epoch"],"fields":{e.name:v.value(e) for e in b["records"] if e.name in BLOCK_ADDRESS_FIELDS}} for b in blocks]
def address_key(address,blocks):
    return address_key_fields(address,compact_blocks(blocks))

def validate_input_addresses(s,scopes,rank):
    providers=matrix(s,'input.provider_key325',4);local=integers(s,'input.local_provider')
    addresses=matrix(s,'input.address',6);blocks=block_storage(scopes,s)
    require(len(addresses)==len(local),'physical.input_address_count')
    for p,a in zip(local,addresses):
        require(1<=p<=len(providers),'physical.input_provider_slot')
        key=address_key(a,blocks);require(tuple(providers[p-1])==(*key[:3],key[3]+1),'physical.input_provider_address')
    ap=integers(s,'input.adapter_provider');ad=integers(s,'input.adapter_domain');an=integers(s,'input.adapter_node')
    require(len(ap)==len(ad)==len(an),'physical.input_adapter_count')
    for p,d,node in zip(ap,ad,an):
        require(rank<160 and d==1 and tuple(providers[p-1])==(641+rank,4,node,1),'physical.input_fixed_provider')

def route_map(q,s,scopes,rank,arm):
    """Actual slot identities, including mixed fixed Domain/moving-owner plans."""
    if s['name']=='rhs.prepare':
        keys=matrix(s,'placement.key',4)
        m=rhs_route_layout(q,len(keys),rank)
        return {'rows':physical_rows(keys,m),'role':'numerical','group':group(rank,arm)}
    if q['route.allocated']==[0]:return {'rows':[],'role':'none'}
    phase,nslot=q['route.header'][:2]
    if s['name'] in ('input.scalar_route','input.vector_route','input.mean_route'):
        p=ancestor(scopes,s,'backend.rk_input' if any(x[0]=='backend.rk_input' for x in s['parent']) else 'backend.input')
        validate_input_addresses(p,scopes,rank)
        providers=matrix(p,'input.provider_key325',4);targets=matrix(p,'input.physical_key325',4)
        m=3 if s['name']=='input.vector_route' else 1
        rows=[('state',*key,c) for key in providers for c in range(1,m+1)]+physical_rows(targets,m)
        require(len(rows)==nslot,'physical.input_slots')
        active=set(integers(p,'input.local_provider')+integers(p,'input.adapter_provider'))
        require(all((slot-1)//m+1 in active for slot in q['route.source']),'physical.input_unfilled_source')
        for slot in q['route.source']:
            source_key=providers[(slot-1)//m]
            expected=cov.owner325(source_key[0]-1,arm) if source_key[0]<=640 else source_key[0]-641
            require(rank==expected,'physical.input_source_owner')
        return {'rows':rows,'role':'fixed' if p['name']=='backend.rk_input' else 'numerical','group':rank if p['name']=='backend.rk_input' else group(rank,arm)}
    if s['name']=='backend.adapter_route':
        p=ancestor(scopes,s,'rhs.prepare');keys=matrix(p,'placement.key',4)
        domain=integers(p,'backend.adapter_domain');nodes=integers(p,'backend.adapter_node')
        require(all(d==1 and rank<160 for d in domain),'physical.adapter_domain')
        rows=physical_rows(keys)+[('physical',rank+1,node,1) for node in nodes]
        require(len(rows)==nslot,'physical.adapter_slots')
        owned=integers(p,'placement.owned')
        require(len(owned)==len(keys) and all(1<=slot<=len(keys) and owned[slot-1]==1 for slot in q['route.source']),
                'physical.adapter_numerical_source',{'rank':rank,'numerical_rows':len(keys),'sources':q['route.source'][:8]})
        return {'rows':rows,'role':'fixed','group':rank}
    if s['name'].startswith('transform.'):
        p=ancestor(scopes,s,'transform.plan');keys=matrix(p,'placement%key',4)
        nd,nrow,nfixed=integers(p,'extent')[:3]
        m=3 if ('vector' in s['name'] or 'outer' in s['name']) else 1
        if s['name'] in ('transform.fixed_scalar_route','transform.fixed_vector_route'):
            domains=integers(p,'fixed_domain');nodes=integers(p,'fixed_physical_node')
            require(all(d==1 and rank<160 for d in domains),'physical.transform_fixed_domain')
            rows=[('physical',rank+1,node,c) for node in nodes for c in range(1,m+1)]
            rows+=physical_rows([keys[i-1] for i in integers(p,'fixed_node')],m)
            require(len(rows)==nslot,'physical.transform_fixed_slots')
            return {'rows':rows,'role':'numerical','group':group(rank,arm),'fixed_prefix':m*nfixed}
        rows=physical_rows(keys,m)
        if nslot>len(rows):
            require(nslot==m*nrow,'physical.transform_nrow')
            extra=[None]*(nrow-nd)
            for fields in (('save_node','save_index'),('fixed_node','fixed_index')):
                indexes=integers(p,fields[1])
                for node in integers(p,fields[0]):
                    j=indexes[node-1];require(nd<j<=nrow,'physical.transform_saved_row')
                    extra[j-nd-1]=keys[node-1]
            pub=matrix(p,'publication.key325',4)
            for row,key in zip(integers(p,'publication_row'),pub):
                if row<=nd:
                    require(1<=row<=nd and key[:2]==keys[row-1][:2] and integers(p,'placement%owned')[row-1]==1,'physical.local_publication_binding')
                else:
                    old=extra[row-nd-1];require(old is None or old[:2]==key[:2],'physical.imported_row_alias')
                    extra[row-nd-1]=key
            require(all(x is not None for x in extra),'physical.unidentified_transform_row')
            rows+=physical_rows(extra,m)
        require(len(rows)==nslot,'physical.transform_slots')

        if 'publication' in s['name']:
            expected=[];seen=set()
            for row in integers(p,'publication_row'):
                if row<=nd or row in seen:continue
                seen.add(row);expected.extend(range(m*(row-1)+1,m*row+1))
            require(q['route.target']==expected and all(x==1 for x in q['route.orientation']),'physical.publication_complete_targets')
        return {'rows':rows,'role':'publication' if 'publication' in s['name'] else 'numerical','group':group(rank,arm),'canonical_rows':m*nd,'owned':integers(p,'placement%owned')}
    raise v.AuditError('physical.unhandled_route_scope: '+s['name'])

# Fields that have independent fixed-Domain or movable numerical identity.
NUMERIC_SCOPES={'rhs.prepared_state','rhs.output322','rhs.evaluate','rhs.retain'}
FIXED_SCOPES={'domain.fields322','physics.column322','cfl322','inverse.values322','inverse.scalar','inverse.inner','inverse.outer','inverse.local_alias','inverse.node_map','inverse.level','inverse.outer_level','inverse.gather','inverse.scatter','inverse.begin','inverse.prepare','rhs.rk_prepare','rhs.rk_install','adapt322'}
NUMERIC_RHS_FIELDS=cov.RHS_PROGRAM_FIELDS|{'placement.key.shape','placement.key','placement.base','placement.owned','backend.surface','backend.owned','executor.scalar_rhs.shape','executor.velocity_rhs.shape','transport.route.bounds','transport.bound','transport.route_phase','transport.route_level','backend.scalar.shape','backend.mean.shape','backend.velocity.shape'}
NUMERIC_TRANSFORM_FIELDS={'masks','placement%key.shape','placement%key','placement%state_key.shape','placement%state_key','placement%level.shape','placement%level','placement%mask.shape','placement%mask','placement%owned.shape','placement%owned','placement%scalar_order.shape','placement%scalar_order','placement%inner_order.shape','placement%inner_order','placement%outer_order.shape','placement%outer_order','info.shape','info','save_node.shape','save_node','save_index.shape','save_index','fixed_node.shape','fixed_node','fixed_index.shape','fixed_index','interior_node.shape','interior_node','forward_inner.shape','forward_inner','scalar.program.shape','scalar.program','scalar.coefficients','inner.program.shape','inner.program','inner.coefficients','outer.program.shape','outer.program','outer.coefficients','scalar_wavelet.shape','vector_wavelet.shape'}

def normalize_block(s,arm):
    ans=[]
    for e in s['records']:
        item=descriptor(e)
        if e.name in ('block_bdry','ghost_storage'):
            width,col=(17,15) if e.name=='block_bdry' else (10,8)
            rows=[list(x) for x in cols(v.value(e),width)]
            for row in rows:
                if 0<=row[col]<N:row[col]=group(row[col],arm)
            item[3]=[x for row in rows for x in row]
        ans.append(item)
    return ans

def transform_owner_identity(s,rank,arm):
    """Normalize canonical owners while preserving dependency-row sentinels.

    distribute_transform_geometry initializes every footprint row to the full
    sentinel tuple. Only a tag1 Node_Record changes that tuple, and it is sent
    exclusively to node%owner. A nonowned operator dependency is not a numerical
    rank and must never enter the rank permutation (Python's -1 indexes a real
    element). Its actual providers are resolved separately from wire plans.
    """
    owners=integers(s,'placement%owner');owned=integers(s,'placement%owned')
    physical_keys=matrix(s,'placement%key',4);state_keys=matrix(s,'placement%state_key',4)
    levels=integers(s,'placement%level');masks=matrix(s,'placement%mask',4)
    nd=len(physical_keys)
    require(all(len(x)==nd for x in (owners,owned,state_keys,levels,masks)),
            'physical.transform_placement_extent',{'rank':rank,'nodes':nd})
    normalized=[]
    for slot,(key,flag,owner,state,level,mask) in enumerate(zip(physical_keys,owned,owners,state_keys,levels,masks),1):
        detail={'rank':rank,'epoch':s['epoch'],'step':s['step'],'slot':slot,'physical_key':key,
                'owned':flag,'owner':owner,'state_key':state,'level':level,'mask':mask}
        require(flag in (0,1),'physical.transform_owned_flag',detail)
        if flag:
            require(owner==rank and 0<=owner<N,'physical.transform_canonical_owner',detail)
            normalized.append(group(owner,arm))
        else:
            require(owner==-1 and state==(0,0,0,0) and level==-1 and mask==(0,0,0,0),
                    'physical.transform_dependency_sentinels',detail)
            normalized.append(-1)
    return normalized

def feature_rows(scopes,rank,arm):
    """Ordered within physical program/consumer, not blanket sorted records."""
    out=[];seen=Counter()
    def emit(role,key,payload):
        full=(role,key);seen[str(full)]+=1
        out.append((repr((full,seen[str(full)])),sha(payload)))
    for s in scopes:
        path=scope_key(s);name=s['name'];step=s['step']
        if name in NUMERIC_SCOPES:
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records']])
        elif name in FIXED_SCOPES:
            emit('fixed',(rank,path),[descriptor(e) for e in s['records']])
        elif name=='rhs.prepare':
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records'] if e.name in NUMERIC_RHS_FIELDS or e.name.startswith('workspace.')])
            # Retention inventory is per physical block. Preserve each block's
            # node/column order; only normalize the storage index.
            groups=[];current=None
            for e in s['records']:
                if e.name=='retention.block':
                    current=[v.value(e)[1],[]];groups.append(current)
                if e.name.startswith('retention.') and current is not None and e.name!='retention.block':current[1].append(descriptor(e))
            for cat,records in groups:emit('block-retention-plan',(cat,step),records)
        elif name=='transform.plan':
            emit('numeric-owners',(group(rank,arm),path),transform_owner_identity(s,rank,arm))
            keys=matrix(s,'placement%state_key',4);blocks=block_storage(scopes,s)
            addresses=matrix(s,'interior_address',6);nodes=integers(s,'interior_node')
            require(len(addresses)==len(nodes),'physical.transform_interior_addresses')
            for address,node in zip(addresses,nodes):
                require(address_key(address,blocks)==keys[node-1],'physical.transform_interior_binding')
            publications=matrix(s,'publication.key325',4);domains=integers(s,'publication_domain');nodes=integers(s,'publication_node')
            require(len(domains)==len(nodes),'physical.publication_domain_count')
            for d,node,key in zip(domains,nodes,publications):
                require(rank<160 and d==1 and key[:2]==(rank+1,node),'physical.publication_domain_key')
                emit('publication-demand',(step,s['epoch'],'domain',rank,node),key)
            boundary=matrix(s,'boundary_address',6)
            require(len(publications)==len(domains)+len(boundary),'physical.publication_boundary_count')
            for address,key in zip(boundary,publications[len(domains):]):
                emit('publication-demand',(step,s['epoch'],'block',address_key(address,blocks)),key)
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records'] if e.name in NUMERIC_TRANSFORM_FIELDS])
        elif name=='transform.level':
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records']])
        elif name=='transform.values322':
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records'] if e.name.startswith('transform.numeric.') or e.name.startswith('transform.scalar_wavelet') or e.name.startswith('transform.vector_wavelet') or e.name=='transform.state'])
        elif name=='rhs.rk_values322':
            emit('numerical',(group(rank,arm),path),[descriptor(e) for e in s['records'] if e.name.startswith('rk.adapter.numeric')])
            emit('fixed',(rank,path),[descriptor(e) for e in s['records'] if e.name.startswith('rk.adapter.domain') or e.name.startswith(('rk.value','rk.scalar','rk.mean','rk.velocity'))])
        elif name=='rhs.rk_installed322':
            emit('fixed',(rank,path),[descriptor(e) for e in s['records']])
        elif name=='block':emit('block',(s['epoch'],s['parent']),normalize_block(s,arm))
        elif name=='rhs.retained322':
            ident=integers(s,'retained.block')
            emit('retained',(ident[0],path),[ident[0],ident[2:],*[descriptor(e) for e in s['records'] if e.name!='retained.block']])
        elif name=='inverse.boundary' and 'family_and_range' in s['values']:
            emit('fixed',(rank,path),[descriptor(e) for e in s['records']])
        elif name=='inverse.plan':
            exclude={'transfer_send.shape','transfer_recv.shape'}
            emit('fixed',(rank,path),[descriptor(e) for e in s['records'] if not e.name.startswith('plan%') and e.name not in exclude])
    require(bool(out),'physical.no_features')
    return out


def _create_cache(path):
    c=v._cache_open(path)
    page_size=c.execute('pragma page_size').fetchone()[0]
    c.execute('pragma max_page_count='+str(SCRATCH_LIMIT//page_size))
    c.execute('create table context(arm text,rank int,kind text,key text,payload text,primary key(arm,rank,kind,key))')
    c.execute('create table feature(arm text,key text,digest text,primary key(arm,key))')
    return c

def expected_storage_order(rank,arm):
    owned=[i+1 for i in range(640) if cov.owner325(i,arm)==rank]
    return [i for i in owned if (i-1)//4==rank]+[i for i in owned if (i-1)//4!=rank]

def validate_storage_orders(scopes,rank,arm):
    stores=[s for s in scopes if s['name']=='storage.installed']
    require(len(stores)==5,'physical.store_frame_count')
    orders=[]
    for store in stores:
        order=[s['epoch'] for s in scopes if s['name']=='block' and store['record']<s['record']<store['end_record']]
        require(order==expected_storage_order(rank,arm),'physical.actual_store_order',rank)
        orders.append(order)
    return orders

def scratch_usage(directory,stats):
    size=sum(p.stat().st_size for p in Path(directory).rglob("*") if p.is_file())
    stats["sampled_peak_scratch_bytes"]=max(stats.get("sampled_peak_scratch_bytes",0),size)
    require(size<=SCRATCH_LIMIT,"compact.scratch_bound",size)
    return size

def cache_rank325(c,path,arm,rank):
    """Run the exact per-rank native gates, also usable for saved-capture preflight."""
    coverage='checked' if arm.startswith('checked') else 'candidate'
    obj=v.parse_stream(path,rank,N,arm,coverage);scopes=obj['scopes']
    orders=validate_storage_orders(scopes,rank,arm)
    obj['reader']['storage_orders']=orders
    v._cache_rank(c,arm,rank,obj)
    lookup=record_scope_routes(scopes,rank)
    for key,q in obj['routes'].items():
        context=route_map(q,lookup[key],scopes,rank,arm)
        c.execute('insert into context values(?,?,?,?,?)',(arm,rank,'routes',key,json.dumps(context,separators=(',',':'))))
    for kind in ('adapters','static'):
        for key,q in obj[kind].items():
            skey,occurrence=ast.literal_eval(key)
            matches=[s for s in scopes if tuple(s['parent'])+((s['name'],s['epoch']),)==skey]
            require(len(matches)>=occurrence,'physical.static_scope_binding');s=matches[occurrence-1]
            if kind=='adapters':
                p=ancestor(scopes,s,'transform.plan');m=q['extent'][0]
                context={'rows':physical_rows(matrix(p,'placement%key',4),m),'group':group(rank,arm),'name':s['name']}
            else:
                if s['name']=='ghost.plan':ghost_local_order325(q,rank)
                elif s['name'] in INVERSE_TRANSFERS325:inverse_transfer_contract325(q)
                blocks=block_storage(scopes,s)
                context={'blocks':compact_blocks(blocks), 'name':s['name']}
            c.execute('insert into context values(?,?,?,?,?)',(arm,rank,kind,key,json.dumps(context,separators=(',',':'))))
    for key,digest in feature_rows(scopes,rank,arm):c.execute('insert into feature values(?,?,?)',(arm,key,digest))
    c.commit()
    return obj['reader']

def cache_run(c,directory,arm,scratch,started,stats,summaries=None,progress=None):
    if summaries is None:summaries=[]
    for rank,path in enumerate(v._stream_set(directory,N)):
        if progress:progress('streams',arm,rank,N,{'rank':rank,'file':path.name},rank%16==0)
        summaries.append(cache_rank325(c,path,arm,rank))
        size=scratch_usage(scratch,stats)
        if progress:progress('streams',arm,rank+1,N,None,(rank+1)%16==0 or rank==N-1)
        if (rank+1)%16==0 or rank==N-1:
            print(f'Stage325 {arm}: validated {rank+1}/{N} ranks; scratch {size/2**20:.1f} MiB; elapsed {time.monotonic()-started:.1f}s',flush=True)
    return summaries


def _contexts(c,arm,kind,key):
    rows=c.execute('select rank,payload from context where arm=? and kind=? and key=? order by rank',(arm,kind,key)).fetchall()
    require([x[0] for x in rows]==list(range(N)),'physical.missing_context')
    return [json.loads(x[1]) for x in rows]

def _provenance(context,rank,slot,arm,source=True):
    rows=context['rows'];require(1<=slot<=len(rows),'physical.source_or_target_slot')
    row=rows[slot-1]
    if source:
        # Catalogue keys identify the authoritative block; uncatalogued adapter
        # state keys identify the fixed Domain. Numerical duplicate providers
        # retain their logical owner so equal values cannot hide wrong writers.
        if row[0]=='state':return tuple(row)
        numerical=slot>context.get('fixed_prefix',0)
        return ('provider',group(rank,arm) if numerical else rank,*row)
    return tuple(row)

INVERSE_TRANSFERS325=frozenset(('inverse.interior','inverse.level_interior','inverse.boundary',
                              'inverse.level_boundary','inverse.outer_install'))

def inverse_transfer_contract325(plan):
    """The sealed sparse-outer=0 configuration uses six-component addresses."""
    require(plan.get('edgewise')==[0],'physical.inverse_edgewise_disabled')
    if plan.get('plan%sc.shape')==[-1]:
        require(plan.get('plan%address.shape')==[-1] and 'plan%address' not in plan,
                'physical.inverse_unallocated_address')
        return
    require('plan%sc' in plan and 'plan%rc' in plan,'physical.inverse_allocated_counts')
    count=sum(plan['plan%rc'])
    require(plan.get('plan%address.shape')==[6,count] and len(plan.get('plan%address',[]))==6*count,
            'physical.inverse_address_shape')
    a=plan.get('plan%source_info.shape');b=plan.get('plan%target_info.shape')
    if a==[-1] or b==[-1]:
        require(a==b==[-1] and 'plan%source_info' not in plan and 'plan%target_info' not in plan,
                'physical.inverse_optional_metadata')
    else:
        require(a==[2,sum(plan['plan%sc'])] and b==[2,count] and
                len(plan.get('plan%source_info',[]))==2*sum(plan['plan%sc']) and
                len(plan.get('plan%target_info',[]))==2*count,'physical.inverse_metadata_extent')

def ghost_local_order325(plan,rank,shared_source=False):
    """Source loop order, with allocation tied to the explicit execution option.

    The sealed Stage325 runs set WAVETRISK_SHARED_GHOST_SOURCE=0. The scalar
    loop then visits request indices in ascending order, skipping remote owners;
    vector execution always does so. Grouped scalar reuse is a distinct mode,
    never inferred from a missing/available array. Shapes must agree with it.
    """
    pre='ghost_exchange_plan%'
    owners=plan[pre+'source_owner']
    expected=[i+1 for i,owner in enumerate(owners) if owner==rank]
    fields=('local_order297','local_start297')
    if not shared_source:
        require(all(plan.get(pre+k+'.shape')==[-1] and pre+k not in plan for k in fields),
                'physical.ghost_ungrouped_allocation',rank)
        return expected
    require(all(pre+k in plan and plan.get(pre+k+'.shape')==[len(plan[pre+k])] for k in fields),
            'physical.ghost_grouped_allocation',rank)
    order=plan[pre+fields[0]];starts=plan[pre+fields[1]]
    require(sorted(order)==expected,'physical.ghost_local_order_coverage',rank)
    require(starts and starts[0]==1 and starts[-1]==len(order)+1 and
            all(a<b for a,b in zip(starts,starts[1:])), 'physical.ghost_local_group_bounds',rank)
    sources=plan[pre+'source_block'];patches=plan[pre+'source_local_patch']
    require(len(sources)==len(patches)==len(owners),'physical.ghost_group_source_extent',rank)
    for a,b in zip(starts,starts[1:]):
        require(len({(sources[i-1],patches[i-1]) for i in order[a-1:b-1]})==1,
                'physical.ghost_group_source',rank)
    # group_local_297 walks global block/patch offsets in increasing order,
    # and retains original request order inside each group.
    expected_order=sorted(expected,key=lambda i:(sources[i-1],patches[i-1]))
    expected_starts=[1]
    for j in range(1,len(expected_order)):
        a,b=expected_order[j-1],expected_order[j]
        if (sources[a-1],patches[a-1])!=(sources[b-1],patches[b-1]):expected_starts.append(j+1)
    if expected_order:expected_starts.append(len(expected_order)+1)
    require(order==expected_order and starts==expected_starts,'physical.ghost_source_group_order',rank)
    return order

def route_certificates(c,arm,progress=None):
    families=c.execute('select kind,key,count(*) from route where arm=? group by kind,key order by kind,key',(arm,)).fetchall()
    result={};total=0;publication=0
    for family_index,(kind,key,num) in enumerate(families,1):
        if progress:progress('providers',arm,family_index-1,len(families),{'kind':kind,'key':key},(family_index-1)%32==0)
        if family_index%32==0:
            print(f'Stage325 {arm}: physical provider family {family_index}/{len(families)}',flush=True)
        require(num==N,'physical.route_family_ranks',key)
        plans=[json.loads(x[0]) for x in c.execute('select payload from route where arm=? and kind=? and key=? order by rank',(arm,kind,key))]
        context=_contexts(c,arm,kind,key)
        if kind=='routes':
            resolved=v.reconcile_routes(plans);rows=[]
            for receiver,entries in enumerate(resolved):
                ctx=context[receiver]
                for target,producer,source,orientation,ordinal in entries:
                    target_key=_provenance(ctx,receiver,target,arm,False)
                    src_ctx=context[producer];source_key=_provenance(src_ctx,producer,source,arm)
                    if ctx['role']=='publication':
                        # Local publication is satisfied by direct canonical
                        # reuse; only imported rows produce a route. Both must
                        # name the exact same physical source key.
                        require(source<=src_ctx['canonical_rows'] and tuple(src_ctx['rows'][source-1])==target_key,'physical.publication_provider')
                        m=3 if 'vector' in key else 1
                        require(src_ctx['owned'][(source-1)//m]==1,'physical.publication_authoritative_owner')
                        publication+=1;continue
                    logical=receiver if ctx['role']=='fixed' else group(receiver,arm)
                    rows.append((logical,target_key,source_key,orientation))
                    total+=1
            # Exchanges snapshot before any destination writes. Preserve the
            # relative sequence for duplicate destinations; independent target
            # locations may arrive in a different peer order.
            result[(kind,key)]=canonical_writes(rows)
        elif kind=='adapters':
            domains=[json.loads(x[0]) for x in c.execute('select payload from domain where arm=? order by rank',(arm,))]
            resolved=v.reconcile_adapters(plans,domains);rows=[]
            for receiver,entries in enumerate(resolved):
                for target,producer,domain,node,component,sign,ordinal in entries:
                    require(1<=target<=len(context[receiver]['rows']),'physical.transform_adapter_target')
                    rows.append((group(receiver,arm),tuple(context[receiver]['rows'][target-1]),('domain',domain,node,component),sign));total+=1
            result[(kind,key)]=canonical_writes(rows)
        else:
            if context[0]['name'] in INVERSE_TRANSFERS325:
                for plan in plans:inverse_transfer_contract325(plan)
            resolved=v.reconcile_static(plans)
            name=context[0]['name'];rows=[]
            if name=='ghost.plan':
                for receiver,entries in enumerate(resolved):
                    for item in entries:
                        _,producer,ident,*rest=item
                        # Manifest includes physical source and destination
                        # block identities; local/remote is a changed route.
                        rows.append((group(receiver,arm),tuple(ident),group(producer,arm)));total+=1

                local_writes=[]
                for receiver,p in enumerate(plans):
                    order=ghost_local_order325(p,receiver)
                    entries={entry[0]:entry for entry in resolved[receiver]}
                    expected={i for i,entry in entries.items() if entry[1]==receiver}
                    require(sorted(i-1 for i in order)==sorted(expected),'physical.ghost_local_order_coverage')
                    for i in order:
                        ident=entries[i-1][2]
                        local_writes.append((group(receiver,arm),tuple(ident[2:4]),tuple(ident[:2]),ident[4]))
                result[(kind,key)]=sha([sorted(rows,key=repr),canonical_writes(local_writes)])
            elif name=='writeback.plan':
                # Check every active wire manifest, then compare canonical
                # complete coverage including retained home blocks separately.
                for receiver,entries in enumerate(resolved):
                    p=plans[receiver];pre='block_writeback_plan%'
                    owned=[i+1 for i in range(640) if cov.owner325(i,arm)==receiver]
                    sends=[cat for cat in owned if (cat-1)//4!=receiver]
                    receives=[i+1 for i in range(640) if i//4==receiver and cov.owner325(i,arm)!=receiver]
                    require(sorted(p[pre+'send_block'][:p['n_send'][0]])==sends and sorted(p[pre+'recv_block'][:p['n_recv'][0]])==receives,'physical.writeback_complete_blocks')
                    require(p['n_retained']==[len(owned)-len(sends)] and p['installed_block_count']==[len(owned)],'physical.writeback_retained_coverage')
                    for producer,ordinal,ident in entries:
                        cat=ident[0];require(1<=cat<=640 and receiver==(cat-1)//4,'physical.writeback_home')
                        expected=cov.owner325(cat-1,arm);require(producer==expected,'physical.writeback_owner');total+=1
                result[(kind,key)]=sha(reconcile_boundary325(plans,context,lambda cat:cov.owner325(cat-1,arm)))
            elif name in ('inverse.alias','inverse.range_alias'):
                # Fixed Domain-to-Domain alias install order is retained.
                result[(kind,key)]=sha(resolved);total+=sum(map(len,resolved))
            else:
                for receiver,entries in enumerate(resolved):
                    for target,producer,key4,address,ordinal in entries:
                        cat=key4[0];require(1<=cat<=640 and producer==cov.owner325(cat-1,arm),'physical.inverse_owner')
                        # Address is checked against the actual block-store
                        # order, not a presumed catalogue-order installation.
                        blocks=context[producer]['blocks'];require(1<=address[0]<=len(blocks),'physical.inverse_storage')
                        require(address_key_fields(address,blocks)==tuple(key4),'physical.inverse_catalogue_address')
                        rows.append((receiver,target,tuple(key4)));total+=1
                result[(kind,key)]=canonical_writes(rows)
    if progress:progress('providers',arm,len(families),len(families),None,True)
    return {repr(k):value for k,value in result.items()},{'families':len(families),'resolved_consumers':total,'publication_imported_consumers':publication}

def canonical_writes(rows):
    grouped=defaultdict(list)
    for row in rows:grouped[repr(row[:2])].append(row[2:])
    return sha(sorted(grouped.items()))

def compare_features(c,left,right):
    a=dict(c.execute('select key,digest from feature where arm=?',(left,)))
    b=dict(c.execute('select key,digest from feature where arm=?',(right,)))
    require(a.keys()==b.keys(),'physical.feature_coverage',{'missing':list(a.keys()-b.keys())[:3],'extra':list(b.keys()-a.keys())[:3]})
    bad=[k for k in a if a[k]!=b[k]]
    require(not bad,'physical.ordered_program_or_value_difference',bad[:3])
    return {'passed':True,'features':len(a),'digest':sha(a)}

def write_report325(path,report):
    """Atomic, bounded progress output; never writes inside a capture stream."""
    if path is None:return
    path=Path(path);data=json.dumps(report,sort_keys=True,indent=2)+'\n'
    require(len(data.encode())<=16*1024*1024,'compact.report_bound')
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(data);temporary.replace(path)

def qualify(results,kit=None,placements=None,*,report_path=None,scratch_root=None):
    started=time.monotonic();stats={}
    results=Path(results);report={'passed':False,'status':'running','stage':325,
        'contract':'stage325-closed-physical-v4',
        'claim':'configured level5 four-step physical program/value/provider qualification; performance held',
        'scratch_limit_bytes':SCRATCH_LIMIT,'runs':{}}
    def persist():
        report.update(stats);report['elapsed_seconds']=time.monotonic()-started
        write_report325(report_path,report)
    def progress(stage,arm,completed,total,detail=None,save=False):
        report['progress']={'stage':stage,'arm':arm,'completed':completed,'total':total,'detail':detail}
        if save:persist()
    try:
        persist()
        with tempfile.TemporaryDirectory(prefix='rk325-semantic-',dir=scratch_root) as scratch:
            c=_create_cache(Path(scratch)/'routes.sqlite')
            try:
                proofs={}
                for arm in ('checked-L256','qual-C0','qual-L256'):
                    report['runs'][arm]={'streams':[],'status':'streams'}
                    cache_run(c,results/arm,arm,scratch,started,stats,report['runs'][arm]['streams'],progress)
                    report['runs'][arm]['status']='providers'
                    proofs[arm],counts=route_certificates(c,arm,progress)
                    report['runs'][arm].update(counts);report['runs'][arm]['status']='complete';persist()
                progress('cross-arm','qual-C0/qual-L256',0,1,None,True)
                changes=[]
                for rank in range(N):
                    a=report['runs']['qual-C0']['streams'][rank]['storage_orders']
                    b=report['runs']['qual-L256']['streams'][PERM[rank]]['storage_orders']
                    if a!=b:changes.append(rank)
                require(len(changes)==48,'physical.changed_storage_group_count')
                report['changed_storage_groups']=changes
                report['storage_frames_checked']=5
                report['optimized_pair']=compare_features(c,'qual-C0','qual-L256')
                require(proofs['qual-C0']==proofs['qual-L256'],'physical.provider_pair_difference',
                        [k for k in proofs['qual-C0'] if proofs['qual-C0'].get(k)!=proofs['qual-L256'].get(k)][:3])
                report['provider_identity_digest']=sha(proofs['qual-C0'])
                report['checked_scope']='independent source-enabled oracle and full coverage/provider closure; checkpoint science handled by campaign'
                progress('complete',None,1,1)
            finally:c.close()
        report['passed']=True;report['status']='passed';persist()
    except BaseException as error:
        report['passed']=False;report['status']='failed'
        report['failure']={'type':type(error).__name__,'message':str(error)}
        try:persist()
        except BaseException as report_error:
            # A full disk or failed final save must not replace the original
            # qualification/write error, or leave an in-memory success flag.
            detail={'type':type(report_error).__name__,'message':str(report_error)}
            report['failure']['failure_report_write']=detail
            if hasattr(error,'add_note'):
                error.add_note('Stage325 failure report could not be saved: '+str(detail))
        raise
    return report
