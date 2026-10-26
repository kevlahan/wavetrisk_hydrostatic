"""End-to-end provider certificates: real SQLite and reconcilers, no mocks.

Small hand-constructed worlds retain all 256 rank slots. They are checker tests,
not native correctness evidence. Expected writes are specified independently.
"""
import copy, gzip, hashlib, json, unittest
from pathlib import Path
import verify325 as q
import verify322 as v
from test_boundary325 import fixture as boundary_fixture, P as WP

N=256

def certify(kind,name,plans,contexts,arm='qual-C0',domains=None):
    c=q._create_cache(':memory:')
    try:
        for rank,(plan,context) in enumerate(zip(plans,contexts)):
            c.execute('insert into route values(?,?,?,?,?)',(arm,rank,kind,name,json.dumps(plan)))
            c.execute('insert into context values(?,?,?,?,?)',(arm,rank,kind,name,json.dumps(context)))
            c.execute('insert into domain values(?,?,?)',(arm,rank,json.dumps((domains or [[r] if r<160 else [] for r in range(N)])[rank])))
        return q.route_certificates(c,arm)
    finally:c.close()


def routes(role='numerical',publication_vector=False,sender=0,receiver=1):
    rows=[['physical',11,21,1],['physical',12,22,1]]
    if publication_vector:rows=[['physical',11,21,i] for i in (1,2,3)]
    plans=[];contexts=[]
    for r in range(N):
        sc=[0]*N;rc=[0]*N
        if r==sender:sc[receiver]=2
        if r==receiver:rc[sender]=2
        sources=[1,2] if r==sender else []
        targets=[1,1] if r==receiver else []
        p={'route.allocated':[1],'route.epoch':[1,1,0],
           'route.header':[501,len(rows),1,r,N,501,0],'route.flags':[0,0],
           'route.sc':sc,'route.sd':v.prefix(sc),'route.rc':rc,'route.rd':v.prefix(rc),
           'route.source':sources,'route.target':targets,'route.address':[0,1] if r==receiver else [],
           'route.orientation':[1,-1] if r==receiver else [],'route.send_peer':[receiver+1] if r==sender else [],
           'route.receive_peer':[sender+1] if r==receiver else [],'route.target_peer':[sender+1,sender+1] if r==receiver else [],
           'route.allocations':[0,0,int(r in (sender,receiver))]}
        ctx={'rows':copy.deepcopy(rows),'role':role,'canonical_rows':len(rows),'owned':[1]*len(rows)}
        if role=='publication' and r==receiver:
            p['route.target']=[1,2];p['route.orientation']=[1,1]
        plans.append(p);contexts.append(ctx)
    return plans,contexts


def ghosts(arm='qual-C0',active=(0,1)):
    physical=lambda r:q.PERM[r] if q.is_locality(arm) else r
    # request tuple: provider, source block, source patch, destination, ghost, level
    requests=[[] for _ in range(N)]
    requests[physical(active[0])]=[(physical(active[0]),101,1,301,1,5),(physical(active[0]),102,2,301,1,5)]
    requests[physical(active[1])]=[(physical(active[0]),101,1,302,1,5)]
    plans=[]
    pre='ghost_exchange_plan%'
    for rank,req in enumerate(requests):
        sc=[sum(x[0]==peer and peer!=rank for x in req) for peer in range(N)]
        incoming=[[x for x in requests[peer] if x[0]==rank and peer!=rank] for peer in range(N)]
        rc=[len(x) for x in incoming]
        indices=[i+1 for peer in range(N) for i,x in enumerate(req) if x[0]==peer and peer!=rank]
        p={'__scope':'ghost.plan',pre+'local_order297.shape':[-1],pre+'local_start297.shape':[-1],
           pre+'source_owner':[x[0] for x in req],pre+'request_index':indices,
           pre+'send_data':[k for i in indices for k in req[i-1][1:]],
           pre+'recv_data':[k for group in incoming for x in group for k in x[1:]]}
        for j,field in enumerate(('source_block','source_local_patch','destination_block','destination_ghost','source_patch_level'),1):
            p[pre+field]=[x[j] for x in req]
        for direction,counts in [('send',sc),('recv',rc)]:
            p[pre+direction+'_record_count']=counts;p[pre+direction+'_record_displ']=v.prefix(counts)
        for field,width in [('scalar',2),('vector',3)]:
            p[pre+field+'_send_count']=[n*width for n in rc]
            p[pre+field+'_recv_count']=[n*width for n in sc]
        plans.append(p)
    return plans,[{'name':'ghost.plan','blocks':[]} for _ in range(N)]


def inverses(name,unallocated=False,ordered=False):
    alias=name in ('inverse.alias','inverse.range_alias')
    plans=[];contexts=[]
    cat=next(c for c in range(1,641) if q.cov.owner325(c-1,'qual-C0')==0)
    block={'catalog':cat,'fields':{'identity':[cat-1,(cat-1)//4,2,5,1,2,1,1,30,1,3,1,0],
            'patch.shape':[5,1],'patch':[0,0,5,1,0],'node.shape':[3,16],
            'bdry_storage':[0,5,0,4,4,16,0],'bdry_node.shape':[3,16]}}
    for rank in range(N):
        p={'__scope':name,'edgewise':[0]};ctx={'name':name,'blocks':[copy.deepcopy(block)] if rank==0 else [],'inverse_nodes':[(rank+1,100+i) for i in range(8)],'alias_component':1}
        if unallocated:p.update({'plan%sc.shape':[-1],'plan%address.shape':[-1]})
        else:
            sc=[0]*N;rc=[0]*N
            if alias:
                if rank==0:sc[1]=2
                if rank==1:rc[0]=2
            else:
                if rank==1:sc[0]=1
                if rank==0:rc[1]=1
            p.update({'plan%sc':sc,'plan%sd':v.prefix(sc),'plan%rc':rc,'plan%rd':v.prefix(rc)})
            if alias:
                p.update({'plan%source':[7,8] if rank==0 else [],'plan%dest':[4,4] if rank==1 else [],'plan%install_order.shape':[-1]})
                if ordered:p.update({'plan%install_order.shape':[sum(rc)],'plan%install_order':[2,1] if rank==1 else []})
            else:
                p.update({'plan%slot':[3] if rank==1 else [],'plan%key.shape':[4,sum(rc)],
                    'plan%key':[cat,1,1,0] if rank==0 else [],'plan%address.shape':[6,sum(rc)],
                    'plan%address':[1,1,0,16,0,30] if rank==0 else [],
                    'plan%source_info.shape':[-1],'plan%target_info.shape':[-1]})
        plans.append(p);contexts.append(ctx)
    return plans,contexts,cat


def adapters():
    plans=[];contexts=[]
    for r in range(N):
        sc=[0]*N;rc=[0]*N
        if r==1:sc[0]=1
        if r==0:rc[1]=1
        plans.append({'extent':[1,1,30],'plan%sc':sc,'plan%sd':v.prefix(sc),'plan%rc':rc,'plan%rd':v.prefix(rc),
                      'plan%source.shape':[3,sum(rc)],'plan%source':[1,7,1] if r==0 else [],
                      'plan%sent.shape':[30,sum(rc)],'plan%received.shape':[30,sum(sc)],
                      'plan%send_peer':[2] if r==0 else [],'plan%receive_peer':[1] if r==1 else [],
                      'plan%target':[1] if r==1 else [],'plan%sign':[-1] if r==1 else [],'plan%index':[1] if r==1 else []})
        contexts.append({'rows':[['physical',2,7,1]],'name':'transform.scalar_adapter'})
    return plans,contexts


def writebacks(arm):
    owners=[q.cov.owner325(i,arm) for i in range(640)]
    plans,contexts=boundary_fixture(owners)
    for rank,(p,ctx) in enumerate(zip(plans,contexts)):
        ctx['name']='writeback.plan';p['__scope']='writeback.plan'
        p['installed_block_count']=[len(ctx['blocks'])]
        for direction in ('send','recv'):
            cats=p[WP+direction+'_block'][:p['n_'+direction][0]]
            counts=[0]*N
            for cat in cats:counts[(cat-1)//4 if direction=='send' else owners[cat-1]]+=1
            p[WP+direction+'_count']=counts;p[WP+direction+'_displ']=v.prefix(counts)
            for field in ('patch_count','scalar_nvalue','vector_nvalue'):
                p[WP+direction+'_'+field]=[cat+{'patch_count':1,'scalar_nvalue':2,'vector_nvalue':3}[field] for cat in cats] or [0]
            for field,width in [('scalar',2),('vector',3)]:p[WP+field+'_'+direction+'_count']=[width*x for x in counts]
    return plans,contexts

class ProviderPathsTests(unittest.TestCase):
    def test_routes_numerical_fixed_state_and_fixed_prefix(self):
        for role in ('numerical','fixed'):
            for provenance in ('numerical','fixed-prefix','state'):
                with self.subTest(role=role,provenance=provenance):
                    p,c=routes(role,sender=2,receiver=69)
                    self.assertNotEqual(q.INV[2],2);self.assertNotEqual(q.INV[69],69)
                    if provenance=='fixed-prefix':c[2]['fixed_prefix']=2
                    if provenance=='state':c[2]['rows']=[['state',1,1,1,1,1],['state',2,1,1,1,1]]
                    proofs,count=certify('routes','example',p,c,'qual-L256')
                    receiver=69 if role=='fixed' else q.INV[69]
                    providers=[tuple(x) if provenance=='state' else ('provider',2 if provenance=='fixed-prefix' else q.INV[2],*x) for x in c[2]['rows']]
                    expected=q.canonical_writes([(receiver,tuple(c[69]['rows'][0]),providers[0],1),(receiver,tuple(c[69]['rows'][0]),providers[1],-1)])
                    self.assertEqual(proofs[repr(('routes','example'))],expected);self.assertEqual(count['resolved_consumers'],2)
    def test_route_duplicate_destination_order_is_significant(self):
        p,c=routes();before=certify('routes','example',p,c)[0]
        p[1]['route.address'].reverse()
        self.assertNotEqual(before,certify('routes','example',p,c)[0])
    def test_unallocated_routes_and_collective_mismatch(self):
        p=[{'route.allocated':[0]} for _ in range(N)];c=[{'rows':[],'role':'none'} for _ in range(N)]
        self.assertEqual(certify('routes','absent',p,c)[1]['resolved_consumers'],0)
        p[1]={'route.allocated':[1]}
        with self.assertRaisesRegex(v.AuditError,'collective_allocation'):certify('routes','absent',p,c)
    def test_scalar_vector_publication_and_wrong_provider(self):
        for vector in (False,True):
            p,c=routes('publication',vector);key='transform.'+('vector' if vector else 'scalar')+'_publication'
            self.assertEqual(certify('routes',key,p,c)[1]['publication_imported_consumers'],2)
            c[0]['owned'][0]=0
            with self.assertRaisesRegex(v.AuditError,'authoritative_owner'):certify('routes',key,p,c)
    def test_adapter_actual_fixed_domain_and_mutation(self):
        p,c=adapters();key='transform.scalar_adapter';proof,count=certify('adapters',key,p,c)
        self.assertEqual(proof[repr(('adapters',key))],q.canonical_writes([(1,('physical',2,7,1),('domain',1,7,1),-1)]))
        self.assertEqual(count['resolved_consumers'],1)
        p[0]['plan%source'][0]=2
        with self.assertRaisesRegex(v.AuditError,'provider_domain'):certify('adapters',key,p,c)
    def test_real_ghost_provider_branch_and_permutation(self):
        a,ca=ghosts(active=(2,69));b,cb=ghosts('qual-L256',active=(2,69))
        pa,counts=certify('static','ghost.plan',a,ca);pb,_=certify('static','ghost.plan',b,cb,'qual-L256')
        self.assertEqual(pa,pb);self.assertEqual(counts['resolved_consumers'],3)
        pre='ghost_exchange_plan%'
        # Two ordered local writes to one destination must not be sorted away.
        for field in ('source_block','source_local_patch'):a[2][pre+field].reverse()
        self.assertNotEqual(pa,certify('static','ghost.plan',a,ca)[0])
    def test_ghost_wire_corruption_rejected(self):
        p,c=ghosts();p[0]['ghost_exchange_plan%recv_data'][0]+=1
        with self.assertRaisesRegex(v.AuditError,'actual_provider_manifest'):certify('static','ghost.plan',p,c)
    def test_ghost_disabled_mode_does_not_accept_payload(self):
        p,c=ghosts();p[0]['ghost_exchange_plan%local_order297']=[]
        with self.assertRaisesRegex(v.AuditError,'ungrouped_allocation'):certify('static','ghost.plan',p,c)
    def test_each_inverse_transfer_scope_and_source_metadata(self):
        for name in ('inverse.interior','inverse.level_interior','inverse.boundary','inverse.level_boundary','inverse.outer_install'):
            for metadata in (False,True):
                with self.subTest(name=name,metadata=metadata):
                    p,c,cat=inverses(name)
                    if metadata:
                        for rank,x in enumerate(p):
                            x['plan%source_info.shape']=[2,sum(x['plan%sc'])];x['plan%target_info.shape']=[2,sum(x['plan%rc'])]
                            x['plan%source_info']=[8,9] if rank==1 else [];x['plan%target_info']=[8,9] if rank==0 else []
                    proof,_=certify('static',name,p,c)
                    expected={'gather':q.canonical_writes([(1,(2,102),(cat,1,1,0))] if name in ('inverse.interior','inverse.level_interior') else []),
                              'scatter':q.canonical_writes([('block',(cat,1,1,0),('fixed-domain',1,(2,102)))])}
                    self.assertEqual(proof[repr(('static',name))],q.sha(expected))
                    if metadata:
                        p[0]['plan%target_info'][0]+=1
                        with self.assertRaisesRegex(v.AuditError,'source_metadata'):certify('static',name,p,c)
                        p[0]['plan%target_info'][0]-=1
                    p[0]['plan%address'][0]=2
                    with self.assertRaisesRegex(v.AuditError,'inverse_storage'):certify('static',name,p,c)
    def test_inverse_unallocated_and_mixed_allocation(self):
        p,c,_=inverses('inverse.level_boundary',True)
        self.assertEqual(certify('static','inverse.level_boundary',p,c)[1]['resolved_consumers'],0)
        p[1]=inverses('inverse.level_boundary')[0][1]
        with self.assertRaisesRegex(v.AuditError,'collective_allocation'):certify('static','inverse.level_boundary',p,c)
    def test_inverse_wrong_option_or_metadata_rejected(self):
        for mutation in ('edgewise','address_width','metadata_asymmetry','metadata_hidden_payload'):
            p,c,_=inverses('inverse.outer_install')
            if mutation=='edgewise':p[0]['edgewise']=[1]
            elif mutation=='address_width':p[0]['plan%address.shape'][0]=7
            elif mutation=='metadata_asymmetry':p[0]['plan%source_info.shape']=[2,0]
            else:p[0]['plan%source_info']=[]
            with self.subTest(mutation=mutation),self.assertRaises(v.AuditError):certify('static','inverse.outer_install',p,c)
    def test_missing_rank_and_unknown_static_family_rejected(self):
        p,c=ghosts()
        with self.assertRaisesRegex(v.AuditError,'route_family_ranks'):certify('static','ghost.plan',p[:-1],c[:-1])
        for row in p:row['__scope']='unknown'
        with self.assertRaisesRegex(v.AuditError,'unknown_family'):certify('static','unknown',p,c)

    def test_alias_default_and_explicit_install_order(self):
        for name in ('inverse.alias','inverse.range_alias'):
            results=[]
            for ordered in (False,True):
                p,c,_=inverses(name,ordered=ordered);proof,_=certify('static',name,p,c);results.append(proof)
                resolved=[(1,('node',2,103),('fixed-domain',0,('node',1,106)),1),
                          (1,('node',2,103),('fixed-domain',0,('node',1,107)),1)]
                if ordered:resolved.reverse()
                self.assertEqual(proof[repr(('static',name))],q.canonical_writes(resolved))
            self.assertNotEqual(*results)
            p[1]['plan%install_order']=[1,1]
            with self.assertRaisesRegex(v.AuditError,'install_permutation'):certify('static',name,p,c)
    def test_signed_vector_alias_repeated_target_order(self):
        results=[]
        for ordered in (False,True):
            p,c,_=inverses('inverse.alias',ordered=ordered)
            for ctx in c:ctx['alias_component']=2
            p[1]['plan%dest']=[4,-4]
            proof,_=certify('static','inverse.alias',p,c);results.append(proof)
            rows=[(1,('edge',2,101,1),('fixed-domain',0,('edge',1,102,1,1)),1),
                  (1,('edge',2,101,1),('fixed-domain',0,('edge',1,102,2,1)),-1)]
            if ordered:rows.reverse()
            self.assertEqual(proof[repr(('static','inverse.alias'))],q.canonical_writes(rows))
        self.assertNotEqual(*results)

    def test_complete_writeback_boundary_path_invariant(self):
        a,ca=writebacks('qual-C0');b,cb=writebacks('qual-L256')
        self.assertEqual(certify('static','writeback.plan',a,ca)[0],certify('static','writeback.plan',b,cb,'qual-L256')[0])
        b[0][WP+'boundary_scalar_domain_send_count'][0]+=1
        with self.assertRaisesRegex(ValueError,'domain_send_count'):certify('static','writeback.plan',b,cb,'qual-L256')

class GhostModeTests(unittest.TestCase):
    def plan(self):return ghosts()[0][0]
    def test_disabled_order_including_empty(self):
        p=self.plan();self.assertEqual(q.ghost_local_order325(p,0),[1,2]);self.assertEqual(q.ghost_local_order325(p,255),[])
    def test_both_unallocated_descriptors_required(self):
        for field in ('local_order297','local_start297'):
            for shape in (None,[],[0],[2]):
                p=self.plan();p['ghost_exchange_plan%'+field+'.shape']=shape
                with self.assertRaisesRegex(v.AuditError,'ungrouped_allocation'):q.ghost_local_order325(p,0)
    def test_grouped_mode_is_explicit_and_validates_groups(self):
        p=self.plan();pre='ghost_exchange_plan%'
        p[pre+'source_block']=[102,101]
        p.update({pre+'local_order297':[2,1],pre+'local_order297.shape':[2],pre+'local_start297':[1,2,3],pre+'local_start297.shape':[3]})
        self.assertEqual(q.ghost_local_order325(p,0,True),[2,1])
        with self.assertRaisesRegex(v.AuditError,'ungrouped_allocation'):q.ghost_local_order325(p,0)
        for order,starts in [([1,1],[1,2,3]),([2,1],[1,3]),([2,1],[0,2,3]),([2,1],[1,2,2,3])]:
            x=copy.deepcopy(p);x[pre+'local_order297']=order;x[pre+'local_start297']=starts;x[pre+'local_start297.shape']=[len(starts)]
            with self.assertRaises(v.AuditError):q.ghost_local_order325(x,0,True)
    def test_saved_native_missing_payload_is_valid_disabled_mode(self):
        path=Path(__file__).parent/'fixture/native-ghost55569-rank0.json.gz'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),'2ff0dadecc9064a2cc4b70b2e49945e6f634a5eee419eaaca8047f827c93e714')
        record=json.loads(gzip.decompress(path.read_bytes()))
        self.assertEqual(q.ghost_local_order325(record['plan'],record['rank']),record['expected_local_request_order'])
        self.assertGreater(len(record['expected_local_request_order']),0)
        with self.assertRaisesRegex(v.AuditError,'grouped_allocation'):q.ghost_local_order325(record['plan'],record['rank'],True)

    def test_grouped_empty_has_one_sentinel(self):
        p=ghosts()[0][255];pre='ghost_exchange_plan%'
        p.update({pre+'local_order297':[],pre+'local_order297.shape':[0],pre+'local_start297':[1],pre+'local_start297.shape':[1]})
        self.assertEqual(q.ghost_local_order325(p,255,True),[])

if __name__=='__main__':unittest.main()
