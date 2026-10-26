import gzip, hashlib, json, struct, tempfile, unittest
from pathlib import Path
from collections import Counter
import verify325 as q
import verify322 as v


def entry(name,values,kind=1):
    data=(bytes(values) if kind==4 else struct.pack('<'+str(len(values))+('i' if kind==1 else 'q'),*values))
    return v.Entry(0,0,name,kind,len(values),data)
def scope(name,fields,epoch=1,record=1,end=9,parent=()):
    records=[entry(k,x) for k,x in fields.items()]
    return dict(name=name,epoch=epoch,record=record,end_record=end,parent=parent,step=1,
                records=records,values={e.name:e for e in records},tags=set(fields),
                record_counts=Counter(e.name for e in records),record_descriptors={e.name:(e.kind,e.count) for e in records})
def block(cat,patchstart=0):
    return scope('block',{'identity':[cat-1,(cat-1)//4,2+(cat-1)%4,5,1,2,1,1,30,1,3,1,0],
                         'patch.shape':[5,1],'patch':[0,patchstart,5,1,0],
                         'node.shape':[3,patchstart+16],
                         'bdry_storage.shape':[7,1],'bdry_storage':[0,5,0,4,4,16,0],
                         'bdry_node.shape':[3,16]},epoch=cat)
def wire(kind=13,count=2,data=None):
    head=v.HEADER.pack(b'RK325V01',1,0,4,8,0x01020304)
    name=b'value';data=hashlib.sha256(struct.pack('<2d',1.,2.)).digest() if data is None else data
    body=v.RECORD.pack(len(name),kind,count)+name+data
    footer=v.RECORD.pack(8,12,2)+b'complete'+struct.pack('<2q',1,len(head)+len(body))
    return head+body+footer

class CodecTests(unittest.TestCase):
    def read(self,data):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'execution316-rank-0.bin.gz'
            with gzip.open(path,'wb') as f:f.write(data)
            reader=q.Reader325(path,0);entries=list(reader)
            return entries,reader.summary
    def test_finite_source_digest_exact(self):
        rows,summary=self.read(wire())
        self.assertEqual((rows[0].kind,rows[0].count),(3,2))
        self.assertEqual(rows[0].data,hashlib.sha256(struct.pack('<2d',1.,2.)).digest())
        self.assertEqual(summary['physical_bytes'],len(wire()))
    def test_raw_real_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'descriptor'):self.read(wire(3,2,struct.pack('<2d',1.,2.)))
    def test_bad_empty_digest_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'empty_digest'):self.read(wire(13,0,b'0'*32))
    def test_good_empty_digest(self):self.read(wire(13,0,hashlib.sha256(b'').digest()))
    def test_truncation_rejected(self):
        with self.assertRaises(ValueError):self.read(wire()[:-1])
    def test_missing_footer_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'complete'):self.read(wire()[:-40])
    def test_after_footer_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'after_footer'):self.read(wire()+b'x'*16)
    def test_original_magic_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'header'):self.read(b'RK316V01'+wire()[8:])
    def test_count_width_bound(self):
        data=bytearray(wire());struct.pack_into('<q',data,36,2**31)
        with self.assertRaisesRegex(v.AuditError,'size_bound'):self.read(bytes(data))

class PhysicalTests(unittest.TestCase):
    def test_map_has_same_owner_groups(self):
        self.assertEqual(sorted(q.PERM),list(range(256)))
        current=[q.cov.owner325(i,'qual-C0') for i in range(640)]
        new=[q.cov.owner325(i,'qual-L256') for i in range(640)]
        self.assertEqual(new,[q.PERM[x] for x in current])
        self.assertEqual(sum(o==i//4 for i,o in enumerate(current)),6)
        self.assertEqual(sum(o==i//4 for i,o in enumerate(new)),400)
        self.assertEqual(sorted(Counter(current).values()),sorted(Counter(new).values()))
    def test_exactly_48_storage_groups_change(self):
        changes=[r for r in range(256) if q.expected_storage_order(r,'qual-C0')!=q.expected_storage_order(q.PERM[r],'qual-L256')]
        self.assertEqual(len(changes),48)
    def test_storage_order_is_physical_not_rank(self):
        old=[block(9),block(7),block(8)]
        new=[block(7),block(8),block(9)]
        self.assertEqual(q.address_key([2,1,0,16,0,30],old),q.address_key([1,1,0,16,0,30],new))
        self.assertNotEqual(q.address_key([2,1,0,16,0,30],old),q.address_key([2,1,0,16,0,30],new))
    def test_wrong_storage_address_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'patch_address'):q.address_key([1,1,16,16,0,30],[block(7)])
    def test_field_layer_address_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'address_layers'):q.address_key([1,1,0,16,1,30],[block(7)])
    def test_boundary_address(self):self.assertEqual(q.address_key([1,2,3,16,0,30],[block(7)]),(7,2,1,3))
    def test_duplicate_alias_final_writer_preserved(self):
        a=[(0,('target',1),('source',4),1),(0,('target',1),('source',5),-1)]
        self.assertNotEqual(q.canonical_writes(a),q.canonical_writes(list(reversed(a))))
    def test_independent_target_order_is_irrelevant(self):
        a=[(0,('target',1),('source',4),1),(0,('target',2),('source',5),-1)]
        self.assertEqual(q.canonical_writes(a),q.canonical_writes(list(reversed(a))))
    def test_orientation_mutation_detected(self):
        self.assertNotEqual(q.canonical_writes([(0,1,2,1)]),q.canonical_writes([(0,1,2,-1)]))
    def test_equal_physical_key_wrong_numerical_writer_is_not_hidden(self):
        ctx={'rows':[('physical',1,2,1)],'role':'numerical'}
        self.assertNotEqual(q._provenance(ctx,0,1,'qual-C0'),q._provenance(ctx,1,1,'qual-C0'))
    def test_fixed_provider_does_not_permute(self):
        ctx={'rows':[('physical',3,2,1)],'role':'numerical','fixed_prefix':1}
        self.assertEqual(q._provenance(ctx,2,1,'qual-C0'),q._provenance(ctx,2,1,'qual-L256'))
    def test_numerical_provider_permutation_preserved(self):
        ctx={'rows':[('physical',3,2,1)],'role':'numerical'}
        self.assertEqual(q._provenance(ctx,2,1,'qual-C0'),q._provenance(ctx,q.PERM[2],1,'qual-L256'))
    def test_raw_provider_slot_out_of_bounds_rejected(self):
        with self.assertRaisesRegex(v.AuditError,'source_or_target_slot'):q._provenance({'rows':[]},0,1,'qual-C0')
    def test_extra_real_replaced_by_integer_is_rejected(self):
        s=scope('rhs.rk_values322',{'rk.adapter.shape':[0,150], 'rk.value.shape':[0,150],
            'rk.adapter.numeric.shape':[0,150], 'rk.adapter.numeric':[],
            'rk.adapter.domain.shape':[0,150], 'rk.adapter.domain':[]})
        with self.assertRaisesRegex(v.AuditError,'extra_real_kind'):
            q.validate_coverage([s],0,256,4,'candidate','qual-C0')
    def test_extra_provider_key_wrong_kind_is_rejected(self):
        s=scope('backend.input',{'input.extent':[0,0],
            'input.provider_key325.shape':[4,0],'input.provider_key325':[],
            'input.physical_key325.shape':[4,0],'input.physical_key325':[]})
        bad=entry('input.provider_key325',[],2)
        s['records']=[bad if e.name==bad.name else e for e in s['records']];s['values'][bad.name]=bad
        with self.assertRaisesRegex(v.AuditError,'extra_key_kind'):
            q.validate_coverage([s],0,256,4,'candidate','qual-C0')
    def test_missing_new_record_rejected_even_in_both_arms(self):
        s=scope('backend.input',{'input.extent':[0,0]})
        with self.assertRaisesRegex(v.AuditError,'extra_coverage'):q.validate_coverage([s],0,256,4,'candidate','qual-C0')


class WorldRouteTests(unittest.TestCase):
    def test_actual_wire_provider_not_equal_value_inference(self):
        # Two actual wire plans. The checker follows the requested source slot;
        # no numerical values participate in provider reconstruction.
        def plan(rank,sc,rc,source,target,address,peer):
            return {'route.allocated':[1],'route.epoch':[1,1,0],
                    'route.header':[501,2,1,rank,2,501,0],'route.flags':[0,0],
                    'route.sc':sc,'route.sd':v.prefix(sc),'route.rc':rc,'route.rd':v.prefix(rc),
                    'route.source':source,'route.target':target,'route.address':address,
                    'route.orientation':[1]*len(target),'route.send_peer':[i+1 for i,c in enumerate(sc) if c],
                    'route.receive_peer':[i+1 for i,c in enumerate(rc) if c],
                    'route.target_peer':peer,'route.allocations':[0,0,sum(c>0 for c in sc)+sum(c>0 for c in rc)]}
        a=plan(0,[0,1],[0,0],[2],[],[],[])
        b=plan(1,[0,0],[1,0],[],[1],[0],[1])
        self.assertEqual(v.reconcile_routes([a,b])[1],[(1,0,2,1,0)])
        a['route.source']=[1]
        self.assertEqual(v.reconcile_routes([a,b])[1],[(1,0,1,1,0)])
    def test_wrong_wire_count_rejected(self):
        # Preserve total receive size but corrupt one sender/receiver agreement.
        def p(rank,sc,rc):
            return {'route.allocated':[1],'route.epoch':[1,1,0],'route.header':[1,2,1,rank,2,1,0],
                    'route.flags':[0,0],'route.sc':sc,'route.sd':v.prefix(sc),'route.rc':rc,'route.rd':v.prefix(rc),
                    'route.source':[1]*sum(sc),'route.target':[],'route.address':[],'route.orientation':[],
                    'route.send_peer':[i+1 for i,c in enumerate(sc) if c],'route.receive_peer':[i+1 for i,c in enumerate(rc) if c],
                    'route.target_peer':[],'route.allocations':[0,0,sum(c>0 for c in sc)+sum(c>0 for c in rc)]}
        with self.assertRaisesRegex(v.AuditError,'wire_count'):v.reconcile_routes([p(0,[0,1],[0,0]),p(1,[0,0],[0,1])])
    def test_transform_publication_cannot_skip_a_consumer(self):
        # Imported row 2 must appear in the actual route target list even when
        # all numerical data would happen to have equal values.
        plan=scope('transform.plan',{
          'extent':[1,2,0,1,0,5,5,30],
          'placement%key.shape':[4,1],'placement%key':[1,7,1,1],
          'placement%owned':[1], 'save_node':[],'save_index':[0],'fixed_node':[],'fixed_index':[0],
          'publication.key325.shape':[4,1],'publication.key325':[2,8,1,1],'publication_row':[2]
        },record=1,end=10)
        child=scope('transform.scalar_publication',{},record=2,end=3,parent=(('transform.plan',1),))
        route={'route.allocated':[1],'route.header':[9301,2,60,0,256,9301,0], 'route.target':[],'route.orientation':[]}
        with self.assertRaisesRegex(v.AuditError,'complete_targets'):q.route_map(route,child,[plan,child],0,'qual-C0')
        route['route.target']=[2];route['route.orientation']=[1]
        result=q.route_map(route,child,[plan,child],0,'qual-C0')
        self.assertEqual(result['rows'][1],('physical',2,8,1))

if __name__=="__main__":unittest.main()
