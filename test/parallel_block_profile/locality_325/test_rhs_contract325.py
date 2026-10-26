"""Closed live RHS contracts, including the phase-4 failure in job55569."""
import unittest
import verify325 as q
import verify322 as v
from test_verify325 import entry, scope

# Derived independently from routes.f90 (components) and distributed.f90
# (capacity). Logical phase18 has no physical route; it reuses physical16.
TABLE={1:(3,2),2:(1,2),3:(3,1),4:(2,1),5:(3,2),6:(2,1),7:(1,1),
       11:(2,1),12:(3,1),13:(3,2),14:(1,2),15:(3,1),16:(1,2),20:(3,2)}

def route(phase,nd=2,rank=0):
    m,nc=TABLE[phase];encoded=500+phase
    return {'route.allocated':[1],'route.header':[encoded,m*nd,nc,rank,256,3*encoded+nc,nc],'route.flags':[0,1]}

def prepared():
    logical={**{p:p for p in TABLE},18:16}
    s=scope('rhs.prepare',{'placement.key.shape':[4,2],'placement.key':[1,7,1,1,2,8,1,1],
        'transport.route.bounds':[5,5],'transport.bound':[int(p in TABLE) for p in range(1,21)],
        'transport.route_phase':[logical.get(p,0) for p in range(1,21)],
        'transport.route_level':[5 if p in logical else 0 for p in range(1,21)]})
    for p in TABLE:
        s['records'].append(entry('transport.route.coordinate',[p,5]))
        s['records'] += [entry(k,x) for k,x in route(p).items()]
    return s

def replace(s,name,values):
    e=entry(name,values);s['values'][name]=e
    s['records']=[e if old.name==name else old for old in s['records']]

def mutate(s,phase,name,values):
    current=None
    for i,e in enumerate(s['records']):
        if e.name=='transport.route.coordinate':current=v.value(e)[0]
        if current==phase and e.name==name:s['records'][i]=entry(name,values);return
    raise AssertionError((phase,name))

class RhsContractTests(unittest.TestCase):
    def test_all_fourteen_source_routes(self):
        self.assertEqual(q.RHS_PHASE325,TABLE);s=prepared();q.validate_rhs_routes(s,0)
        self.assertEqual(len(q.record_scope_routes([s],0)),14)
        for p,(m,nc) in TABLE.items():
            with self.subTest(phase=p):
                ctx=q.route_map(route(p),s,[s],0,'qual-C0')
                self.assertEqual(ctx['rows'],[('physical',d,n,c) for d,n in [(1,7),(2,8)] for c in range(1,m+1)])
    def test_second_triangle_component_preserves_provider_identity(self):
        s=prepared()
        for p in (4,6,11):
            with self.subTest(phase=p):
                ctx=q.route_map(route(p),s,[s],0,'qual-C0')
                self.assertEqual(ctx['rows'][1],('physical',1,7,2))
                actual=(0,ctx['rows'][1],q._provenance(ctx,0,2,'qual-C0'),1)
                wrong=(0,ctx['rows'][1],q._provenance(ctx,0,1,'qual-C0'),1)
                self.assertNotEqual(q.canonical_writes([actual]),q.canonical_writes([wrong]))
    def test_all_wrong_component_multipliers(self):
        for p,(m,nc) in TABLE.items():
            for bad in set((0,1,2,3,4))-{m}:
                with self.subTest(phase=p,components=bad):
                    r=route(p);r['route.header'][1]=2*bad
                    with self.assertRaisesRegex(v.AuditError,'rhs_components'):q.rhs_route_layout(r,2,0,[p,5])
    def test_all_wrong_channel_capacities_and_last_columns(self):
        for p,(m,nc) in TABLE.items():
            for column in (2,6):
                with self.subTest(phase=p,column=column):
                    r=route(p);r['route.header'][column]=3-nc
                    with self.assertRaisesRegex(v.AuditError,'rhs_columns'):q.rhs_route_layout(r,2,0,[p,5])
    def test_all_inactive_phases_rejected(self):
        for p in set(range(22))-TABLE.keys():
            with self.subTest(phase=p):
                r=route(1);r['route.header'][0]=500+p
                with self.assertRaisesRegex(v.AuditError,'rhs_phase'):q.rhs_route_layout(r,2,0,[p,5])
    def test_other_levels_rejected(self):
        for level in (0,4,6):
            r=route(4);r['route.header'][0]=100*level+4
            with self.assertRaisesRegex(v.AuditError,'rhs_phase'):q.rhs_route_layout(r,2,0,[4,level])
    def test_wrong_rank_world_tag_rejected(self):
        for column in (3,4,5):
            r=route(4);r['route.header'][column]+=1
            with self.assertRaisesRegex(v.AuditError,'rhs_header_contract'):q.rhs_route_layout(r,2,0,[4,5])
    def test_wrong_flags_rejected(self):
        for flags in ([1,1],[0,0],[1,0]):
            r=route(4);r['route.flags']=flags
            with self.assertRaisesRegex(v.AuditError,'shared_variable'):q.rhs_route_layout(r,2,0,[4,5])
    def test_zero_nodes_not_in_sealed_map(self):
        with self.assertRaisesRegex(v.AuditError,'rhs_components'):q.rhs_route_layout(route(4,0),0,0,[4,5])
    def test_unallocated_compiled_route_rejected(self):
        s=prepared()
        with self.assertRaisesRegex(v.AuditError,'rhs_allocated'):q.route_map({'route.allocated':[0]},s,[s],0,'qual-C0')
    def test_same_shape_header_swap_rejected(self):
        s=prepared();mutate(s,4,'route.header',route(6)['route.header']);mutate(s,6,'route.header',route(4)['route.header'])
        with self.assertRaisesRegex(v.AuditError,'coordinate_header'):q.validate_rhs_routes(s,0)
    def test_missing_compiled_phase_rejected(self):
        s=prepared();current=None;kept=[]
        for e in s['records']:
            if e.name=='transport.route.coordinate':current=v.value(e)[0]
            if current!=4:kept.append(e)
        s['records']=kept
        with self.assertRaisesRegex(v.AuditError,'physical_census'):q.validate_rhs_routes(s,0)
    def test_duplicate_phase_rejected(self):
        s=prepared();s['records']+=s['records'][-4:]
        with self.assertRaisesRegex(v.AuditError,'physical_census'):q.validate_rhs_routes(s,0)
    def test_coordinate_without_route_rejected(self):
        s=prepared();i=next(i for i,e in enumerate(s['records']) if e.name=='transport.route.coordinate')
        s['records'].insert(i,entry('transport.route.coordinate',[4,5]))
        with self.assertRaisesRegex(v.AuditError,'coordinate_without_route'):q.validate_rhs_routes(s,0)
    def test_route_without_coordinate_rejected(self):
        s=prepared();i=next(i for i,e in enumerate(s['records']) if e.name=='transport.route.coordinate');del s['records'][i]
        with self.assertRaisesRegex(v.AuditError,'route_without_coordinate'):q.validate_rhs_routes(s,0)
    def test_bound_census_rejects_missing_and_extra_routes(self):
        for p in (4,6,11,18):
            s=prepared();b=q.integers(s,'transport.bound');b[p-1]=1-b[p-1];replace(s,'transport.bound',b)
            with self.assertRaisesRegex(v.AuditError,'bound_census'):q.validate_rhs_routes(s,0)
    def test_injected_bernoulli_reuses_sixteen(self):
        s=prepared();p=q.integers(s,'transport.route_phase');p[17]=18;replace(s,'transport.route_phase',p)
        with self.assertRaisesRegex(v.AuditError,'logical_phase_binding'):q.validate_rhs_routes(s,0)
    def test_unused_and_active_logical_level_bindings(self):
        for name,slot,value in [('transport.route_phase',16,17),('transport.route_level',16,5),('transport.route_level',17,0)]:
            s=prepared();a=q.integers(s,name);a[slot]=value;replace(s,name,a)
            with self.assertRaisesRegex(v.AuditError,'logical_.*_binding'):q.validate_rhs_routes(s,0)
    def test_failure_details(self):
        r=route(4);r['route.header'][1]=6
        with self.assertRaises(v.AuditError) as error:q.rhs_route_layout(r,2,0,[4,5])
        for text in ("'phase': 4","'nodes': 2","'slots': 6","'components': 3.0","'expected_components': 2"):
            self.assertIn(text,str(error.exception))

if __name__=='__main__':unittest.main()
