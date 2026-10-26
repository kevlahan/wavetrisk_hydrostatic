"""Source-derived transform ownership and untouched dependency-row sentinels."""
import unittest
import verify325 as q
import verify322 as v
from test_verify325 import entry,scope

def placement(rank):
    return scope('transform.plan',{
        'placement%key.shape':[4,2],'placement%key':[1,7,1,1,1,8,1,1],
        'placement%owner':[rank,-1],'placement%owned':[1,0],
        'placement%state_key.shape':[4,2],'placement%state_key':[1,1,1,0,0,0,0,0],
        'placement%level':[5,-1],
        'placement%mask.shape':[4,2],'placement%mask':[1,2,3,4,0,0,0,0],
    })

def replace(s,name,values):
    e=entry(name,values);s['values'][name]=e
    s['records']=[e if old.name==name else old for old in s['records']]

class TransformOwnerTests(unittest.TestCase):
    def adapter_fixture(self):
        parent=scope('rhs.prepare',{'placement.key.shape':[4,2],'placement.key':[1,7,1,1,1,8,1,1],
            'placement.owned':[1,0],'backend.adapter_domain':[1],'backend.adapter_node':[9]},record=1,end=10)
        child=scope('backend.adapter_route',{},record=2,end=3,parent=(('rhs.prepare',1),))
        route={'route.allocated':[1],'route.header':[8001,3,150,2,256,8001,150],'route.source':[1]}
        return parent,child,route
    def test_backend_adapter_moves_numerical_source_to_fixed_consumer(self):
        p,s,route=self.adapter_fixture();ctx=q.route_map(route,s,[p,s],2,'qual-C0')
        self.assertEqual(ctx['role'],'fixed')
        self.assertEqual(q._provenance(ctx,2,1,'qual-C0'),('provider',2,'physical',1,7,1))
        self.assertEqual(q._provenance(ctx,2,3,'qual-C0',False),('physical',3,9,1))
    def test_backend_adapter_cannot_source_fixed_tail_or_unowned_prefix(self):
        for slot in (0,2,3,4):
            p,s,route=self.adapter_fixture();route['route.source']=[slot]
            with self.assertRaisesRegex(v.AuditError,'adapter_numerical_source'):q.route_map(route,s,[p,s],2,'qual-C0')
    def test_owned_and_dependency_rows_are_distinct(self):
        self.assertEqual(q.transform_owner_identity(placement(2),2,'qual-C0'),[2,-1])
    def test_exact_rank_permutation_preserves_sentinel_all_ranks(self):
        for rank in range(256):
            with self.subTest(rank=rank):
                a=q.transform_owner_identity(placement(rank),rank,'qual-C0')
                b=q.transform_owner_identity(placement(q.PERM[rank]),q.PERM[rank],'qual-L256')
                self.assertEqual(a,b);self.assertEqual(b[1],-1)
    def test_owned_row_cannot_name_different_valid_rank(self):
        s=placement(2);replace(s,'placement%owner',[3,-1])
        with self.assertRaisesRegex(v.AuditError,'canonical_owner'):q.transform_owner_identity(s,2,'qual-C0')
    def test_owned_row_cannot_use_unassigned_sentinel(self):
        s=placement(2);replace(s,'placement%owner',[-1,-1])
        with self.assertRaisesRegex(v.AuditError,'canonical_owner'):q.transform_owner_identity(s,2,'qual-C0')
    def test_dependency_owner_must_be_exactly_minus_one(self):
        for owner in (-2,0,2,255,256):
            s=placement(2);replace(s,'placement%owner',[2,owner])
            with self.assertRaisesRegex(v.AuditError,'dependency_sentinels'):q.transform_owner_identity(s,2,'qual-C0')
    def test_dependency_state_key_all_four_components_stay_zero(self):
        for column in range(4):
            s=placement(2);state=q.integers(s,'placement%state_key');state[4+column]=1;replace(s,'placement%state_key',state)
            with self.assertRaisesRegex(v.AuditError,'dependency_sentinels'):q.transform_owner_identity(s,2,'qual-C0')
    def test_dependency_level_stays_minus_one(self):
        for level in (-2,0,5):
            s=placement(2);replace(s,'placement%level',[5,level])
            with self.assertRaisesRegex(v.AuditError,'dependency_sentinels'):q.transform_owner_identity(s,2,'qual-C0')
    def test_dependency_all_four_masks_stay_zero(self):
        for column in range(4):
            s=placement(2);mask=q.integers(s,'placement%mask');mask[4+column]=1;replace(s,'placement%mask',mask)
            with self.assertRaisesRegex(v.AuditError,'dependency_sentinels'):q.transform_owner_identity(s,2,'qual-C0')
    def test_wrong_owned_flag_rejected(self):
        s=placement(2);replace(s,'placement%owned',[1,2])
        with self.assertRaisesRegex(v.AuditError,'owned_flag'):q.transform_owner_identity(s,2,'qual-C0')
    def test_extent_mismatch_rejected(self):
        s=placement(2);replace(s,'placement%owner',[2])
        with self.assertRaisesRegex(v.AuditError,'placement_extent'):q.transform_owner_identity(s,2,'qual-C0')
    def test_scaffold_owned_state_key_may_be_zero(self):
        # build_rhs_owner_map assigns scaffold nodes to the Domain's first
        # block owner and leaves their state_key zero; owned does not mean
        # this row must have an interior block address.
        s=placement(2);replace(s,'placement%state_key',[0]*8)
        self.assertEqual(q.transform_owner_identity(s,2,'qual-C0'),[2,-1])
    def test_sentinel_cannot_index_rank_permutation(self):
        for arm in ('qual-C0','qual-L256'):
            for rank in (-1,256):
                with self.assertRaisesRegex(v.AuditError,'numerical_rank_range'):q.group(rank,arm)
    def test_failure_contains_row_identity(self):
        s=placement(2);replace(s,'placement%level',[5,4])
        with self.assertRaises(v.AuditError) as error:q.transform_owner_identity(s,2,'qual-C0')
        for text in ("'rank': 2","'slot': 2","'physical_key': (1, 8, 1, 1)","'owner': -1","'level': 4"):
            self.assertIn(text,str(error.exception))

if __name__=='__main__':unittest.main()
