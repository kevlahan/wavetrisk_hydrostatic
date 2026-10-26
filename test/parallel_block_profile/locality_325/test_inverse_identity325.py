"""Adversarial offline checks: a storage bijection is not a value proof."""
import hashlib
import struct
import unittest
import inverse_identity325 as inv
import verify322 as v
from test_verify325 import scope, entry


def real(s, name, vals):
    e = v.Entry(0, 0, name, 3, len(vals), hashlib.sha256(struct.pack('<' + 'd' * len(vals), *vals)).digest())
    s['records'].append(e); s['values'][name] = e
    s['tags'].add(name); s['record_counts'][name] += 1
    return s


def world(order=(0, 1, 2)):
    p = scope('inverse.plan', {'extent': [3, 0, 2, 30, 5, 5], 'options': [1, 0],
        'node_domain': [1, 1, 1], 'node_id': list(order), 'node_flags': [0, 1, 2],
        'scaffold_slot': [], 'consumer_alias_slot': [1, 2, 3],
        'restore203_index': [0, 0, 0], 'restore203_slot': [],
        'scalar.program.shape': [30, 0], 'scalar.program': [],
        'inner.program.shape': [104, 0], 'inner.program': [],
        'outer.program.shape': [18, 0], 'outer.program': [],
        'operation_first': [1, 1], 'operation_last': [0, 0],
        'publication.program.shape': [18, 0], 'publication.program': []}, record=1, end=100)
    real(p, 'outer.coefficients', [])
    m = scope('inverse.node_map', {'slot': [order.index(i) + 1 for i in range(3)],
        'key.shape': [4, 3], 'key': [0] * 12}, epoch=102, record=2, end=3,
        parent=[('inverse.plan', 1)])
    return [p, m], p


def payload(vals):
    return real(scope('inverse.values322', {'inverse.value.shape': [1, 1, len(vals)]}), 'inverse.value', vals)


class IdentityTests(unittest.TestCase):
    def test_bijection_has_physical_global_domain(self):
        s, p = world((2, 0, 1)); x = inv.InverseSlots325(s, p)
        self.assertEqual(x.nodes, [(102, 2), (102, 0), (102, 1)])
        self.assertEqual(x.node(2), (102, 0))

    def test_slot_zero_and_range_rejected(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        for code in (0, -1, 4):
            with self.subTest(code=code), self.assertRaisesRegex(v.AuditError, 'node_slot'): x.node(code)

    def test_broken_reciprocal_rejected(self):
        s, p = world(); s[1]['values']['slot'] = entry('slot', [2, 1, 3])
        with self.assertRaisesRegex(v.AuditError, 'reciprocal'): inv.InverseSlots325(s, p)

    def test_hidden_duplicate_reverse_slot_rejected(self):
        s, p = world(); s[1]['values']['slot'] = entry('slot', [1, 2, 2])
        with self.assertRaisesRegex(v.AuditError, 'reciprocal'): inv.InverseSlots325(s, p)

    def test_edge_component_and_sign_are_preserved(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        self.assertEqual(x.edge(-5), ('edge', 102, 1, 2, -1))
        self.assertNotEqual(x.edge(5), x.edge(-5)); self.assertNotEqual(x.edge(4), x.edge(5))

    def test_edge_zero_only_explicit_unused_operand(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        self.assertEqual(x.edge(0, allow_zero=True), ('unused',))
        with self.assertRaisesRegex(v.AuditError, 'zero_edge'): x.edge(0)

    def test_negative_publication_is_restore_index_not_edge_sign(self):
        s, p = world()
        p['values']['restore203_slot'] = entry('restore203_slot', [3])
        p['values']['restore203_index'] = entry('restore203_index', [0, 0, 1])
        p['values']['consumer_alias_slot'] = entry('consumer_alias_slot', [1, 2])
        s[1]['values']['key'] = entry('key', [0] * 8 + [7, 1, 1, 1])
        x = inv.InverseSlots325(s, p)
        self.assertEqual(x.publication(-1), ('restore', 102, 2))
        self.assertNotEqual(x.publication(-1), x.edge(-1))
        with self.assertRaisesRegex(v.AuditError, 'publication_restore'): x.publication(-2)

    def test_negative_scalar_alias_rejected(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        with self.assertRaisesRegex(v.AuditError, 'scalar_alias_sign'): x.alias(-1, 1, destination=True)

    def test_negative_vector_source_rejected_but_destination_preserved(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        with self.assertRaisesRegex(v.AuditError, 'source_alias_sign'): x.alias(-1, 2)
        self.assertEqual(x.alias(-1, 2, destination=True), ('edge', 102, 0, 1, -1))

    def test_nonempty_scalar_tape_rejected_in_sealed_level5_configuration(self):
        s, p = world()
        p['values']['scalar.program.shape'] = entry('scalar.program.shape', [30, 1])
        p['values']['scalar.program'] = entry('scalar.program', [0] * 30)
        with self.assertRaisesRegex(v.AuditError, 'inactive_tapes'): inv.InverseSlots325(s, p)

    def test_physical_topology_preserves_pernode_flags(self):
        s, p = world(); x = inv.InverseSlots325(s, p)
        t, q = world((2, 0, 1)); q['values']['node_flags'] = entry('node_flags', [2, 0, 1])
        self.assertEqual(x.topology(), inv.InverseSlots325(t, q).topology())
        q['values']['node_flags'] = entry('node_flags', [0, 1, 2])
        self.assertNotEqual(x.topology(), inv.InverseSlots325(t, q).topology())

    def test_consumer_value_axis_retained_in_storage_order(self):
        s, p = world((2, 0, 1))
        self.assertEqual(inv.InverseSlots325(s, p).bindings()['inverse.alias'], [(102, 2), (102, 0), (102, 1)])


    def test_local_alias_chain_operation_order_is_preserved(self):
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})
        s, p = world(); x = inv.InverseSlots325(s, p)
        left = scope('inverse.local_alias', {'source': [1, 2], 'dest': [2, 3], 'range': [1]}, epoch=1)
        right = scope('inverse.local_alias', {'source': [2, 1], 'dest': [3, 2], 'range': [1]}, epoch=1)
        self.assertNotEqual(inv.normalize_inverse_local_alias325(left, x),
                            inv.normalize_inverse_local_alias325(right, x))

    def test_changed_physical_program_operand_rejected(self):
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})
        s, p = world(); x = inv.InverseSlots325(s, p)
        left = scope('inverse.local_alias', {'source': [1], 'dest': [2], 'range': [1]}, epoch=1)
        right = scope('inverse.local_alias', {'source': [3], 'dest': [2], 'range': [1]}, epoch=1)
        self.assertNotEqual(inv.normalize_inverse_local_alias325(left, x),
                            inv.normalize_inverse_local_alias325(right, x))

    def test_unbound_source_is_failclosed_after_bad_manifest(self):
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})
        with self.assertRaisesRegex(v.AuditError, 'source_binding'):
            inv.bind_inverse_bridge_sources325({'candidate': {}})
        s, p = world(); x = inv.InverseSlots325(s, p)
        with self.assertRaisesRegex(v.AuditError, 'unbound_source'):
            inv.normalize_inverse_plan_records325(x)


class ValueEvidenceTests(unittest.TestCase):
    def test_identical_raw_digest_permuted_physical_axis_must_fail(self):
        # Bytes [1,2] in both arms with reversed physical node identities means
        # node10 has value1 versus2 despite equal SHA256.
        with self.assertRaisesRegex(v.AuditError, 'value_digest_not_permutable'):
            inv.compare_value_binding325(payload([1., 2.]), payload([1., 2.]),
                [(1, 10), (1, 20)], [(1, 20), (1, 10)], 'inverse.value')

    def test_equal_physical_axis_and_exact_digest_passes(self):
        self.assertTrue(inv.compare_value_binding325(payload([1., 2.]), payload([1., 2.]),
            [(1, 10), (1, 20)], [(1, 10), (1, 20)], 'inverse.value')['passed'])

    def test_exact_zero_bytes_prove_uniform_values_under_permutation(self):
        self.assertTrue(inv.compare_value_binding325(payload([0., 0.]), payload([0., 0.]),
            [(1, 10), (1, 20)], [(1, 20), (1, 10)], 'inverse.value')['passed'])

    def test_uniform_negative_zero_is_not_inferred(self):
        with self.assertRaisesRegex(v.AuditError, 'value_digest_not_permutable'):
            inv.compare_value_binding325(payload([-0., -0.]), payload([-0., -0.]),
                [(1, 10), (1, 20)], [(1, 20), (1, 10)], 'inverse.value')

    def test_changed_physical_value_rejected(self):
        with self.assertRaisesRegex(v.AuditError, 'raw_value_difference'):
            inv.compare_value_binding325(payload([1., 2.]), payload([1., 3.]),
                [(1, 10), (1, 20)], [(1, 10), (1, 20)], 'inverse.value')

    def test_changed_coverage_rejected_even_for_zero_payload(self):
        with self.assertRaisesRegex(v.AuditError, 'value_physical_coverage'):
            inv.compare_value_binding325(payload([0., 0.]), payload([0., 0.]),
                [(1, 10), (1, 20)], [(1, 10), (1, 30)], 'inverse.value')

    def test_duplicate_physical_axis_rejected_even_with_equal_digest(self):
        with self.assertRaisesRegex(v.AuditError, 'value_binding_extent'):
            inv.compare_value_binding325(payload([1., 1.]), payload([1., 1.]),
                [(1, 10), (1, 10)], [(1, 10), (1, 10)], 'inverse.value')


if __name__ == '__main__': unittest.main()
