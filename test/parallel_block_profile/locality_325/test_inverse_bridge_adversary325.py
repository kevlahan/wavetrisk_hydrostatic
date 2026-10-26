"""Independent counterexamples for the source-bound physical value bridge.

These small complete schedules exercise the actual context builder and pair
checker. The fixtures supply real-byte hashes and fixed physical publication
coordinates; they do not grant native qualification.
"""
import hashlib
import struct
import unittest

import inverse_identity325 as inv
import verify322 as v
from test_inverse_identity325 import world
from test_verify325 import entry, scope


def append_record(s, e):
    s['records'].append(e)
    s['values'][e.name] = e
    s['tags'].add(e.name)
    s['record_counts'][e.name] += 1


def digest_entry(name, count, number=1.):
    return v.Entry(0, 0, name, 3, count,
                   hashlib.sha256(struct.pack('<d', number) * count).digest())


def schedule(order=(0, 1, 2)):
    scopes, plan = world(order)
    plan['step'] = 1
    cursor = 101
    parent = [('rk.final322', 1)]

    def add(name, fields, owner=parent):
        nonlocal cursor
        s = scope(name, fields, record=cursor, end=cursor + len(fields) + 1,
                  parent=owner)
        cursor = s['end_record'] + 1
        scopes.append(s)
        return s

    def payload(phase, owner=parent):
        nonlocal cursor
        s = add('inverse.values322', {'inverse.phase': [phase, 3, 2, 30]}, owner)
        for name, shape in (('inverse.value', [5, 30, 3, 2]),
                            ('inverse.alias', [5, 30, 3]),
                            ('inverse.scaffold', [5, 30, 0]),
                            ('inverse.restore', [5, 30, 0])):
            append_record(s, entry(name + '.shape', shape))
            count = 1
            for extent in shape:
                count *= extent
            append_record(s, digest_entry(name, count))
        s['end_record'] = s['record'] + len(s['records']) + 1
        cursor = s['end_record'] + 1
        return s

    def publication(mode):
        nonlocal cursor
        p = add('inverse.publish322', {'publication.mode': list(mode)})
        f = add('domain.fields322', {'fields.header': [31, 3, 30, 1]},
                owner=[('inverse.publish322', 1)])
        append_record(f, entry('fields.domain', [1, 101, 3, 0, 0]))
        for k in range(1, 31):
            for variable in range(1, 4):
                count = 9 if variable == 1 else 3
                append_record(f, entry('fields.coordinate',
                              [1, k, variable, 2 if variable == 1 else 1, 0, count]))
                append_record(f, digest_entry('fields.values', count))
        f['end_record'] = f['record'] + len(f['records']) + 1
        p['end_record'] = f['end_record'] + 1
        cursor = p['end_record'] + 1
        return p, f

    add('inverse.prepare', {'first_level': [5]})
    for component in (1, 2):
        for family in (1, 2):
            add('inverse.gather', {'family': [component, family], 'flags': [0, 0, 0]})
    for component in (1, 2):
        add('inverse.boundary', {'family_and_range': [component, 2, 5, 5]})
        add('inverse.boundary', {'family_and_range': [component, 1, 5, 5]})
        for family in (2, 1):
            add('inverse.scatter', {'family': [component, family], 'flags': [0, 0, 0]})
    payload(1)
    begin = add('inverse.begin', {'options': [1, 0, 0]})
    cursor = begin['record'] + 2
    phase2 = payload(2, [('inverse.begin', 1)])
    begin['end_record'] = phase2['end_record'] + 1
    cursor = begin['end_record'] + 1
    add('inverse.gather', {'family': [1, 1], 'flags': [1, 1, 1], 'level': [5]})
    add('inverse.boundary', {'family_and_range': [1, 1, 5, 5]})
    add('inverse.scatter', {'family': [1, 1], 'flags': [0, 1, 1], 'level': [5]})
    for mode in ((1, 1, 0), (1, 2, 0), (2, 1, 0)):
        publication(mode)
    end = add('inverse.end', {'work': [0] * 6})
    cursor = end['record'] + 2
    phase3 = payload(3, [('inverse.end', 1)])
    end['end_record'] = phase3['end_record'] + 1
    cursor = end['end_record'] + 1
    publication((1, 1, 1))
    return scopes, plan


class IndependentBridgeTests(unittest.TestCase):
    def setUp(self):
        inv.bind_inverse_bridge_sources325({'candidate': dict(inv.EXPECTED_SOURCE_HASHES325)})

    def contexts(self):
        a, ap = schedule()
        b, bp = schedule((0, 2, 1))
        return (inv.inverse_bridge_context325(a, ap),
                inv.inverse_bridge_context325(b, bp))

    def test_uniform_values_with_genuine_storage_permutation_pass(self):
        a, b = self.contexts()
        result = inv.compare_inverse_bridge325(a, b)
        self.assertTrue(result['passed'])
        self.assertEqual(result['moved_physical_nodes'], 2)

    def test_equal_raw_hashes_cannot_hide_wrong_physical_publication(self):
        a, b = self.contexts()
        self.assertEqual(a['phases'], b['phases'])
        b['publication']['(1, 1, 1)'][-1][-1] = 'different physical bytes'
        with self.assertRaisesRegex(v.AuditError, 'physical_publication_difference'):
            inv.compare_inverse_bridge325(a, b)

    def test_every_bridge_publication_is_required(self):
        for mode in ('(1, 1, 0)', '(1, 2, 0)', '(1, 1, 1)'):
            with self.subTest(mode=mode):
                a, b = self.contexts()
                del b['publication'][mode]
                with self.assertRaisesRegex(v.AuditError, 'physical_publication_difference'):
                    inv.compare_inverse_bridge325(a, b)

    def test_raw_hash_guard_retained_for_every_phase(self):
        for phase in ('1', '2', '3'):
            with self.subTest(phase=phase):
                a, b = self.contexts()
                b['phases'][phase]['inverse.value']['digest'] = 'changed raw bits'
                with self.assertRaisesRegex(v.AuditError, 'raw_value_difference'):
                    inv.compare_inverse_bridge325(a, b)

    def test_moved_catalogued_node_rejected(self):
        a, b = self.contexts()
        b['keys'][1] = [17, 1, 1, 0]
        with self.assertRaisesRegex(v.AuditError, 'moved_catalogued'):
            inv.compare_inverse_bridge325(a, b)

    def test_missing_consumer_snapshot_node_rejected(self):
        a, b = self.contexts()
        b['bindings']['inverse.alias'] = b['bindings']['inverse.alias'][:1]
        with self.assertRaisesRegex(v.AuditError, 'moved_auxiliary'):
            inv.compare_inverse_bridge325(a, b)

    def test_moved_scaffold_and_restore_nodes_rejected(self):
        for field in ('inverse.scaffold', 'inverse.restore'):
            with self.subTest(field=field):
                a, b = self.contexts()
                b['bindings'][field] = [b['nodes'][1]]
                with self.assertRaisesRegex(v.AuditError, 'moved_auxiliary'):
                    inv.compare_inverse_bridge325(a, b)

    def test_changed_step_and_source_contract_rejected(self):
        for field, value, code in (('step', 2, 'schedule_binding'),
                                   ('source_contract', {}, 'source_contract'),
                                   ('schedule_verified', False, 'schedule_binding')):
            with self.subTest(field=field):
                a, b = self.contexts()
                b[field] = value
                with self.assertRaisesRegex(v.AuditError, code):
                    inv.compare_inverse_bridge325(a, b)

    def test_inserted_inverse_arithmetic_rejected(self):
        scopes, plan = schedule()
        scopes.insert(-1, scope('inverse.scalar', {'level_phase': [4, 1]},
                               record=1000, end=1003))
        with self.assertRaisesRegex(v.AuditError, 'optimized_final_schedule'):
            inv.inverse_bridge_context325(scopes, plan)

    def test_family2_write_after_phase1_rejected(self):
        scopes, plan = schedule()
        late = [s for s in scopes if s['name'] == 'inverse.gather'][-1]
        late['values']['family'] = entry('family', [1, 2])
        with self.assertRaisesRegex(v.AuditError, 'optimized_scope_binding'):
            inv.inverse_bridge_context325(scopes, plan)

    def test_intervening_phase2_record_rejected(self):
        scopes, plan = schedule()
        phase2 = [s for s in scopes if s['name'] == 'inverse.values322'][1]
        phase2['record'] += 1
        with self.assertRaisesRegex(v.AuditError, 'phase2_adjacency'):
            inv.inverse_bridge_context325(scopes, plan)

    def test_wrong_publication_family_and_missing_publication_rejected(self):
        for remove in (False, True):
            with self.subTest(remove=remove):
                scopes, plan = schedule()
                publication = [s for s in scopes if s['name'] == 'inverse.publish322'][1]
                if remove:
                    scopes.remove(publication)
                else:
                    publication['values']['publication.mode'] = entry('publication.mode', [1, 1, 0])
                with self.assertRaisesRegex(v.AuditError, 'optimized_'):
                    inv.inverse_bridge_context325(scopes, plan)

    def test_missing_numeric_payload_rejected_at_each_phase(self):
        for phase in range(3):
            with self.subTest(phase=phase + 1):
                scopes, plan = schedule()
                s = [s for s in scopes if s['name'] == 'inverse.values322'][phase]
                del s['values']['inverse.value']
                with self.assertRaisesRegex(v.AuditError, 'bridge_numeric_extent'):
                    inv.inverse_bridge_context325(scopes, plan)

    def test_wrong_numeric_shape_or_count_rejected(self):
        for wrong_count in (False, True):
            with self.subTest(wrong_count=wrong_count):
                scopes, plan = schedule()
                s = [s for s in scopes if s['name'] == 'inverse.values322'][0]
                if wrong_count:
                    s['values']['inverse.value'] = digest_entry('inverse.value', 899)
                else:
                    s['values']['inverse.value.shape'] = entry('inverse.value.shape', [5, 30, 3, 1])
                with self.assertRaisesRegex(v.AuditError, 'bridge_numeric_extent'):
                    inv.inverse_bridge_context325(scopes, plan)

    def test_wrong_fixed_domain_identity_and_component_extent_rejected(self):
        for field, replacement, code in (
                ('fields.domain', [1, 99, 3, 0, 0], 'domain_identity'),
                ('fields.coordinate', [1, 1, 1, 1, 0, 9], 'field_extent')):
            with self.subTest(field=field):
                scopes, plan = schedule()
                pub = [s for s in scopes if s['name'] == 'domain.fields322'][0]
                index = next(i for i, e in enumerate(pub['records']) if e.name == field)
                pub['records'][index] = entry(field, replacement)
                with self.assertRaisesRegex(v.AuditError, code):
                    inv.inverse_bridge_context325(scopes, plan)


class AliasOverwriteAdversaryTests(unittest.TestCase):
    def test_opposite_sign_writes_share_one_destination_and_preserve_order(self):
        # +4 and -4 are the same workspace edge, with opposite write signs.
        # Grouping by the signed code would wrongly erase the last-writer
        # distinction. Use the real SQLite provider path and reconciler.
        from test_provider_paths325 import inverses, certify
        plans, contexts, unused_catalogue = inverses('inverse.alias', ordered=True)
        for rank, context in enumerate(contexts):
            context['alias_component'] = 2
            context['inverse_nodes'] = [[rank + 1, node] for node in range(3)]
        plans[1]['plan%dest'] = [4, -4]
        plans[1]['plan%install_order'] = [1, 2]
        forward = certify('static', 'inverse.alias', plans, contexts)[0]
        plans[1]['plan%install_order'] = [2, 1]
        reversed_writes = certify('static', 'inverse.alias', plans, contexts)[0]
        self.assertNotEqual(forward, reversed_writes)


if __name__ == '__main__':
    unittest.main()
