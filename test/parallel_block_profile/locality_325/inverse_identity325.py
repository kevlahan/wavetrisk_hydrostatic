"""Offline inverse-slot identity audit; this module never grants qualification.

The Stage325 compact codec hashes each complete real array. A node-slot
bijection can normalize integer addresses, but does not permute such a digest.
This module exposes that distinction and fails closed on unproved value rows.
Source contract: parallel_block_inverse.f90 node_slot/edge_slot, preparation,
audit_inverse_values322, transfer_nodes and publish_native_inverse.
"""
import hashlib
import math
import verify322 as v


def insist(ok, code, detail=''):
    v.require(ok, 'inverse.identity.' + code, detail)


def values(scope, field):
    insist(field in scope['values'], 'missing_field', field)
    return v.value(scope['values'][field])


def columns(scope, field, width):
    shape = values(scope, field + '.shape')
    data = values(scope, field)
    insist(len(shape) == 2 and shape[0] == width and shape[1] >= 0 and
           len(data) == width * shape[1], 'matrix_extent', field)
    return [tuple(data[i:i + width]) for i in range(0, len(data), width)]


class InverseSlots325:
    """Checked, reciprocal mapping from buffer nodes to fixed physical nodes.

    Physical node ids are zero-based, as node_id and Domain Node_Map are in
    Fortran. Global Domain ids follow inverse.node_map epochs and are one-based.
    Vector addresses are signed, one-based EDGE=3 codes. Publication negative
    ids are restore indices, a different namespace from vector sign encoding.
    """
    def __init__(self, scopes, plan):
        insist(plan['name'] == 'inverse.plan', 'plan_scope')
        extent = values(plan, 'extent')
        insist(len(extent) == 6 and extent[0] >= 0 and extent[2:] == [2, 30, 5, 5],
               'sealed_extent', extent)
        insist(values(plan, 'options') == [1, 0], 'sealed_options')
        self.plan = plan
        self.nnode = extent[0]
        maps = [s for s in scopes if s['name'] == 'inverse.node_map' and
                plan['record'] < s['record'] < s['end_record'] < plan['end_record']]
        self.domains = [s['epoch'] for s in maps]
        insist(len(set(self.domains)) == len(self.domains), 'duplicate_domain')
        node_domain = values(plan, 'node_domain')
        node_id = values(plan, 'node_id')
        insist(len(node_domain) == len(node_id) == self.nnode, 'node_extent')
        self.nodes = []
        self.state_keys = []
        for slot, (domain, node) in enumerate(zip(node_domain, node_id), 1):
            insist(1 <= domain <= len(maps), 'local_domain', (slot, domain))
            m = maps[domain - 1]
            entries = values(m, 'slot')
            keys = columns(m, 'key', 4)
            insist(len(entries) == len(keys) and 0 <= node < len(entries),
                   'node_map_extent', (slot, domain, node))
            insist(entries[node] == slot, 'reciprocal', (slot, domain, node))
            self.nodes.append((m['epoch'], node))
            self.state_keys.append(keys[node])
        insist(len(set(self.nodes)) == self.nnode, 'duplicate_physical_node')
        for domain, m in enumerate(maps, 1):
            entries = values(m, 'slot')
            keys = columns(m, 'key', 4)
            insist(len(entries) == len(keys), 'node_map_extent')
            for node, slot in enumerate(entries):
                insist(0 <= slot <= self.nnode, 'slot_range', (domain, node, slot))
                if slot:
                    insist((node_domain[slot - 1], node_id[slot - 1]) == (domain, node),
                           'reverse_reciprocal', (domain, node, slot))
        self.restore_slots = values(plan, 'restore203_slot')
        restore_index = values(plan, 'restore203_index')
        insist(len(restore_index) == self.nnode and
               all(0 <= q <= len(self.restore_slots) for q in restore_index),
               'restore_index_extent')
        insist(len(set(self.restore_slots)) == len(self.restore_slots), 'duplicate_restore_slot')
        for index, slot in enumerate(self.restore_slots, 1):
            self.node(slot)
            insist(restore_index[slot - 1] == index, 'restore_reciprocal')
        insist(sum(q > 0 for q in restore_index) == len(self.restore_slots),
               'restore_reverse_reciprocal')
        flags = values(plan, 'node_flags')
        insist(len(flags) == self.nnode and all(0 <= flag <= 511 for flag in flags), 'node_flags_extent')
        insist(values(plan, 'consumer_alias_slot') == [i for i, key in enumerate(self.state_keys, 1) if key[0] == 0],
               'consumer_alias_inventory')
        insist(values(plan, 'scaffold_slot') == [i for i, key in enumerate(self.state_keys, 1) if key[1] == -1],
               'scaffold_inventory')
        insist(self.restore_slots == [i for i, (key, flag) in enumerate(zip(self.state_keys, flags), 1)
                                      if key[0] > 0 and flag & 3], 'restore_inventory')
        # The sealed level_start=level_end means these tapes are not merely
        # assumed irrelevant: nonempty tapes violate the audited configuration.
        insist(columns(plan, 'scalar.program', 30) == [] and
               columns(plan, 'inner.program', 104) == [], 'inactive_tapes')

    def node(self, slot):
        insist(1 <= slot <= self.nnode, 'node_slot', slot)
        return self.nodes[slot - 1]

    def edge(self, code, allow_zero=False):
        if code == 0:
            insist(allow_zero, 'zero_edge')
            return ('unused',)
        node = (abs(code) - 1) // 3 + 1
        component = (abs(code) - 1) % 3 + 1
        return ('edge', *self.node(node), component, 1 if code > 0 else -1)

    def alias(self, code, component, destination=False):
        insist(component in (1, 2), 'alias_component')
        if component == 1:
            insist(code > 0, 'scalar_alias_sign')
            return ('node', *self.node(code))
        if not destination:
            insist(code > 0, 'source_alias_sign')
        return self.edge(code)

    def publication(self, code):
        if code > 0:
            return ('workspace', *self.node(code))
        insist(code < 0 and -code <= len(self.restore_slots), 'publication_restore', code)
        return ('restore', *self.node(self.restore_slots[-code - 1]))

    def bindings(self):
        result = {'inverse.value': list(self.nodes)}
        for payload, slots in (('inverse.scaffold', 'scaffold_slot'),
                               ('inverse.alias', 'consumer_alias_slot'),
                               ('inverse.restore', 'restore203_slot')):
            indices = values(self.plan, slots)
            insist(len(set(indices)) == len(indices), 'duplicate_value_slot', slots)
            result[payload] = [self.node(slot) for slot in indices]
        return result

    def topology(self):
        flags = values(self.plan, 'node_flags')
        insist(len(flags) == self.nnode, 'flags_extent')
        return sorted((node, key, flag) for node, key, flag in
                      zip(self.nodes, self.state_keys, flags))


def all_zero_digest(count):
    """Known exact constant proof; no hash inversion or value inference."""
    insist(count >= 0, 'negative_count')
    digest = hashlib.sha256()
    chunk = bytes(8192)
    remaining = count * 8
    while remaining:
        n = min(remaining, len(chunk))
        digest.update(chunk[:n])
        remaining -= n
    return digest.digest()


def compare_value_binding325(left_scope, right_scope, left_nodes, right_nodes, field):
    """Prove physical equality only when the captured digest permits it.

    Equal byte digests with different physical bindings prove array equality,
    not per-physical-node equality. A uniform +0 buffer is independently
    reconstructible and invariant under node permutation. Other constants
    require an explicit separate source proof and are deliberately not inferred.
    """
    a = left_scope['values'][field]
    b = right_scope['values'][field]
    shape_a = values(left_scope, field + '.shape')
    shape_b = values(right_scope, field + '.shape')
    insist(a.kind == b.kind == 3 and shape_a == shape_b and
           a.count == b.count == math.prod(shape_a) and len(shape_a) in (3, 4),
           'value_extent', field)
    insist(shape_a[2] == len(left_nodes) == len(right_nodes) and
           len(set(left_nodes)) == len(left_nodes) and
           len(set(right_nodes)) == len(right_nodes), 'value_binding_extent', field)
    insist(set(left_nodes) == set(right_nodes), 'value_physical_coverage', field)
    insist(a.data == b.data, 'raw_value_difference', field)
    if list(left_nodes) == list(right_nodes):
        return {'passed': True, 'proof': 'same physical axis and exact byte digest'}
    if a.data == all_zero_digest(a.count):
        return {'passed': True, 'proof': 'exact all-positive-zero bytes'}
    raise v.AuditError('inverse.identity.value_digest_not_permutable: ' + field)


def scope_identity(scope):
    return (tuple(tuple(q) for q in scope['parent']), scope['name'], scope['epoch'])


def component_for_alias325(scopes, scope):
    local = [epoch for name, epoch in scope['parent'] if name == 'inverse.local_alias']
    if local:
        insist(len(local) == 1 and local[0] in (1, 2), 'local_alias_component')
        return local[0]
    parents = [s for s in scopes if s['name'] == 'inverse.outer_level' and
               s['record'] < scope['record'] < scope['end_record'] < s['end_record']]
    insist(len(parents) == 1, 'remote_alias_parent')
    components = [v.value(e) for e in parents[0]['records']
                  if e.name == 'component' and e.index < scope['record']]
    insist(components and components[-1] in ([1], [2]), 'remote_alias_component')
    return components[-1][0]


def integer_certificate325(slots):
    """Normalize address namespaces, preserving every arithmetic tape order.

    Sorting applies only to the physical topology census. Program rows, their
    operands, publication patch and element order, and per-node flags remain
    associated explicitly. The original value storage axes remain in bindings.
    This is an offline structural certificate, never a real-value certificate.
    """
    plan = slots.plan
    outer = []
    for op in columns(plan, 'outer.program', 18):
        level, kind = op[:2]
        insist(level == 4 and 0 <= kind <= 3, 'outer_kind_level')
        source = op[2:15]
        want = 13 if kind == 0 else 2 if kind == 1 else 4
        insist(all(source[:want]) and not any(source[want:]), 'outer_source_prefix')
        target = op[15:]
        insist(all(target) if kind == 0 else target[0] == 0 and all(target[1:]),
               'outer_targets')
        outer.append((level, kind,
                      [slots.edge(q, allow_zero=True) for q in source],
                      [slots.edge(q, allow_zero=True) for q in target]))
    publication = []
    for row in columns(plan, 'publication.program', 18):
        domain, start = row[:2]
        insist(1 <= domain <= len(slots.domains) and start >= 0, 'publication_geometry')
        publication.append((slots.domains[domain - 1], start,
                            [slots.publication(q) for q in row[2:]]))
    coefficients = plan['values']['outer.coefficients']
    insist(coefficients.kind == 3 and coefficients.count == 9 * len(outer),
           'outer_coefficients_extent')
    return {
        'physical_topology': slots.topology(),
        'operation_first': values(plan, 'operation_first'),
        'operation_last': values(plan, 'operation_last'),
        'outer_program_ordered': outer,
        'outer_coefficients': coefficients.data.hex(),
        'scalar_program_empty': True,
        'inner_program_empty': True,
        'publication_program_ordered': publication,
        'storage_bindings': slots.bindings(),
    }


def slot_consumers325(scopes, slots, selected):
    """List compiled and live buffer users of selected node slots.

    The live runtime scope census is separate from the dormant outer tape.
    Transfers are bidirectional: gather writes workspace while scatter reads it.
    Alias sources snapshot before local in-place copies; sequence is retained.
    """
    selected = set(selected)
    plan = slots.plan
    result = []
    def save(scope, role, items):
        hits = [(index, slot, slots.node(slot)) for index, slot in items if slot in selected]
        if hits:
            result.append({'scope': scope_identity(scope), 'role': role, 'uses': hits})
    def node_code(code, component):
        return abs(code) if component == 1 else (abs(code) - 1) // 3 + 1
    for scope in scopes:
        if not (plan['record'] <= scope['record'] < scope['end_record'] <= plan['end_record']):
            continue
        name = scope['name']
        if name in ('inverse.interior', 'inverse.boundary', 'inverse.level_interior',
                    'inverse.level_boundary', 'inverse.outer_install'):
            if 'plan%slot' in scope['values']:
                insist(values(scope, 'edgewise') == [0], 'sparse_outer_disabled')
                save(scope, 'gather_write_and_scatter_read_workspace',
                     list(enumerate(values(scope, 'plan%slot'), 1)))
            if 'plan%self_slot' in scope['values']:
                save(scope, 'direct_self_gather_write_and_scatter_read_workspace',
                     list(enumerate(values(scope, 'plan%self_slot'), 1)))
        elif name in ('inverse.alias', 'inverse.range_alias', 'inverse.local_alias'):
            component = scope['epoch'] if name == 'inverse.local_alias' else component_for_alias325(scopes, scope)
            for field, role in (('source', 'ordered_in_place_alias_read' if name == 'inverse.local_alias' else 'snapshot_alias_read'),
                                ('dest', 'ordered_in_place_alias_write' if name == 'inverse.local_alias' else 'ordered_alias_install_write')):
                tag = field if name == 'inverse.local_alias' else 'plan%' + field
                if tag in scope['values']:
                    save(scope, role, [(i, node_code(q, component)) for i, q in enumerate(values(scope, tag), 1)])
    for field, role in (('scaffold_slot', 'fixed_sol_seed_and_scaffold_reset'),
                        ('consumer_alias_slot', 'consumer_snapshot_and_physical_domain_publication'),
                        ('restore203_slot', 'canonical_restore_save_and_install')):
        save(plan, role, list(enumerate(values(plan, field), 1)))
    # Seed is a direct lookup of fixed Domain scaling/wavelet for uncatalogued
    # nodes; this is a source-derived relation, not observed per-node values.
    save(plan, 'fixed_domain_scaling_and_wavelet_seed',
         [(i, i) for i, key in enumerate(slots.state_keys, 1) if key[0] == 0])
    for index, op in enumerate(columns(plan, 'outer.program', 18), 1):
        source = [(i, node_code(q, 2)) for i, q in enumerate(op[2:15], 1) if q]
        target = [(i, node_code(q, 2)) for i, q in enumerate(op[15:], 1) if q]
        save(plan, 'dormant_outer_program_read_' + str(index), source)
        save(plan, 'dormant_outer_program_write_' + str(index), target[1:])
    return result


def audit_inverse_pair325(left_scopes, right_scopes, rank):
    """Bounded offline audit with precise unresolved value evidence per phase."""
    left_plans = {scope_identity(s): s for s in left_scopes if s['name'] == 'inverse.plan'}
    right_plans = {scope_identity(s): s for s in right_scopes if s['name'] == 'inverse.plan'}
    insist(left_plans.keys() == right_plans.keys(), 'plan_coverage')
    report = {'rank': rank, 'claim': 'offline inverse identity and evidence sufficiency audit',
              'qualification': False, 'plans': [], 'value_evidence': [], 'missing_evidence': []}
    def owner(scopes, scope):
        plans = [s for s in scopes if s['name'] == 'inverse.plan' and
                 s['step'] == scope['step'] and s['record'] < scope['record']]
        insist(len(plans) == 1, 'value_plan_binding')
        return plans[0]
    slot_pairs = {}
    for ident, lp in left_plans.items():
        rp = right_plans[ident]
        a, b = InverseSlots325(left_scopes, lp), InverseSlots325(right_scopes, rp)
        ca, cb = integer_certificate325(a), integer_certificate325(b)
        insist(a.topology() == b.topology(), 'physical_topology_difference')
        moved_a = [i for i, node in enumerate(a.nodes, 1) if b.nodes[i - 1] != node]
        moved_b = [i for i, node in enumerate(b.nodes, 1) if a.nodes[i - 1] != node]
        # Binding equality is reported independently, never dropped to make a
        # normalized structural certificate pass the real payload claim.
        ka = {k: x for k, x in ca.items() if k != 'storage_bindings'}
        kb = {k: x for k, x in cb.items() if k != 'storage_bindings'}
        report['plans'].append({'scope': ident, 'node_count': a.nnode,
            'ordered_integer_programs_equal': ka == kb,
            'moved_storage_slots': [{'slot': i, 'left': a.node(i), 'right': b.node(i)} for i in moved_a],
            'left_consumers': slot_consumers325(left_scopes, a, moved_a),
            'right_consumers': slot_consumers325(right_scopes, b, moved_b)})
        slot_pairs[ident] = (a, b)
    def value_table(scopes):
        seen = {}
        out = {}
        for scope in scopes:
            if scope['name'] != 'inverse.values322':
                continue
            ident = scope_identity(scope)
            seen[ident] = seen.get(ident, 0) + 1
            out[(ident, seen[ident])] = scope
        return out
    av, bv = value_table(left_scopes), value_table(right_scopes)
    insist(av.keys() == bv.keys(), 'value_phase_coverage')
    for ident, left in av.items():
        right = bv[ident]
        plan = scope_identity(owner(left_scopes, left))
        a, b = slot_pairs[plan]
        ab, bb = a.bindings(), b.bindings()
        for field in ab:
            if field not in left['values'] and field not in right['values']:
                continue
            detail = {'scope': ident, 'phase': values(left, 'inverse.phase'), 'field': field,
                      'shape': values(left, field + '.shape'),
                      'same_raw_digest': left['values'][field].data == right['values'][field].data,
                      'same_physical_axis': ab[field] == bb[field],
                      'physical_axis_differences': [{'axis_index': i, 'left': x, 'right': y}
                          for i, (x, y) in enumerate(zip(ab[field], bb[field]), 1) if x != y]}
            try:
                detail.update(compare_value_binding325(left, right, ab[field], bb[field], field))
            except v.AuditError as error:
                detail['passed'] = False
                detail['reason'] = str(error)
                report['missing_evidence'].append(detail)
            report['value_evidence'].append(detail)
    report['direct_physical_value_comparison_passed'] = not report['missing_evidence']
    if _SOURCE_BOUND325:
        bridges = [compare_inverse_bridge325(
            inverse_bridge_context325(left_scopes, lp),
            inverse_bridge_context325(right_scopes, right_plans[ident]))
            for ident, lp in left_plans.items()]
        report['source_physical_publication_bridges'] = bridges
        report['source_derived_physical_value_proof_passed'] = all(q['passed'] for q in bridges)
    report['runtime_inverse_kernel_calls'] = {
        label: {name: sum(s['name'] == name for s in scopes)
                for name in ('inverse.scalar', 'inverse.inner', 'inverse.outer')}
        for label, scopes in (('left', left_scopes), ('right', right_scopes))}
    return report


EXPECTED_SOURCE_HASHES325 = {
    'src/parallel_block_inverse.f90': '65b936677196c97fabf0f4ec8363132c3e48e9e4ff1f17ba4e789d7781d60db0',
    'src/parallel_block_mpi.f90': 'c915cc2b6529b38ec7a74d41afb74bdb5caddd8c7731cb15ed76c8ea89755027',
    'src/domain.f90': '54a09b07a979ba54d9695ea23d1b6ea76096893cf1106f9c73acc49f05f27ecd',
    'src/shared.f90': 'bd203af9a687cc2831e4df61fd802211e37b60a2517b8b91bb933ef479b7ca84',
}
_SOURCE_BOUND325 = False


def bind_inverse_bridge_sources325(manifest, source_root=None):
    """Pin the exact audited producer/reader and fixed-Domain capture code.

    The campaign separately binds this manifest to its original compiled inputs.
    Passing source_root also verifies actual bytes; worker processes must call
    this function themselves rather than assuming parent process global state.
    """
    from pathlib import Path
    global _SOURCE_BOUND325
    _SOURCE_BOUND325 = False
    candidate = manifest.get('candidate', {})
    for relative, digest in EXPECTED_SOURCE_HASHES325.items():
        insist(candidate.get(relative) == digest, 'source_binding', relative)
        if source_root is not None:
            insist(hashlib.sha256((Path(source_root) / relative).read_bytes()).hexdigest() == digest,
                   'source_bytes', relative)
    _SOURCE_BOUND325 = True
    return {'passed': True, 'sources': dict(EXPECTED_SOURCE_HASHES325)}


def _record_identity325(e):
    return [e.name, e.kind, e.count, e.data.hex() if e.kind == 3 else v.value(e)]


def _publication_physical325(scopes, publication, slots):
    matches = [s for s in scopes if s['name'] == 'domain.fields322' and
               publication['record'] < s['record'] < s['end_record'] < publication['end_record']]
    insist(len(matches) == 1, 'publication_fields')
    fields = matches[0]
    insist(values(fields, 'fields.header') == [31, 3, 30, len(slots.domains)],
           'publication_field_header')
    domains = [v.value(e) for e in fields['records'] if e.name == 'fields.domain']
    insist(len(domains) == len(slots.domains) and
           [(row[0], row[1] + 1) for row in domains] == list(enumerate(slots.domains, 1)),
           'publication_domain_identity')
    lengths = {row[0]: row[2] for row in domains}
    coordinates = [v.value(e) for e in fields['records'] if e.name == 'fields.coordinate']
    expected = [(d, k, variable) for d in range(1, len(slots.domains) + 1)
                for k in range(1, 31) for variable in range(1, 4)]
    insist([tuple(row[:3]) for row in coordinates] == expected, 'publication_field_coordinates')
    for d, k, variable, position, tag, count in coordinates:
        insist(position == (2 if variable == 1 else 1) and
               count == lengths[d] * (3 if variable == 1 else 1), 'publication_field_extent')
    for domain, node in slots.nodes:
        d = slots.domains.index(domain) + 1
        insist(0 <= node < lengths[d], 'publication_node_coverage')
    real_records = [e for e in fields['records'] if e.name == 'fields.values']
    insist(len(real_records) == len(coordinates) and all(e.kind == 3 and e.count == row[-1]
            for e, row in zip(real_records, coordinates)), 'publication_real_coverage')
    return [_record_identity325(e) for e in fields['records']]


def inverse_bridge_context325(scopes, plan):
    """Prove the exact optimized final-inverse schedule used by the bridge.

    This deliberately rejects checked/oracle trajectories and other level
    ranges. All real hashes remain required. It supplies an additional direct
    physical publication proof for changed storage axes, not a replacement for
    the captured numeric comparison or ordered provider proof.
    """
    insist(_SOURCE_BOUND325, 'unbound_source')
    slots = InverseSlots325(scopes, plan)
    step = plan['step']
    sequence = [s for s in scopes if s['step'] == step and s['name'].startswith('inverse.') and
                s['name'] != 'inverse.plan' and not any(x[0] == 'inverse.plan' for x in s['parent'])]
    expected = [('inverse.prepare', 'first_level', [5])]
    expected += [('inverse.gather', 'family', [component, family])
                 for component in (1, 2) for family in (1, 2)]
    for component in (1, 2):
        expected += [('inverse.boundary', 'family_and_range', [component, 2, 5, 5]),
                     ('inverse.boundary', 'family_and_range', [component, 1, 5, 5]),
                     ('inverse.scatter', 'family', [component, 2]),
                     ('inverse.scatter', 'family', [component, 1])]
    expected += [('inverse.values322', 'inverse.phase', [1, slots.nnode, 2, 30]),
                 ('inverse.begin', 'options', [1, 0, 0]),
                 ('inverse.values322', 'inverse.phase', [2, slots.nnode, 2, 30]),
                 ('inverse.gather', 'family', [1, 1]),
                 ('inverse.boundary', 'family_and_range', [1, 1, 5, 5]),
                 ('inverse.scatter', 'family', [1, 1]),
                 ('inverse.publish322', 'publication.mode', [1, 1, 0]),
                 ('inverse.publish322', 'publication.mode', [1, 2, 0]),
                 ('inverse.publish322', 'publication.mode', [2, 1, 0]),
                 ('inverse.end', 'work', [0, 0, 0, 0, 0, 0]),
                 ('inverse.values322', 'inverse.phase', [3, slots.nnode, 2, 30]),
                 ('inverse.publish322', 'publication.mode', [1, 1, 1])]
    insist([s['name'] for s in sequence] == [name for name, field, want in expected],
           'optimized_final_schedule', [s['name'] for s in sequence])
    for s, (name, field, want) in zip(sequence, expected):
        insist(s['epoch'] == plan['epoch'] and values(s, field) == want, 'optimized_scope_binding', name)
        if name in ('inverse.gather', 'inverse.scatter'):
            before_phase1 = s['record'] < sequence[13]['record']
            flags = [0, 0, 0] if before_phase1 else [1, 1, 1] if name == 'inverse.gather' else [0, 1, 1]
            insist(values(s, 'flags') == flags, 'optimized_transfer_flags')
            insist(('level' not in s['values']) if before_phase1 else values(s, 'level') == [5],
                   'optimized_transfer_level')
    phase1 = sequence[13]
    begin, phase2 = sequence[14:16]
    # The two scopes immediately following phase1 have no intervening records
    # except begin's audited options. This closes the no-mutation phase1->2 edge.
    insist(begin['record'] == phase1['end_record'] + 1 and
           phase2['record'] == begin['record'] + 2 and phase2['end_record'] + 1 == begin['end_record'],
           'phase2_adjacency')
    pubs = {tuple(values(s, 'publication.mode')): s for s in sequence if s['name'] == 'inverse.publish322'}
    publication = {str(mode): _publication_physical325(scopes, pubs[mode], slots)
                   for mode in ((1, 1, 0), (1, 2, 0), (1, 1, 1))}
    phases = {}
    for s in sequence:
        if s['name'] == 'inverse.values322':
            number = values(s, 'inverse.phase')[0]
            bindings = slots.bindings()
            for field, nodes in bindings.items():
                shape = [5, 30, len(nodes)] + ([2] if field == 'inverse.value' else [])
                insist(field in s['values'] and values(s, field + '.shape') == shape and
                       s['values'][field].kind == 3 and s['values'][field].count == math.prod(shape),
                       'bridge_numeric_extent', field)
            phases[str(number)] = {field: {'shape': values(s, field + '.shape'),
                'kind': s['values'][field].kind, 'count': s['values'][field].count,
                'digest': s['values'][field].data.hex()} for field in slots.bindings()
                if field in s['values']}
    return {'nodes': slots.nodes, 'keys': slots.state_keys, 'bindings': slots.bindings(),
            'publication': publication, 'phases': phases,
            'source_contract': dict(EXPECTED_SOURCE_HASHES325),
            'step': step, 'epoch': plan['epoch'], 'schedule_verified': True}


def compare_inverse_bridge325(a, b):
    """Physical proof for moved rows, with exact raw numeric guards retained."""
    insist(_SOURCE_BOUND325 and a['source_contract'] == b['source_contract'] == EXPECTED_SOURCE_HASHES325,
           'bridge_source_contract')
    insist(a['schedule_verified'] is True and b['schedule_verified'] is True and
           (a['step'], a['epoch']) == (b['step'], b['epoch']), 'bridge_schedule_binding')
    an, bn = [tuple(q) for q in a['nodes']], [tuple(q) for q in b['nodes']]
    insist(len(an) == len(bn) and len(set(an)) == len(an) and set(an) == set(bn), 'bridge_node_coverage')
    changed = {node for i, node in enumerate(an) if bn[i] != node}
    changed |= {node for i, node in enumerate(bn) if an[i] != node}
    for context in (a, b):
        nodes = [tuple(q) for q in context['nodes']]
        insist(len(context['keys']) == len(nodes), 'bridge_key_extent')
        alias = {tuple(q) for q in context['bindings']['inverse.alias']}
        scaffold = {tuple(q) for q in context['bindings']['inverse.scaffold']}
        restore = {tuple(q) for q in context['bindings']['inverse.restore']}
        insist(all(context['keys'][nodes.index(node)][0] == 0 for node in changed), 'bridge_moved_catalogued')
        insist(changed <= alias and not changed & (scaffold | restore), 'bridge_moved_auxiliary')
        for field in ('inverse.scaffold', 'inverse.restore'):
            insist(context['bindings'][field] == a['bindings'][field], 'bridge_auxiliary_axis')
    # Each audited publication serializes all initialized fixed-Domain fields
    # including boundaries. The reciprocal Node_Map makes the projection of
    # every moved physical node injective. This is direct physical byte proof.
    insist(a['publication'].keys() == b['publication'].keys() and
           a['publication'] == b['publication'], 'bridge_physical_publication_difference')
    insist(a['phases'].keys() == b['phases'].keys() == {'1', '2', '3'}, 'bridge_phase_coverage')
    for phase in ('1', '2', '3'):
        insist(set(a['phases'][phase]) == set(b['phases'][phase]) ==
               {'inverse.value', 'inverse.alias', 'inverse.scaffold', 'inverse.restore'},
               'bridge_value_field_coverage', phase)
        insist(a['phases'][phase] == b['phases'][phase], 'bridge_raw_value_difference', phase)
    return {'passed': True, 'moved_physical_nodes': len(changed), 'step': a['step'],
            'proof': 'exact fixed-physical publication bridge plus unchanged raw numeric hashes',
            'phase1_and_2_family1': 'immutable consumer alias snapshot -> physical consumer publication',
            'family2_all_phases': 'no family2 writer after phase1 -> physical wavelet publication',
            'phase3_family1': 'physical solution publication -> end without uncatalogued writes'}


def require_inverse_source_binding325():
    insist(_SOURCE_BOUND325, 'unbound_source')


def normalize_inverse_plan_records325(slots):
    """Physical integer metadata; all numeric hashes and shape guards retained."""
    require_inverse_source_binding325()
    plan = slots.plan
    certificate = integer_certificate325(slots)
    bindings = slots.bindings()
    restore_index = values(plan, 'restore203_index')
    normalized = {
        'node_domain': sorted(slots.nodes),
        'node_id': sorted(slots.nodes),
        'node_flags': [(node, flag) for node, key, flag in slots.topology()],
        'restore203_index': sorted((node, None if index == 0 else slots.node(slots.restore_slots[index - 1]))
                                    for node, index in zip(slots.nodes, restore_index)),
        'outer.program': certificate['outer_program_ordered'],
        'publication.program': certificate['publication_program_ordered'],
        # These loops copy each separate physical target independently. Their
        # storage order is recorded in the separate strict value bridge context.
        'scaffold_slot': sorted(bindings['inverse.scaffold']),
        'consumer_alias_slot': sorted(bindings['inverse.alias']),
        'restore203_slot': sorted(bindings['inverse.restore']),
    }
    exclude = {'transfer_send.shape', 'transfer_recv.shape'}
    result = []
    for e in plan['records']:
        if e.name.startswith('plan%') or e.name in exclude:
            continue
        item = _record_identity325(e)
        if e.name in normalized:
            item[3] = normalized[e.name]
        result.append(item)
    return result


def normalize_inverse_node_map325(scope, slots):
    require_inverse_source_binding325()
    result = []
    for e in scope['records']:
        item = _record_identity325(e)
        if e.name == 'slot':
            item[3] = [None if index == 0 else slots.node(index) for index in v.value(e)]
        result.append(item)
    return result


def normalize_inverse_local_alias325(scope, slots):
    require_inverse_source_binding325()
    component = scope['epoch']
    result = []
    for e in scope['records']:
        item = _record_identity325(e)
        if e.name in ('source', 'dest'):
            item[3] = [slots.alias(q, component, destination=e.name == 'dest') for q in v.value(e)]
        result.append(item)
    # No grouping or sorting: local alias copies can read earlier writes to a
    # different target, and the exact complete operation sequence is semantic.
    return result


def physical_inverse_node325(context, slot):
    nodes = context.get('inverse_nodes')
    insist(nodes is not None and 1 <= slot <= len(nodes), 'provider_node_slot', slot)
    return tuple(nodes[slot - 1])


def physical_inverse_alias325(context, code, destination=False):
    component = context.get('alias_component')
    insist(component in (1, 2), 'provider_alias_component')
    if component == 1:
        insist(code > 0, 'provider_scalar_sign')
        return ('node', *physical_inverse_node325(context, code))
    insist(code != 0 and (destination or code > 0), 'provider_vector_sign')
    return ('edge', *physical_inverse_node325(context, (abs(code) - 1) // 3 + 1),
            (abs(code) - 1) % 3 + 1, 1 if code > 0 else -1)
