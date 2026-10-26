#!/usr/bin/env python3
"""Independent bounded-evidence review. Never imports or executes captured code.

The compact review cannot independently rerun the 1024 physical comparisons:
their source-bound checker and hash-bound native databases are explicit premises.
This script verifies their evidence chain, not merely their status fields.
"""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tarfile

PACKAGE = '1a47b001cad28d82de231623f3c85e7b271e744369f9549d8a7b6b30d56dd196'
APPROVED_SOURCE_ARCHIVE = '99fe77316934f2991f65ecf2e0c22ee15051f490ccb8f38a27c9fb78bc9f86b7'
PACKAGE_FILES = 197
EXPECTED_TESTS = 292
FROZEN = 'ff430542d8ab81078d522a2468a1ce1baa22bac501f922e534c839c7c47ee1bb'
FAMILIES = '25d18c76eec08586b05ef3d7238e062e51f66bedc46d1347c3885d5bad327b8c'
ROWS = dict(route=65024, domain=256, context=66048, feature=508840)
PHASES = {'cache': (), 'providers-C0': ('cache',), 'providers-L256': ('cache',),
          'pair': ('cache', 'providers-C0', 'providers-L256'),
          'work': ('pair',), 'preservation': ('work',)}
ARMS = ('qual-C0', 'qual-L256')
LABELS = ('checked-L256', *ARMS)
SOURCE = {
 'src/parallel_block_inverse.f90': '65b936677196c97fabf0f4ec8363132c3e48e9e4ff1f17ba4e789d7781d60db0',
 'src/parallel_block_mpi.f90': 'c915cc2b6529b38ec7a74d41afb74bdb5caddd8c7731cb15ed76c8ea89755027',
 'src/domain.f90': '54a09b07a979ba54d9695ea23d1b6ea76096893cf1106f9c73acc49f05f27ecd',
 'src/shared.f90': 'bd203af9a687cc2831e4df61fd802211e37b60a2517b8b91bb933ef479b7ca84'}
READER = '5feaf0dca377947d78e7840dfb3b167eb95e03234ba59cb9175ba0a503179ae1'
FROZEN_NAMES = ('recover325.py', 'recover_parallel325.py', 'binding325.json',
                'parallel-history325.json', 'recovery-history325.json',
                'prior-native-input-anchors.json')


class ReviewError(ValueError):
    pass


def insist(ok, message):
    if not ok:
        raise ReviewError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False).encode('ascii') + b'\n'


def digest(value):
    return sha(canonical(value))


def decode(data):
    def pairs(items):
        result = {}
        for key, value in items:
            insist(key not in result, 'duplicate JSON key: ' + key)
            result[key] = value
        return result
    return json.loads(data, object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ReviewError('nonfinite JSON')))


def hash_value(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None


def read_review(path):
    data, total = {}, 0
    with tarfile.open(path, 'r|gz') as archive:
        for member in archive:
            insist(member.isfile() and member.name.startswith('review/') and
                   '..' not in Path(member.name).parts and not Path(member.name).is_absolute(),
                   'unsafe review member')
            name = member.name[len('review/'):]
            insist(name not in data and member.size <= 2*1024**2, 'duplicate/oversized review member')
            total += member.size
            insist(total <= 32*1024**2 and len(data) < 10000, 'review cap exceeded')
            data[name] = archive.extractfile(member).read()
    manifest = decode(data['files-sha256.json'])
    insist(set(data) == set(manifest) | {'files-sha256.json'}, 'review inventory differs')
    insist(all(hash_value(h) and sha(data[n]) == h for n, h in manifest.items()), 'review file hash differs')
    return data, total


def empty_failures(value):
    insist(value['count'] == 0 and value['counts'] == {} and value['examples'] == []
           and value['sha256'] == sha(b''), 'nonempty mismatch census')


def scope_zero(value):
    for key in ('solver_runs', 'performance_runs'):
        if key in value:
            insist(value[key] == 0, 'unexpected run scope: ' + key)
    for key in ('compiler_runs', 'builds'):
        if key in value:
            insist(value[key] == 0, 'unexpected compilation scope: ' + key)


def check_source(value):
    insist(value == {'passed': True, 'sources': SOURCE}, 'inverse source proof differs')


def check_reuse(value, kit):
    h = decode(kit['recovery_v4p1/parallel-history325.json'])
    semantic_data = kit['recovery_v4p1/predecessor-v3/semantic-replay.json']
    insist(sha(semantic_data) == h['certificate_sha256']['predecessor-v3/semantic-replay.json'],
           'checked reuse report hash differs')
    earlier = decode(semantic_data)['runs']['checked-L256']
    insist(earlier['status'] == 'complete' and earlier['families'] == 254 and
           [x['rank'] for x in earlier['streams']] == list(range(256)), 'checked closure incomplete')
    insist(value == dict(passed=True, source_job=55573, streams=256, families=254,
           semantic_sha256=sha(semantic_data), actual_read_hashes=256,
           scope='Earlier checked closure only; optimized proof recomputed.'), 'checked reuse widened or changed')


def original_science(root, binding):
    manifest_data = (root/'files-sha256.json').read_bytes()
    insist(sha(manifest_data) == binding['metadata_return_manifest_sha256'], 'original review manifest differs')
    manifest = decode(manifest_data)
    def original(name):
        raw = (root/name).read_bytes()
        insist(sha(raw) == manifest[name], 'original review file differs: ' + name)
        return decode(raw)
    result = original('results/results.json')
    builds = original('results/builds.json')
    source_data = (root/'kit/source-files.json').read_bytes()
    insist(sha(source_data) == binding['source_manifest_sha256'] == builds['source_manifest_sha256'],
           'original compiled source manifest differs')
    insist(sum(map(len, decode(source_data).values())) == 258, 'original source census')
    insist(builds['status'] == 'passed' and builds['toolchain']['status'] == 'passed'
           and builds['toolchain']['expected']['gcc_version'] == '13.2.0'
           and builds['toolchain']['expected']['mpi_version'] == '4.0', 'original exact toolchain gate')
    insist(result['job_id'] == '55569' or result['job_id'] == 55569, 'original science job differs')
    insist(result['status'] == 'failed' and result['error'] == binding['expected_failure']
           and result['schedule'] == list(LABELS) and result['performance_measurement'] is False,
           'original failure ledger/scope changed')
    for label in LABELS:
        run = result['runs'][label]
        insist(run['status'] == 'passed' and run['ranks'] == 256 and run['steps'] == 4,
               'original native solver record incomplete')
        insist(run['solver_placement']['passed'] and run['solver_placement']['observed_ranks'] == 256,
               'original placement incomplete')
        for name, gate in run['science'].items():
            insist(gate['passed'] is True and gate['failures'] == [], 'original science gate failed')
            if name.endswith('exact'):
                insist(gate['sha256_exact'] is True, 'original exact science gate not exact')
            report = gate.get('report')
            if report:
                insist(report['header_differences'] == report['topology_differences'] ==
                       report['threshold_nonfinite'] == 0, 'original science topology/nonfinite difference')
                insist(all(x['nonfinite'] == 0 and math.isfinite(x['max_abs'])
                           for x in report['fields'].values()), 'original science nonfinite fields')
        insist(run['binary_sha256'] == builds['builds'][run['kind']]['binary_sha256'],
               'original run/build binary mismatch')
    insist(result['runs']['qual-C0']['checkpoints'] == result['runs']['qual-L256']['checkpoints'],
           'optimized checkpoints not exact')
    return result, builds


def review(archive, job, post_path, original_root, expected_package=PACKAGE, expected_tests=EXPECTED_TESTS):
    data, size = read_review(archive)
    kit = {name[4:]: raw for name, raw in data.items() if name.startswith('kit/')}
    insist(hash_value(expected_package) and sha(kit['native-files-sha256.json']) == expected_package,
           'wrong sealed native package')
    package = decode(kit['native-files-sha256.json'])
    insist(len(kit) == PACKAGE_FILES and set(kit) == set(package) | {'native-files-sha256.json'}, 'sealed package inventory')
    insist(all(sha(kit[n]) == h for n, h in package.items()), 'sealed package bytes differ')
    insist(sha(kit['recovery_v4p1/files-sha256.json']) == FROZEN, 'wrong frozen historical payload')
    frozen = decode(kit['recovery_v4p1/files-sha256.json'])
    insist(all(sha(kit['recovery_v4p1/'+n]) == h for n, h in frozen.items()), 'frozen payload changed')
    submission = decode(data['submission.json'])
    insist(submission['job'] == job and submission['package_manifest_sha256'] == expected_package,
           'submission identity differs')
    scope_zero(submission)
    preflight_log = data['work/preflight-tests.log'].decode()
    insist(re.search(r'Ran '+str(expected_tests)+r' tests in [0-9.]+s', preflight_log) and
           preflight_log.rstrip().endswith('OK') and '... skipped ' not in preflight_log,
           'native preflight count/zero-skip condition not met')
    insist(sha(data['work/preflight-tests.log']) == submission['preflight']['log_sha256'],
           'native preflight receipt differs')
    account = decode(data['accounting.json'])
    rows = [line.split('|') for line in account['stdout'].splitlines() if line.split('|')[0] == str(job)]
    insist(account['returncode'] == 0 and len(rows) == 1, 'scheduler accounting missing/ambiguous')
    row = rows[0]
    insist(row[1] == 'COMPLETED' and row[2] == '0:0' and row[7] == 'bb13', 'native job did not complete cleanly')
    elapsed = int(row[6]); insist(0 < elapsed <= 1800, 'new allocation exceeded its declared budget')
    insist(post_path is not None, 'independent post-terminal native hashes are required')
    post = decode(Path(post_path).read_bytes())
    insist(post['passed'] is True, 'post-terminal read-only audit did not pass')
    output_identity = post['output_identity']
    canonical_out = Path(output_identity['canonical_path'])
    insist(output_identity['lexical_path'] == submission['out'] and
           canonical_out.is_absolute() and '..' not in canonical_out.parts and
           output_identity['samefile'] is True and
           output_identity['lexical_stat'] == output_identity['canonical_stat'] and
           set(output_identity['canonical_stat']) == {'device', 'inode'} and
           all(type(x) is int and x >= 0 for x in output_identity['canonical_stat'].values()) and
           output_identity['canonical_stat']['inode'] > 0,
           'native output lexical/canonical identity lacks fresh samefile/stat proof')
    get = lambda name: decode(data['work/out/'+name])
    supervisor, attempt = get('supervisor-report.json'), get('state.attempt.json')
    insist(supervisor['status'] == 'passed' and supervisor['gate_passed'] is True and
           supervisor['child_returncode'] == 0 and supervisor['signals'] == [] and
           supervisor['remaining_live_child_group_pids'] == [] and
           supervisor.get('cleanup_wait_timed_out') is False, 'supervisor/cleanup did not pass')
    insist(supervisor['prior_elapsed_seconds'] == 2656 and supervisor['new_allocation_cap_seconds'] == 1800
           and supervisor['cumulative_maximum_seconds'] == 4456, 'budget history was rewritten')
    insist(supervisor['sampled_peak_derived_bytes'] <= 4*1024**3 and
           supervisor['sampled_peak_virtual_bytes'] <= 16*1024**3, 'recorded resource guard exceeded')
    insist(supervisor['package_manifest_sha256'] == expected_package, 'native executed package differs')
    environment = supervisor['environment']
    insist(environment['job'] == job and environment['hostname'].split('.')[0] == 'bb13' and
           len(set(environment['cpus'])) == len(environment['cpus']) == 16 and
           environment['address_space_limits'] == [16*1024**3]*2, 'native allocation differs')
    limit_proof = supervisor['address_limit_installation']
    insist(limit_proof['after'] == [16*1024**3]*2 and
           supervisor['address_limit_before_installation'] == limit_proof['before'] and
           len(limit_proof['before']) == 2 and all(type(x) is int for x in limit_proof['before']) and
           (limit_proof['before'][1] == -1 or limit_proof['before'][1] >= 16*1024**3) and
           limit_proof['installation'] == 'inside Python verification step before allocation/capture checks',
           'supervisor did not prove exact in-step address-limit installation')
    worker_environment = get('worker-environment.json')
    insist(worker_environment['environment'] == environment and
           worker_environment['address_limit']['before'] == [16*1024**3]*2 and
           worker_environment['address_limit']['after'] == [16*1024**3]*2 and
           worker_environment['address_limit']['installation'] == limit_proof['installation'],
           'worker allocation/inherited hard limit differs')
    insist(attempt['status'] == 'passed' and attempt['phase'] == 'complete' and
           attempt['native_gate_passed'] is True and attempt['native_qualification'] is False
           and attempt['signal'] is None and not attempt['attempt_history'], 'workflow did not complete one fresh attempt')
    insist([x['phase'] for x in attempt['phases']] == list(PHASES) and
           all(not x['reused'] for x in attempt['phases']), 'unexpected reuse or phase schedule')
    scope_zero(supervisor);scope_zero(attempt)
    identity = get('state/identity.json'); binding = identity['binding']
    insist(identity['binding_sha256'] == digest(binding) and
           identity['root'] == str(canonical_out/'state') and
           identity['cap_bytes'] == 4*1024**3-18*1024**2 and identity['metadata_limit'] == 1024**2,
           'state binding, native root or storage cap differs')
    completion = get('worker-complete.json')
    certs, artifact_rows = {}, {}
    for phase, dependencies in PHASES.items():
        raw = data['work/out/state/phases/'+phase+'/certificate.json']; cert = decode(raw)
        insist(set(cert) == {'schema','phase','binding_sha256','dependencies','artifacts','census','notes'}
               and cert['schema'] == 1 and cert['phase'] == phase and
               cert['binding_sha256'] == digest(binding) and cert['census']['passed'] is True,
               'phase certificate structure/binding differs: '+phase)
        insist(cert['dependencies'] == {d:certs[d] for d in dependencies}, 'phase dependency differs: '+phase)
        certs[phase] = sha(raw)
        for name, artifact in cert['artifacts'].items():
            insist(set(artifact) == {'path','size','sha256'} and hash_value(artifact['sha256']), 'artifact certificate structure')
            rel = Path(artifact['path'])
            insist(not rel.is_absolute() and '..' not in rel.parts, 'artifact path escape')
            full = 'work/out/state/phases/'+phase+'/'+str(rel)
            artifact_rows[(phase,name)] = (full,artifact)
            if phase == 'cache' and name in ('left_db','right_db'):
                insist(full not in data and artifact['size'] > 0, 'native database should remain excluded')
            else:
                insist(full in data and len(data[full]) == artifact['size'] and sha(data[full]) == artifact['sha256'],
                       'phase artifact differs: '+full)
    insist(certs == completion['phase_certificates'] == attempt['certificate_sha256'] ==
           supervisor['phase_certificates'], 'final certificate maps disagree')
    def artifact(phase,name):return decode(data[artifact_rows[(phase,name)][0]])
    anchors = decode(kit['recovery_v4p1/prior-native-input-anchors.json'])
    native_binding = decode(kit['recovery_v4p1/binding325.json'])
    insist(len(anchors) == 1536 and all(hash_value(x['sha256']) for x in anchors.values()), 'historical anchor census')
    cache = artifact('cache','receipt')
    expected_capture = {f'{arm}/execution316-rank-{rank}.bin.gz' for arm in ARMS for rank in range(256)}
    reads = cache['actual_capture_reads']
    insist(set(reads) == expected_capture and cache['optimized_streams'] == 512, 'new optimized capture census')
    for name, value in reads.items():
        insist(value == dict(compressed_sha256=anchors[name]['sha256'],stored_bytes=anchors[name]['stored_bytes'],complete=True),
               'new actual capture hash differs: '+name)
    check_source(cache['source_binding']);check_reuse(cache['checked_reuse'],kit)
    insist(cache['expected_provider_keys_sha256'] == FAMILIES, 'cache provider inventory digest differs')
    for arm, db_name in zip(ARMS, ('left_db', 'right_db')):
        database_cert = artifact_rows[('cache', db_name)][1]
        census = cache['cache_census'][arm]
        insist(census['row_counts'] == ROWS and census['database_sha256'] == database_cert['sha256'] and
               census['database_bytes'] == database_cert['size'] and census['provider_families'] == 254 and
               census['inverse_contexts'] == 1024, 'cache exact row census/database differs')
        insist(cache['merge_manifests'][arm]['database_sha256'] == database_cert['sha256'] and
               cache['merge_manifests'][arm]['row_counts'] == ROWS, 'merged/cache artifact receipt differs')
        workers = cache['worker_source_bindings'][arm]
        insist(len(workers) == 8, 'worker source binding census differs')
        for worker in workers:
            insist(worker['passed'] is True and worker['sources'] == SOURCE and
                   worker['source_manifest_sha256'] == native_binding['source_manifest_sha256'] and
                   worker['native_reader']['sha256'] == READER and
                   worker['native_reader']['path'] == str(Path(submission['kit'])/'recovery_v4p1/recover325.py'),
                   'worker source/reader binding differs')
    # Python Pool.__exit__ terminates already-consumed workers on normal exit;
    # -SIGTERM is expected cleanup, while successful shard/read/source receipts
    # above (not process exits) establish actual capture processing.
    insist(len(cache['child_exit_statuses']) == 16 and
           all(x['status_observed'] and type(x['exitcode']) is int and x['exitcode'] in (0, -15) and
               x['cleanup_policy'] == 'legacy Pool context terminates workers on exit' and
               x['scope'] == 'cache shard Pool worker' for x in cache['child_exit_statuses']),
           'worker cleanup status/policy census incomplete')
    providers = [artifact(phase,'proof') for phase in ('providers-C0','providers-L256')]
    for arm, provider in zip(ARMS,providers):
        insist(provider['arm'] == arm and provider['passed'] is True and provider['rank_count'] == 256
               and provider['counts']['families'] == len(provider['proofs']) == 254, 'provider completeness differs')
        insist(all(hash_value(x) for x in provider['proofs'].values()) and
               digest(provider['proofs']) == provider['proofs_sha256'], 'provider proof digest differs')
        keys = [ast.literal_eval(key) for key in provider['proofs']]
        insist(all(type(x) is tuple and len(x) == 2 and x[0] in ('routes','adapters','static')
                   and type(x[1]) is str for x in keys), 'provider key structure differs')
        insist(digest(sorted(keys)) == provider['families_sha256'] == FAMILIES,
               'provider identities differ from authenticated optimized family inventory')
    insist(providers[0]['proofs'] == providers[1]['proofs'], 'physical provider programs differ')
    pair = artifact('pair','proof')
    insist(pair['passed'] and pair['providers']['passed'] and pair['physical']['passed'], 'physical pair did not pass')
    empty_failures(pair['providers']['failures']);empty_failures(pair['physical']['failures'])
    physical = pair['physical']
    insist(physical['arms'] == list(ARMS) and physical['expected_ranks'] == 256 and physical['steps'] == 4,
           'physical pair scope changed')
    insist(physical['features']['common'] == physical['features']['equal'] == 508840 and
           all(x['rows'] == x['unique'] == 508840 for x in physical['features']['sides']), 'feature coverage incomplete')
    bridge = physical['inverse_bridge']
    insist(bridge['expected'] == bridge['paired'] == bridge['proved'] == 1024 and
           all(x['rows'] == x['unique'] == 1024 for x in bridge['sides']), 'physical bridge coverage incomplete')
    insist(type(bridge['moved_physical_nodes']) is int and bridge['moved_physical_nodes'] >= 24,
           'known inverse storage renumbering not represented')
    work = artifact('work','report'); work_reads = artifact('work','actual_reads')
    expected_work = {f'{arm}/work325-rank-{rank}.tsv' for arm in LABELS for rank in range(256)}
    insist(work['passed'] is True and work['actual_read_hashes_matched'] == 768 and
           set(work_reads) == expected_work, 'work actual-read census incomplete')
    for name,value in work_reads.items():
        insist(value == dict(sha256=anchors[name]['sha256'],bytes=anchors[name]['stored_bytes'],complete=True),
               'work read hash differs: '+name)
    original, builds = original_science(original_root,native_binding)
    for label in LABELS:
        old = original['runs'][label]['work_ledger']; new = work['runs'][label]
        insist(new['summary'] == old['summary'], 'replayed work counters differ from original: '+label)
        insist(new['per_rank_artifact']['sha256'] == old['per_rank_artifact']['sha256'] and
               new['per_rank_artifact']['bytes'] == old['per_rank_artifact']['bytes'],
               'replayed vectors differ from exact original serializer: '+label)
        insist(new['per_rank_artifact']['retained'] is False and 'path' not in new['per_rank_artifact'],
               'work report contains dangling vector path')
    preservation = artifact('preservation','report'); raw_binding = preservation['binding']
    insist(preservation['passed'] is True and preservation['full_input_hashes_matched'] == 1536 and
           preservation['original_ledger_preserved'] is True, 'preservation incomplete')
    check_source(preservation['source_binding']);check_reuse(preservation['checked_reuse'],kit)
    metadata = preservation['metadata']
    insist(metadata['passed'] is True and metadata['source_files'] == 258 and
           metadata['original_kit_files'] == 79 and metadata['metadata_files'] == 3889,
           'full original preservation census differs')
    for binary in metadata['binaries']:
        insist(binary['sha256'] == builds['builds'][binary['kind']]['binary_sha256'], 'preserved binary hash differs')
    insist(len(metadata['binaries']) == 2 and len(metadata['checkpoints']) == 6, 'binary/checkpoint preservation incomplete')
    for checkpoint in metadata['checkpoints']:
        path = Path(checkpoint['path']);label = path.parent.name
        index = path.name.split('_')[2].split('.')[0]
        insist(checkpoint['sha256'] == original['runs'][label]['checkpoints'][index], 'preserved checkpoint differs')
    insist(raw_binding['captures'] == digest(anchors), 'raw input binding differs from1536anchors')
    source_value = dict(original_binding_sha256=sha(kit['recovery_v4p1/binding325.json']), metadata=metadata,
        source_manifest_sha256=native_binding['source_manifest_sha256'],
        inverse_sources=preservation['source_binding'], checked_reuse=preservation['checked_reuse'])
    insist(raw_binding['source'] == digest(source_value), 'raw source binding cannot be reconstructed')
    checker = {name: sha(raw) for name,raw in kit.items() if '/' not in name and name.endswith('.py') and not name.startswith('test_')}
    checker.update({name:sha(kit[name]) for name in ('map325.json','work325-schema.json')})
    checker.update({'frozen/'+name:sha(kit['recovery_v4p1/'+name]) for name in FROZEN_NAMES})
    insist(raw_binding['checker'] == digest(checker), 'executed checker binding cannot be reconstructed')
    config = dict(workers=16,expected_ranks=256,steps=4,expected_families=254,
        arm_names=list(ARMS),native_table_counts=ROWS,contract='stage325-native-resume-v1')
    insist(raw_binding['config'] == digest(config), 'native runtime config binding differs')
    limits={'providers-C0':2*1024**2,'providers-L256':2*1024**2,'pair':2*1024**2,'work':192*1024**2,'preservation':32*1024**2}
    family_keys=sorted(ast.literal_eval(x) for x in providers[0]['proofs'])
    contract=dict(ranks=256,steps=4,families=254,arms=list(ARMS),rows=ROWS,native=True,limits=limits,family_keys=family_keys)
    wrapped=dict(raw_binding, config=digest(dict(backend=raw_binding['config'],contract=contract)),
        checker=digest(dict(backend=raw_binding['checker'],files={name:sha(kit[name]) for name in
            ('resume_workflow325.py','resume_state325.py','stream_compare325.py')})))
    insist(binding == wrapped == attempt['binding'], 'complete state is not bound to reconstructed raw evidence')
    insist(post['job'] == job and post['package_manifest_sha256'] == expected_package and
           post['fresh_binding'] == raw_binding and post['metadata'] == metadata,
           'independent post-terminal source/input binding differs')
    insist(post['full_input_hashes_matched'] == 1536 and post['phase_certificate_sha256'] == certs and
           post['state_file_sha256'] == sha(data['submission.json']) == sha(data['work/submission.json']),
           'post-terminal complete input/certificate/submission identity differs')
    expected_account = [str(job), 'COMPLETED', '0:0', str(elapsed), row[3], row[4], 'bb13']
    for receipt in (post['scheduler'], post['scheduler_after']):
        insist(receipt['returncode'] == 0 and receipt['terminal'] is True and
               receipt['clean_completion'] is True and receipt['row'] == expected_account,
               'post-terminal scheduler identity differs from archived accounting')
    check_source(post['source_binding']);check_reuse(post['checked_reuse'],kit)
    for name in ('left_db','right_db'):
        full,cert=artifact_rows[('cache',name)];actual=post['databases'][name]
        insist(actual['sha256'] == cert['sha256'] and actual['bytes'] == cert['size'] and
               actual['path'] == str(canonical_out/'state/phases/cache'/cert['path']) and
               actual['canonical_path'] == actual['path'] and actual['samefile'] is True and
               actual['stat_before'] == actual['stat_after'] and
               actual['stat_before']['bytes'] == cert['size'] and
               actual['expected_sha256'] == cert['sha256'] and actual['expected_bytes'] == cert['size'],
               'independent retained database hash differs: '+name)
    return dict(passed=True,job=job,review_sha256=sha(Path(archive).read_bytes()),
        review_files=len(data),review_logical_bytes=size,package_manifest_sha256=expected_package,
        approved_source_archive_sha256=APPROVED_SOURCE_ARCHIVE,package_files=len(kit),
        native_preflight_tests=expected_tests,native_preflight_skipped=0,in_step_resource_installation_proved=True,
        checked_reuse_job=55573,optimized_actual_reads=512,provider_families_each=254,
        features_each=508840,physical_bridges=1024,moved_physical_node_instances=bridge['moved_physical_nodes'],
        work_actual_reads=768,preserved_input_hashes=1536,preserved_source_files=258,
        source_job=55569,scheduler_elapsed_seconds=elapsed,cumulative_target_seconds=2656+elapsed,
        post_terminal_report_sha256=sha(Path(post_path).read_bytes()),
        native_qualification_supported=True,performance_qualified=False,
        limitations=['No new solver or performance run; original science gates are reused through exact source/binary/checkpoint/input hashes.',
          'Compact evidence does not permit an independent second execution of all semantic comparisons. Source-reviewed checker completion and exact native cache hashes are explicit premises.',
          'Communication/work changes are counters, not elapsed-time benefit or a cause of negative scaling.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path);parser.add_argument('--job',type=int,required=True)
    parser.add_argument('--post-terminal',type=Path,required=True)
    parser.add_argument('--original-review',type=Path,default=Path(__file__).parent/'native-review-55569/review')
    parser.add_argument('--manifest-sha256',default=PACKAGE,
                        help='Trusted separately reviewed sealed package manifest; never inferred from the archive')
    parser.add_argument('--expected-tests',type=int,default=EXPECTED_TESTS)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:result=review(args.archive,args.job,args.post_terminal,args.original_review,
                      args.manifest_sha256,args.expected_tests)
    except (ReviewError,KeyError,ValueError,TypeError,OSError) as error:
        result=dict(passed=False,job=args.job,review_error_type=type(error).__name__,review_error=str(error),
                    native_qualification_supported=False,performance_qualified=False)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));return 0 if result['passed'] else 1


if __name__=='__main__':
    sys.exit(main())
