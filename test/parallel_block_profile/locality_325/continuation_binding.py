"""Reuse the completed job55567 gate without reading any execution captures.

This module only reads existing files. Its authority is the pinned, completed
review archive, not a mutable passed flag in the original failed campaign.
Current compiler/library/placement checks remain the continuation runner's duty.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import re
import tarfile

REVIEW_NAME = 'rk322-semantic-review-55567.tar.gz'
REVIEW_SHA256 = '44fbd01b84dee5afd4cdcc26e9ea764050b41c6d04a5baca769021ce698ba874'
BINDING_SHA256 = 'a384b5b664ac29a54e9a2fd606a95eeac115572c773f291d39dfb1d84dfefbfc'
KIT_MANIFEST_SHA256 = 'caf586098c44cb2b9fdd339d9944c63715b82678ee950f76be3706a655abe428'
SOURCE_MANIFEST_SHA256 = '8d5882a1bb24e2c74bc5ca5ca4f04d387ad3636a49c4208e6f3160b63404fbc5'
REVIEW_ROOT = 'rk322-semantic-replay-55567/'
ARMS = {'checked-B': 256, 'A0': 160, 'B': 256, 'C0': 256, 'Cs': 256}
HEX = re.compile(r'[0-9a-f]{64}\Z')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def decode(data):
    return json.loads(data, object_pairs_hook=_unique_object)


def sha(path):
    """Match the original source/metadata manifest's symlink representation."""
    path = Path(path)
    if path.is_symlink():
        return hashlib.sha256(str(path.readlink()).encode()).hexdigest()
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def under(root, relative):
    """Manifest names must stay under their declared root, including parents."""
    root = Path(root).resolve(strict=True)
    rel = PurePosixPath(relative)
    require(isinstance(relative, str) and relative and not rel.is_absolute()
            and '..' not in rel.parts and str(rel) == relative,
            'Unsafe relative evidence path: ' + str(relative))
    path = root.joinpath(*rel.parts)
    require(path.parent.resolve(strict=True).is_relative_to(root),
            'Evidence parent escapes root: ' + relative)
    return path


def match(path, expected, *, regular=False):
    path = Path(path)
    require(isinstance(expected, str) and HEX.fullmatch(expected),
            'Invalid expected SHA256: ' + str(path))
    if regular:
        require(path.is_file() and not path.is_symlink(),
                'Expected regular evidence file: ' + str(path))
    require(sha(path) == expected, 'Evidence hash mismatch: ' + str(path))


def _read_archive(path):
    """No extraction, links, ambiguous members, or unbounded decompression."""
    payload = {}
    with tarfile.open(path, 'r:gz') as archive:
        total = 0
        for entry in archive:
            name = entry.name
            rel = PurePosixPath(name)
            require(entry.isfile() and not rel.is_absolute()
                    and '..' not in rel.parts and str(rel) == name,
                    'Unsafe review member: ' + name)
            require(name not in payload, 'Duplicate review member: ' + name)
            require(0 <= entry.size <= 16 * 1024**2, 'Review member too large')
            total += entry.size
            require(total <= 64 * 1024**2, 'Review archive too large')
            payload[name] = archive.extractfile(entry).read()
    manifest_name = REVIEW_ROOT + 'review-manifest.json'
    require(manifest_name in payload, 'Missing review manifest')
    manifest = decode(payload[manifest_name])
    require(len(manifest) == 9 and set(manifest) == set(payload) - {manifest_name},
            'Incomplete or extra review members')
    for name, digest in manifest.items():
        require(hashlib.sha256(payload[name]).hexdigest() == digest,
                'Review member hash mismatch: ' + name)
    return payload


def _validate_review(payload):
    def read(name):
        return decode(payload[REVIEW_ROOT + name])
    report = read('report.json')
    require(report['schema'] == 'stage322-corrected-replay-v1'
            and report['source_job'] == 55561 and report['verification_job'] == '55567'
            and report['passed'] is True and report['performance_qualified'] is False
            and report['all_streams_verified'] == 1184
            and report['metadata_files_matched'] == 6793
            and all(report[k] == 0 for k in ('solver_runs', 'builds', 'source_writes')),
            'Completed semantic review gate differs')
    rows = payload[REVIEW_ROOT + 'accounting.txt'].decode().splitlines()
    require(len(rows) == 2, 'Unexpected verification accounting rows')
    account = dict(zip(rows[0].split('|'), rows[1].split('|')))
    require(account['JobID'] == '55567' and account['State'] == 'COMPLETED'
            and account['ExitCode'] == '0:0', 'Verification accounting did not pass')
    qualification = read('qualification.json')
    require(qualification['passed'] is True
            and qualification['performance_qualified'] is False
            and qualification['source_binding'] == report['schema_binding'],
            'Qualification status/source binding differs')
    checks = qualification['checks']
    require(len(checks) == 3 and all(c['passed'] is True for c in checks),
            'Incomplete qualification checks')
    checked, pair_ab, pair_cs = checks
    require(checked['arm'] == 'checked-B' and checked['ranks'] == 256,
            'Checked arm differs')
    require((pair_ab['left'], pair_ab['right'], pair_ab['compared_ranks']) == ('A0', 'B', 160)
            and (pair_cs['left'], pair_cs['right'], pair_cs['compared_ranks']) == ('C0', 'Cs', 256),
            'Exact rank-pair qualification differs')
    streams = {'checked-B': checked['streams'], 'A0': pair_ab['left_streams'],
               'B': pair_ab['right_streams'], 'C0': pair_cs['left_streams'],
               'Cs': pair_cs['right_streams']}
    for arm, n in ARMS.items():
        require([s['rank'] for s in streams[arm]] == list(range(n)),
                'Rank-contiguous stream inventory differs: ' + arm)
        require(all(HEX.fullmatch(s['stream_sha256']) and s['records'] > 0
                    and s['uncompressed_bytes'] > 0 for s in streams[arm]),
                'Incomplete stream summary: ' + arm)
    families = qualification['crossrank_reconstructed_families']
    require(set(families) == set(ARMS)
            and all(v == {'adapters': 8, 'routes': 140, 'static': 106}
                    for v in families.values()), 'Incomplete provider reconstruction')
    digests = qualification['route_identity_digests']
    require(set(digests) == set(ARMS) and all(HEX.fullmatch(v) for v in digests.values())
            and digests['A0'] == digests['B'] and digests['C0'] == digests['Cs'],
            'Paired provider identity differs')
    checkpoints = read('checkpoint-checks.json')
    expected_checks = {'checked-B', 'qual-A0', 'qual-B', 'qual-C0', 'qual-Cs',
                       'qual-A0/qual-B', 'qual-C0/qual-Cs'}
    require(set(checkpoints) == expected_checks
            and all(c['passed'] is True for c in checkpoints.values()),
            'Checkpoint comparisons incomplete')
    for key in ('qual-A0/qual-B', 'qual-C0/qual-Cs'):
        require(checkpoints[key]['sha256_exact'] is True
                and checkpoints[key]['exact_fields_required'] is True,
                'Exact paired science endpoint differs')
    native = read('native-gates.json')
    require(native['passed'] is True and native['callback_cases'] == 17
            and native['serial_full_builds'] == 3,
            'Recorded exact-native build gates incomplete')
    return report, qualification, native


def load_qualification(continuation_kit):
    """Read and validate only the small sealed completed-review artifact."""
    kit = Path(continuation_kit).resolve(strict=True)
    archive = kit / REVIEW_NAME
    match(archive, REVIEW_SHA256, regular=True)
    payload = _read_archive(archive)
    report, qualification, native = _validate_review(payload)
    binding_path = kit / 'replay-binding.json'
    match(binding_path, BINDING_SHA256, regular=True)
    binding = decode(binding_path.read_bytes())
    public_binding = {k: v for k, v in binding.items()
                      if k not in ('metadata_hashes', 'raw_inventory')}
    require(public_binding == report['binding'], 'Reviewed replay binding differs')
    return report, qualification, native, binding


def metadata_entries(binding):
    """Select the same three evidence roots as the completed V4 replay.

    Collection-level files such as ``slurm.log`` have no family separator and
    are deliberately outside this native-input inventory.
    """
    for name, digest in binding['metadata_hashes'].items():
        if not name.startswith(('kit/', 'results/', 'preflight/')):
            continue
        family, relative = name.split('/', 1)
        yield family, relative, digest


def verify_binding(continuation_kit, state_path=None):
    """Fail closed unless the completed qualification still binds every input.

    Read small metadata, sealed sources, executables and compressed checkpoints.
    Never open raw execution captures or decompress their contents.
    """
    report, qualification, native, binding = load_qualification(continuation_kit)
    state_path = Path(state_path) if state_path else Path.home() / 'rk322-v1-job-55561.json'
    state = decode(state_path.read_bytes())
    require(int(state['job']) == 55561, 'Wrong source campaign job')
    kit = Path(state['kit']).resolve(strict=True)
    results = Path(state['out']).resolve(strict=True)
    require(kit == Path(report['original_kit']).resolve(strict=True)
            and results == Path(report['original_results']).resolve(strict=True),
            'Source campaign paths differ from completed qualification')
    match(kit / 'files-sha256.json', KIT_MANIFEST_SHA256, regular=True)
    original = decode((kit / 'files-sha256.json').read_bytes())
    for name, digest in original.items():
        match(under(kit, name), digest)
    roots = {'kit': kit, 'results': results, 'preflight': results.parent / 'preflight'}
    count = 0
    for family, rel, digest in metadata_entries(binding):
        require('execution316-rank-' not in rel, 'Raw capture present in metadata binding')
        match(under(roots[family], rel), digest)
        count += 1
    require(count == 6793, 'Reviewed metadata inventory differs')
    match(kit / 'source-files.json', SOURCE_MANIFEST_SHA256, regular=True)
    sources = decode((kit / 'source-files.json').read_bytes())
    require(set(sources) == {'candidate', 'main'}, 'Unexpected source families')
    source_count = 0
    for family, files in sources.items():
        root = results / 'sources' / family
        for name, digest in files.items():
            match(under(root, name), digest)
            source_count += 1
    ledger = decode((results / 'results.json').read_bytes())
    require(ledger['status'] == 'failed' and ledger['qualification'] is None,
            'Original failed ledger must remain unchanged')
    expected_runs = {'bootstrap-fresh', 'legacy-four', 'checked-B',
                     'qual-A0', 'qual-B', 'qual-C0', 'qual-Cs'}
    require(set(ledger['runs']) == expected_runs
            and all(not r['measurement'] and r['status'] == 'passed'
                    for r in ledger['runs'].values()),
            'Original seven-launch qualification ledger differs')
    builds = decode((results / 'builds.json').read_bytes())
    toolchain = decode((results / 'toolchain313.json').read_bytes())
    require(builds == native['builds'] == ledger['builds'], 'Recorded native builds differ')
    require(toolchain == builds['toolchain'] and toolchain['status'] == 'passed'
            and toolchain['expected'] == decode((kit / 'expected-toolchain.json').read_bytes()),
            'Recorded compiler/runtime gate differs')
    bound = []
    def artifact(path, digest):
        match(path, digest, regular=True)
        row = {'path': str(path), 'sha256': digest}
        if row not in bound:
            bound.append(row)
        return str(path)
    binaries = {}
    for kind, row in builds['builds'].items():
        binaries[kind] = artifact(under(results, row['binary']), row['binary_sha256'])
    artifact(results / 'callback-gate/callback-check', builds['callback_gate']['binary_sha256'])
    probe_hash = binding['metadata_hashes']['results/mpi-identity']
    probe = artifact(results / 'mpi-identity', probe_hash)
    seed = ledger['fresh_fixture']
    seed_path = artifact(under(results, seed['path']), seed['sha256'])
    references = {}
    for label, row in ledger['runs'].items():
        require(row['steps'] == (0 if label == 'bootstrap-fresh' else 4),
                'Qualified step count differs: ' + label)
        for index, digest in row['checkpoints'].items():
            require(index in ('9900', '9901'), 'Unexpected qualified checkpoint index')
            artifact(results / label / ('test_checkpoint_' + index + '.bin.zst'), digest)
        if row['steps']:
            require(row['seed_sha256'] == seed['sha256'], 'Per-run seed differs')
            artifact(results / label / 'test_checkpoint_0003.bin.zst', row['seed_sha256'])
            require(set(row['checkpoints']) == {'9900', '9901'}, 'Missing reference endpoint')
            references[label] = {
                'initial': str(results / label / 'test_checkpoint_9900.bin.zst'),
                'initial_sha256': row['checkpoints']['9900'],
                'final': str(results / label / 'test_checkpoint_9901.bin.zst'),
                'final_sha256': row['checkpoints']['9901']}
    for left, right in (('qual-A0', 'qual-B'), ('qual-C0', 'qual-Cs')):
        require(references[left]['final_sha256'] == references[right]['final_sha256'],
                'Qualified paired endpoint hash differs')
    return dict(passed=True, source_job=55561, qualification_job=55567,
                performance_qualified=False, original_kit=str(kit), original_results=str(results),
                qualification_archive_sha256=REVIEW_SHA256, replay_binding_sha256=BINDING_SHA256,
                original_kit_manifest_sha256=KIT_MANIFEST_SHA256,
                source_manifest_sha256=SOURCE_MANIFEST_SHA256,
                reviewed_metadata_files=count, source_files_checked=source_count,
                candidate_binary=binaries['candidate'],
                candidate_sha256=builds['builds']['candidate']['binary_sha256'],
                binaries=binaries, private_seed=seed_path, private_seed_sha256=seed['sha256'],
                mpi_identity_probe=probe, mpi_identity_probe_sha256=probe_hash,
                reference_checkpoints=references, original_builds=builds,
                original_toolchain=toolchain, bound_artifacts=bound,
                qualification={'passed': True, 'all_streams_verified': 1184,
                               'crossrank_reconstructed_families': qualification['crossrank_reconstructed_families'],
                               'route_identity_digests': qualification['route_identity_digests']},
                raw_capture_reads=0,
                limitation='Reuses completed captured-state/science qualification; current runtime and placement must be checked separately. No timing result is inferred.')


def recheck_artifacts(binding_result):
    """Recheck the native executables and checkpoints after continuation."""
    require(binding_result['passed'] is True, 'Binding was not successful')
    for row in binding_result['bound_artifacts']:
        match(row['path'], row['sha256'], regular=True)
    return {'passed': True, 'files': len(binding_result['bound_artifacts'])}
