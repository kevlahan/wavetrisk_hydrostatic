"""Read-only native adapter for resumable verification of job55569 captures.

This adapter has no scheduler, compiler or solver entry point. Its constructor
does not run a replay. ``binding`` freshly reads and hashes protected inputs;
only explicit workflow phase calls construct derived state. The old recovery
helpers remain byte-for-byte frozen. Their execution/cleanup entry points are
never called here.
"""
import ast
from contextlib import contextmanager
import functools
import hashlib
import importlib
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import sys
import types

from resume_state325 import canonical325, sha256325


FROZEN_FILES325 = {
    'recover325.py': '5feaf0dca377947d78e7840dfb3b167eb95e03234ba59cb9175ba0a503179ae1',
    'recover_parallel325.py': '85133b6ae716598fc2e879346e8224ce602ebbb3485c2e71512cac5bfa9f6ea0',
    'binding325.json': 'f76036ce9035f64c2c60f221fc0e9049a675f5a12b4744dd55dab937a09d98b5',
    'parallel-history325.json': '47ea2f4df9d39d3cbcada411bbb07614f6ba897d351e519b8940071cfc4dee43',
    'recovery-history325.json': 'e7d2d5671dd361cbc71c7520e73e4570c40b495ed2ef3c3e3b6819d8ff6274d3',
    'prior-native-input-anchors.json': 'f1a3b67af02f8a70c25e420027b9c3b3a4c2527278940cfc7301ff655593deb7',
}
NATIVE_TABLE_COUNTS325 = {'route': 65024, 'domain': 256, 'context': 66048,
                         'feature': 508840}
LABELS325 = ('checked-L256', 'qual-C0', 'qual-L256')
ARMS325 = ('qual-C0', 'qual-L256')
COMPACT_LIMIT325 = 1024**2
CACHE_RECEIPT_RESERVE325 = 2 * COMPACT_LIMIT325
EXECUTED_MODULES325 = ('verify325', 'verify322', 'coverage322', 'inverse_identity325',
                       'parallel_replay325', 'account325', 'resume_state325',
                       'resume_native325', 'stream_compare325', 'reconcile_boundary325')


def require325(ok, reason):
    if not ok:
        if any(word in reason for word in ('budget', 'preallocation', 'receipt_reserve')):
            raise RuntimeError('resume.native325.resource.' + reason)
        raise ValueError('resume.native325.' + reason)


def semantic_phase325(method):
    @functools.wraps(method)
    def wrapped(*args, **kwargs):
        try:
            return method(*args, **kwargs)
        except ValueError as error:
            # Account/frozen helper guards use plain ValueError. Preserve the
            # workflow's failure latch instead of treating a failed invariant
            # as an interrupt that permits an unreviewed restart.
            from resume_workflow325 import SemanticFailure325
            raise SemanticFailure325(str(error)) from error
    return wrapped


def digest325(value):
    return hashlib.sha256(canonical325(value)).hexdigest()


def actual_checker_paths325():
    root = Path(__file__).resolve().parent
    paths = {}
    for name in EXECUTED_MODULES325:
        module = importlib.import_module(name)
        actual = Path(module.__file__).resolve(strict=True)
        require325(actual == root / (name + '.py'), 'executed_checker_path:' + name)
        paths[name + '.py'] = actual
    return paths


def frozen_reader_path325(name, expected_package=None):
    """Resolve a namespace reader without importing any package initializer."""
    require325(type(name) is str and len(name.split('.')) == 2 and
               name.endswith('.recover325'), 'worker_reader_name')
    package_name, leaf = name.split('.')
    specification = importlib.machinery.PathFinder.find_spec(package_name, sys.path)
    require325(specification is not None and specification.origin is None and
               specification.submodule_search_locations is not None, 'reader_namespace_required')
    roots = {Path(path).resolve(strict=True) for path in specification.submodule_search_locations}
    require325(len(roots) == 1, 'reader_namespace_ambiguous')
    package = next(iter(roots))
    require325(expected_package is None or package == Path(expected_package).resolve(strict=True),
               'conflicting_frozen_reader')
    loaded = sys.modules.get(package_name)
    if loaded is not None:
        require325(getattr(loaded, '__file__', None) is None and
                   {Path(path).resolve(strict=True) for path in loaded.__path__} == roots,
                   'loaded_reader_namespace_changed')
    child = importlib.machinery.PathFinder.find_spec(leaf, [str(package)])
    require325(child is not None and child.origin is not None and
               Path(child.origin).resolve(strict=True) == package / 'recover325.py',
               'worker_reader_missing')
    path = Path(child.origin)
    require325(sha256325(_regular325(path)) == FROZEN_FILES325['recover325.py'] and
               sha256325(_regular325(package / 'binding325.json')) == FROZEN_FILES325['binding325.json'],
               'worker_reader_or_binding_bytes')
    return path


def native_cache_worker325(job):
    """Check fresh-spawn imports and reader bytes before executing that reader."""
    actual_checker_paths325()
    reader_name = job[-1]
    path = frozen_reader_path325(reader_name)
    import parallel_replay325 as parallel
    result = parallel._worker325(job)
    result['source_binding']['native_reader'] = dict(path=str(path.resolve()), sha256=sha256325(path))
    return result


def _regular325(path):
    path = Path(path)
    require325(path.is_file() and not path.is_symlink(), 'regular_file_required')
    return path


def tree_bytes325(root):
    root = Path(root)
    total = 0
    for path in root.rglob('*'):
        require325(not path.is_symlink(), 'derived_symlink')
        require325(path.is_dir() or path.is_file(), 'derived_file_type')
        if path.is_file():
            require325(path.stat().st_nlink == 1, 'derived_hardlink')
            total += path.stat().st_size
    return total


def bounded_json325(root, path, value, allowance, *, compact_limit=COMPACT_LIMIT325):
    """Budget actual encoded bytes *before* the only temporary-file write."""
    root, path = Path(root), Path(path)
    require325(path.parent.resolve() == root.resolve() and not path.exists(), 'fresh_compact_path')
    encoded = canonical325(value)
    require325(len(encoded) <= compact_limit, 'compact_receipt_bound')
    require325(tree_bytes325(root) + len(encoded) <= allowance, 'phase_write_budget')
    temp = path.with_suffix(path.suffix + '.tmp')
    require325(not temp.exists() and not temp.is_symlink(), 'fresh_compact_temp')
    descriptor = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'wb') as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


def remove_owned_tree325(root, tree):
    """Only a checked, non-symlink private phase subtree may be discarded."""
    root, tree = Path(root).resolve(strict=True), Path(tree)
    require325(tree.parent.resolve(strict=True) == root and tree.name == 'work-input',
               'owned_work_tree_required')
    tree_bytes325(tree)
    for base, directories, files in os.walk(tree, topdown=False, followlinks=False):
        for name in files:
            (Path(base) / name).unlink()
        for name in directories:
            (Path(base) / name).rmdir()
    tree.rmdir()


@contextmanager
def tracked_pool_context325(parallel, statuses):
    """Observe real exit codes without changing the old pool/merge algorithm."""
    original = parallel.multiprocessing

    class Context:
        def __init__(self, context):
            self.context = context

        def Pool(self, *args, **kwargs):
            pool = self.context.Pool(*args, **kwargs)

            class TrackedPool:
                def __enter__(self):
                    self.workers = tuple(pool._pool)
                    return pool.__enter__()

                def __exit__(self, kind, value, traceback):
                    try:
                        return pool.__exit__(kind, value, traceback)
                    finally:
                        statuses.extend(dict(pid=worker.pid, exitcode=worker.exitcode,
                                             status_observed=worker.exitcode is not None,
                                             scope='cache shard Pool worker',
                                             cleanup_policy='legacy Pool context terminates workers on exit')
                                        for worker in self.workers)

            return TrackedPool()

    parallel.multiprocessing = types.SimpleNamespace(
        get_context=lambda method: Context(original.get_context(method)))
    try:
        yield
    finally:
        parallel.multiprocessing = original


@contextmanager
def bounded_account_writes325(module, root, allowance):
    """Intercept the frozen account module's seven write_text destinations.

    Reads can address the protected schema. Every write must be below this one
    owned tree and fit the already-reserved allowance before opening a writer.
    This changes no parse or accounting operation.
    """
    root = Path(root).resolve(strict=True)
    original_path = module.Path
    concrete = type(Path())

    class AccountPath325(concrete):
        def write_text(self, data, encoding=None, errors=None, newline=None):
            path = Path(self)
            require325(path.parent.resolve(strict=True).is_relative_to(root) and
                       not path.is_symlink(), 'account_write_outside_owned_stage')
            for component in (path.parent, *path.parent.parents):
                if component == root:
                    break
                require325(not component.is_symlink(), 'account_symlink_parent')
            require325(not path.exists(), 'account_overwrite')
            encoded = data.encode(encoding or 'utf-8', errors or 'strict')
            require325(tree_bytes325(root) + len(encoded) <= allowance, 'account_write_budget')
            descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
            with os.fdopen(descriptor, 'wb') as stream:
                stream.write(encoded)
                stream.flush()
                os.fsync(stream.fileno())
            return len(data)

    module.Path = AccountPath325
    try:
        yield
    finally:
        module.Path = original_path


class NativeBackend325:
    native = True
    expected_ranks = 256
    steps = 4
    expected_families = 254
    arm_names = ARMS325
    expected_table_rows = dict(NATIVE_TABLE_COUNTS325)

    def __init__(self, package, *, binding_path=None, runtime_config=None, checker_files=None):
        require325(checker_files is None, 'native_checker_override_forbidden')
        self.package = Path(package).resolve(strict=True)
        require325('.' not in self.package.name, 'top_level_frozen_namespace_required')
        self.binding_path = self.package / 'binding325.json' if binding_path is None else Path(binding_path)
        self._check_frozen()
        require325(sha256325(_regular325(self.binding_path)) == FROZEN_FILES325['binding325.json'],
                   'original_binding_changed')
        self.original_binding = json.loads(self.binding_path.read_text())
        self.original = Path(self.original_binding['original_results']).resolve()
        self.original_kit = Path(self.original_binding['original_kit']).resolve()
        config = {} if runtime_config is None else dict(runtime_config)
        allowed = {'workers', 'expected_ranks', 'steps', 'expected_families', 'arm_names'}
        require325(set(config) <= allowed, 'unknown_runtime_configuration')
        self.workers = config.get('workers', 16)
        require325(type(self.workers) is int and 2 <= self.workers <= 16 and self.workers % 2 == 0,
                   'worker_count')
        for key, expected in (('expected_ranks', 256), ('steps', 4), ('expected_families', 254)):
            require325(type(config.get(key, expected)) is int and config.get(key, expected) == expected,
                       'native_scope_' + key)
        require325(tuple(config.get('arm_names', ARMS325)) == ARMS325, 'native_arms')
        self.runtime_config = dict(workers=self.workers, expected_ranks=256, steps=4,
                                   expected_families=254, arm_names=list(ARMS325),
                                   native_table_counts=NATIVE_TABLE_COUNTS325,
                                   contract='stage325-native-resume-v1')
        self.checker_files = checker_files
        self.helper = self._load_helper()
        h = json.loads((self.package / 'parallel-history325.json').read_text())
        self.protected_roots = (self.original, self.original_kit, self.package,
                                self.binding_path.resolve(strict=True),
                                self.original.parent / 'preflight',
                                *tuple(Path(x).parent.resolve() for x in h['remote_predecessor_files']))
        self.expected_family_keys = None
        self.anchors = None
        self._last_binding = None
        self._bound_metadata = None
        self._checked_reuse = None
        self.child_exit_statuses = []

    def _check_frozen(self):
        require325(not (self.package / '__init__.py').exists(), 'frozen_namespace_required')
        for name, expected in FROZEN_FILES325.items():
            require325(sha256325(_regular325(self.package / name)) == expected,
                       'frozen_helper_changed:' + name)

    def _load_helper(self):
        # Compile only this pinned helper module, without creating __pycache__
        # anywhere inside an old, protected package.
        name = '_stage325_frozen_native_reader'
        module = types.ModuleType(name)
        module.__file__ = str(self.package / 'recover325.py')
        exec(compile(Path(module.__file__).read_bytes(), module.__file__, 'exec'), module.__dict__)
        return module

    def _history(self):
        # Reuse the exact previously reviewed checked-closure guard while never
        # importing its allocation, process-launch or cleanup entry points.
        path = self.package / 'recover_parallel325.py'
        parsed = ast.parse(path.read_bytes(), filename=str(path))
        definitions = [node for node in parsed.body if isinstance(node, ast.FunctionDef)
                       and node.name == 'history_check']
        require325(len(definitions) == 1, 'frozen_history_function')
        scope = dict(Path=Path, json=json, LABELS=self.helper.LABELS,
                     **{name: getattr(self.helper, name)
                        for name in ('require', 'load_history', 'sha', 'safe')})
        exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), 'exec'), scope)
        history, anchors, semantic = scope['history_check'](self.package, self.original_binding)
        checked = semantic['runs']['checked-L256']
        require325(checked['status'] == 'complete' and checked['families'] == 254 and
                   [x['rank'] for x in checked['streams']] == list(range(256)),
                   'checked_closure_incomplete')
        return anchors, dict(passed=True, source_job=55573, streams=256, families=254,
                             semantic_sha256=history['certificate_sha256'][
                                 'predecessor-v3/semantic-replay.json'],
                             actual_read_hashes=256,
                             scope='Earlier checked closure only; optimized proof recomputed.')

    def _fresh_inputs(self, historical):
        anchors = {}
        for inventory, stem in (('capture_inventory', 'execution316-rank-*.bin.gz'),
                                ('work_inventory', 'work325-rank-*.tsv')):
            expected = {row['relative'] for row in self.original_binding[inventory]}
            actual = {str(path.relative_to(self.original)) for arm in LABELS325
                      for path in (self.original / arm).glob(stem)}
            require325(expected == actual and len(expected) == 768, 'native_input_inventory')
            for row in self.original_binding[inventory]:
                path = self.helper.safe(self.original, row['relative'])
                sha = self.helper.hash_bound_input(path, row, historical)
                anchors[row['relative']] = dict(row, sha256=sha, family=inventory)
        require325(len(anchors) == 1536, 'native_anchor_census')
        return anchors

    def bind_sources(self):
        import inverse_identity325 as inv
        path = self.original_kit / 'source-files.json'
        require325(sha256325(_regular325(path)) == self.original_binding['source_manifest_sha256'],
                   'original_source_manifest_changed')
        proof = inv.bind_inverse_bridge_sources325(json.loads(path.read_text()),
                                                  source_root=self.original / 'sources' / 'candidate')
        require325(proof.get('passed') is True, 'inverse_sources_unbound')
        return proof

    def _family_keys(self, anchors):
        import verify325 as q
        observed = {}
        original_reader = q.v.open_stream
        q.v.open_stream = lambda p: self.helper.verified_stream(p, self.original, anchors, observed)
        expected = None
        try:
            # The checked and optimized saved schedules have different
            # transform epochs, so their family keys cannot be interchanged.
            # Bind both optimized inventories to their actual saved bytes;
            # preserve every key verbatim rather than normalizing epochs.
            for arm in ARMS325:
                path = self.original / arm / 'execution316-rank-0.bin.gz'
                obj = q.v.parse_stream(path, 0, 256, arm, 'candidate')
                keys = sorted((kind, key) for kind in ('routes', 'adapters', 'static') for key in obj[kind])
                require325(len(keys) == len(set(keys)) == 254, 'expected_provider_key_census')
                for kind, key in keys:
                    require325(kind in ('routes', 'adapters', 'static') and type(key) is str and
                               0 < len(key) <= 4096, 'expected_provider_key_schema')
                    # Route keys append '/route/N' to a captured tuple path.
                    # Adapter/static keys themselves are Python tuple literals.
                    if kind != 'routes':
                        require325(type(ast.literal_eval(key)) is tuple, 'expected_provider_key_schema')
                if expected is None:
                    expected = keys
                else:
                    require325(keys == expected, 'optimized_expected_provider_key_difference')
        finally:
            q.v.open_stream = original_reader
        require325(set(observed) == {arm + '/execution316-rank-0.bin.gz' for arm in ARMS325},
                   'expected_provider_actual_read_census')
        return expected

    def _worker_bindings(self, report, source_proof):
        bindings = report['source_bindings']
        require325(type(bindings) is dict and set(bindings) == set(ARMS325),
                   'worker_binding_arms')
        for arm in ARMS325:
            rows = bindings[arm]
            require325(type(rows) is list and len(rows) == self.workers // 2,
                       'worker_binding_census')
            for row in rows:
                require325(row.get('passed') is True and row.get('sources') == source_proof['sources']
                           and row.get('source_manifest_sha256') == self.original_binding['source_manifest_sha256'],
                           'worker_source_binding')
                reader = row.get('native_reader', {})
                require325(reader.get('sha256') == FROZEN_FILES325['recover325.py'] and
                           Path(reader.get('path', '')).resolve() ==
                           (self.package / 'recover325.py').resolve(), 'worker_reader_attestation')
        return bindings

    def _checker_hashes(self):
        root = Path(__file__).resolve().parent
        if self.checker_files is None:
            paths = {path.name: path for path in root.glob('*.py')
                     if not path.name.startswith('test_')}
            for name in ('map325.json', 'work325-schema.json'):
                paths[name] = root / name
        else:
            paths = {str(name): Path(path) for name, path in self.checker_files.items()}
        required = ('verify325.py', 'verify322.py', 'coverage322.py', 'inverse_identity325.py',
                    'parallel_replay325.py', 'account325.py', 'resume_state325.py',
                    'resume_native325.py', 'stream_compare325.py', 'reconcile_boundary325.py',
                    'map325.json', 'work325-schema.json')
        require325(set(required) <= set(paths), 'checker_binding_incomplete')
        executed = actual_checker_paths325()
        require325(all(Path(paths[name]).resolve(strict=True) == path
                       for name, path in executed.items()), 'checker_binding_not_executed_bytes')
        hashes = {name: sha256325(_regular325(path)) for name, path in sorted(paths.items())}
        for name, expected in FROZEN_FILES325.items():
            hashes['frozen/' + name] = sha256325(self.package / name)
        # Work-output bounds use the unchanged deterministic serializer.
        original_manifest = json.loads((self.original_kit / 'files-sha256.json').read_text())
        require325(hashes['account325.py'] == original_manifest['account325.py'],
                   'account_serializer_changed')
        return hashes

    def binding(self):
        self._check_frozen()
        metadata = self.helper.metadata_check(self.original_binding)
        historical, checked = self._history()
        anchors = self._fresh_inputs(historical)
        sources = self.bind_sources()
        keys = self._family_keys(anchors)
        checker = self._checker_hashes()
        source = dict(original_binding_sha256=sha256325(self.binding_path), metadata=metadata,
                      source_manifest_sha256=self.original_binding['source_manifest_sha256'],
                      inverse_sources=sources, checked_reuse=checked)
        value = {'source': digest325(source), 'checker': digest325(checker),
                 'config': digest325(self.runtime_config), 'captures': digest325(anchors)}
        self._last_binding = value
        self.anchors = anchors
        self.expected_family_keys = keys
        self._bound_metadata = metadata
        self._checked_reuse = checked
        return dict(value)

    def _ready(self):
        require325(self._last_binding is not None and self.anchors is not None and
                   self.expected_family_keys is not None, 'fresh_binding_required')

    def _stage(self, stage, allowance):
        self._ready()
        stage = Path(stage)
        require325(stage.is_dir() and not stage.is_symlink() and not any(stage.iterdir()),
                   'fresh_owned_stage')
        require325(type(allowance) is int and allowance > 0, 'phase_allowance')
        actual = stage.resolve(strict=True)
        for root in self.protected_roots:
            require325(not actual.is_relative_to(root) and not root.is_relative_to(actual),
                       'stage_overlaps_protected_input')
        return stage

    def _activate_reader(self):
        # Spawned shard workers import the pinned historical reader by its old
        # name, whose __file__ also supplies the exact original binding JSON.
        # -B propagates into spawn and keeps all protected package paths read-only.
        sys.dont_write_bytecode = True
        self._check_frozen()
        actual_checker_paths325()
        # A qualified namespace import prevents a frozen package from shadowing
        # current verify/parallel/compare modules in a fresh spawned interpreter.
        if str(self.package.parent) not in sys.path:
            sys.path.append(str(self.package.parent))
        name = self.package.name + '.recover325'
        frozen_reader_path325(name, self.package)
        module = importlib.import_module(name)
        require325(Path(module.__file__).resolve() == (self.package / 'recover325.py').resolve() and
                   sha256325(module.__file__) == FROZEN_FILES325['recover325.py'], 'worker_reader_bytes')
        return name

    def validate_cache(self, database, arm):
        import stream_compare325 as compare
        with compare.readonly_database325(database) as (connection, identity):
            tables = tuple(NATIVE_TABLE_COUNTS325)
            compare._schema325(connection, tables)
            compare._arm325(connection, tables, arm)
            require325(connection.execute('pragma integrity_check').fetchall() == [('ok',)],
                       'cache_integrity')
            actual = {table: connection.execute('select count(*) from ' + table).fetchone()[0]
                      for table in tables}
            require325(actual == NATIVE_TABLE_COUNTS325, 'cache_native_table_census')
            compare._ranks325(connection, 'domain', arm, 256)
            family_rows = connection.execute('select kind,key,count(*) from route where arm=? '
                                             'group by kind,key order by kind,key', (arm,)).fetchall()
            require325([(kind, key) for kind, key, n in family_rows] == self.expected_family_keys
                       and all(n == 256 for kind, key, n in family_rows), 'cache_exact_provider_families')
            context_rows = connection.execute("select kind,key,count(*) from context where arm=? "
                                              "and kind!='inverse-bridge' group by kind,key order by kind,key",
                                              (arm,)).fetchall()
            require325(context_rows == family_rows, 'cache_context_family_census')
            bridge_rows = connection.execute("select key,count(*) from context where arm=? and "
                                             "kind='inverse-bridge' group by key order by key", (arm,)).fetchall()
            require325(bridge_rows == [(str(step), 256) for step in range(1, 5)], 'cache_bridge_census')
            for table in ('route', 'context'):
                require325(connection.execute('select count(*) from ' + table +
                                               ' where typeof(rank)!=\'integer\' or rank<0 or rank>=256').fetchone()[0] == 0,
                           'cache_rank_range')
            bad = connection.execute("select count(*) from feature where typeof(key)!='text' or "
                                     "typeof(digest)!='text' or length(digest)!=64 or "
                                     "digest glob '*[^0-9a-f]*'").fetchone()[0]
            require325(bad == 0, 'cache_feature_schema')
            return dict(passed=True, arm=arm, row_counts=actual, provider_families=len(family_rows),
                        inverse_contexts=1024, database_bytes=identity['bytes'],
                        database_sha256=sha256325(database))

    @semantic_phase325
    def cache(self, stage, allowance):
        stage = self._stage(stage, allowance)
        require325(allowance > CACHE_RECEIPT_RESERVE325, 'cache_receipt_reserve')
        source_proof = self.bind_sources()
        import parallel_replay325 as parallel
        reader_module = self._activate_reader()
        scratch = stage / 'cache-build'
        scratch.mkdir(mode=0o700)
        original_worker = parallel._worker325
        parallel._worker325 = native_cache_worker325
        try:
            with tracked_pool_context325(parallel, self.child_exit_statuses):
                report = parallel.cache_shards325(parallel.tasks_from_results325(self.original), scratch,
                                                 workers=self.workers, original_root=self.original,
                                                 anchors=self.anchors,
                                                 limit=allowance - CACHE_RECEIPT_RESERVE325,
                                                 reader_module=reader_module)
        finally:
            parallel._worker325 = original_worker
        observed = report['capture_reads']
        expected = {name for name, row in self.anchors.items()
                    if row['family'] == 'capture_inventory' and name.split('/', 1)[0] in ARMS325}
        require325(set(observed) == expected and len(observed) == 512, 'cache_actual_read_census')
        for name, row in observed.items():
            anchor = self.anchors[name]
            require325(row == dict(compressed_sha256=anchor['sha256'],
                                   stored_bytes=anchor['stored_bytes'], complete=True),
                       'cache_actual_read_hash')
        worker_bindings = self._worker_bindings(report, source_proof)
        proofs, artifacts = {}, {}
        for arm, artifact, filename in zip(ARMS325, ('left_db', 'right_db'), ('left.sqlite', 'right.sqlite')):
            path = Path(report['databases'][arm])
            require325(path.parent.resolve() == scratch.resolve() and not path.is_symlink(),
                       'cache_database_path')
            require325(report['merges'][arm]['row_counts'] == NATIVE_TABLE_COUNTS325,
                       'cache_merge_census')
            require325([r['rank'] for r in report['streams'][arm]] == list(range(256)),
                       'cache_stream_rank_census')
            proofs[arm] = self.validate_cache(path, arm)
            require325(proofs[arm]['database_sha256'] == report['merges'][arm]['database_sha256'],
                       'cache_merged_bytes_changed')
            os.rename(path, stage / filename)
            artifacts[artifact] = filename
        require325(not any(scratch.iterdir()), 'cache_unexpected_sidecars_or_shards')
        scratch.rmdir()
        receipt = dict(passed=True, source_job=55569, optimized_streams=512,
                       actual_capture_reads=observed, cache_census=proofs,
                       expected_provider_keys_sha256=digest325(self.expected_family_keys),
                       source_binding=self.bind_sources(), worker_source_bindings=worker_bindings,
                       checked_reuse=self._checked_reuse,
                       child_exit_statuses=self.child_exit_statuses,
                       merge_manifests={arm: {key: value for key, value in report['merges'][arm].items()
                                             if key != 'database'} for arm in ARMS325})
        bounded_json325(stage, stage / 'cache-receipt.json', receipt, allowance)
        artifacts['receipt'] = 'cache-receipt.json'
        require325(tree_bytes325(stage) <= allowance, 'cache_final_budget')
        return dict(artifacts=artifacts, census=dict(passed=True, streams=512, capture_reads=512,
                                                   arms=proofs))

    @semantic_phase325
    def work(self, stage, allowance):
        stage = self._stage(stage, allowance)
        source_rows = self.original_binding['work_inventory']
        raw_bytes = sum(row['stored_bytes'] for row in source_rows)
        ledger = json.loads((self.original / 'results.json').read_text())
        vector_bytes = sum(ledger['runs'][arm]['work_ledger']['per_rank_artifact']['bytes']
                           for arm in LABELS325)
        # All actual writes are additionally checked against allowance below.
        # The deterministic serializer is byte-bound to the original source.
        require325(raw_bytes + 2 * vector_bytes + 8 * COMPACT_LIMIT325 <= allowance,
                   'work_phase_preallocation')
        shadow = stage / 'work-input'
        shadow.mkdir(mode=0o700)
        for arm in LABELS325:
            (shadow / arm).mkdir(mode=0o700)
        copied_bytes = 0
        for row in source_rows:
            original = self.helper.safe(self.original, row['relative'])
            self.helper.stat_match(original, row)
            target = shadow / row['relative']
            require325(target.parent in tuple(shadow / arm for arm in LABELS325) and
                       copied_bytes + row['stored_bytes'] <= allowance, 'work_copy_budget')
            digest, count = hashlib.sha256(), 0
            with original.open('rb') as source, target.open('xb') as destination:
                while data := source.read(1024**2):
                    count += len(data)
                    require325(count <= row['stored_bytes'], 'work_copy_extent')
                    destination.write(data)
                    digest.update(data)
                destination.flush()
                os.fsync(destination.fileno())
            self.helper.stat_match(original, row)
            require325(count == row['stored_bytes'] and
                       digest.hexdigest() == self.anchors[row['relative']]['sha256'], 'work_copy_hash')
            copied_bytes += count
        import account325 as account
        with bounded_account_writes325(account, stage, allowance):
            report = account.qualify(shadow, Path(__file__).resolve().parent)
        require325(report.get('passed') is True, 'work_replay_failed')
        observed, run_reports = {}, {}
        for arm in LABELS325:
            run_path = shadow / (arm + '-work325.json')
            run = json.loads(run_path.read_text())
            vector_path = shadow / arm / 'work325-rank-vectors.json'
            require325(Path(run['per_rank_artifact']['path']).resolve() == vector_path.resolve()
                       and sha256325(vector_path) == run['per_rank_artifact']['sha256'],
                       'work_vector_certificate')
            vectors = json.loads(vector_path.read_text())['per_rank']
            require325([row['rank'] for row in vectors] == list(range(256)), 'work_rank_census')
            for rank, row in enumerate(vectors):
                name = f'{arm}/work325-rank-{rank}.tsv'
                anchor = self.anchors[name]
                require325(row['sha256'] == anchor['sha256'] and row['bytes'] == anchor['stored_bytes'],
                           'work_actual_read_hash')
                observed[name] = dict(sha256=row['sha256'], bytes=row['bytes'], complete=True)
            compact = dict(run)
            compact['per_rank_artifact'] = {k: v for k, v in run['per_rank_artifact'].items() if k != 'path'}
            compact['per_rank_artifact'].update(retained=False, role='verified derived vectors discarded')
            run_reports[arm] = compact
        require325(len(observed) == 768, 'work_actual_read_census')
        # Compact reports preserve all accounting summaries and vector hashes;
        # no path points to a transient directory after its explicit deletion.
        result = {k: v for k, v in report.items() if k != 'runs'}
        result['runs'] = run_reports
        result['actual_read_hashes_matched'] = 768
        bounded_json325(stage, stage / 'work-report.json', result, allowance)
        bounded_json325(stage, stage / 'work-read-checks.json', observed, allowance)
        remove_owned_tree325(stage, shadow)
        return dict(artifacts={'report': 'work-report.json', 'actual_reads': 'work-read-checks.json'},
                    census=dict(passed=True, work_reads=768, copied_tsv_bytes=raw_bytes,
                                derived_vector_bytes=vector_bytes, performance_qualified=False))

    @semantic_phase325
    def preservation(self, stage, allowance):
        stage = self._stage(stage, allowance)
        previous = dict(self._last_binding)
        current = self.binding()
        require325(current == previous, 'final_source_checker_config_or_capture_changed')
        report = dict(passed=True, original_source_job=55569, binding=current,
                      full_input_hashes_matched=1536, metadata=self._bound_metadata,
                      checked_reuse=self._checked_reuse, source_binding=self.bind_sources(),
                      original_ledger_preserved=True, solver_runs=0, compiler_runs=0,
                      performance_runs=0, performance_qualified=False)
        bounded_json325(stage, stage / 'preservation.json', report, allowance)
        return dict(artifacts={'report': 'preservation.json'},
                    census=dict(passed=True, input_hashes=1536, original_ledger_preserved=True))
