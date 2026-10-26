# Stage325 verification workflow repair, 2026-10-26

The subsequent native job 55576 and independent review passed all six phases for
the saved fixed level-5, tol=0, 256-rank, four-step case. The authorized conditional
verification commit is supported. See
[native qualification](parallel-blocks-stage325-native-v5-qualification.md).
The following account records the earlier offline repair; it does not establish
a performance benefit or the cause of the original negative scaling.

At the offline milestone, this change repaired the interrupted verification
workflow without establishing the locality map's scientific qualification,
performance benefit, or the cause of the original 160-to-256-task slowdown. No
new native job, solver, compiler or performance run had then been launched. The
original campaign had used2656of2700allocationseconds; its remaining44seconds
were unchanged. Conditional commit and performance evaluation were held pending
the later native qualification reported above.

The adversarial design review preceded implementation. It rejected reuse of
partially merged SQLite files: the existing builders disable rollback journals,
and a completed progress count does not prove a complete database or provider
certificate. A second review challenged the composed implementation. In
particular, equality can pass when both caches lose the same feature; the new
workflow checks the exact native census before accepting cache completion.

## Completion and reuse boundaries

`test/parallel_block_profile/locality_325/resume_workflow325.py` orchestrates six
phases: the pair of caches, C0 providers, L256 providers, pair comparison, work
replay, and final preservation. Each phase writes in a private staging directory.
Only closed files with an exact inventory, full SHA256, bound input/checker/config
identity and dependency certificates become reusable. Certificate publication uses
an atomic directory rename after file and directory fsync. Progress records and
the attempt report have no qualification authority.

Resume is explicit. It freshly hashes protected inputs, then rehashes every
completed artifact. An uncommitted phase requires a separate explicit partial
discard; no shard or interrupted merge is resumed in place. Unknown files,
symlink ancestors or leaves, hardlinks, changed bytes with preserved size/mtime,
incompatible configurations, concurrent writers and semantic-failure latches are
rejected. A semantic failure creates an immutable fail-closed latch before its
diagnostic write. Even an interruption while writing that latch blocks resume.

Native cache completion requires each arm to contain exactly65024route,
256Domain,66048context and508840feature rows, integrity/schema/rank checks, and
512complete source-bound stream/actual compressed-read receipts. Provider family
identities must match both historical byte-bound optimized rank0 inventories,
not merely a count or agreement between unbound candidate worlds. The native
preflight corrected the earlier checked-inventory premise:170familykeys differ
between checked and optimized transform epochs. Exact optimized keys are retained
without normalization. Native ranks/steps/families remain256/4/254.

## Bounded finalization

`stream_compare325.py` reuses the existing provider reconciler and physical inverse
bridge checker. It persists each compact254-family proof immediately after that
provider closure. The subsequent pair phase reads sorted feature cursors and one
bridge pair at a time; it does not export or reload all feature rows and bridge
contexts as JSON. It checks the complete1024bridge census and retains complete
mismatch counts, rolling digests and at most32diagnostic examples. Row-stream
digests have an explicit new contract and are not equated with the old JSON
dictionary digest.

SQLite reads explicitly enable URI interpretation, use read-only immutable mode
and query_only, and reject transaction sidecars and WAL headers before opening.
The phase store supplies the full-file byte binding. The native adapter reuses
the existing exact stream gate, worker/shard reconciliation and work accounter;
original captures, source, executable and checkpoint data remain read-only.
Work staging uses owned copies, never symlinks into the phase store.

The default4GiBcap covers logical regular-file bytes, retained payloads and
in-progress derived files. Phase reservations leave metadata headroom; twoMiB
are withheld for the bounded current attempt report and its atomic temporary
file. The report retains bounded prior attempt summaries so resumption preserves
earlier signals. Reservations are accounting, not an operating-system disk
quota. SQLite page ceilings and preflight work/report write bounds enforce
external-writer allowances; a future native wrapper must retain aggregate memory
and allocation enforcement. Filesystem allocated blocks are not claimed as the
logical-byte measure.

SIGTERM, SIGINT and the explicit runtime alarm record numeric signal, name, UTC,
monotonic time and active phase. Signal handlers do no filesystem work; unwinding
records the event while holding the state lock. Before ownership, a bounded
stderr event records termination without overwriting another writer's attempt
report. SIGKILL cannot be logged by a handler: a later resume retains the last
unclean attempt as unknown termination, validates certificates and discards only
the explicitly authorized partial phase. Signals never trigger automatic retry.
Phase timings include commit and full artifact rehash cost.

## Validation and remaining gate

Tests use real SQLite databases, the unchanged route reconciler and strict inverse
checker. They include full256x4bridge census, exact old/new verdict agreement,
254provider families, bounded cursor memory, same-feature deletion on both arms,
real subprocess SIGKILL at publication boundaries, SIGTERM/alarm identification,
explicit resume/discard, skipped completed phase callbacks, altered source,
checker, config, capture and database bytes, malformed certificates, concurrent
access and storage reservations. These are checker software tests using
constructed fixtures, not new solver or native campaign results.

The adversarial findings, final test counts and timings are recorded alongside
the code under `validation/` and `RESUME_VALIDATION.json`. Fixture timings do not
justify a bbserv runtime estimate. Original55574derived databases were deleted
by the old cleanup and cannot be recovered by this repair. A future complete
native verification therefore needs a fresh derived-cache build from preserved
captures, a separately agreed time budget, runtime preflight, aggregate resource
supervision, bounded collection and independent review. Only that full pass can
release the conditional commit. No performance campaign is implied.

Final local validation passed258tests in34.739seconds, including the183existing
tests and75newtests. An independent agent reran all75newtests in24.263seconds,
confirmed unchanged before/after code hashes, and recorded no open offline
findings in `validation/adversarial-resume325.json`. The saved checked-L256rank0
probe extracts all254families with its historical compressed-byte hash and four
actual source pins. Its initial route-key-suffix parsing failure, correction and
passing probe remain recorded. All62original v4p1recipe source identities and
the sealed manifest are preserved. These local checks do not release the native
qualification or conditional commit holds.
