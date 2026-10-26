# Stage325: exact locality-map implementation and qualification gate

2026-10-26. The audited 256-owner permutation is implemented as the separately
guarded `WAVETRISK_OWNER_MAP322=L256` selector. Local compilation, initialization,
codec and checker tests pass. **Native qualification passed in job 55576 and
independent review for the saved fixed level-5, tol=0, 256-rank, four-step case.**
The original job 55569 supplied science, placement and exact compiler gates;
the corrected replay completed semantic, provider, physical inverse and work
verification without another solver run. No performance improvement, primary
slowdown cause or accepted-solver promotion is established. See
[full native review](parallel-blocks-stage325-native-v5-qualification.md) and the
[original campaign review](parallel-blocks-stage325-job55569-review.md).

The isolated branch is `codex/stage325-locality-qualification`. Commit
`462600109b281af30e34e5f71832fe61cc2a51fb` first records the exact 258-file candidate
source used in Stage322. Those bytes match the frozen source manifest qualified
by job55567 and measured by job55568. The subsequent Stage325 change is separate
from that snapshot and from the accumulated changes in the research worktree.

## Implemented mechanism and guards

L256 composes the current 640-block map with the exact audited 256-label bijection.
It preserves all 256 nonempty numerical groups and coownership relations, and
changes Domain-home blocks from 6 to 400. It targets migration serialization and
Domain-facing transfers; it does not reduce peak numerical work or the internal
numerical halo graph. Its competing hardware exposures and unchanged geometry
preparation remain part of the assessment.

The existing off/A160/B256/C256 selectors remain intact. The native planner's
exact map, 640 physical catalogue entries, 160 fixed Domain homes, 256 ranks,
level 5 / tol 0, 30 layers, Simple physics, RK4, MODE1 and subdivision 2 are checked before
L256 substitution. Committed/group maps remain checked at every construction.
The new permutation guard and selector operations are included in work accounting.

## Qualification evidence produced by the campaign

The same optimized executable runs the off control and L256. A separate checked
build runs L256 with the existing numerical/reference checks. Each of the three
runs includes initialization and four steps. There are five catalogue/migration
frames and seven adaptation calls; their validators deliberately distinguish them.

The semantic gate retains the corrected Stage322 closed coverage rules and adds
physical-key provider, address and ordered-write verification. It verifies all 48
changed storage orders across all five frames, ordered numerical programs and
consumed state under the rank permutation, fixed Domain programs/state, numerical
backup values, local/imported publication, retained blocks and live boundary
imports. Source-side finite-value checks precede SHA256 encoding of real arrays;
integer metadata remains exact. This is digest-based bit-identity evidence, not
retention of complete raw floating arrays. Unknown scopes, missing fields and
unsupported providers fail closed.

The ledger records actual full serialized block sizes, local/remote migration,
storage moves, unchanged geometry preparation, metadata and numerical routes,
mixed publication imports, scratch, writeback and boundary-import replacement
work. Directed peer payloads reconcile across ranks. Compact review files retain
separate phase/rank counters and actual-placement cross-socket/shared-L3 endpoint
burdens. No independent phase maxima are added into an elapsed critical path.
Full native peer vectors remain hash-bound and available for targeted review.

Science/checkpoint, actual solver placement and exact GCC13.2.0 C/Fortran plus
MVAPICH4.0 runtime/linkage checks are separate required gates. The new control
must match its qualified predecessor exactly; optimized control/candidate initial
and final checkpoints must match exactly. All runs retain the independent legacy
trajectory and established science tolerance checks.

## Local validation and its limits

- Complete checked J5 solver compilation passes with runtime checks and warnings
  as errors, using local GCC15/OpenMPI. This does not substitute for the native
  compiler gate. One missing-symbol error in an early accounting hook was fixed;
  subsequent compilation and all final hooks pass.
- Two private four-rank, zero-step initialization smoke runs passed during
  integration; the final binary completed in 21.08 seconds. No local performance
  campaign or positive-step solver trajectory was run. The L256 selector requires
  256 ranks and was not activated by these initialization tests.
- The compiled selector reproduces all five selector outputs, rejects an invalid
  selector, and proves the sealed permutation/home counts. The C/Fortran codec
  passes 20 SHA256 vectors, nonfinite rejection, empty/noncontiguous input and exact
  physical footer-length checks.
- All 60 offline tests pass, including malformed capture/footer/shape tests,
  physical-provider and alias-order mutations, boundary-import corruption,
  work/peer conservation, hardware endpoint projection, scheduler failure and
  compact-collection recovery. Source patch reproduction matches all 258 files.
- The final four initialization streams pass framing and peer-ledger checks:
  3,556,108 physical capture bytes in total. This does not bound four-step native
  capture volume or qualify the new owner map on target hardware.

The source/checker reviews found and closed local-publication ownership,
writeback/retained coverage and live boundary-import gaps before sealing. Native
execution can still reveal a contract defect or a correctness failure. Stop and
review such a failure; do not weaken a gate to obtain a passing result.

## Bounded bbserv campaign

The kit is `rk325-locality-qualification-kit-v1`. It schedules two serial full
solver builds and three 256-rank qualification launches on exclusive bb13/nk,
plus the existing small compiler/callback/runtime and placement gates. There are
zero performance runs. Planning estimate 15–30 minutes is uncertain; the fixed
allocation cap is 45 minutes. There is no automatic repeat or budget extension.

The original private resume=-1 seed and independent legacy reference are reused
through the pinned completed 55567 evidence. New runs receive private read-only
copies. Missing/changed retained evidence stops preflight. Original user
checkpoints and earlier source, binaries and captures are preserved.

Per-rank compact streams have a 128 MiB physical hard limit. The runner monitors
64 MiB compressed per rank, 2 GiB per run and 3 GiB across runs every 0.5 seconds;
these are stop thresholds, not atomic aggregate quotas. SQLite verification
scratch is limited to 4 GiB and progress/observed usage is reported. The collector
streams at most64 MiB of selected metadata without a staging copy. Raw captures,
checkpoints, work TSV and full peer vectors stay native; compact review counters
and hash-bound reports return. Timeout preserves evidence and does not require
rerunning the solvers merely to inspect it.

The campaign must complete the semantic, science, placement and work gates and
receive independent review before any performance comparison is considered.
Passing would qualify this fixed-case candidate, not demonstrate a speedup.

## Reproduction

Tooling, tests, map certificates and the source baseline manifest are in
`test/parallel_block_profile/locality_325`. `prepare_kit.py CHECKOUT OUTPUT`
packages only the manifest-listed source inputs into a fresh sealed kit. The
patch generators independently reproduce the candidate from the frozen baseline.
`LOCAL_VALIDATION.json` binds the local evidence and preserves the native-pending
status. See the kit README for launch/status/collection instructions.
