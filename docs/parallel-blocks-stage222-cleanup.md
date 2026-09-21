# Stage222 — remove dead Domain publications and prepare empty-Domain ranks

This stage implements architectural cleanup. It has no minimum speedup target
and makes no cluster performance claim. Stage215 remains the frozen cluster
reference; these changes are on the investigation branch, not main.

**Local validation passed.** Optimized, checked cleanup, oracle-enabled and
checked empty-Domain ownership trajectories all completed 12 steps, checkpoint
reloads 2 and 3, and two vertical remaps. All four saved checkpoints match
Stage215 exactly, including the changed ownership case. The latter needed no
numerical tolerance relaxation. The Git `main` branch remains at
`b82647b36eb825aed0c0aae4e92adcf4ab46d608`.

## Implemented changes

1. **Remove the early production-provisional inverse publication.** Native
   boundary work reads compact/workspace state; the later publication installs
   the version actually consumed. Final transforms and oracle boundary refresh
   retain their earlier publication. No inverse arithmetic, alias ordering or
   final-owner routing changes.
2. **Pack and stage only physical temperature for temperature-only writeback.**
   The compact scalar getter supports a bounded temperature slice. Remote sends
   no longer pack a whole scalar patch before slicing. Local reconstruction
   skips unused scalar/vector preparation and whole-array clearing; coarse
   temperature still comes from the timestep-base scaffold. MPI already sent
   only temperature, so this change does not claim a new wire-volume reduction.
   Full staging consumers continue to rebuild their own complete images.
3. **Reuse the final native publication at grid change.** The completed native
   transform has installed SOL and WAV_COEFF; its final refresh changes ghost
   and derived state. Production skips the two repeat full-family writebacks.
   The saved physics-consumer alias version is still installed afterward. The
   existing pending-alias/completion state identifies this interval; no new
   persistent cache or validity flag is introduced. Dynamics/adaptation oracle
   paths and other callers retain the full synchronization. Writeback counters
   count actual operations, including zero operations in the reuse case.
4. **Allow an empty source store when a rank owns no Domains.**
   `build_source_blocks` previously stopped before migration on such a rank.
   It now allocates empty source/index arrays and participates normally in
   migration. Missing source blocks on a rank that does own Domains still fail.
   Empty-source summary handling also covers rank zero.
5. **Load global restart metadata on ranks with no Domains.** Previously the
   simulation clock, step/output counters, thresholds and case-specific restart
   prefix were read only inside the local Domain loop. The ownership probe
   reached seven steps but then stalled: the empty-Domain rank entered timestep
   reduction while Domain owners entered checkpoint writeback. The new path
   broadcasts the existing serialized global prefix from the owner of Domain
   zero and decodes it on empty-Domain ranks. It does not replicate Domain
   fields or change the checkpoint format. The extra broadcast is used only
   when the ownership map contains an empty rank.
6. **Run shared restart cleanup even with an empty Domain array.** After the
   metadata fix, the checked probe wrote an exact checkpoint but stopped during
   restart because vertical-coordinate arrays were still allocated. The early
   `size(grid)<=0` return in `deallocate_structures` skipped all shared cleanup.
   Removing that return lets empty Domain loops do no work while shared arrays
   and zero-sized field containers are freed normally.

## Evidence and validation design

The preliminary optimized and checked runs containing the first two deletions
completed the 12-step four-task J4/J6, 30-layer adaptive/restart fixture. Both
matched the accepted Stage215 checkpoint exactly. A separate production-mode
probe checked existing Domain interiors against reconstructed compact SOL and
WAV_COEFF **before** the old grid-change writeback: all 48 rank-step checks
passed across the trajectory. This motivated the third deletion.

The final combined candidate is validated using:

- An uninstrumented optimized trajectory and a checked trajectory, each compared
  with frozen Stage215 fields, thresholds, topology and printed adaptive states.
- Poisoning the omitted noncanonical inverse stores and unused temperature
  staging lanes on every step. The checked build enables bounds checks and
  floating-point traps. It also checks exact Domain/compact agreement whenever
  the two grid-change publications are omitted.
- Byte-exact finite intermediate streams for the first complete RK4 step:
  inverse Domain output, native mass and velocity tendencies, Domain RK state,
  compact scalar/vector RK state, and compressed forward output. These are
  bounded correctness observations, not timers or a broad profiling campaign.
- An oracle-enabled optimized run, retaining the independent comparison paths.
- A four-task ownership probe with Domains assigned to ranks 0–2 and the real
  block distributor free to use all four ranks. It requires rank 3 to receive
  nonzero owned blocks, complete the adaptive trajectory, and survive checkpoint
  reload/remap. This is distinct from testing an idle communication participant.
  The checked probe also poisons discarded grid-writeback staging and transport
  buffers without reconstructing the omitted writes, detecting hidden consumers
  that the separate Domain/compact comparison could otherwise mask.

The checked cleanup trajectory predates the restart-metadata fix; its 28 finite
intermediate streams (1,636,639,904 bytes) and final checkpoint match Stage215
exactly. The final optimized/oracle trajectories and the checked ownership probe
test the additional restart metadata fix. The checked ownership probe also
tests the subsequent empty-rank cleanup correction; normal-layout runs predate
that correction, which affects only the empty-rank early return. The validator
verifies this exact source difference. The case callback is empty for this climate
fixture; other case-specific payloads are preserved by the serialization path
but have not been exercised by these tests.

Final outcomes and source/binary identities are recorded by
`test/parallel_block_profile/validate_cleanup_222.py` in
`parallel-blocks-stage222-cleanup-results.json`. It requires 12 steps and reloads
2 and 3. Same-layout cleanup comparisons require exact fields. The changed
geometry-owner probe uses the existing conservative numerical policy if exact
equality is lost, with identical topology, thresholds and adaptive trajectory;
its actual exactness is reported rather than assumed.

The checked cleanup run exercised 144 poisoned provisional publications, 192
poisoned temperature transactions and 48 exact checks before skipped grid
writebacks. The final ownership run passed 48 independent discarded-buffer
poison checks. Rank 3 owned zero Domains and 36–38 blocks across 14 installation
events, including restart and adaptation. Its solver exited normally; the
launcher's already-imported old checkpoint reader failed postprocessing, which
was then rerun with the corrected reader. The original record and recovery
provenance are retained. The checkpoint-reader regression suite passed all four
tests. No local elapsed times are used to claim speedup.

Audit and probe code is injected only into disposable test trees by
`audit_cleanup_222.py`, `probe_grid_overlap_222.py` and
`probe_domainless_222.py`. Production has no new trace/poison switches. Local
tests include saving and restart intentionally for correctness; their elapsed
times are not performance evidence.

The semantic checkpoint reader now validates physical record coverage after
sorting directory offsets, while continuing to compare payloads by Domain ID.
Rank-ordered checkpoint output can legitimately permute records when ownership
changes. A regression test proves exact comparison under record permutation and
rejection of overlapping records; corruption checks are retained.

## Capability boundary and next gate

The startup `n_process > N_GLO_DOMAIN` guard is unchanged. This stage does not
claim a production run above 160 tasks. The ownership probe isolates one
architectural requirement without changing that guard or fabricating final
block ownership.

The next capability campaign should test an explicitly bounded startup candidate
above the fixture's Domain count: 41 tasks for a 40-Domain fixture, then 161 or
167 for the 160-Domain production fixture. Require useful owned-block execution
on every newly available rank, exact structural invariants, the established
numerical comparison policy, adaptation, remap and restart. Include ranks with
zero Domains and, where applicable, zero blocks without skipping collectives.
Only then remove or replace the production startup restriction. Speedup is
secondary to correct execution for this capability gate.

The first above-limit prototype should change the startup guard only in a
disposable candidate and retain normal Domain and block distribution. Include
both an empty last rank and a deliberately empty rank zero in the bounded
ownership tests, since the metadata broadcaster need not be MPI rank zero.
Record actual Domain and block counts per rank after each redistribution.

Remaining Domain representations are not certified minimal. In particular,
uncompressed vector wavelets still pass through Domain storage before compact
import; primitive/tendency producers and external physics/adaptation consumers
still require Domain-addressed values. Removing those representations needs a
direct producer-to-consumer replacement with explicit alias/version and geometry
contracts. The present cleanup does not establish that all such work is dead.
