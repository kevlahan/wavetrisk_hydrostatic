Final status: qualification PASSED. The authoritative current result is
`parallel-blocks-stage236-local-results.json` and the completed handoff is
`parallel-blocks-handoff-236.md`. The notes below preserve development history;
their pending/failed statuses describe earlier source versions.

# Stage236 development history — superseded by final acceptance

Authoritative worktree remains 7173 on codex/performance-investigation, HEAD
4ea640128a71a857f9b8c6b9411dc34fc0aeb641. Main unchanged. No commits or staging.
The earlier Stage236 preflight handoff is now historical: implementation has begun.

BLOCK_RK defaults off. Current source edits join RK closure to sparse inverse
aliases, bypass production q1 arithmetic/publication, remove q1 numerical storage,
retain mass tendencies on owners, and bypass full mass/velocity RK returns.
The current development tree is .stage236-validation/candidate; it is mutable,
not qualified frozen evidence. Stage235 qualified-candidate is untouched.

Corrections discovered by checked runs:
- Fixed coarse scaffold consumers need explicit registration in inverse metadata.
- Batched STAGE203 mass consumption must read retained owner tendencies too.
- Canonical velocity aliases require exact base+h*gradient from existing block
  checkpoint; recomposed block candidate has different last-bit rounding.
- Provisional low scaffold retains q1 initialization/previous provisional epochs,
  distinct from global sol fixed scaffold reset. Sparse keyed low-scaffold state
  now preserves this lifetime and resets with init_RK_mem.
- Initial adaptation oracles run outside RK; retained scaffold reads must be
  limited to retained_rk_cycle_active, not the build selector alone.

Further corrections:
- Scaffold packing must honor zmin:zmax, including soil/surface offsets; only
  physical layers use the retained RK scaffold. This fixed the first inverse
  mismatch. All four stages then passed bitwise alias comparisons.
- Final coarse scaffold must capture the completed signed boundary epoch. The
  earlier pre-boundary capture failed a final full field comparison in the
  runtime oracle. oracle-eighth-4 is checking this correction.

The development production-first-4 trajectory is still the pre-boundary version;
it has matched at least eight step states but is not final source evidence.
pre-boundary-candidate is a stopped, superseded frozen-build attempt.
qualified-candidate is the new corrected source freeze; builds are in progress.
No complete acceptance pass exists yet. No remote jobs or commits were made.

Next: finish the corrected oracle check, full frozen 4/41 trajectories, native
and feature-off controls, poison and fault rejection; then run review_rk.py.
Do not claim completion or speedup until all required gates pass.

Latest qualification progress:
- oracle-eighth-4 passed all four exact alias boundaries and final transform.
- production-first-4 passed 12 steps/2 remaps/exact restart, checkpoint and 62
  VTK members, but is the prior source version and excluded from final evidence.
- Corrected qualified-candidate has all four binaries and build-236.json; all
  69 inputs, 115 resolved sources and four binary identities were verified.
- qualified-native-4 and poison-4 passed; both step states match the reference.
  Poison leaves original registered atmospheric sol inputs NaN throughout all
  four RK closures, with no Domain backup/restoration. q1 and native RK arrays
  are absent in the checked native route.
- reject-source, reject-missing and reject-epoch passed specific diagnostics.
- qualified-production-4, qualified-all-oracles-4 and reject-base are running.
  The live command groups then run qualified-production-41, control-4 and
  reject-fallback. No final acceptance report exists yet.

Final gates now passed: frozen qualified-production-4 (12 steps, two remaps,
exact checkpoint/restart fields, exact step states and 62 VTK members), frozen
qualified-all-oracles-4 (four exact alias stages plus final transform), and all
five fault cases including base and fallback. qualified-production-41 and
control-4 remain active. Rank40 reports zero Domains/aliases and 304 owned RK
operands; its first completed step is exact. After these finish, run review_rk.py,
replace the preflight handoff with the final handoff, update the implementation
document's qualification status, and run git diff --check. Sources must remain
identical to qualified-candidate; defaults off; no commits/staging or remote jobs.

Completed: qualified-production-41 and control-4 passed, and review_rk.py
accepted all eleven cases against frozen sources and binaries. Both full runs
have 12 exact steps, two remaps, exact restart/checkpoint fields and 62 exact VTK
members. No commits, staging, main changes or remote jobs were made.
