# RHS phase exchange development — campaign remains unfinished

This development step addresses the missing distributed final-writer mechanism.
It does **not** complete the requested connected RHS/RK cutover or bbserv kit.
Stage266 remains accepted; production is inactive and nothing is promoted.

`parallel_block_rhs_exchange` compiles physical storage identities, accepted
writer priorities and consumer orientations into deduplicated peer routes.
The distributed directory exists only during compilation. Ordinary exchanges
visit participating peers, send one message per peer, and perform direct local
copies after snapshotting all source values. They use neither global collectives
nor arithmetic dependency-depth dispatch. Every peer message carries a state
stamp; topology, mask, phase and state mismatches are rejected. Plans own a
communicator and must be explicitly released before rebuilding.

`rhs_schedule_definitions` connects existing flux/PV/vorticity writer records to
that interface. Direct writes precede their immediate SW writes; pentagons follow
all accepted Domain-level patches, including patches assigned to other owners.
The caller must still distribute real physical alias identities and accepted
patch ordinals. The implementation cannot substitute coordinates or local maxima.
Actual production reader binding and transformed inter-Domain halo semantics are
not qualified by these synthetic tests.

## Evidence

[Source-bound results](../.stage278-validation/exchange-validation.json) preserve
303 files: current and earlier builds/runs, the initial strict-compiler failure,
a partial Makefile setup failure, corrected results and development scripts.
The preceding 26-file draft manifest and its exact files are preserved separately.

- Optimized and checked MPI tests pass at one and four ranks. Each normal run
  checks three nonempty epochs and 15 transactions, then a completely empty
  workspace. Cases include an empty rank, a consumer-only rank, component-specific
  final writers, reordered definitions, repeated aliases, negative orientation,
  local overwrite cycles and rebuild/release.
- Fourteen negative cases per build reject missing/ambiguous writers, bad
  destinations/identities/extents, invalid recompilation, stale local epochs and
  mismatched distributed mask/state epochs.
- Large-message tests pass at four ranks for three transactions and 30 layers:
  8,192 values per peer, 1,966,080 payload bytes plus an 8-byte state stamp per
  message; 2,949,120 received values checked exactly per build. This exercises
  payloads beyond the ordinary small-message case; it is not a speedup measurement.
- The schedule tests pass in optimized/checked builds, including all four pentagon
  writer exports, direct/SW precedence and noncontiguous global patch ordinals.
- Both changed modules compile through the actual Makefile in optimized and
  strict checked modes. The corrected partial-build command explicitly builds
  `param_J4` and uses the existing generated physics include directory. This is
  a module build, not a new full climate replay.

The earlier three exact J4 replays still describe their archived schedule source.
They are **not** retroactively attributed to this changed draft or to an active
connected production RHS.

## Cost and dual-goal decision

Ordinary outgoing bytes are `8 * (layers * unique_remote_values + send_peers)`;
the second term is the state header. Multiple aliases on one consumer rank share
one transmitted value. Local aliases use no MPI payload. These are structural
counts, not J5 traffic measurements. Compilation includes metadata collectives,
a distributed writer directory and communicator creation. Persistent index/buffer
payload is reported; transient compilation storage, MPI internals and peak RSS
remain to be measured on the connected candidate.

This mechanism supports block-owned derived fields without introducing the
Stage251–255 per-value/depth interpreter. It earns **zero** production ownership,
retirement, performance or P>=4D credit. There is no new defensible runtime gain.
The previous quantitative sensitivity case and matched legacy gaps are unchanged.

Still required: bind real topology/masks/metrics and producer/reader identities;
connect complete pressure, primitive, diffusion, restriction and physics phases;
feed authoritative mass/velocity RHS to existing resident RK consumers; retire
matched Domain computation/publication; then run active optimized/checked/oracle
qualification and seal the source-bound campaign. The entry gate still rejects
the draft. No launchable kit has been produced and no cluster job was launched.
