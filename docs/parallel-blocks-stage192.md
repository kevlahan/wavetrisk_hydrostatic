# Stage 192: send written outer-vector edges

Stage 191's matched whole WT transaction was 4.33 times the legacy cost in job
54561. Gather/alias/scatter helpers consumed 57% of native inverse inclusive time;
scatter accounted for 25.35 GB of the 42.26 GB of remote inverse payloads over the
profile window. These are overlapping scopes, not additive costs or guaranteed
recoverable savings.

This experiment narrows the outer-vector scatter from geometric route owners back
to final block owners. Regular scalar and inner-vector kernels already run on
final owners. Sending unchanged companion vector edges back is unnecessary.

`parallel_block_inverse.f90` marks the exact edges written by outer operations
(targets 2 and 3) and signed local/remote aliases. It sends this conservative union
of edge masks in the existing topology-key exchange. Masks include all levels;
existing fine/coarse selection still chooses each epoch. The outer plan
expands only written interior/boundary edges into wire records, preserves complete
scaffold resets, and retains one exchange/completion per phase. No vector kernel
intervenes since the previous vector synchronization or initial seed; other edges
still hold the gathered value. Arithmetic,
alias snapshot order, local/remote alias precedence, edge signs and final
consumer publication are unchanged. The existing maximum-sized buffers suffice
because the new wire payload cannot exceed the old three-edge payload.

`parallel_block.f90` accepts an optional selected edge per compiled route; old
callers retain the original packing/install behavior. Integer routes are compiled
per generation, without persistent pointers into swappable numerical buffers.
The edge schedule increases integer routing metadata and can increase dispatch
work; cluster RSS and complete-step timings must determine whether this pays off.
It does not eliminate gathers, geometric aliases or Domain consumer publication.

`WAVETRISK_STAGE192_SPARSE_OUTER=1` selects this experiment (candidate default).
0 restores Stage 191's full outer installation. Settings must agree on every rank
and are logged once. Stage 191 installation remains on and mass overlap off.
The sparse route also runs when phase oracles are enabled.

An initial node-only restriction passed local exact restart comparison but did
not reduce remote payloads: staged nodes were already relevant. The final
candidate restricts individual edges within those nodes. The interior-edge-only prototype reduced local scatter by 4.61%; applying the
same mask to boundary records reduces it by 7.57%. Prototype runs are not evidence
for the final edge schedule.

## Validation and cluster decision

Use `test/parallel_block_profile/local_inverse_192.py` with the same `--build`,
`--fixture`, `--out`, `--binary` arguments as local_j4; add `--control`, `--profile`
or `--oracles` explicitly. It records selector values in the run identity. Local
J4/J6 uses four ranks, 12 steps, resume 2, checkpoint 3/reload and post-restart
remapping. Source/build/run identities and comparison results are recorded in
`parallel-blocks-stage192-local-results.json` after validation.

Focused tests exercise actual Fortran route selection, sparse peer ordering,
selected-edge writes at atmospheric offsets in both storage families, and MPI
selector agreement/failure. The campaign tests cover build hashes, common physics,
large provenance through a file argument, scheduling, workload rejection and
fixed legacy numerical screening.

The cluster kit compares pinned Stage 191 `30c83114`, this candidate and pinned
legacy with the common I/O fix. See `parallel-blocks-cluster-192.md`. Require exact
block/control/checked checkpoints, the unchanged legacy screen, matched workload,
restart/remap coverage and six balanced timing triples. Complete WT transaction
and transfer evidence are collected separately from headline timing. Traffic
reduction alone does not justify accepting this optimization. Local elapsed
seconds are not performance evidence on the memory-compressed workstation.

Local validation passed: the final optimized and checked builds complete all 12
steps with checkpoint 3/reload and two remaps. Checkpoint fields match Stage 191
and the saved independent RK4 oracle exactly; the unchanged legacy screen also
passes. Fresh full phase-oracle checking remains a required cluster run.

The local profile records 564,272,640 -> 521,582,400 remote scatter bytes (7.57%
less). Total inverse remote payload falls 1.45%, with the same 12,628 messages;
gathers and aliases are unchanged. This is a modest traffic reduction, not proof
that the 2.4x ordinary-step gap has closed. Cluster timing must include the added
edge-route metadata and dispatch cost. The full Python/Fortran suite passes 107
of 110 tests with three optional MPI skips; the four focused Fortran/MPI tests
pass separately, including rank disagreement and strided edge installation.
