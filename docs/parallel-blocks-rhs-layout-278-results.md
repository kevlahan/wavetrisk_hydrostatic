# Stage278 resident input-layout development result

The actual resident input-address compiler is implemented in
`src/parallel_block_rhs_layout.f90`, with its Makefile dependency. It resolves
owned, boundary and ghost storage from `Block_Data` topology and stores one
seven-node footprint per owned node, shared across layers and input components.
It does not duplicate geometry or state arrays. Exact structural snapshots detect
stale topology reuse; field extents, layer bounds and neighbor identities are
checked. Topology changes rebuild the addresses. Numerical mask plans and
state/halo epochs remain separate obligations and are not qualified by this module.

The source is compiled in the development tree, but remains inactive in normal
production. Its gather routine is a validation consumer, not the final patch-loop
implementation. It performs no new production exchange and retires no Domain work.

## Local verification

Two full climate builds pass: optimized and strict checked. Three serial J4
replays pass with four MPI ranks: optimized, independent-oracle and checked.
They use an isolated Stage276 validation overlay on accepted266, replacing its
resident input gather with the new compiler and checking it against both the old
resident gather and independent Domain reference packets. These packets are
validation-only; this is not the proposed production transport or cluster source.
The normal current solver is not presented as a connected Stage278 replacement.

Each replay checks600,736 sampled node observations,9,611,776 inputs and5,406,624
primitive outputs, with zero differences or skipped mapped nodes. Input counts
are8,774,496 owned,467,212 boundary and370,068 ghost observations. Optimized,
checked and oracle science checks pass against the accepted-family local fixture,
including the adaptive/reload/remap trajectory. These are repeated observations
at the existing Stage276 selected calls/layers, not all-call/layer RHS proof.
There are no zero-Domain ranks in the four-rank J4 fixture.

Optimized and checked unit tests also pass: owned/internal-neighbor/ghost
addressing, layer layout, reuse and topology invalidation; five expected failures
reject stale generation, changed topology without an epoch, missing neighbor,
invalid layer and invalid field extent. Mask-only changes deliberately leave
this input footprint reusable; they do not authorize reuse of numerical masks.

Evidence: `.stage278-validation/layout-validation.json`, with complete sources,
binaries, logs and per-rank records in `.stage278-validation/layout-replay`.
Earlier setup/build failures and the initial sandbox socket failure are preserved.
The latter was resolved by an approved local four-rank execution; no bbserv job ran.
No local timing is claimed as a cluster performance estimate.

## Campaign status and remaining work

`connected_rhs_278/campaign-spec.json` records the proposed four-build/sixteen-run
nk/bb13 matrix and immutable J5 fixture hashes. `campaign_gate.py` checks the
source-bound local production qualification before sealing and rejects the
current draft. It does not submit jobs or constitute an approval mechanism.
No launchable kit, launch command or candidate performance result is supplied.

**The user's requested connected production changes and bbserv campaign remain
unfinished.** Still required: actual mask-aware numerical plans and geometry
placement; complete pressure/diffusion/boundary phases and derived-field routes;
authoritative resident RHS and real RK integration; paired Domain work/transport
retirement; replacement-cost evidence; full connected optimized/checked/oracle
qualification. The input-layout result earns only development integration credit.
Stage266 stays accepted and both project goals remain open.
