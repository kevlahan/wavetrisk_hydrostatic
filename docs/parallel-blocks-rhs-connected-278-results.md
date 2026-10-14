# Connected RHS arithmetic and placement metadata — Stage278 development

Two four-rank checked J4 replays pass against accepted266: the first complete
connected-chain shadow and a revised chain with independent oracle and real
block-catalogue plan selection/serialization checks. Each covers 12 steps,
48 RHS transactions, 30 layers, 12 topology generations, reloads 2 and 3,
two remaps, an exact final checkpoint and 62 exact VTK members. Each rank checks
1,440 connected layer results.

The connected numerical executor computes pressure, fused mass primitives,
temperature flux, scalar/divu/rotational diffusion, ordered B/Exner injection,
scalar restriction/divergence, velocity physics/source and final gradient.
Mass/velocity tendencies and source, B, Exner, PV and vorticity are compared
bit-for-bit with the accepted path. The oracle run additionally compares the
coarsest-level mass/temperature tendency writes. This does not yet establish
all-level temperature ownership or a new production-science result: the accepted
solver remains authoritative in these overlays.

The revised replay selects operations by actual block owner, preserves original
producer ordinals, compacts physical addresses, verifies complete producer/
restriction/source/injection coverage, and round-trips every selected plan
through a topology-only byte format. It does not execute these split plans.

Evidence is archived in `.stage278-validation/connected-replay`, with source,
binary, per-rank result and file hashes in `connected-validation.json`. Earlier
source and logs are preserved separately. Local wall/CPU time is not cluster
performance evidence.

Subsequent draft work adds actual plan distribution, phase routing, raw producer
delivery before cross-owner restrictions, resident input binding, and identity
copy elimination. Those later source bytes are not qualified by these two runs.
Their distributed replay and resident/RK integration are in progress.

The campaign remains unsealed. No production Domain pathway or transfer has
been retired, and no speedup, legacy-gap reduction, above-Domain scaling or
>=4D credit is awarded. Accepted266 remains unchanged. Remaining work is the
production resident RHS/RK selector, boundary/inverse reader closure, paired
retirement, exact active-candidate qualification, replacement-cost accounting,
and sealed bbserv delivery. No further user guidance is needed.
