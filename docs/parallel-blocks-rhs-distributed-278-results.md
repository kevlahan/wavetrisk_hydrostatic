# Distributed connected RHS — Stage278 development

Current status: final active checked/optimized/oracle qualification and kit sealing
are complete. See [the production result and launch commands](parallel-blocks-rhs-production-278-results.md).
The development chronology below preserves the earlier evidence scopes.


The four-rank checked J4 distributed replay passes 48 complete RHS transactions,
30 layers each, all 12 timesteps, reloads 2 and 3, two remaps, an exact final
checkpoint and 62 exact VTK members against accepted266. Each rank records
1,440 connected layer comparisons. The actual block-catalogue ownership splits
the numerical plans; pressure, shared primitives, scalar and velocity diffusion,
restriction, physics, sources and both scalar/velocity tendencies execute there.

The first placement replay exposed a missing raw-divergence delivery before
coarse restriction and Laplacian(divu). The repaired phase delivers that completed
field by physical identity without applying a legacy boundary alias. Original
failed runs, the successful run, exact source and binary are preserved in
`.stage278-validation/distributed-replay`; `distributed-validation.json` hashes
370 archived files. This replay uses validation-only state/scratch input and RHS
return streams, so it earns distributed arithmetic credit, not production
retirement or speedup credit.

A subsequent short checked resident-input replay passes all four substages of
one timestep: resident catalogue state/means plus requested uncatalogued aliases
match accepted post-halo operands bit-for-bit and drive the distributed chain.
The production selector now connects that input path to the real resident RK
accumulator and a bounded inverse/boundary reader adapter. It bypasses the old
Domain RHS, mass/velocity interior RK and full temperature-stage publication.
Active qualification is in progress. Its first oracle replay reached the
wavelet consumer, exposing geometry initialization formerly coupled to the
scalar-tendency capture. A geometry-only adapter is being qualified; no old
numerical capture is restored. No accepted solver or performance claim changes.

The full campaign request remains active. Required completion is active checked,
optimized and independent-oracle science qualification, retirement/cost evidence,
source-bound campaign gate and a verified sealed bbserv kit. No further user
guidance is needed and no Slurm job has been submitted.

## Active production correction and complete checked-oracle result

`checked-active-footprint-v1` completed the trajectory but failed exact checkpoint
comparison (maximum temperature-wavelet difference about1.47e-9). The accepted
velocity closure deliberately adds and subtracts a reconstructed thermodynamic
term; omitting that cancellation changed rounding. The resident retention now
preserves it using block state, and the bounded raw mass/velocity RK reader aliases
retain their accepted values. Complete-node edge lengths are delivered as topology
metadata. No numerical Domain capture was restored and no tolerance was relaxed.

`checked-active-roundtrip-oracle-v2` passes12 steps,48 actual RHS/RK stages,
all30 layers, both reload/remap events, exact checkpoint and62 exact VTK members.
This is active production numerical evidence, unlike the reference-fed distributed
replay above. Final source-bound checked/optimized/optimized-oracle replays include
the subsequently completed replacement-cost ledger and are underway. Cluster
performance and actual above-Domain ownership remain campaign gates.
