# Stage325 provider-checker audit and offline preflight

2026-10-26. The checker defect from job 55572 is repaired and the offline
provider preflight passes. Full native qualification remains incomplete. This
work changes the verifier and its tests, not the solver, owner map, numerical
work or saved captures. No new target-hardware job was submitted in this step.

## Repair and source contract

The failed checker unconditionally read `local_order297` and `local_start297`.
All three source-bound environments set `WAVETRISK_SHARED_GHOST_SOURCE=0`.
`build_block_ghost_exchange_plan` therefore leaves those arrays unallocated,
and `audit_ghost316` emits shape `[-1]` without a payload. Scalar execution
visits request indices in ascending order and skips remote source owners;
vector execution also uses that request order.

The corrected gate requires the disabled-option allocation state explicitly,
then derives the exact local request sequence from `source_owner`. It rejects
missing or inconsistent shape descriptors and unexpected grouped payloads.
It does not replace missing arrays with an empty list or omit local writes.
The same check now runs during each rank's normal capture validation, so the
selected-capture preflight reaches the contract that previously failed late.

The helper has an explicitly selected grouped-mode check for source-derived
regression tests: complete permutation, group boundaries, common source per
group, increasing global block/patch groups, and stable request order within
one group. The sealed Stage325 verification path always selects disabled mode;
these tests do not qualify a new grouped-mode solver configuration.

Inverse transfer contracts were also reviewed against their actual producer
and consumer. The recorded sparse-outer option is disabled, so addresses have
six components and `edgewise` is false. Optional source/target metadata are
paired and have the exact producer/request extents. These conditions are now
checked before provider reconciliation. An unsupported edgewise layout is
rejected explicitly rather than being accepted or inferred from array length.

## Coverage of the remaining provider paths

A complete saved rank stream contains 254 provider families in 26 scope names.
The audit maps all of those scopes to the following reconciliation paths:

| Path | Recorded families per arm | Contract exercised |
|---|---:|---|
| Snapshot routes | 140 | Actual wire source slots, signs and ordered repeated destinations; numerical, fixed, mixed-prefix and publication roles |
| Transform adapters | 8 | Reversed count convention, packed receive index, fixed Domain source identity and signed consumer |
| Ghost plans | 5 | Local request order, remote manifests and numerical-owner permutation |
| Writeback and boundary import | 5 | All 640 catalogue identities, retained blocks, reverse Domain sources and exact replacement extents |
| Inverse transfers | 64 | Each of the five scope names, address-to-catalogue binding, optional metadata and empty plans |
| Inverse aliases | 32 | Fixed Domain identity, default install order and explicit install permutations |

The source audit used the exact job55569-bound versions of
`parallel_block_mpi.f90`, `parallel_block_inverse.f90`,
`parallel_block_rhs_exchange.f90`, `parallel_block_transform_adapter.f90`, and
`parallel_block_transform_backend.f90`. Their hashes match the original native
source manifest. In particular, route sends snapshot before ordered writes,
adapter source indices resolve through fixed Domain IDs, inverse aliases retain
explicit install order, and boundary import runs in the reverse direction to
bulk block writeback.

Nineteen new tests execute the real provider function, SQLite cache and
reconcilers on constructed 256-rank worlds. No reconciler is mocked in those
tests. The cases exercise all dispatch branches, including empty participants,
actual moved labels, repeated destinations and retained/remote boundary work.
Mutations cover missing ranks, wrong wire data, unowned publication providers,
invalid Domain indices, wrong storage addresses, inconsistent metadata,
unsupported options, bad group bounds and corrupt install permutations.

The previous reporting tests still use mocks to test report persistence; they
are not counted as provider-path evidence. Instrumented execution of the full
suite reached 74 of 75 provider-function statement lines. The only unexecuted
line prints periodic progress. All statement lines in the two new contract
helpers were exercised. This is line-execution evidence, not proof of every
boolean combination or complete solver correctness.

## Results and reproducibility

- All 128 offline tests pass, including the 19 new provider tests.
- The unchanged original provider implementation reproduces the exact missing-
  field exception on a complete constructed world. The repaired implementation
  passes on the same input database, resolving both local and remote consumers.
- A 961-byte fixture extracted from job55569's rank0 capture preserves the exact
  unallocated ghost-plan state. Its capture and fixture hashes are recorded.
- Nine selected actual saved captures pass the normal per-rank cache path,
  including the new ghost/inverse contracts. All 5,450 available cross-arm
  feature comparisons remain exact.
- The sealed failed v2 package is still reproducible byte-for-byte at the payload
  level. Its verifier is frozen separately so this repair cannot alter it. The
  normalized reproduction archive retains SHA256
  `52ea34e3abd82eb9942ee8ce589904ece03061ac6eed8348b320b185be67f489`.

The repaired verifier SHA256 is
`7adce8a70ebf6b923f06414a3650dac5344337d5511bc2bd48c4f0dc8411314d`.
`test/parallel_block_profile/locality_325/PROVIDER_PREFLIGHT.json` records the
source bindings, test results, source-derived scope census, selected-capture
results and explicit limitations. Reproduce the unit tests from the worktree:

```bash
python3 -m unittest discover -s test/parallel_block_profile/locality_325 -p 'test_*.py'
```

## Decision

The offline prerequisite for another full verification replay is satisfied.
This is stronger preflight evidence than the earlier selected-stream checks:
the provider-consumption branches themselves now execute in tests, and the
failed option state is checked while parsing real captures. It still cannot
establish the untested full-world identities in the original native data.

The remaining gate is a complete verification-only replay of those existing
captures, followed by independent result review. Solver reruns and performance
runs are unnecessary for that gate. A new replay must use a separately sealed
package and remain inside the original budget: 881 seconds have been consumed,
leaving 1,819 seconds. Reusing the former 35-minute allocation would violate
that remaining budget. No claim of speedup or primary slowdown attribution,
accepted-solver promotion, or validation-conditional commit follows from this
offline result alone.
