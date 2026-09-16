# Stage202: scalar adapters and inverse staging prototype

This prototype addresses Stage201's evidence of excessive representation
maintenance. Cluster job **54693 passes** the predeclared target of at least
**0.10 seconds per ordinary time step in both reversed bare pairs**: measured
recoveries are 0.127 and 0.133 s/step, with all numerical, placement and observer
quality gates passing. The production default remains `STAGE202=off`; use
`STAGE202=on` in a fresh build to select the validated combined implementation.
See [the cluster review](parallel-blocks-cluster-202-54693.md) for the evidence,
remaining 0.3835 s/step legacy gap and limits of the timing comparison.

## Changes and exact scope

* Mixed scalar range reads, writes and fills resolve sample/node/geometry class
  once per record segment, instead of calling the old single-slot adapter for
  every value. This affects capture, restriction, transport and execution callers
  of these shared APIs. The 17 field slots and 33 shared geometry slots retain
  their existing layout and canonical geometry-write ownership.
* Production direct scalar flux no longer reconstructs a 50-slot record: the mass
  shortcut reads one slot; physical non-mass flux reads two edge triplets;
  nonphysical non-mass levels read none. Validation still reads the independently
  captured oracle record. The mass shortcut remains before level rejection.
* Local-owner inverse gather/install routes transfer directly between block fields
  and the inverse workspace. They bypass pack-to-wire, self-buffer copy and
  unpack-from-wire, replacing three copies with one. Self route metadata is
  prepared once per topology generation. Remote routes retain the original wire
  layout, snapshot epoch and MPI completion ordering; there is no message-volume
  reduction. The same mechanism covers scalar/vector and solution/wavelet families,
  full and level-restricted routes, and selected-edge installs.
* Inverse alias operations resolve addresses once per route rather than per level
  and field. Remote snapshots precede local copies; local p/k/v ordering, signed
  destinations, receive-install order and consumer alias epochs remain intact.

The inverse workspace, scaffold state and consumer alias snapshot remain necessary
in this prototype. It does **not** remove an entire second inverse computation or
claim that every compatibility boundary update is redundant. Full scalar records
still exist at interfaces requiring them, but reconstruction uses bulk addressing.
This is a bounded test of removing repeated representation work across the two
measured source families, not a completed storage redesign.

`STAGE202=records`, `inverse`, `on`, `off` compile the respective variants. Always
use fresh `BUILD_DIR`/`BIN_DIR` when changing these values. The supplied cluster
baseline is reconstructed from the exact Stage201 source patch, independently of
the new switches. The candidate and checked/observer candidates use `on`.

## Expected mechanism and decision

Scalar adapter instruction counts should fall in remaining dynamics and adaptation.
Inverse staging/alias instructions and copying should fall in wavelet dynamics.
A reduction in subsequent MPI progress is an indirect effect, not an independent
saving to add to the producer's improvement. Local J4 results validate behavior;
bb13 J5 results determine performance.

Six fresh builds: Stage201 baseline, Stage202 candidate, checked candidate, matched
baseline/candidate coarse observers, and the same legacy Domain source with common
checkpoint-I/O hardening. Identical 23 external physics objects and contraction off
are required. Observers use the same Stage201 whole-phase clock/PMU module and C
helper, with only coarse hooks: no fine operation hooks or stack sampling.

Seventeen runs, all fixed 83 ranks split 42/41 across sockets on exclusive bb13:

1. Checked candidate, oracle off and on.
2. Baseline/candidate/legacy warmups.
3. Baseline/candidate/legacy bare timings, then the reverse order.
4. Baseline/candidate coarse controls.
5. Baseline/candidate core-PMU pair, then reverse order.

Retain exact checkpoint comparisons across baseline/candidate/checked/observers,
the existing tolerance-based legacy screen, printed adaptive-workload equality,
checkpoint reload and post-restart remap, CPU placement and memory evidence.
Raw counters must cover every rank/step and run for at least 95% of enabled time.
The core PMU group includes user-space MPI progress. Instructions/CPI are evidence
of mechanism, not a causal conversion to removable seconds. Memory-event and
sampled-stack passes are omitted because this prototype tests the work-removal
hypothesis established in Stage201.

`prototype-summary.json` reports seconds recovered, remaining legacy gap, both
pair differences, per-phase instruction/cycle changes, and observer calibration.
An insufficient timing gain produces `do_not_promote_insufficient_gain` even if
all numerical tests pass. A substantial gain with calibration issues requires
observer review. Two repeats are a prototype screen, not confidence intervals.

## Local checks

The production scalar-storage module is tested against the retained path for
mixed ranges crossing record boundaries, canonical/noncanonical geometry writes,
inactive/physical classes, empty ranges, and signed zero. Tests execute the actual
alias routine against the retained implementation for snapshot ordering, chained
local aliases, signed destinations, both families and receive-install permutations.
The routed-copy tests exercise the actual block module at O0 and O2, including
local-owner bypass, untouched wire slots, selected edges and empty plans.

```bash
python3 -m unittest discover -s test/parallel_block_profile -p 'test_*202.py' -v
sh test/parallel_block_inverse/run.sh
```

`local_batching_202.py` freezes the previous algorithm selectors for local J4/J6
runs. `validate_local_202.py` compares candidate, runtime-check, optimized-oracle
and observer runs with the prior Stage201 checkpoint, checks workload/reload/remap,
checks the fixed legacy screen, and validates the observer's disjoint clock ledger.
The combined runtime-check-plus-oracle run is retained in the cluster protocol.
See `parallel-blocks-stage202-local-results.json` for completed local evidence and
any deliberately excluded diagnostic runs.
