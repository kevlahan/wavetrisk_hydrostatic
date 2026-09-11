# Stage 182: numerical diagnosis before shared horizontal geometry

Stage 181 was committed as `78ba3574`. The next step is to resolve the
independent legacy comparison before changing computational storage. No
production numerical source, authoritative `main` source, or acceptance
threshold has been changed by this diagnostic work.

## Experiment and instrumentation

The evolved J6 checkpoint-4 fixture and four-rank RK4 interval through
checkpoints 5 and 6 are unchanged from Stage 180. All three block oracles are
disabled. Legacy and block runs execute sequentially, not concurrently.
These are diagnostic runs, not timing measurements.

`prepare_restart_probe.py` prepares external source copies. Its optional
physics probe records the **existing** packed single-precision arrays before
and after the physics solver, together with the double-precision raw column
inputs. It does not repack, refresh aliases, or invoke numerical producers.
Phase probes observe existing Domain interiors after dynamics, physics,
adaptation, and at step end. A Domain filter limits phase output; an independent
option can record packed physics on all Domains.

The first, read-only paired experiment is in
`/private/tmp/wavetrisk-stage182-probe-runs-01`. Both instrumented executables
produce **exactly the same compared checkpoint fields**, metadata, thresholds,
and topology as their respective uninstrumented references at checkpoints 5
and 6. This is instrumentation transparency, not legacy/block equivalence.
The complete phase comparison is `trace-results.json` in that directory.

## First amplification event

Through step 259, the packed physics inputs on Domain 77 are identical.
At step 260, global Domain 77, node address 123, atmospheric layer 26,
the input U velocity straddles a float32 rounding midpoint:

| Quantity | Value |
| --- | ---: |
| Legacy raw double U | 0.06222941167660471 |
| Block raw double U | 0.062229411676523724 |
| Absolute double difference | 8.098383075250126e-14 |
| Float32 midpoint | 0.06222941167652607 |
| Legacy packed U | 0.06222941353917122 |
| Block packed U | 0.06222940981388092 |
| One float32 ULP at this magnitude | 3.725290298461914e-9 |

The packed input bit patterns are 1031726149 and 1031726148. The physics
outputs are two ULPs apart (7.450580596923828e-9). The phase probes locate
the jump **inside physics**, before the subsequent remap. Temperature and
mass then diverge during later dynamics steps. The worst original checkpoint-6
differences are on this same Domain and atmospheric layer.

This identifies amplification by mixed-precision conversion. It does not, by
itself, prove that every remaining numerical difference is acceptable.

## Controlled causal experiment — not a production correction

An explicitly marked isolated legacy executable replaces only this one packed
U input at step 260 with the block bit pattern. No other input, timestep,
checkpoint, or numerical operation is changed. The replacement is logged once.
The corresponding preparation option is:

```sh
--physics --physics-replay-u 260 77 123 26 1031726148
```

The experiment is in `/private/tmp/wavetrisk-stage182-replay-runs-01`.
After the substitution, all compared packed physics inputs on Domain 77 match
the block inputs through step 266. Global checkpoint-6 differences shrink:

| Atmospheric wavelets | Original legacy vs block | One-input replay vs block |
| --- | ---: | ---: |
| Velocity | 5.997233210042197e-9 | 1.657954606937295e-10 |
| Mass | 5.704197870992633e-9 | 2.874203882912547e-11 |
| Mass-weighted temperature | 2.787588158525978e-6 | 1.887845790804477e-8 |

The remaining worst differences are on Domain 8, atmospheric layer 28.
Velocity still exceeds the unchanged Stage 180 screen of 1e-10. Therefore this
experiment explains the dominant error but is **not a passing global numerical
gate**. A second trace examines the residual before any geometry optimization.

## Second amplification event

The all-Domain packed-physics trace is in
`/private/tmp/wavetrisk-stage182-global-probe-runs-01`. Its legacy executable
includes the first causal substitution; this is explicitly not an unmodified
legacy reference. Phase snapshots focus on Domain 8.

At step 262, Domain 8, node 123, layer 28, W crosses another float32 midpoint:

| Quantity | Value |
| --- | ---: |
| Legacy raw double W | 0.001824796025367625 |
| Block raw double W | 0.0018247960250686623 |
| Absolute double difference | 2.989626775307119e-13 |
| Float32 midpoint | 0.0018247960251756012 |
| Legacy packed W | 0.0018247960833832622 |
| Block packed W | 0.0018247959669679403 |

The additional isolated intervention is
`--physics-replay-velocity W 262 8 123 28 988753406`.
The two-event causal run in `/private/tmp/wavetrisk-stage182-two-replay-runs-01`
completed both checkpoint/reload cycles and the remaps. Each declared
substitution occurred exactly once. Both global checkpoint comparisons pass
the original, unchanged screen. Checkpoint-6 atmospheric wavelet maxima are
velocity **1.100630697692395e-11**, mass **5.552364124028486e-12**, and
mass-weighted temperature **5.165159189246538e-9**. Topology, thresholds and
soil remain exact; no nonfinite values occur.

Thus the two float32 midpoint crossings causally account for the errors above
the screening bounds on this fixture. This is not a universal error bound or
a claim of bitwise legacy equivalence. The original unmodified comparison
still fails its fixed screen and remains reported as such.

`summarize_rounding_diagnosis.py` records original, one-event and two-event
comparisons separately, verifies both observed midpoint crossings, checks
the two intervention log entries, and verifies instrumentation transparency
against the correct reference for each executable. Its complete evidence is
`rounding-results.json` in the two-event run directory. The broader block
trace also matches the original block checkpoint outputs exactly, and the
broader one-event legacy trace matches the earlier one-event legacy outputs
exactly. The authoritative legacy archive still verifies against `main`.

The all-Domain trace also includes inactive/scaffold column values. Their
raw packed differences must not be classified as authoritative solution
errors without examining masks and downstream consumption. This is why
global checkpoint comparison, phase location, and controlled replay are all
required; a count of differing packed numbers is not a correctness verdict.

The user authorized proceeding to horizontal-geometry sharing **once the
rounding diagnosis is complete**, while preserving the block baseline and
explicitly retaining the original legacy-screen failure in reporting. This
does not authorize adding replay to production or claiming the original
screen has passed.

## Safety and acceptance

- Replay is an intervention to test causation, not a way to manufacture a
  passing production comparison. Never transfer a replay executable as a
  production solver.
- Do not claim the unmodified legacy comparison passes because a perturbed
  diagnostic comparison improves.
- Keep the original comparison and fixed screening policy visible.
- Shared horizontal geometry remains a separate numerical-equivalence and
  capacity-reduction change; no speedup is inferred from this diagnostic run.

## Stage 188 arithmetic control

Compiling both unmodified legacy and current block solvers with floating-point
contraction disabled passes the unchanged numerical screen on this original
J5 fixture at both checkpoints, without replay. It does not give bitwise legacy
equality or change the original default-build failure above. The smaller J4
fixture passes the same screen under this control. See the
[Stage 188 diagnosis](parallel-blocks-arithmetic-diagnosis-188.md) for measured
residuals, exact checked/optimized comparisons and the opt-in build setting.
