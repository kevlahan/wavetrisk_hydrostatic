# Stage241 — native physics input ownership

All fourteen qualification cases and `physics_241/review_physics.py` passed.
The source-bound evidence is in `parallel-blocks-stage241-local-results.json`.

`PHYSICS_NATIVE_INPUT=on/check` requires the batch physics path and retained RK.
The production queue skips the old atmospheric packing routine. Each column and
its six velocity operands carry native storage keys. Execution owners fetch
canonical atmospheric values, means and soil directly from final-owner Block_Data.
Zero-key atmospheric aliases come from the saved final-inverse consumer epoch;
readiness, registration, storage class and owner checks reject invalid access.
The tendency epoch is never substituted for that saved consumer version.

The original serial packing remains in the independent check build. Every owner
packet and double-precision returned state is compared exactly. Native alias
inputs are also compared against the reference source values before execution.

Passed: four-rank replay (35,812 columns), full 4/41 trajectories (12
states, two remaps and all 62 VTK members each), serial control, inherited DAG
faults, destructive publication poison and two-step package-soil poison. Each full
trajectory used 3,122,100 canonical native requests and 102,156 alias requests.
Native storage-owner, epoch and input-corruption checks were rejected at their
intended guards. Entry poison passed and matched the reference timestep exactly.
It exercised 20 primary alias columns across the four source ranks.

The entry-poison test poisons all old Domain atmospheric solution values and
canonical mean/soil inputs during gathering. It deliberately restores those test
buffers after gathering because the downstream compatibility lifecycle still
uses them. The separate publication-poison test has no backup or restoration.
These are distinct claims; neither proves the whole lifecycle has migrated.

Zero-key means and soil still use a bounded source adapter. Physics result
publication, soil transforms, adaptation reconstruction, initial alias imports,
restart and output still require further migration. Defaults remain unchanged.
No bbserv job or performance claim is part of this checkpoint.
