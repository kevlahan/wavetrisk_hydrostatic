# Stage242 — retained native physics state

All sixteen qualification cases and `physics_242/review_physics.py` passed.
Evidence is recorded in `parallel-blocks-stage242-local-results.json`.

`PHYSICS_NATIVE_STATE=on/check` requires native physics inputs. Means and soil
references for aliases are captured at RK entry, where those fields are unchanged
through the RK stages. Atmospheric aliases retain their distinct final-inverse
consumer epoch. Production physics gathering now reads no Domain prognostic,
mean or soil values. Geometry and fixed reference inputs remain geometry metadata.

Completed temperature, velocity and soil values are retained on final block
owners; primary aliases retain their result in the native alias store. The 27th
column-header field carries the pole flag. Pole velocity, mass, means and the
pre-physics atmospheric wavelet decision image remain unchanged. Installation
happens after the complete DAG executes, preserving version-zero input values.
The pending grid-change phase and writeback accounting are retained explicitly.

The check build reads stored values back and compares them bit for bit, including
unchanged mass/mean and pole-velocity invariants. The original serial physics
replay independently verifies complete package packets and returned state.

Passed: oracle-4 with 35,812 storage readbacks; full 4/41 trajectories with
12 exact states, two remaps and 62 exact VTK members each; control; DAG route,
dependency, version and output rejection; destructive publication poison. Each
full run retains 460,368 canonical and 240 primary alias columns. Package-soil poison and stronger entry poison of all old mean/soil inputs also
matched exactly. Native input faults and corruptions of retained canonical/alias
storage were rejected at their intended guards. The alias corruption marker was
required to prove that injection was reached before the shared mismatch guard.

The entry-poison test restores test buffers after gathering so remaining legacy
consumers can run. The separate publication-poison test uses no restoration.
Native physical state is retained, but the current downstream soil transform and
adaptation still use compatibility publication and rebuild native storage afterward.
They, the initial alias import, restart/output and obsolete persistent mirrors
remain migration work. No complete-lifecycle or performance claim is made here.

Local solver cases run sequentially on the four-core host with library threads
one. Forty-one ranks is oversubscribed correctness coverage. No bbserv job runs.
