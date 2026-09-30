# Stage243: native soil transform and inverse

SOIL_NATIVE=on/check retains soil wavelets and reconstructed values in native
storage. Selected-range transfer plans preserve atmospheric decision wavelets.
Native inverse preserves the full legacy geometric alias traversal and arithmetic.

All ten qualification cases pass: exact oracle4, production4/41, control4,
epoch/range/alias rejection, full Domain poison, owner readback4/41. Full cases
match12 timestep states, two remaps and62 VTK members; checkpoint numerical
checks pass. The poison harness saves every allocation independently and prohibits
implicit resizing on restoration. The failed earlier harness is preserved.

Evidence is source/binary bound in parallel-blocks-stage243-results.json;
15 kernels pass the normalized arithmetic audit. Frozen source and binaries are
in .stage243b-validation/candidate. Builds use -j1 and solver runs are sequential.
41 ranks provide oversubscribed correctness coverage on the four-core machine.

Adaptation still consumes compatibility publication and rebuilds blocks from
Domain. Native adaptation and bounded lifecycle adapters remain required before
performance qualification. Main is unchanged; no bbserv job has run.
