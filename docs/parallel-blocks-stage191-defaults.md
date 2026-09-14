# Accepted Stage 191 defaults

Enable installation transfer/proof gating by default. Disable mass overlap by
default; it remains available only with `WAVETRISK_STAGE191_MASS_OVERLAP=1`.
`WAVETRISK_STAGE191_INSTALL=0` remains the installation control. Explicit runtime
settings still undergo the same rank-consistency checks, with requested values
printed once.

The exact original cluster-tested source is committed at
`37ac536405afbb8b2f8d59037db126950af574cd`; it defaulted both selectors on. Job
54561's three installation-only runs explicitly selected `1 0`, passed all
numerical/workload gates and each recorded a 1.34 s ordinary-step median.
Mass overlap alone provided no repeatable gain, so the accepted default selects
that installation-only configuration. The historical kit and results retain
their original source hashes and explicit selector values.

The only production-source delta from that tested snapshot is the initial
selector value and its explanatory comment. No MPI schedule, numerical kernel,
buffer lifetime or ownership-transfer implementation changes. Validation of
the default selection includes actual two-rank selector tests for unset values
and explicit overrides, alongside the ownership-transfer and MPI buffer tests
(three tests pass). A fresh O2 contraction-off J4 build completes 12 local steps,
checkpoint 2 -> 3/reload and two remaps, with exact equality to the tested
candidate and saved independent RK4 oracle. Local timing is excluded.

Details are in `parallel-blocks-stage191-default-validation.json`. The cluster
kit continues to set selectors explicitly for the historical combined candidate
and ablations; its results are not relabeled as default-off runs.
