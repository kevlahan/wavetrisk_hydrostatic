# Stage275: minimal inputs and resident-state epoch diagnostic

Starts from accepted Stage266 and changes only a disposable `multi_level.f90`.
Stage273 remains excluded. This is the input gate for the compact block execution
prototype, not a production ownership change or performance comparison.

At every live source and gradient call, a shadow execution uses only the required
input components, with unused input locations filled with quiet NaNs. Results
must match the unchanged production kernel bit for bit. Source uses 19 of the
27 declared components in each flux/PV stencil. Ordered restriction is retained.
The reviewer rebuilds per-field budgets from the sampled actual plans.

At gradient call ordinals 1, 2, 3, 9, 17, 33, 41 and 48, first and last atmospheric
layers, the probe records actual resident block interior scalars and the matching
mapped Domain gradient operands. Offline joining uses block ID, mesh level and
exact coordinate bits within the same call/layer. State, mean and their sum are
compared separately. Mismatches are reported as a failed reuse hypothesis; they
are not hidden or turned into science failures. Unmapped aliases, halo freshness,
intervening layers/calls and upstream derived fields remain unproven.

Four serial builds and eight full J5 validations run on exclusive nk/bb13:
16-rank main/Stage266/candidate/checked/oracle, then 83/160 candidate and 256 oracle.
Allow **20–30 minutes after allocation**. Stage274's combined runner time was
14:10; the new shadow computation and field export add uncertain cost. Queue and
collection/transfer are additional. No automatic repeats or promotion.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage275-validation/rk275-input-contract-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk275-input-contract-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk275-input-contract-kit-v1/launch.sh"
```

Check progress:

```bash
bash "$HOME/rk275-input-contract-kit-v1/status.sh"
```

After completion, collect and transfer the review to the local Desktop:

```bash
bash "$HOME/rk275-input-contract-kit-v1/collect.sh"
```

Only the established hash-checked J5 fixture is required on bbserv; source and
helpers are self-contained. The review includes accounting and failed-job logs.
