# Stage277: cached fused primitive production integration

Stage276 proved sampled primitive arithmetic on block owners. This stage places
that arithmetic in the production pathway for every positive patch interior,
all RHS calls and all layers, preserving the existing consumers and final shared
write order. It is deliberately **not yet an ownership migration**. Execution
stays on Domain owners because physics, Bernoulli and other live consumers still
need those fields there. Returning block-produced intermediates would introduce
a new bulk transfer, which this stage avoids.

The isolated Stage266 overlay changes only `multi_level.f90` and `ops.f90`.
It caches 37 geometry values and 16 operand addresses per interior node at the
writeback topology/ownership generation. The fused array kernel replaces direct
interior mass-advection flux, PV, kinetic energy and vorticity arithmetic. Physics
residuals are added in the original order. The original traversal installs all
four shared vorticity writes, with boundary, pentagon and restriction processing
unchanged. Cached geometry contains no mask-dependent choices; mask predicates
continue to use current metadata. No new MPI route or collective is introduced.

The independent legacy arithmetic is retained for boundary work, the full oracle
and explicit assignment checks. Production timing runs disable assignment checks
and do not execute the displaced interior calculations. Validation checks 19
values per interior call, including final mass flux after physics. Rank logs must
cover all 48 RHS transactions and all 30 layers. Ranks without Domains correctly
report zero primitive work: no improvement to ownership distribution is claimed.

Four serial builds; eight validations at 16/83/160/256 ranks; eight single timing
runs: main/Stage266/candidate at 83 and 160, Stage266/candidate at 256. Startup,
checkpoint and reload steps use the established timing exclusions. No automatic
promotion and no repeated identical timing configurations. If a result regresses
or is marginal, resolve it before stacking the next removal. An accepted kernel
can later be reused in a complete block-owned producer/consumer chain; that larger
migration and >=640-rank qualification remain outstanding. MPI shared memory is
kept separate.

Budget roughly 20–30 minutes after allocation, with uncertainty from full-call
assignment checks and three extra timing runs versus Stage276's measured 16:11.
Queue and collection/transfer are extra. The kit records exact Slurm accounting.

On the local machine:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage277-validation/rk277-primitive-production-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk277-primitive-production-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk277-primitive-production-kit-v1/launch.sh"
```

Status:

```bash
bash "$HOME/rk277-primitive-production-kit-v1/status.sh"
```

After completion, collect and transfer the review to the local Desktop:

```bash
bash "$HOME/rk277-primitive-production-kit-v1/collect.sh"
```

The kit contains its sources and helpers. It requires only the unchanged,
hash-checked established J5 fixture on bbserv. Partition nk contains bb13 and
needs no bb04 exclusion. Launch guards against duplicate submissions; collection
also preserves failed-job evidence. The accepted Stage266 source and separate
Stage273 candidate remain untouched.
