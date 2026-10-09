# Stage276: fused block primitive prototype

A bounded executable prototype on sealed accepted Stage266. The source overlay
adds the prototype to multi_level.f90 and an optional assignment-point reference
output in ops.f90. Stage273 partitioning is excluded.

One fused loop computes positive-interior mass flux, PV, kinetic energy and
centre-triangle vorticity at its initial assignment on actual block owners. Resident mass (state + mean)
and velocity are gathered through addresses compiled at the velocity-plan
generation. Direct neighbour records resolve patch, boundary and ghost coordinates without
depending on downstream restriction-field metadata.

The original Domain step1 remains the reference. At the original positive-interior assignment, before later shared vorticity
writes and restriction/pentagon postprocessing, records contain exact geometry, mass-physics residuals, input
operands and reference outputs. A bulk validation exchange sends these to block
owners. Reference-input execution must match step1 bitwise. Resident-input
execution separately reports differences for patch, boundary and ghost storage.
A failed residency hypothesis is retained as evidence, never silently repaired
by substituting the reference input. Output stays in a private prototype buffer.

Scope: positive patch interiors, RHS ordinals 1,2,3,9,17,33,41,48 and atmospheric
layers 1/30. Domain scaffold/boundary production, physics, B/Exner, restriction,
output halo publication and downstream production consumers remain unchanged.
The supplied reference geometry/physics and validation packets are not a proposed
production communication path. This stage cannot qualify a production speedup.
No displaced work is claimed as removed. A subsequent complete replacement must
connect consumers, retire displaced Domain work, and pass a full-step gate.

Four serial builds; eight validations at 16/83/160/256 ranks (including checked
and oracle); five separate uninstrumented timing controls: main/Stage266 at 83,
main/Stage266 at 160, Stage266 at 256. No repeated identical timing configuration.
The 256-rank prototype must report work on ranks with zero owned Domains.
J5 has 160 Domains, so the >=640-rank strategic target remains a later multinode
qualification, not a claim from this intermediate test.

Allow approximately 18–25 minutes after allocation, with uncertainty for the new
prototype. Stage275's four builds/eight validations took 14:35; this adds five
controls and bounded prototype overhead. Queue and collection time are extra.
The review preserves Slurm accounting, source seals and failed-job diagnostics.

Local transfer:

```bash
scp /Users/kevlahan/.codex/worktrees/26d8/wavetrisk_hydrostatic/.stage276-validation/rk276-primitive-prototype-kit-v1.tar.gz kevlahan@bbserv.math.mcmaster.ca:~/
```

On bbserv:

```bash
tar -xzf "$HOME/rk276-primitive-prototype-kit-v1.tar.gz" -C "$HOME"
bash "$HOME/rk276-primitive-prototype-kit-v1/launch.sh"
```

Status:

```bash
bash "$HOME/rk276-primitive-prototype-kit-v1/status.sh"
```

After completion, collect and transfer to the local Desktop:

```bash
bash "$HOME/rk276-primitive-prototype-kit-v1/collect.sh"
```

Only the established hash-checked J5 fixture is required on bbserv. Source and
helpers are self-contained. Partition nk contains bb13; no exclusion is needed.
