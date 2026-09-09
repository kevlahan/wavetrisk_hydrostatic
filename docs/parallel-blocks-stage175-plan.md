# Proposed Stage 175: complete native velocity tendency

Placement amendment: the user subsequently approved native geometry-owner
execution first, with final-owner kernel placement deferred. See
`parallel-blocks-stage175.md` for the implementation and current acceptance
status. The original final-owner/shared-producer migration below remains the
longer-term architectural plan, not a claim about this delivery.

Base: validated Stage 174, commit `d75e4634`.
Status: design and acceptance contract. A numerical prototype and independent
source oracle now exist; see `parallel-blocks-stage175-work-in-progress.md`.
The complete production replacement is not implemented. Stage 174 is the
rollback checkpoint.

## Readiness assessment

The project is ready to tackle the complete velocity tendency as a coordinated
stage. It is not ready to disable the Domain producer using the current code.
Stage 174's native residual starts from an already completed Domain velocity
tendency; exact residual agreement does not prove native Qperp, source
restriction, or restricted Bernoulli/Exner production.

The current oracle measurement in `multi_level.f90` records direct Qperp,
physics, edge length, integrated source, and activity/coverage. Its `direct`
flag is true only on `level_end`; the coarse integrated source is captured
from `dvelo`. Thus its direct-source identity is not an independent test of
the multilevel velocity-source restriction implementation.

The existing full Domain tendency oracle, component measurements, native
prognostic storage, thermodynamic workspaces, boundary machinery, and validated
restart test provide the foundation for the replacement. They do not remove
the need to implement the missing producers and their exact phase dependencies.

## Delivery boundary

One substantial delivery should compute the velocity tendency from block
prognostic fields and required primitive inputs, publish it to native RK
accumulation, and complete its boundary dependencies. It must not consume a
completed Domain velocity tendency or a Domain-produced Qperp/gradient/source
as its production answer.

Mass integration and external physics interfaces can remain compatible for
this stage, provided shared producers run once and any remaining handoff is
explicit and measured. Domain topology may be used when compiling plans, and
Domain state can remain for output, checkpoint, and external interfaces. A
mass-only compatibility consumer must not cause velocity kernels to run again.

## Implementation workstreams (internal gates, not separate cluster releases)

1. **Compile and close the velocity input dependencies.** Build persistent
   per-level final-owner stencils/routes for edge mass flux, edge potential
   vorticity, geometry weights, pressure/geopotential, kinetic energy,
   Bernoulli, Exner, theta and physics inputs. Resolve ownership, pentagons,
   scaffold nodes and coarse/fine dependencies explicitly. Rebuild metadata
   only on the existing topology/ownership generation. Reset dynamic validity
   each RK stage. Preserve the ordering of boundary completion and restriction.

2. **Produce the complete velocity source natively.** Implement the original
   Qperp stencil and geometry-weight arithmetic, apply the once-evaluated
   physics source, and reproduce `cpt_or_restr_u_source`: direct evaluation
   and child-edge summation depend on masks, absent children and level order.
   Compare direct values and every restricted component before proceeding.
   Captured Domain integrated sources are oracle data, not production inputs.

3. **Produce and apply the complete gradient natively.** Integrate pressure
   and geopotential in the original physical-layer order; construct kinetic
   energy/Bernoulli and preserve their fine-to-coarse restriction. Carry the
   correctly restricted Exner values into the gradient. Do not substitute the
   Stage 173/174 dynamically reconstructed Exner cache for that field. Apply
   metric scaling, edge theta averaging and masks with the original operation
   order, including the first physical layer and pentagon cases.

4. **Cut the production dependency, including RK closure.** Stop calling
   Domain velocity-source and velocity-gradient routines on non-oracle runs.
   Remove raw compatibility velocity-tendency transport and its residual
   reference dependency from production. Publish native velocity tendency and
   stage values directly; complete boundary/scaffold values needed by the next
   RK stage. `RK_sub_step_compatibility` currently still integrates Domain
   velocity, so replacing only the final tendency kernel is insufficient.
   Preserve any accepted rounding operations needed for numerical continuity
   without retaining a Domain-produced completed velocity field.

Audit shared basic operators carefully: native velocity must not trigger a
second mass-flux/PV/pressure pass while Domain compatibility evaluates the
same fields. Where mass remains compatible, share a single authoritative
producer through a bounded, measured handoff, rather than duplicate evaluation.

## Expected code areas

- New native velocity workspace/kernel module, with topology-specific plans
  kept separate from numeric stage state.
- `parallel_block_mpi.f90`: primitive routing, component oracles, native
  tendency publication, removal of raw velocity-tendency import.
- `multi_level.f90`: independent full oracle and genuinely mass-only production
  compatibility path, without Domain velocity-source/gradient execution.
- `time_integr.f90`: native velocity RK and boundary-stage closure.
- `parallel_block.f90`: native velocity cache/tendency interfaces as required.
- `ops.f90` and `diagnostics.f90`: original reference formulas remain intact;
  they define arithmetic, masks and ordering for independent comparison.

## Required evidence before delivery

- Fresh checked RK3/RK4 and optimized RK4 builds, normal line limit and checks.
- Component checks for Qperp, physics contribution, direct and restricted
  sources, pressure/geopotential/Bernoulli/Exner and final velocity tendency.
- No tolerance relaxation to conceal a phase mismatch. Exact checks where
  arithmetic is unchanged; investigate rather than dismiss reference changes.
- Poison or deny access to completed Domain velocity tendency on the production
  path, proving the native result does not depend on it. Diagnostic reference
  buffers must be separately owned and available only with the oracle.
- Assert zero non-oracle Domain velocity-source/gradient evaluations and zero
  completed Domain velocity-tendency imports, including RK boundary handling.
- Full local optimized restart through adaptation/remap and post-restart
  steps; compare checkpoint/numerical results with the accepted baseline.
- Then full 83-rank RK3/RK4 oracles, optimized timing and profiling. No series
  of cluster deliveries that merely transport already computed Domain answers.

## Performance judgement

Stage 174's production time is 22.418 s on the current cluster test. Its first
eight profiled steps include 3.8775 s rank-average Domain tendency compatibility,
3.8182 s scalar restriction/divergence and 3.6940 s native inverse transform.
These inclusive timers overlap; they are not additive removable budgets.
In particular the 3.8775 s pass includes mass/shared work, not just velocity.

The delivery must show actual removal of velocity computation/traversals and
imports, plus account for replacement kernel, route-setup and compatibility
costs. Report repeated paired unprofiled Stage 174/175 runs on the same CPU
types, rank count, input/checkpoint, compiler flags and allocation. A useful
performance target is at least 10% lower total runtime on this test, but it is
an acceptance target to investigate, not a speedup guarantee. If complete
native velocity is correct but does not meet it, report that result honestly
and use the new profile to select the next shared-work/communication target.
