# Performance reset after Stage255

The current path is not closing the J5 regression quickly enough. Stage249's
matched direct-main comparison measured 0.540682 s/ordinary step for main
and 1.677614 s for the migrated executable, a 1.136932 s gap. Stage250's
logging calibration left a 1.123003 s main-relative gap. The later Stage251
feature ladder found a 0.386426 s matched c4-to-c5 jump when enabling
`BLOCK_VELOCITY_OWNER`; Stage252 reproduced a 0.405565 s jump and located
0.343379 s/step of the all-step gap in native velocity source.

Stage253's empty-collective change saved about 15 ms/ordinary step. Stage255's
depth index saved 70.408 ms in c5, but only 27.963 ms per profiled all-step
inside native velocity source; the index also speeds owner programs present
in c4. Thus the full 70 ms must not be counted against the velocity-owner
feature gap. The depth-indexed c5 run is a partial feature stage, not the
full migrated executable from the direct-main campaign. There is no new
matched main-relative result after Stage255; the last confirmed full
main-relative gap remains Stage250's 1.123003 s/step.

The owner velocity path currently compiles a per-operation dependency graph.
For every physical layer and adaptive level, `execute_owner_velocity` runs
that graph; `program_execute` reads each dependency depth with
`MPI_Alltoallv`, then reads the rotation delivery. The gradient phase later
performs another source-value exchange. Stage254 measured 0.152 s per
all-step inside owner-program read MPI and 0.111 s in producer loops; MPI
time can include waiting for uneven producer work. Removing one scan or one
empty collective cannot plausibly eliminate a ~0.4 s feature gap by itself.

The next implementation should replace the hot velocity owner graph with a
compact block-resident stencil plan: retain owner-owned operands and results,
compile the existing direct/restriction action order into contiguous local
work, and exchange only values actually needed across owners at causal
boundaries. The pure `execute_velocity_sources` kernel already expresses the
required direct and restriction arithmetic without a Domain or MPI
dependency. Its current caller still gathers from Domain-backed geometry;
the new caller must obtain the same inputs from block-owned primitives,
mass flux and diffusion values. This is a design direction, not an assertion
that the input routing is already solved. First audit ownership and
cross-block dependencies, then implement a narrow locally testable prototype
and compare exact J4 state/checkpoint/VTK outputs. A matched bb13 c5
comparison should follow only after the prototype is correct.

Stage256's prepared depth-indexed c4/c5 campaign measures the residual
feature gap but cannot improve it. Defer that full campaign until the
measurement is needed as a baseline for the larger design. After a material
velocity-source improvement, repeat a direct main comparison: the common
Stage202/203, primitive, diffusion and initial-tendency costs remain even
if velocity ownership is repaired. Keep Domain removal and promotion paused
until the main-relative regression is addressed.
