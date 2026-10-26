# IEEE118 EMT `perf` findings

## Scope and method

These profiles compare the conventional composite EMT network with the
`allssn` EMT network in `IEEE118_Ph3`. In this example, `allssn` switches the
PiLines, impedance loads, and transformers to SSN implementations; the
generators remain conventional in both profiles. Both runs used a 1.0 s simulated interval,
a 0.5 ms timestep, and disabled data logging (`nolog`). The executable was
sampled at 997 Hz on user-space CPU cycles with DWARF call chains. The captured
profiles had no lost samples: 9,582 conventional samples and 750 `allssn`
samples. The reported percentages are shares of sampled user-space cycles over
the whole process, including power-flow/setup work; they are not percentages of
simulation-step time alone.

For profiling, `kernel.perf_event_paranoid` was temporarily changed from 2 to 1
and restored to 2 when sampling finished. Raw profiles were written outside the
repository to `/tmp/dpsim-perf/emt-conventional-sample.data` and
`/tmp/dpsim-perf/emt-ssn-sample.data`.

## Conventional EMT hot spots

| Function / work | Sampled cycles |
| --- | ---: |
| `CompositePowerComp::mnaCompApplyRightSideVectorStamp` | 13.77% |
| `MnaSolverDirect::addRightVectorStamps` | 7.09% |
| `memset` (mostly called under composite RHS aggregation) | 5.10% |
| `CompositePowerComp::mnaCompPostStep` | 4.94% |
| `CompositePowerComp::mnaCompPreStep` | 3.46% |
| EMT resistor current update | 2.99% |
| EMT resistor voltage update | 2.59% |
| `Math::invertMatrix` | 2.15% |
| `klu_refactor` | 1.90% |

The call chain places the 13.77% RHS aggregation cost under
`Simulation::step` → `SequentialScheduler::step` →
`CompositePowerComp::mnaCompPreStep`. The accompanying `memset` samples also
come mainly from this path. This is direct evidence that composite components
clear and aggregate dense, network-sized RHS vectors every step. The solver's
separate `addRightVectorStamps` pass costs another 7.09% and remains a distinct
optimization opportunity.

Composite pre/post dispatch accounts for a further 8.40% combined. Individual
R/L/C state updates and matrix inversion are smaller hotspots. The profile also
contains formatting and logger work from setup, so those samples should not be
read as recurring timestep cost.

## `allssn` EMT profile

The largest SSN symbols were `klu_refactor` (10.00%), `SSN::PiLine::mnaCompPostStep`
(4.98%), `MnaSolverDirect::addRightVectorStamps` (4.14%), and
`SSN::PiLine::mnaCompPreStep` (2.02%). There is no
`CompositePowerComp::mnaCompApplyRightSideVectorStamp` hotspot in this run.
This supports the conclusion that dense RHS aggregation in conventional
composites is one important source of the gap; it does not explain all of the
SSN speedup, since solver and component work also differ.

## Composite optimizations and comparison

`CompositePowerComp` now caches the matrix-vector indices touched by its parent
and RHS-contributing subcomponents during MNA initialization. Each step it
clears and combines only those local entries rather than clearing and adding
full network-sized vectors. It also exposes a compact RHS contribution to the
solver, so the global accumulation adds only the composite's touched rows.
Nested composites contribute their cached internal rows as well. Child task
lists are traversed by reference, avoiding shared-pointer copies in the loops.
These changes keep the same component stamps, topology, and timestep.

The updated `IEEE118_Ph3` target built successfully. For a correctness check,
the conventional EMT case with compact accumulation was compared against the
same case before compact accumulation. The CSVs contain 1,001 time points and
28 matching columns; every logged value is identical at the recorded precision
(maximum absolute difference: 0).

The established per-task profile was rerun for both the optimized conventional
case and its all-SSN counterpart. The historical conventional profile completed
in 10.348 s; the optimized conventional profile completed in 7.300 s (29.5%
lower). The matching SSN profile completed in 3.877 s. Average wall time per
profiled step was 5.172 ms, 3.648 ms, and 1.937 ms respectively. All runs used
1 s of EMT simulation at 0.5 ms, 2,000 steps, KLU, task/step instrumentation,
and no data logger. The generators remained conventional in all cases.

For conventional lines, summed pre-step time fell from 1.336 ms to 0.582 ms per
step; post-step time fell from 1.103 ms to 1.013 ms. The solver task fell from
0.511 ms to 0.215 ms. Post-step remains comparatively costly because composite
lines still dispatch the child component state updates. The full comparison,
including loads, transformers, generators, and the estimate after pre/post, is
in [IEEE118_EMT_Composite_Optimization_Comparison.ipynb](IEEE118_EMT_Composite_Optimization_Comparison.ipynb).

These are single profiled runs. Task timing instrumentation adds overhead, and
the “after pre/post” value includes solver and scheduler overhead; it is not a
pure integration-kernel timer.

## Compact RHS follow-up for the all-SSN network

The all-SSN case already accumulates SSN PiLine and impedance-load RHS values
only at their touched terminal rows. The remaining full-vector additions came
from its conventional VBR generators and SSN transformers. I tried giving
those two component types compact row-wise accumulation as well. Five 1 s
profile runs averaged 3.790 s (1.894 ms/step), compared with 3.767 s
(1.882 ms/step) for five matched runs with that change disabled. The difference
was within run-to-run variation, so the extra hooks were removed. Clearing the
single global RHS vector remains, but at this case size it is not a demonstrated
hotspot. The next optimization should target measured KLU refactorization or
SSN component post-step work instead of adding more RHS compaction.

## Standalone EMT PiLine without child components

The separate `EMT::Ph3::PiLine_NoSub` implementation directly stamps the same
three-phase Pi R/L/C/G model and updates its state without constructing R/L/C
children or a virtual node. In a matching 1 s, 0.5 ms IEEE118 EMT profile with
conventional loads, shunts, transformers, and generators, it completed in
4.810 s (2.403 ms/step). The optimized composite run completed in 7.300 s
(3.648 ms/step), and the all-SSN network in 3.877 s (1.937 ms/step).

Line pre-step time was 0.272 ms/step for `PiLine_NoSub` versus 0.582 ms/step for
the optimized composite; line post-step time was 0.301 ms/step versus 1.013
ms/step. This is consistent with removing child-task dispatch and repeated
child voltage/current updates. The profile is a single run; electrical output
comparison is still needed for validation. The plots and per-task breakdown are
in [IEEE118_EMT_Composite_Optimization_Comparison.ipynb](IEEE118_EMT_Composite_Optimization_Comparison.ipynb).

## Reproduction

From the build directory, run:

```sh
perf record -F 997 -e cycles:u --call-graph dwarf,8192 \
  -o /tmp/emt-conventional.data -- ./dpsim/examples/cxx/IEEE118_Ph3 emt short nolog

perf record -F 997 -e cycles:u --call-graph dwarf,8192 \
  -o /tmp/emt-ssn.data -- ./dpsim/examples/cxx/IEEE118_Ph3 emt allssn short nolog

perf report --stdio --no-children --call-graph none --sort=symbol \
  -i /tmp/emt-conventional.data
```

The local system may require temporary authorization to lower
`kernel.perf_event_paranoid` for sampling. Restore its original value after the
profile is complete.
