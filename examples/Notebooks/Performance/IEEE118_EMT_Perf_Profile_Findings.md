# IEEE118 EMT `perf` findings

Historical note: older benchmark artifacts in this folder refer to a separate
`PiLine_NoSub` class. That class has been removed; the optional direct MNA path
now lives in `EMT::Ph3::PiLine`. The old standalone timing is not representative
of the integrated constructor, and its claim of removing the internal midpoint
was incorrect for the validated implementation. See
`examples/Notebooks/Benchmarks/IEEE300_EMT_Publishable_Performance_Benchmark.ipynb`
for current trace-checked results.

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

## Compact RHS follow-up for reduced-order VBR generators

The IEEE118 benchmark uses `ReducedOrderSynchronGeneratorVBR` (not the
separate full `SynchronGeneratorVBR` class). Its source stamp has only three
active rows, but previously each of 54 generators added a full network-sized
vector into the global RHS on every step. The reduced-order base now advertises
and implements the existing compact RHS interface, selecting the three bus
rows for Norton mode or the three internal voltage-source rows otherwise.

On matched 10 s, 50 µs profiles (200,000 steps), generator dense stamp calls
fell from 10.8 million to zero. The resulting RHS-assembly totals were:

| Case | Before (s) | After (s) | Change |
| --- | ---: | ---: | ---: |
| Conventional | 4.315 | 2.682 | −37.8% |
| All-SSN | 2.925 | 2.492 | −14.8% |
| Aggregate | 1.387 | 0.950 | −31.5% |
| Grouped shared-map | 2.676 | 2.639 | −1.4% |

The grouped case gains little because the group callback and state updates
dominate its RHS phase. In five clean 1 s, 50 µs aggregate runs, the wall-time
mean was 2.258 s after the change versus 2.435 s in the matched five-run
baseline; the medians were 2.249 s versus 2.367 s. Full-process timings include
power-flow setup and vary more than the isolated RHS phase.

The change preserves summation order at each RHS row. A same-revision disabled
versus enabled trace comparison matched all 10,001 logged samples across 28
channels bit-for-bit. Detailed 10 s per-task and solver-phase data, call counts,
and plots are in
[IEEE118_EMT_Deep_Solver_Profile.ipynb](IEEE118_EMT_Deep_Solver_Profile.ipynb).

## Half-decoupling line constant-matrix work

The half-decoupling line history-current update had been inverting its fixed
terminating impedance and its square inside every timestep. Those matrices are
constant after `setParameters()`, so their inverses are now calculated once
during parameter setup and reused by both the direct stamp and history update.
The timestep expression and multiplication order are otherwise unchanged.

Five post-change runs of the EMT direct-component comparison executable were
compared with ten pre-change runs. For the legacy and direct half-line cases,
the median full-case timings changed from 0.10950 s to 0.09177 s (−16.2%) and
from 0.11130 s to 0.09474 s (−14.9%), respectively. The pre-change direct mean
was affected by one 0.189 s outlier, so medians are the more representative
comparison. All four old/new legacy/direct trace files had 10,001 samples and
13 columns; all 260,026 compared numeric values matched exactly at the logged
precision.

The DP `PQLoadCS` direct path was also retested in five runs against its legacy
component path. It was about 2% slower in that small fixture, with both paths
under 1.4 ms per run and identical full traces. It is therefore not retained
as a claimed optimization target; the fixture is too small for a useful speed
gain and the observed difference is within timing noise.

## EMT PiLine non-subcomponent constructor validation

The optional no-subcomponent path now lives in the existing `EMT::Ph3::PiLine`
class; the historical standalone `PiLine_NoSub` class was removed. Its direct
constructor stamps and updates the three-phase Pi R/L/C/G model without
constructing the child components. Five matched IEEE118 EMT runs at 1 s and
50 µs compared the conventional composite topology with `nosubline`; all runs
used the same conventional loads, shunts, transformers, and generators.

Both variants produced 10,002×28 finite traces. The traces were numerically
very close but not byte-identical: maximum absolute logged difference was
1.0e-6, peak-normalized maximum difference 3.38e-12, and peak-normalized RMSE
2.13e-13. The simulation-calculation mean was 5.759 s composite versus 4.084 s
non-subcomponent (29.1% faster); the median was 5.772 versus 4.080 s, with the
non-subcomponent path faster in all five pairs. This validates the prior
profiled speed advantage against the electrical output and replaces the older
single-run, 0.5 ms profile-only claim. The per-task plots and breakdown remain
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
