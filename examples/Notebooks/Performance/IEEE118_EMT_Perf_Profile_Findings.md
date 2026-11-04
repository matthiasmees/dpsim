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

## AMD GPU solver fit

An optional HIP dense direct-solver adapter was built and run on the laptop's
AMD Radeon 610M. The benchmark was a short IEEE118 EMT run with logging
disabled. Conventional EMT took 7.23 s with KLU, 100.26 s on the first HIP
run, and 99.19 s on a warm HIP run. In the reduced `allssn` case, KLU took
1.40 s and HIP dense took 22.77 s. The GPU was active during HIP runs, but
the dense conversion and repeated host/device traffic outweighed its compute
throughput. The dense GPU solver is therefore an experimental backend, not a
performance choice for this case.

The existing optimized conventional profile attributes about 0.215 ms of a
3.648 ms step to the solver task (roughly 6%). Even an infinitely fast linear
solve would cap total step speedup near 1.06x under that profile. HIP's legacy
sparse QR entry point also does not provide the factor-once/refactor-many
interface used by DPsim's timestep loop, so replacing dense conversion with
that API is not a promising next implementation.

For this laptop and IEEE118, the useful optimization target remains the CPU
work around the solve: conventional composite component pre/post-step work,
RHS assembly, and the remaining component state updates. The existing compact
RHS work already reduced the optimized conventional profile to 3.648 ms per
step; GPU offload should be reconsidered only for a substantially larger,
regular kernel with enough device-resident work to amortize transfers and
launches.

## HIP device-resident prototype: implementation map

The benchmark's current time-step path is host driven: `Simulation::step()`
calls `SequentialScheduler::step()`, which invokes each component task on the
CPU. Component state is held in Eigen matrices and attributes. IEEE118 EMT
uses 54 conventional four-order VBR generators by default; their model code
updates rotor state, dq/abc transforms, Norton conductance, history voltage,
and RHS contributions in separate CPU callbacks. Keeping only the linear
solve on the GPU cannot keep this data resident.

The HIP dense adapter now exposes its persistent device RHS buffer and an
in-place device solve entry point while retaining its host `solve()` wrapper.
The wrapper build and a short IEEE118 `allssn` HIP run succeeded (21.09 s).
This is only the first data-path seam: the system matrix is still converted
from sparse Eigen storage to dense host storage and copied to the device, and
component callbacks still execute on the CPU.

The build now enables CMake's HIP language when `WITH_HIP=ON` and compiles a
first device kernel for clearing a resident RHS vector and scattering indexed
component contributions with atomic accumulation. A small Radeon execution
probe with duplicate row indices produced the expected sum. This kernel is
not yet connected to the IEEE118 timestep; it is the initial assembly primitive
for the later batched component task.

The installed ROCm libraries provide sparse matrix primitives and a legacy
sparse QR solve call, but no drop-in reusable sparse direct factorization with
the lifecycle DPsim uses. A resident implementation therefore needs either a
GPU iterative solve with convergence/residual checks and a suitable
preconditioner, or a purpose-built reusable sparse factorization. The installed
rocSPARSE headers include CSR matrix-vector multiply, ILU(0), and triangular
solve routines, which are the building blocks for a preconditioned Krylov
solver. The next implementation stages are to validate that sparse solver on
the IEEE118 MNA matrices, pack component parameters and state into flat device
arrays, batch generator and line/load update kernels, assemble the variable CSR
values and RHS on-device, and replace per-component scheduler callbacks with a
device-resident step task. CPU synchronization should happen only for events,
measurements, and output. Numerical traces and full-case timings are required
before treating this as a usable backend.

The first sparse-solver experiment adds a GPU BiCGSTAB path with device-side
CSR SpMV and vector updates. On the one-second `allssn` run it completed with
KLU-comparable traces (1,001 rows, 28 columns; maximum logged absolute
difference 5.66e-4 V, about 2.93e-6 relative at the largest discrepancy), but
took 45.41 s versus 0.98 s with KLU. The conventional EMT matrix has zero
diagonal entries from its MNA constraints; a Jacobi preconditioner is therefore
not valid, and row equilibration still did not converge. This solver is useful
as a device-resident data-path prototype, but it is not yet a robust solver for
the default conventional IEEE118 case. The next solver work needs an
MNA-aware preconditioner or a sparse direct method before extending device
residency to all component states.

A sparse QR reference adapter does solve the one-second conventional IEEE118
case, but it took 223.22 s with logging disabled, compared with 6.44 s for KLU
under the same command. It refactorizes on every sparse QR call. This confirms
that preserving sparsity alone is insufficient on this device; reusable
factorization or a much cheaper iterative path is necessary for a viable
backend.

## GPU-resident port status and next implementation boundary

The HIP path now has persistent device state for two IEEE118 component groups:
the fused Pi-line/load/transformer models in `allssn aggregate`, and the 54
reduced-order SG4 Norton generators used by the conventional topology. The SG4
batch kernel performs the swing update, Park transforms, changing conductance,
Norton source, and post-solve dq current update. It currently accepts only
ungrounded Norton generators without attached controllers; unsupported models
continue through their CPU callbacks. In both paths, model state persists on the
GPU, but the host receives component outputs for the existing Eigen MNA stamps
and RHS accumulation. The host also uploads the changing CSR values and RHS
each step, and copies the solution back for node propagation.

The one-second conventional IEEE118 EMT run produced the same 1,002-by-28
logged trace shape as KLU. HIP dense differed by at most 1.01e-6 V across
28,028 finite values (peak-scaled relative difference 4.53e-10). Its full
process time was 99.81 s versus 6.68 s for KLU; simulation-calculation times
were 96.689 s versus 3.920 s. Thus the SG4 kernel is numerically validated in
the default topology, but the HIP path is about 15 times slower end to end on
this Radeon 610M.

On the one-second `allssn aggregate` run, HIP dense differed from KLU by at
most 1.01e-6 V, while HIP BiCGSTAB differed by 5.64e-4 V due to iterative-solve
error. Full-process times were 0.88 s for KLU, 21.55 s for HIP dense, and
33.29 s for HIP BiCGSTAB. Neither GPU solver is faster. The component batches
work numerically, but they do not yet avoid per-step host assembly and data
transfers.

To remove those transfers, the next changes are to port the conventional
passive component histories and stamps, then let the generator/fused component
kernels update the CSR entries and RHS directly in the solver's device
buffers. The current intermediate result moves component state updates to the
GPU while leaving MNA assembly on the host; it is not yet an all-GPU timestep.

## HIP-resident SG4 assembly prototype

The sparse HIP BiCGSTAB path now keeps the IEEE118 SG4 Norton conductance and
source assembly on the GPU. The host still resets and stamps the remaining
variable matrix entries and assembles RHS contributions from components that
do not have device kernels. It transfers only the variable CSR values and the
remaining host RHS to the device; SG4 conductance values are scattered into
the resident CSR array and SG4 source currents are atomically accumulated into
the resident RHS. The solution is copied to the host for node voltage
propagation, while SG4 post-step state is updated on the device.

For a one-second `allssn aggregate` run, the device-assembly HIP trace had
1,001 rows and 28 columns and differed from KLU by at most 1.87e-3 V across
the logged values (peak-scaled maximum 6.32e-9). The HIP process took 48.1 s;
the matched KLU process took 0.87 s. This prototype is about 55 times slower
than KLU despite moving SG4 matrix and RHS contributions to the device. The
conventional case still does not converge with BiCGSTAB because of its MNA
constraint matrix (reported residual 5.64e4 after 3.7 s). The existing HIP
dense route remains the numerically validated conventional GPU option, at
99.81 s versus 6.68 s for KLU in the prior one-second run.

The CPU-only build and KLU solver remain available; the device assembly is
enabled only for the HIP BiCGSTAB implementation and supported ungrounded SG4
Norton generators. Remaining host assembly, solver convergence cost, and
solution synchronization dominate this early prototype, so it does not yet
deliver an IEEE118 speedup.

## Fully resident IEEE118 aggregate timestep

The device-resident path now covers the `allssn aggregate` IEEE118 model for
the complete electrical timestep. The fused PiLine/load/transformer batch and
SG4 batch keep their histories and states on the GPU; the fused and generator
models restore/update CSR entries and accumulate the RHS directly in device
memory. A single HIP kernel performs the full BiCGSTAB iteration, reductions,
and convergence decision without returning dot products or solver state to the
host. With logging disabled, the solved vector and component state stay on the
device and host node propagation is skipped. Logging mode copies component
outputs back so the existing CSV logger can record them.

The logged one-second run still produced a 1,001-by-28 trace, with maximum
absolute difference 1.87e-3 V from KLU. For matched one-second `nolog` runs,
the resident HIP path took 7.49 s and KLU took 0.875 s; the GPU path is about
8.6 times slower on this Radeon 610M. Moving the computation off the CPU did
not make this sparse iterative solve faster.

`hipbicgstab` selects the device-compatible `allssn aggregate` representation
automatically for the standard IEEE118 command. Pass `conventional` or `paper`
with it and the benchmark rejects the run because those component formulations
still need CPU-side assembly. The DPsim host scheduler advances time and submits
HIP work; electrical state updates, matrix/RHS assembly, and solver iterations
run on the GPU. HIP mode disables the host-side per-step logger so solved
voltages remain resident. DPsim's existing `RealTimeDataLogger` buffers rows in
CPU memory, but it still reads attributes on the CPU every step, so it does not
meet that requirement. The `hipdense` and `hipsparseqr` benchmark flags are
rejected because those earlier prototypes use CPU-side assembly.

The short GPU-only invocation is `emt hipbicgstab short`; `nolog` is no longer
needed. A device-side trace buffer with a single end-of-run transfer would be
needed to preserve a timestep trace without returning values to the CPU each
step.

## rocSOLVER sparse LU refactorization prototype

The `hiplu` mode now seeds rocSOLVER's CSRRF sparse refactorization from one
KLU factorization during solver initialization. It then refactors the fixed
CSR pattern and solves each timestep on the GPU, using the device-resident
matrix and RHS. This removes BiCGSTAB iteration count and single-block
underutilization from the solve. A startup probe compares the rocSOLVER seed
solve with KLU and fails early if the initial factors or permutations disagree.
The one-time seed factorization runs on the CPU; the timestep refactorization
and solves run on the GPU. `hiplu` is available when the build includes both
HIP and KLU.

On matched one-second, 50 µs IEEE118 aggregate runs, `hiplu` took 2.55 s,
BiCGSTAB took 7.46 s, and KLU took 0.83 s. The refactorization prototype is
about 2.9x faster than the current GPU BiCGSTAB path, but remains about 3.1x
slower than KLU on the Radeon 610M. A logged 1,001-row trace stayed finite;
the bus-voltage channels differed from KLU by at most 1 µV, and all logged
channels by at most 1.87e-3 in their native units. `hiplog` enables that
host-copied trace for verification; without it, results stay on the GPU.

The run confirms that fixed-pattern GPU refactorization is numerically
feasible on IEEE118, but KLU remains the faster choice for this system. Larger
networks may provide more GPU parallel work, but the IEEE300 test is not yet
stable enough to confirm an advantage.

The CPU-only configuration was also rebuilt with HIP disabled and ran the
one-second IEEE118 aggregate case through KLU successfully. Its build exposed
an existing signed/unsigned comparison warning promoted to an error in the KLU
adapter; the comparison now explicitly converts Eigen's nonnegative nonzero
count to `std::size_t`.

## IEEE300 GPU validation and profiling follow-up

The initial one-second IEEE300 comparison used the default 500 µs step, and
both CPU and GPU traces became unstable. Repeating it with the benchmark's
50 µs step did not resolve the case's instability. On matched 1 s runs
(20,000 solver steps), both KLU and `hiplu` produced non-finite traces at
21.55 ms in the final HIP-enabled `-O3` build. The timed no-logger runs took
10.24 s for KLU and 37.37 s for `hiplu`; these are not valid performance
results because the simulated trajectory is not finite.

The `hiplog` validation path now synchronizes the SG4 mechanical state
attributes to the host, so logged torque and speed channels reflect device
state. It also checks the relative infinity-norm residual of every rocSOLVER
LU solve when host trace synchronization is enabled. This separates
linear-solver residual accuracy from instability in the model trajectory.

HIP-enabled Release builds use `-O3` instead of `-Ofast`, which permits
unsafe floating-point reassociation. This did not resolve IEEE300's instability;
CPU and GPU traces now become non-finite at the same logged timestep.

A 2,000-step profile of the IEEE300 500 µs run attributed 3.02 s total
(1.51 ms per solve) to GPU numeric solves, compared with 45 ms total
(22.7 µs per solve) for KLU. This confirms the sparse refactorization is the
dominant performance cost on the Radeon 610M. The IEEE300 case needs a stable
CPU reference trajectory before further GPU solver tuning can be evaluated
reliably. The same-build finite IEEE118 run took 2.62 s with `hiplu` and
0.86 s with KLU, with bus-voltage traces within 1 µV.

A single-precision rocSOLVER refactorization with one double-precision
residual correction was also benchmarked as `hiplu32`. It took 3.01 s on
IEEE118 and 44.11 s on IEEE300, compared with 2.59 s and 37.37 s respectively
for double-precision `hiplu`; the mixed-precision prototype was removed.

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
