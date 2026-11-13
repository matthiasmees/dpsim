# FMI 2.0/3.0 import for DPsim EMT

DPsim now supports both Model Exchange and Co-Simulation FMUs. `ModelExchange`
and `NonlinearSSN` expose continuous states to the MNA/Newton solver. The new
`CoSimulation` wrapper keeps the FMU instance alive and drives it through its
standard `doStep` interface; `CoSimulationSSN` connects its scalar or abc
voltage/current port to DPsim with fixed communication steps. This adapter is
explicit: each network solve uses the last completed FMU current, then the FMU
advances with the preceding PCC voltage held over one communication interval.
It does not claim implicit Newton convergence or rollback support.

`FMI_EMT_CoSimulation`, `FMI_EMT_Ph3_CoSimulation`, and
`FMI_EMT_Ph3_TwoCoSimulations` are scalar, three-phase, and parallel two-FMU
examples with finite grid impedance. They accept FMI 2.0 or FMI 3.0
Co-Simulation FMUs with scalar Real/Float64 electrical variables. See
`examples/cosimulation/README.md` for how they map to the existing example
families and for their port signs.
The checked-in generator and RLC archives declare Model Exchange; they cannot
be run as Co-Simulation FMUs. Matching Co-Simulation exports must be generated
by their source tools before using the Co-Simulation example.

Nonlinear three-phase import is available through `ModelExchange` and
`NonlinearSSN`; see [the nonlinear importer guide](NONLINEAR.md). The sections
below describe the fixed LTI extraction path and its limitations.

For a conventional three-phase nonlinear FMU, the short C++ form is
`auto fmu = DPsim::FMI::importFMU("converter.fmu"); fmu->connect(pcc);`.
The equivalent Python entry point is
`dpsimpy.fmi.import_fmu("converter.fmu")`, followed by
`fmu.connect(pcc)`. See the [Python import notebooks](examples/python/README.md)
for executable examples and explicit mapping options.

`dpsim-fmi` imports a **declared linear, time-invariant** Model Exchange FMU
from FMI 2.0 or FMI 3.0
as a two-terminal scalar or three-phase EMT component. Its internal state count and equations
come from the FMU. The component connects to ordinary DPsim nodes and contributes
an SSN Norton equivalent to the existing MNA solver. The Co-Simulation path is
documented separately below; MATLAB runtime and VILLAS are not required by
either importer.

This first implementation supports homogeneous models

```text
dx/dt = A x + B v
    y = C x + D v
```

with one or three scalar voltage inputs, matching current outputs, and optional
scalar diagnostic outputs. Three-phase inputs/currents are ordered abc. DPsim's existing `CPS::EMT::VTypeSSNComp` implements
trapezoidal integration, history terms, scheduling and post-solve state updates.
`LinearSSN` supplies scalar or full 3×3 electrical stamps and metadata-derived
attributes. Cross-phase terms in B/C/D and the resulting Norton matrix are retained.
The initial state is obtained from FMI initialization and is retained, including
nonzero initial states; it is not replaced by an AC steady-state solution.

## Architecture

- `LinearModel::load` / `loadMulti`: owns temporary extraction, parses XML, loads the native
  library, performs FMI Model Exchange initialization, recovers A/B/C/D by
  symmetric basis evaluations, checks mixed-state/amplitude/time probes, and
  releases the FMI instance and library. The returned immutable model can be
  shared by several SSN component instances without instantiating the FMU again.
- `LinearSSN`: a new DPsim electrical component derived from the general SSN
  base. The model's mapped current-output rows determine its MNA conductance and history
  contribution; all output rows remain available as attributes. Terminal order
  is `{negative, positive}`, so voltage is `v1-v0` and current flows `1 -> 0`.
- `VoltageBoundary` / `VoltageBoundary3Ph`: example/test helpers imposing a prescribed voltage with
  a branch-current unknown in MNA. It represents excitation outside the FMU.

The archive contains a black-box electrical branch, not an internal netlist.
The importer therefore does **not** infer resistor/inductor/capacitor elements
or invent internal circuit nodes. The caller explicitly maps FMI variable names
to an electrical port and connects that port to the surrounding DPsim topology.
The current implementation uses SI values without automatic unit conversion.

## Build

The feature is optional and disabled by default. It requires LibXml2 and fetches
the pinned miniz 3.0.2 ZIP library. FMI 2.0 and official FMI 3.0.2 C headers and
their licenses are included under `third_party/fmi2` and `third_party/fmi3`.
No FMI dependencies are searched or fetched
when `WITH_FMI=OFF`.

From the repository root, with normal DPsim dependencies installed:

```sh
cmake -S . -B build-fmi -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DWITH_FMI=ON -DDPSIM_FMI_BUILD_TESTS=ON \
  -DFETCH_SPDLOG=ON -DFETCH_GRID_DATA=OFF \
  -DDPSIM_BUILD_DOC=OFF -DWITH_VILLAS=OFF -DWITH_OPENMP=OFF
cmake --build build-fmi --parallel 3 --target \
  FMI_EMT_LinearSSN FMI_EMT_Ph3_LinearSSN FMI_ImportTest FMI_ImportTest3Ph
ctest --test-dir build-fmi -R '^fmi_' --output-on-failure
```

Build the named targets rather than every upstream example when OpenMP is off;
some unrelated Diakoptics examples assume it is available.

## Run the FMI examples

Build with `WITH_FMI=ON` and `DPSIM_BUILD_EXAMPLES=ON`, then run the scalar or
three-phase linear examples with an FMU archive:

```sh
build-fmi/dpsim-fmi/FMI_EMT_LinearSSN path/to/model.fmu u i
build-fmi/dpsim-fmi/FMI_EMT_Ph3_LinearSSN path/to/model.fmu
```

The checked-in sample archives and additional nonlinear/co-simulation examples
are listed in [`examples/README.md`](examples/README.md). Logs and other generated
files belong in the configured build directory.

For an application, the essential construction is:

```cpp
auto model = DPsim::FMI::LinearModel::load("branch.fmu", {"u", "i"});
auto branch = DPsim::FMI::LinearSSN::make("imported", model);
branch->connect({negativeNode, positiveNode});
// Add branch, the nodes, and external boundary components to SystemTopology.
// Run a DPsim Simulation with Domain::EMT and doSteadyStateInit(false).
```

Available attributes are `x`, `v_intf`, `i_intf`,
`fmi.input.<inputName>`, `fmi.output.<outputName>` and
`fmi.state.<stateName>`. Derived output attributes depend explicitly on `x`
and `v_intf`, so consumers participate in DPsim's scheduler dependency graph.
State names/order follow ModelStructure and the runtime state vector, not
assumptions about the exporter's internal coordinates.

## Three-phase connections

The scalar and three-phase examples use the same `LinearSSN` implementation.
For a three-phase archive, specify the voltage/current port names and initial
terminal voltages when the FMU uses nonzero initialization conditions:

```cpp
DPsim::FMI::VoltagePorts ports{
    {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"},
    {0.0, -281.69132042006549, 281.69132042006549}};
auto model = DPsim::FMI::LinearModel::loadMulti("branch.fmu", ports);
auto branch = DPsim::FMI::LinearSSN::make("imported", model);
branch->connect({negativeABCNode, positiveABCNode}); // ground is also accepted
```

`initialVoltages` are instantaneous SI voltages at t=0, supplied to FMI
initialization and retained as the initial SSN input history. An empty vector
means zeros. They must match the application's initial terminal voltage
condition. Clones reuse the same recovered matrices and initialized state;
they do not reinitialize the FMU for a different starting condition.

When a `LinearSSN` component is connected to DPsim nodes, its initialization
uses the connected nodes' instantaneous initial voltage difference for the SSN
history and electrical input attributes. This applies to scalar, grounded
three-phase, and floating three-phase connections. The recovered linear model's
state remains the state obtained during FMU import; changing the network start
voltage therefore aligns the terminal history but does not reconstruct a new
internal FMU operating point. For a state-dependent operating point, provide
matching `initialVoltages` when importing the model or use the nonlinear
Model Exchange path, which evaluates the FMU at the connected voltage.

Both terminals must be ABC nodes or ground for a three-channel model.
`v_intf` and `i_intf` are 3×1 abc vectors; each `fmi.input.<name>` and
`fmi.output.<name>` remains a scalar attribute. Output declaration order in XML
does not determine the phase mapping: current rows are selected by name.

Build `FMI_EMT_Ph3_LinearSSN` to run the three-phase example. The
[Python import quickstart](examples/python/README.md) demonstrates automatic
port discovery and explicit mapping with the checked-in three-phase FMU.
Simulation outputs are written under the configured build directory.

## Nonlinear example

`examples/synchronous_generator/matlab/create_sg_grid_fmu_3ph.m` builds and simulates a small
three-phase synchronous-generator grid with field/stator flux dynamics, rotor
motion, an AVR, a droop governor, and scheduled load/voltage-reference steps.
The complete grid has sixteen states; the optional generator/control Model
Exchange FMU has ten. See
[the example model documentation](examples/synchronous_generator/matlab/nonlinear_3ph/README.md)
for equations, operating-point initialization, parameters and reference checks.
This is reference material for a future nonlinear adapter and does not extend
`LinearSSN` to nonlinear equations.

## Boundaries of this implementation

- Calling `LinearModel::load` asserts that the model is homogeneous LTI.
  Numerical probes detect many violations but do not prove global linearity.
  This is not a general nonlinear FMU or local-linearization simulator.
- Arrays at the boundary/in state vectors, clocks,
  discrete/tunable variables, affine offsets, state/time events and external
  execution-tool dependencies are rejected. Fixed parameters use FMU defaults.
- FMI 2.0 and FMI 3.0 Model Exchange share the same public `ModelExchange`,
  `LinearModel`, `LinearSSN`, and `NonlinearSSN` readers. FMI 2.0 uses its
  native lifecycle and event API; it is not converted by editing the FMU XML.
- FMI 2.0 and FMI 3.0 Co-Simulation import scalar Real/Float64 variables and
  synchronous `doStep`. Early return, asynchronous/pending steps, intermediate
  updates, clocks and non-scalar numeric types are rejected. FMI 3.0 event mode
  is serviced when the FMU declares it. The electrical adapter supports one or
  three voltage inputs and current outputs and fixed DPsim communication steps.
- The electrical adapter supports one scalar pair or one ordered abc bundle
  at each terminal. Scalar Float64 FMI variables represent the three channels;
  FMI array variables, arbitrary multiport layouts, dynamic-phasor conversion
  and automatic electrical topology discovery are outside this adapter.
- The native binary must match the running architecture and operating system.
  Source-only FMUs are not compiled automatically. The FMU archive must include a binary for the target operating system and
  architecture.
- `canBeInstantiatedOnlyOncePerProcess` is honored even after extraction has
  finished. Reuse the returned model to construct multiple components.
- Runtime events and parameter changes are not forwarded to the FMU after
  extraction. Standard SSN integration owns all subsequent state evolution.
- Archive paths/links, extraction sizes, XML DTDs, missing symbols and FMI errors
  are checked. An FMU binary is native code running inside the process.

The CTest fixtures exercise the importer ABI, matrix recovery, initial-state
preservation, derived attributes, parallel branch KCL, 100 steps against an
independent trapezoidal recurrence, capability/error handling and malformed
archives. Additional three-phase tests cover shuffled XML declarations,
cross-phase coupling, feedthrough, unbalanced excitation, input-dependent
initial states, clone behavior, phase validation, and both grounded and
ungrounded negative terminals. They are minimal test doubles, not certified general-purpose FMUs.

## Nonlinear FMU capability inspection

The example `FMI_Capabilities` inspects a Model Exchange FMU and reports its
continuous interface, event indicators, FMU-state capabilities, process
instantiation limit, and a 50 µs rollback replay probe. The probe distinguishes
rollback advertised by metadata from rollback deterministic enough for Newton
coupling:

```sh
build-fmi/dpsim-fmi/FMI_Capabilities path/to/model.fmu
```

Use `rollback_newton` only when the probe reports deterministic. Otherwise use
the generic `single_pass_explicit` multirate mode.
