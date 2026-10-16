# FMI 3.0 Model Exchange to EMT SSN

Nonlinear three-phase import is available through `ModelExchange` and
`NonlinearSSN`; see [the nonlinear importer guide](NONLINEAR.md). The sections
below describe the fixed LTI extraction path and its limitations.

`dpsim-fmi` imports a **declared linear, time-invariant** Model Exchange FMU
as a two-terminal scalar or three-phase EMT component. Its internal state count and equations
come from the FMU. The component connects to ordinary DPsim nodes and contributes
an SSN Norton equivalent to the existing MNA solver. No dedicated R/L/C classes,
Co-Simulation `doStep()`, MATLAB runtime, or VILLAS interface is involved.

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
the pinned miniz 3.0.2 ZIP library. Official FMI 3.0.2 C headers and their license
are included under `third_party/fmi3`. No FMI dependencies are searched or fetched
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

On macOS, the selected compiler and SDK must match. If `xcode-select` points to
an older Xcode for MATLAB while your existing DPsim build uses newer Command Line
Tools, prefix **both** configure and build with
`env DEVELOPER_DIR=/Library/Developer/CommandLineTools`. This only selects tools
for that command and leaves the system setting used by MATLAB unchanged.

Build the named targets rather than every upstream example when OpenMP is off;
some unrelated Diakoptics examples assume it is available.

## Run the MATLAB RLC benchmark

```sh
mkdir -p outputs/fmi
cd outputs/fmi
../../build-fmi/dpsim-fmi/FMI_EMT_LinearSSN /absolute/path/dpsim_rlc_me.fmu u i
```

This constructs one imported branch in parallel with an ideal 230 V RMS,
50 Hz voltage boundary, using a 10 microsecond step and zero-time FMU
initialization. No R, L or C parameters are supplied by the C++ example.
Outputs are written to `logs/fmi_rlc.csv`.

For a complete runnable comparison, open
[`examples/FMI_EMT_LinearSSN_MATLAB_Comparison.ipynb`](examples/FMI_EMT_LinearSSN_MATLAB_Comparison.ipynb)
in VS Code, select the **DPsim Python** kernel, and run all cells. It runs the
actual C++ example, plots the MATLAB and DPsim waveforms and errors, and saves
maximum/RMS error metrics under `outputs/fmi-notebook`. The first cell contains
the configurable FMU, reference-data and build paths.

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

To compare the example's CSV with the MATLAB script's references:

```sh
python3 dpsim-fmi/tests/compare_reference.py \
  outputs/fmi/logs/fmi_rlc.csv /absolute/path/to/matlab-output
```

The script checks the complete 10,001-sample time grid and agreement with the
trapezoidal reference within DPsim's six-decimal CSV rounding. It separately
reports errors against the tightly toleranced Simulink `ode45` reference.

## Three-phase MATLAB benchmark

The MATLAB exporter is `examples/matlab/create_rlc_fmu_3ph.m`. It generates a
six-state grounded-wye RLC benchmark, an FMI 3.0 Model Exchange FMU, and both
adaptive Simulink and fixed-step trapezoidal reference CSVs. See its adjacent
README for MATLAB instructions.

From the repository root after building the targets above:

```sh
mkdir -p outputs/fmi-3ph
cd outputs/fmi-3ph
../../build-fmi/dpsim-fmi/FMI_EMT_Ph3_LinearSSN /absolute/path/dpsim_rlc_3ph_me.fmu
cd ../..
python3 dpsim-fmi/tests/compare_reference_3ph.py \
  outputs/fmi-3ph/logs/fmi_rlc_3ph.csv /absolute/path/to/matlab-output
```

The example defaults to input names `u_a/u_b/u_c` and current names
`i_a/i_b/i_c`. To override, append all six names in that order. Only the external
source and experiment settings are specified in C++; branch matrices and state
count are recovered from the FMU. The source is 230 V phase-to-neutral RMS,
50 Hz, with phase angles `[0,-120,+120]` degrees.

Applications construct the same `LinearSSN` class with an explicit mapping:

```cpp
DPsim::FMI::VoltagePorts ports{
    {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"},
    {0.0, -281.69132042006549, 281.69132042006549}};
auto model = DPsim::FMI::LinearModel::loadMulti("branch.fmu", ports);
auto branch = DPsim::FMI::LinearSSN::make("imported", model);
branch->connect({negativeABCNode, positiveABCNode}); // ground is also accepted
```

`initialVoltages` are instantaneous SI voltages at t=0, supplied to FMI
initialization **and** retained as the initial SSN input history. An empty vector
means zeros. They must match the application's initial terminal voltage
condition. Clones reuse the same recovered matrices and initialized state;
they do not reinitialize the FMU for a different starting condition. For the
balanced sine benchmark phases B/C start at approximately ±281.691 V; assuming
all three inputs start at zero would produce an incorrect first step.

Both terminals must be ABC nodes or ground for a three-channel model.
`v_intf` and `i_intf` are 3×1 abc vectors; each `fmi.input.<name>` and
`fmi.output.<name>` remains a scalar attribute. Output declaration order in XML
does not determine the phase mapping: current rows are selected by name.

Open
[`examples/FMI_EMT_Ph3_LinearSSN_MATLAB_Comparison.ipynb`](examples/FMI_EMT_Ph3_LinearSSN_MATLAB_Comparison.ipynb)
with the DPsim Python kernel to run and compare all nine traces. Configure the
paths in the first code cell. Saving the executed notebook embeds the plots.
The comparison script checks 10,001 samples over 0.1 s, nine channels within
six-decimal logging precision, each phase's source/branch KCL, and balanced
current sums. The ode45 comparison includes interpolation to DPsim timestamps.

## Nonlinear benchmark for the next importer stage

`examples/matlab/create_sg_grid_fmu_3ph.m` builds and simulates a small
three-phase synchronous-generator grid with field/stator flux dynamics, rotor
motion, an AVR, a droop governor, and scheduled load/voltage-reference steps.
The complete grid has sixteen states; the optional generator/control Model
Exchange FMU has ten. See
[the nonlinear benchmark documentation](examples/matlab/nonlinear_3ph/README.md)
for equations, operating-point initialization, parameters and reference checks.
This is reference material for a future nonlinear adapter and does not extend
`LinearSSN` to nonlinear equations.

## Boundaries of this implementation

- Calling `LinearModel::load` asserts that the model is homogeneous LTI.
  Numerical probes detect many violations but do not prove global linearity.
  This is not a general nonlinear FMU or local-linearization simulator.
- Co-Simulation-only FMUs, arrays at the boundary/in state vectors, clocks,
  discrete/tunable variables, affine offsets, state/time events and external
  execution-tool dependencies are rejected. Fixed parameters use FMU defaults.
- The electrical adapter supports one scalar pair or one ordered abc bundle
  at each terminal. Scalar Float64 FMI variables represent the three channels;
  FMI array variables, arbitrary multiport layouts, dynamic-phasor conversion
  and automatic electrical topology discovery are outside this adapter.
- The native binary must match the running architecture and operating system.
  Source-only FMUs are not compiled automatically. Only macOS ARM64 has been
  exercised here; Linux/Windows loading paths are present but unverified.
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
The actual MATLAB-exported benchmark is a separate end-to-end check and is not
committed to this repository.
