# Nonlinear FMI Model Exchange in EMT three-phase MNA

`NonlinearSSN` keeps an FMI 3.0 Model Exchange instance alive and solves its
nonlinear continuous states together with the surrounding DPsim network.
The FMU is evaluated through its native functions; no fixed A/B/C/D fit or
exporter-specific generator equations are substituted. This implementation
uses one abc electrical port at each terminal and scalar Float64 variables.

The previous `LinearModel` / `LinearSSN` API remains available. Both import paths
share ZIP/XML/native-library handling through `ModelExchange`.

## Coupled numerical step

Write the FMU equations as `x'=f(x,v,c,t)` and passive terminal current as
`i=g(x,v,c,t)`, with control inputs c. At an endpoint trial `(x*,v*)`, form

```text
r = x* - x_old - h/2 (f_old + f*)
M = I - h/2 f_x
q = -M^-1 r
S = M^-1 (h/2 f_v)
W = g_v + g_x S
H = g* + g_x q - W v*
```

The component stamps `i = W v + H` into MNA. After the network solve, the
candidate state is `x* + q + S(v-v*)`. The true FMI equations are evaluated
again. Both the scaled trapezoidal residual and the difference between true
current and the stamped current must pass. Otherwise all opt-in nonlinear
components rebuild their tangents and the MNA solve repeats.

The new `MNANonlinearCompInterface` is opt-in. Its iteration loop is confined to
the direct variable-matrix EMT solver. Components that do not implement it keep
their previous behavior, including existing variable/affine SSN components and
legacy synchronous-generator iteration. Existing SSN classes are unchanged.
A solver path that does not perform the required convergence check is rejected
before the new component commits a state.

Default component settings are absolute/relative state tolerances 1e-11/1e-10,
absolute/relative current tolerances 1e-9/1e-10, a central finite-difference
Jacobian perturbation of 6e-6 times variable scale, and at most 20 iterations.
State scales include FMI state nominals. The aggregate solver has an additional
40-iteration bound. Failed convergence throws; it is not accepted as a step.

Only post-step commits state/history. `fmi3CompletedIntegratorStep` is called
when the model declares `needsCompletedIntegratorStep=true`; the FMI default is
false. FMI trial evaluations explicitly set continuous states, inputs and time
and do not call completed-step or discrete-state updates during Newton probes.
Accepted outputs/derivatives are read once and cached. When a hybrid FMU provides
FMI state snapshots, provisional evaluations save and restore the instance.
Snapshot support must still be validated for the particular exporter/model.
Completed-step and internal scheduled time events are handled at accepted
endpoints. Internal events must coincide with the fixed time grid; crossing one
is rejected. External reference changes coinciding with a controller tick are
applied in the same event iteration, so the controller sees the right-side
reference. Discrete updates never run during Newton probes.

## Public mapping and attributes

Create a persistent `ModelExchange` with all input names and initial values,
then create `NonlinearSSN` with three voltage names, three current names and an
explicit current sign. Terminal order is `{negativeABC, positiveABC}` and
voltage is positive minus negative. Sign +1 means the FMU reports passive
current into the positive terminal; -1 means it reports injection into the bus.
Raw FMI output attributes preserve their original sign.

The model provides `x`, `v_intf`, `i_intf`, `fmi.inputs`, `fmi.outputs`,
`fmi.input.<name>`, `fmi.output.<name>`, `fmi.state.<name>`,
`newton_iterations`, and `newton_residual` Attributes. Non-voltage inputs also
have `control.<name>` Attributes. A control may be driven with
`setInputFunction(name,function,discontinuityTimes)`; explicit event times must
align to the timestep. Input functions must match initial FMI input values.

Each component owns an independent runtime. Sharing one FMU instance between
clones is not supported. The exporter's singleton capability is respected;
load a separate instance only if that capability allows it. The SG MATLAB FMU
permits one instance per process, so run different scenarios in separate
processes.

## Reproduce the generator benchmark on this SSD

All newly generated results belong under `/Volumes/SSDMatthias/development`.
The original exported generator and references were copied to
`outputs/fmi-nonlinear/sg/matlab-origin`. The C++ example takes external-grid
parameters from `benchmark.json`, not from hard-coded machine equations.

Build from the repository root:

```sh
env DEVELOPER_DIR=/Library/Developer/CommandLineTools \
  cmake --build build-macos-check --parallel 3 --target \
  FMI_EMT_Ph3_Nonlinear FMI_NonlinearTest FMI_ImportTest FMI_ImportTest3Ph
```

In MATLAB, generate the independently implemented full-state trapezoidal
references and configuration:

```matlab
addpath('/Volumes/SSDMatthias/development/dpsim/dpsim-fmi/examples/matlab')
create_sg_trapezoidal_reference( ...
  '/Volumes/SSDMatthias/development/dpsim/outputs/fmi-nonlinear/sg/matlab', ...
  [50e-6,25e-6]);
```

From the repository root, a single run is:

```sh
mkdir -p outputs/fmi-nonlinear/sg/dpsim-50us
cd outputs/fmi-nonlinear/sg/dpsim-50us
env TMPDIR=/Volumes/SSDMatthias/development/dpsim/outputs/fmi-nonlinear/tmp \
  /Volumes/SSDMatthias/development/dpsim/build-macos-check/dpsim-fmi/FMI_EMT_Ph3_Nonlinear \
  ../matlab-origin/dpsim_sg_3ph_me.fmu ../matlab/benchmark.json 5e-5
```

Open `examples/FMI_EMT_Ph3_Nonlinear_MATLAB_Comparison.ipynb` with the DPsim
Python kernel and run all cells. It runs 50 us, 25 us, and 50 us with ten-times
tighter tolerances; validates all 34 saved physical/control signals; and plots
the matched waveforms, errors, controls and timestep refinement. Executed copies,
a standalone HTML and figures are in `outputs/fmi-nonlinear/sg/comparison`.
The notebook reuses the saved MATLAB references; it does not rerun MATLAB.

The comparison can also be invoked directly after all three runs:

```sh
dpsim-python/bin/python dpsim-fmi/tests/compare_nonlinear_sg.py \
  outputs/fmi-nonlinear/sg
```

Results from the initial validated run:

| Comparison | Maximum absolute discrepancy |
| --- | --- |
| DPsim vs MATLAB trapezoidal, 50 us, phase current | about 1.4 microampere |
| DPsim vs MATLAB trapezoidal, 50 us, bus voltage | about 2 microvolt |
| MATLAB trapezoidal vs ode15s, phase-A bus voltage, 50 us | 0.17054 V |
| MATLAB trapezoidal vs ode15s, phase-A bus voltage, 25 us | 0.04264 V |

All 34 signals pass peak/RMS normalized errors below 0.1% of their declared
engineering bases. The largest normalized error versus ode15s is the load-power
peak: 0.0671% at 50 us and 0.0168% at 25 us. Halving the step reduces the voltage
error by approximately four. At 50 us the importer uses two Newton iterations
per accepted step. Tighter Newton settings leave six-decimal CSV traces
unchanged. The logger's quantization limits interpretation of tiny differences.
Full absolute/normalized metrics and solver residuals are in `comparison.json`.

## Reference independence and disturbance convention

The MATLAB reference independently wires the sixteen-state equations and
solves one full trapezoidal residual, using complex-step Jacobians. DPsim uses
FMI evaluations and local state elimination with central finite differences.
Both share physical parameter values and initial conditions; the tighter
ode15s trajectory is an additional integration check. This validates numerical
implementation and wiring, not the manufacturer fidelity of the teaching model.

A step ending exactly on a known discontinuity uses the input from the prior
interval. At the boundary the new input is applied and the old derivative for
the next interval is recomputed. No half-step anticipation or trace shifting is
used. The external grid helper recomputes its state-space history for the load
conductance applicable to each interval; capacitor current history is not
silently carried across a changed algebraic load.

## Tests and current limits

`FMI_NonlinearTest` checks a nonlinear cubic/cross-phase fixture with shuffled
variable declarations, explicit injection polarity, two independent instances,
100 implicit steps against an analytical-Jacobian reference and MNA KCL. A
separate case forces non-convergence and checks that no state was committed.
The existing scalar/abc linear import and malformed-FMU tests continue to pass.
Native GFM Jacobian and induction-motor load-switch examples have also been run.

Current supported subset: continuous scalar Float64 boundary variables and
continuous states, internal discrete variables, scheduled external control changes,
and time/completed-step events that preserve continuous states and terminal
current. Internal event times must be strictly advancing and align with the fixed
DPsim timestep. An event that changes terminal current is rejected because it
requires an instantaneous algebraic network iteration. State events,
continuous-state reinitialization, clocks, FMI arrays, structural parameters
and unsupported execution-tool dependencies are rejected explicitly. Parameters
marked tunable retain their exported defaults and cannot be changed by this API.
Parameter values and initial conditions are the exported defaults. Runtime
parameter editing and nonlinear steady-state/load-flow initialization are not
implemented. No dynamic-phasor conversion or native circuit reconstruction is
performed. VILLAS remains separate.

FMI lifecycle/function semantics follow the
[FMI 3 specification](https://fmi-standard.org/docs/3.0.2/). The condensed Newton
method above is implemented in `NonlinearSSN.cpp` and the optional solver hook
in `MNASolverDirect.cpp`; it is not a claim that all FMI 3 FMUs are supported.

## Modelica exporter check

The same importer also runs a six-state nonlinear magnetic circuit exported by
OpenModelica as FMI 3.0 Model Exchange. It adapts three phase-shifted copies of
MSL 4.0.0 `ShowSaturatingInductor`; the native reference uses the original MSL
example instances. No exporter-specific equations are added to C++.
Peak current error falls from 58.1 µA at 10 ms to 14.9 µA at 5 ms, with RMS
refinement ratio 3.99. At 200 µs it reaches about 0.52 µA, near CSV rounding.
See [provenance and reproduction](examples/modelica/README.md) and
`examples/FMI_EMT_Ph3_Modelica_Comparison.ipynb`. Saved outputs and standalone
HTML are in `outputs/fmi-nonlinear/modelica`.

## Sampled-control tests and inspection

`FMI_SampledTest` compares a nonlinear continuous plant with a 1 ms sampled
controller against an independent trapezoidal recurrence. It checks a coincident
external command step at 2 ms, per-phase KCL and all accepted states. Separate
tests reject a timestep that crosses a controller tick and prevent the LTI
extractor from treating sampled dynamics as a frozen linear model.
`FMI_EventCurrentTest` rejects a completed-step event that changes terminal current
without a coupled algebraic network iteration.

`FMI_Inspect model.fmu [input_values.json]` initializes a model and reports state,
input and output names and the next internal time event. Every input must be
provided in the optional JSON object. This instantiates native FMU code.
`FMI_DiscreteProbe model.fmu stop_s step_s` is a standalone lifecycle diagnostic
for closed, zero-continuous-state models. It does not demonstrate DPsim MNA
coupling or continuous-state integration; its purpose is checking sampled event
semantics against an exporter's original simulation.

## Larger Simulink PV compatibility result

The official MathWorks `PVArrayGridAverageModelExample` was downloaded to the SSD
and exported as FMI 3.0 Model Exchange with its MPPT, boost converter, PLL,
DC-link and current-control loops retained. The original is a **fully discrete**
model: its FMU has zero continuous states and schedules internal events every
50 µs. `FMI_DiscreteProbe` replays all 120001 samples over six seconds, matching
native Simulink to about 2.2e-12 A and independent FMPy 0.3.32 to about 1.2e-12 A.
Those are event-lifecycle results, **not coupled DPsim MNA or trapezoidal plant
integration results**.

A provisional output getter before a tick changes this PV export's subsequent
history; even its advertised FMI state restore does not undo that effect in the
isolated protocol test. It is therefore not validated for Newton probing.
A separate continuous-powergui copy retains sampled controls but fails MATLAB's
ME export because `sfun_spid_contc` is a non-inlined S-function. An attempted
replacement with ordinary State-Space blocks is regenerated back to the native
S-function during model update. Also, the network contains a physical switch;
one fixed `power_analyze` matrix must not be substituted for all switch states.
The original model files remain unchanged.

The coupled larger-converter benchmark remains unfinished. Results and diagnostic
scripts are under `outputs/fmi-nonlinear/pv` and related `pv-*.log` files. These
limitations must not be described as a successful general PV/GFL import.
MathWorks documents the [non-inlined S-function restriction](https://www.mathworks.com/help/fmuexport/ug/generate-model-exchange-FMU.html).
