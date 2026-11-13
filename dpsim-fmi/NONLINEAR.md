# Nonlinear FMI Model Exchange in EMT three-phase MNA

`NonlinearSSN` keeps an FMI 2.0 or FMI 3.0 Model Exchange instance alive and solves its
nonlinear continuous states together with the surrounding DPsim network.
The FMU is evaluated through its native functions; no fixed A/B/C/D fit or
exporter-specific generator equations are substituted. This implementation
uses one abc electrical port at each terminal and scalar Float64 variables.

The importer selects the FMI backend from `modelDescription.xml`; the FMU binary
and lifecycle remain native to that FMI version. The previous `LinearModel` /
`LinearSSN` API remains available. Both import paths
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
endpoints. If an internal time event falls inside a network step and the FMU
supports rollback, `NonlinearSSN` automatically replays the interval with
implicit substeps that land on every reported event. Network-endpoint voltage is
linearly interpolated across those substeps. Each outer Newton trial restores
the complete FMI checkpoint; only the converged replay is retained. FMUs without
rollback reject a crossed internal event by default. Such an FMU can instead
use the explicitly selected `MultirateMode::SinglePassExplicit` policy. That
policy solves the network with the previously accepted terminal current, then
advances the FMU exactly once through every internal event using the solved
endpoint-voltage trajectory. It reports the terminal-current discrepancy and
does not claim converged implicit coupling. External reference changes
coinciding with a controller tick are
applied in the same event iteration, so the controller sees the right-side
reference. Discrete updates never run during Newton probes.

## Public mapping and attributes

For the conventional three-phase interface, the short C++ path is:

```cpp
auto fmu = FMI::importFMU("converter.fmu");
fmu->connect(pcc);
```

This reads the component name and scalar input start values from the FMU and
infers `v_a/v_b/v_c` (or `u_a/u_b/u_c`) and `i_a/i_b/i_c`. The one-node
connection grounds the negative terminal. Add the returned component to the
system's component list as usual. Use the named `importFMU(name, path, ...)`
overload when you need a chosen component name, or create a `ModelExchange`
and provide `ElectricalPort3Ph` explicitly for nonstandard ports or current
polarity. The one-argument convenience import enables endpoint handling for
FMUs with state events; pass `ModelExchangeOptions` explicitly to control that
behavior.

The Python binding follows the same pattern:

```python
fmu = dpsimpy.fmi.import_fmu("converter.fmu")
fmu.connect(pcc)
```

It is included in `dpsimpy` when both `WITH_FMI` and `WITH_PYBIND` are enabled.
Use `import_fmu_with_mapping(...)` for nonstandard abc variable names or
current polarity. See [the Python examples](examples/python/README.md).

`ModelExchange(path)` discovers continuous scalar Float64 inputs and uses their
declared start values. `NonlinearSSN(name, model)` recognizes the conventional
`v_a/v_b/v_c` and `i_a/i_b/i_c` electrical interface with passive-positive
current. Nonstandard names or polarity use the explicit mapping overload.
Terminal order is `{negativeABC, positiveABC}` and
voltage is positive minus negative. Sign +1 means the FMU reports passive
current into the positive terminal; -1 means it reports injection into the bus.
Raw FMI output attributes preserve their original sign.

The model provides `x`, `v_intf`, `i_intf`, `fmi.inputs`, `fmi.outputs`,
`fmi.input.<name>`, `fmi.output.<name>`, `fmi.state.<name>`,
`newton_iterations`, `newton_residual`, `fmi_internal_steps`, and
`fmi_coupling_error_A` Attributes. The last value is zero for converged
single-rate steps and records the largest phase-current discrepancy of an
explicit multirate step. Non-voltage inputs also have `control.<name>`
Attributes. A control may be driven with
`setInputFunction(name,function,discontinuityTimes)`; explicit event times must
align to the timestep. Input functions must match initial FMI input values.

During DPsim initialization, `NonlinearSSN` first validates its connected
terminals and then uses their initial voltage difference as the FMU's initial
electrical input whenever it is a finite three-phase value. This automatically
aligns the imported FMU with a native network initial condition instead of
always forcing the voltage declared in `modelDescription.xml`. The FMI initial
continuous state is still retained; this is therefore a PF-like terminal
initialization, not a general hidden-controller steady-state solve. FMUs whose
controller, PLL, or integrator states cannot be set through FMI may still need
a short settling interval or an exporter-specific operating point.

Each component owns an independent runtime. Sharing one FMU instance between
clones is not supported. The exporter's singleton capability is respected;
load a separate instance only if that capability allows it.

## Build and run

Build with `WITH_FMI=ON` and `DPSIM_BUILD_EXAMPLES=ON`. The reusable nonlinear
examples are `FMI_EMT_Ph3_Nonlinear` and `FMI_EMT_Ph3_TwoGenerators`; their
source and sample FMUs are under `dpsim-fmi/examples/`. Generated logs belong
in the configured build directory. See [the FMI examples guide](examples/README.md)
for the available targets and inputs.

## Regression coverage and supported subset

The FMI CTest fixtures cover FMI 2.0 and FMI 3.0 scalar Model Exchange, linear
scalar/ABC models, nonlinear coupled solves, sampled controls, rollback and
single-pass scheduling, singleton handling, malformed imports, nonconvergence,
and unsupported event-current changes. Native SSN solver changes are opt-in through
`MNANonlinearCompInterface`; existing native components retain their previous
path.

Supported FMI 2.0/3.0 boundary/state variables are scalar Real/Float64. Scheduled events must
advance time. Crossed events use automatic substeps with either reproducible
FMI rollback or the explicitly selected single-pass policy. Explicit external input discontinuities must still align with the
network grid. State-indicator endpoint handling and continuous-state resets
require explicit opt-in; root localization between reported internal substeps is
not implemented. In the single-rate path, an event changing terminal current
beyond tolerance is rejected because instantaneous algebraic network re-solving
is unsupported. Clocks, FMI arrays, structural parameters and unsupported
execution-tool dependencies are rejected. Runtime parameter editing, generic
automatic steady-state initialization and dynamic-phasor conversion are not
implemented. A model's singleton declaration is enforced; two-generator
validation uses two genuine exports with distinct identities rather than editing
capability metadata.

`NewtonSettings::maxJacobianReuseSteps` optionally reuses a local Jacobian as a chord predictor. It defaults to zero. Actual state and current residuals still govern convergence; failure refreshes the Jacobian before the next iteration.
