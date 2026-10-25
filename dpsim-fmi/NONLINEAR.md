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
load a separate instance only if that capability allows it. The SG MATLAB FMU
permits one instance per process, so run different scenarios in separate
processes.

## Build and reproduce

Work from the repository root on the SSD. Example sources and retained inputs are under `dpsim-fmi/examples/`; generated outputs belong under build.

```sh
env DEVELOPER_DIR=/Library/Developer/CommandLineTools \
  cmake --build build-macos-check --parallel 3 --target \
  FMI_EMT_LinearSSN FMI_EMT_Ph3_LinearSSN FMI_EMT_Ph3_Nonlinear \
  FMI_EMT_Ph3_Mixed FMI_EMT_Ph3_TwoGenerators
ctest --test-dir build-macos-check -R fmi_ --output-on-failure
```

Execute each case's comparison notebook with the DPsim Python kernel. Each notebook contains its comparison logic and runs the C++ example against retained MATLAB references. The inverter notebook also reproduces direct 50/25 µs capability checks. Open the retained physical reference models in MATLAB to regenerate converter results; they do not invoke the FMU. The generator reference functions independently solve the full model equations with implicit trapezoid and segmented ode15s.

| Case directory | Validated scope |
| --- | --- |
| `linear_scalar` | Existing scalar linear FMI import |
| `linear_three_phase` | Existing ABC linear FMI import |
| `synchronous_generator` | Nonlinear coupled generator, matched full-state MATLAB reference, 50/25 µs refinement and ode15s cross-check |
| `inverter` | Switched converter with native voltage source, R/L, base and switched loads; balanced/unbalanced load, sag, phase and frequency disturbances at supported 1 µs step |
| `two_generators` | Two separately exported generators on one shared bus, independent control changes, model-order isolation, singleton rejection and meaningful nonzero operating point, 50/25 µs refinement |

Each directory keeps one C++ example and one executed comparison notebook, final FMU/XML/source models in `sources/`, and meaningful MATLAB reference results/model definitions in `matlab/`. All work remains local, uncommitted and unpushed on `feature/FMU_import`. MATLAB/Simulink is the benchmark scope.

## Numerical evidence and interpretation

The generator comparison checks all 34 saved signals. The original 0.1% maximum/RMS engineering-base criterion passes; the largest normalized discrepancy against ode15s is approximately 0.0671% at 50 µs and 0.0168% at 25 µs. The two-generator comparison checks 52 physical/control signals; its largest matched normalized discrepancy is approximately 5.23e-7, and the ode15s cross-check approximately 3.44e-6 (fractions of declared bases). Refinement reduces representative RMS errors by about four. Model-supplied initial conditions have nonzero power and network currents, with stationary machine/controller residual about 1.14e-13. This is not arbitrary-FMU automatic load-flow initialization.

The switched converter preserves the IGBT/diode plant, transformer/filter, original loads, PLL, PI controllers, PWM and internal breaker. The independent physical Simscape reference preserves the original grid source and its initial state. A stateless series correction applies grid disturbances; additional reference resistors match the native switch resistances. Both controller subsystem XMLs are byte-identical to the original reference. Metrics include every coincident sample and retain startup/event peaks, raw angle maxima and wrapped angle errors. No trace shifts are applied.

For the five matched 1 µs disturbance runs, the worst converter phase-current peaks are approximately 0.182%, 0.163%, 0.139%, 0.121% and 0.106% of the declared 16.32993 A engineering base. The largest PI/control peaks are approximately 0.834% of their declared pu base. These exceed the original provisional 0.1% threshold; the notebook shows the individual failures explicitly. The user accepted a wider margin, without a new universal numerical guarantee. All five runs satisfy network KCL to logged precision, converged Newton residuals and unchanged physical/control CSV outputs under tenfold tighter Newton tolerances.

The requested practical network step is 50 µs. The switched converter now
runs through 0.2 s at 50 and 25 µs while retaining all 1 µs scheduled events.
Its MathWorks R2025a Model Exchange export advertises state rollback, but an
identical replay changes the continuous state norm by about 0.394 and output
norm by about 9.38e-4. The importer detects and rejects that rollback path.
These two runs therefore use the generic single-pass policy. After the first
10 ms, the 50 µs run has maximum converter-current and PCC-voltage errors of
about 0.095% and 0.028% of their declared bases against the 1 µs import. The
largest controller difference is 0.556% (`pi_q`). At 25 µs these become about
0.086%, 0.020%, and 0.295%. Full-window results retain a very large initial PCC
voltage difference caused by the zero-current line initialization at the coarse
network step. The direct coupling discrepancy peaks at 1.268 A at 50 µs and
0.683 A at 25 µs during startup; after 10 ms the maxima are 0.0523 A and
0.0258 A. These are measured explicit-coupling results, not converged implicit
coupling claims.

The new multirate path implements that rollback/replay interval map. Its Norton
tangent linearizes terminal output current at the converged substep state; the
outer residual is checked against a fresh complete replay. This is a partitioned
multirate approximation rather than a monolithic solve of every internal state
and every network node. Its accuracy must therefore be measured for each model,
especially around ideal switching and network discontinuities.

## Regression coverage and supported subset

The current FMI CTest suite has 27 passing checks, including FMI 2.0 and FMI
3.0 scalar Model Exchange, linear scalar/ABC, nonlinear coupled solves, sampled
controls, rollback multirate scheduling, explicit single-pass scheduling without
rollback, singleton handling, malformed imports, explicit nonconvergence and
unsupported event-current changes. Legacy example notebooks also execute
successfully. Native SSN solver changes are opt-in through
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

## Official PV example limitation

The official MathWorks `PVArrayGridAverageModelExample` was investigated separately. Its original Model Exchange export has zero continuous states and 50 µs internal ticks: successful lifecycle replay does not establish continuous, tightly coupled EMT MNA import. Provisional getter calls altered its later history in the isolated protocol check, even across advertised state restore. A continuous-powergui variant failed MATLAB ME export on the non-inlined `sfun_spid_contc` S-function, and generated network code overwrote an attempted ordinary State-Space replacement. The network also contains switching; one frozen network matrix would remove meaningful behavior. Consequently this PV/MPPT model is not claimed as validated. Development evidence remains outside the examples under the SSD benchmark/build directories.

The successful switched inverter example is a different official MathWorks model and does not contain an MPPT loop. Neither example demonstrates universal FMI support. The implemented lifecycle follows the [FMI 3 specification](https://fmi-standard.org/docs/3.0.2/); MATLAB documents the [non-inlined S-function export restriction](https://www.mathworks.com/help/fmuexport/ug/generate-model-exchange-FMU.html).

`NewtonSettings::maxJacobianReuseSteps` optionally reuses a local Jacobian as a chord predictor. It defaults to zero. Actual state and current residuals still govern convergence; failure refreshes the Jacobian before the next iteration.
