# Inverter

Build the C++ example with CMake in the repository build directory. `sources/` contains the final FMU, extracted `modelDescription.xml`, and Simulink model. `matlab/` contains saved reference results and model definitions where needed. The comparison notebook contains its comparison logic and writes generated runs and plots under `build-macos-check/fmi-examples/inverter/`.

## Multirate timestep support

The requested practical network timestep is 50 µs (`5e-5` s). The switched converter export schedules FMI events every 1 µs; its controller parameter file specifies `Ts = 1e-6`, a 6 kHz switching frequency, and `Tsc = 5e-4`. These are different from the network integration step.

`NonlinearSSN` now has a generic rollback/replay scheduler that lands implicit
substeps on every reported internal event. `fmi_internal_steps` records the
executed count. An independent sampled FMI fixture passes with a network step
twice its internal event period.

This R2025a switched Model Exchange export advertises FMI-state rollback, but an
identical 50 µs replay changes its final continuous state by about 0.394 and its
outputs by about 9.38e-4. The importer detects and rejects that behavior. Trying
the separately advertised serialized state aborts inside the exported binary,
so it is not used.

For such FMUs, `MultirateMode::SinglePassExplicit` is a generic opt-in policy.
It solves the native network with the previous accepted terminal current,
advances the FMU once through all internal events using the solved voltage
trajectory, and logs `fmi_coupling_error_A`. It does not contain converter names,
controller periods, or equations. The default remains rollback/Newton because
that retains converged implicit coupling when state restoration is valid.

Executed 0.2 s runs at 50 and 25 µs cross the internal breaker closure while
executing respectively 50 and 25 one-microsecond internal steps per network
step. Against the validated 1 µs import, maximum errors after 10 ms are:

| network step | converter current | PCC voltage | largest control error |
|---|---:|---:|---:|
| 50 µs | 0.095% of base | 0.028% of base | 0.556% (`pi_q`) |
| 25 µs | 0.086% of base | 0.020% of base | 0.295% (`pi_q`) |

The explicit coupling-current discrepancy peaks during startup at 1.268 A and
0.683 A, then stays below 0.0523 A and 0.0258 A after 10 ms. The full-window
table retains large first-sample voltage errors caused by applying a coarse
trapezoidal step to the zero-current native line initialization.

For a practical loaded start, `sources/loaded_start_50us.json` enables the
optional line-history initialization. It retains the FMU's initialized
terminal current and adds the native load current at the initial PCC voltage;
this is automatically enabled whenever `initial_voltage_scale` is nonzero.
The zero-voltage reference configuration therefore retains its zero-current
startup behavior.

`sources/weak_grid_sag_50us.json` provides a reproducible weak-grid stress
case: doubled series impedance and a 30% three-phase voltage sag from 5 to
10 ms, using the same loaded-start and 50 µs single-pass settings.

Per-phase sag factors are also supported through `grid_sag_scale_abc`; the
`sources/unbalanced_sag_50us.json` case exercises a 50/20/35% unbalanced sag
on the same weak grid.

The longer `sources/long_weak_grid_50us.json` case extends this disturbance to
0.2 s for sustained controller and recovery validation.

For long runs, `log_downsampling` reduces CSV output without changing the
50 µs simulation step. The loaded-start configuration initializes the native
line history from the load and converter current; this is the preferred
starting point for long-term inverter runs.

Open [`Long_WeakGrid_50us_Results.ipynb`](Long_WeakGrid_50us_Results.ipynb)
to plot the generated CSV and inspect voltages, currents, coupling error,
internal FMI steps, controller outputs, and PLL recovery.

Open [`TwoConverters_50us_Results.ipynb`](TwoConverters_50us_Results.ipynb) to
inspect the two identical FMU instances coupled through one native DPsim
three-phase network, including shared PCC voltage, per-converter currents,
residuals, and internal-step activity.

Open [`TenSecond_FMI_Results.ipynb`](TenSecond_FMI_Results.ipynb) for voltage,
current, residual, and internal-step plots from the completed 10-second FMI
examples.

Sub-microsecond agreement is not an acceptance requirement. Loaded-start
validation remains the next accuracy check; existing finer MATLAB runs provide
reference uncertainty information. The executed notebook, plots, metrics,
configs and logs are under `build-macos-check/fmi-examples/inverter/`.
