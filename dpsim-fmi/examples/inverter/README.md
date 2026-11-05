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

## 50 µs gate runs

Loaded-start, balanced weak-grid sag, and unbalanced-sag cases were run for
20 ms at 50 µs with 50 internal 1 µs FMI steps per macro step. The compact
full-resolution metrics are in
`build-fmi/fmi-validation/gate3/gate3_metrics.json`; case logs are alongside
that file. The sag intervals are 5–10 ms. The table reports sag-window values
for the disturbances and 10–20 ms recovery values for all cases.

| Case | PCC voltage RMS [V] | Converter current RMS [A] | Mean P [W] | Mean Q [var] | Peak reported lag error [A] |
|---|---:|---:|---:|---:|---:|
| Loaded-start, 10–20 ms | 14,431 | 1.746 | 75,557 | 2,414 | 0.047 after 1 ms |
| Balanced weak-grid sag, 5–10 ms | 10,105 | 1.238 | 37,446 | 2,067 | 0.281 |
| Unbalanced sag, 5–10 ms | 9,749 | 1.185 | 34,397 | 1,847 | 0.406 |

The stored Gate 3 logs predate the consistent initialization below; their first
sample has zero line and converter currents and is an initialization artifact.
The full-window coupling maximum is 1.169 A; after 1 ms it is at most 0.047 A
in loaded-start. The residual field is a one-macro-step implicit mismatch in
`single_pass_explicit` mode, not a Newton convergence metric; each macro step
performs one network pass. Separately, 1 µs MATLAB/Simulink physical
references were regenerated for matching 0.2 s balanced-load, unbalanced-load,
and voltage-sag cases. These long cases close the converter breaker at 150 ms
and apply their load/grid events at 170 ms; they are distinct from the short
20 ms Gate 3 scenarios. Their aligned comparison metrics are in
`build-fmi/fmi-validation/gate6-weak-grid/paper/matlab_load_step_metrics.csv`
and `matlab_voltage_sag_metrics.csv`. In the 170–199.9 ms load-step windows,
PCC phase-voltage maximum error is at most 0.94 V (0.0046% of phase-peak base),
converter phase-current maximum error is at most 0.030 A (0.18% of rated-current
base), the added load-current maximum error is below 0.08 mA, and calculated
instantaneous P/Q maximum errors are below 0.87 kW/0.91 kvar (0.182% of the
500 kVA base); Simulink converter current for P/Q is reconstructed by PCC KCL.
In the
170–190 ms voltage-sag window, phase-A PCC-voltage maximum/RMS error is
0.724/0.113 V and converter-current error is 0.0140/0.00669 A. The reference
model and DPsim cases use coincident 1 µs timestamps; DPsim CSV logging is at
100 µs intervals. Startup is reported separately: before the 150 ms breaker
closure, phase-A voltage/current errors stay below 2.1 µV/0.53 µA; from 150 to
160 ms they stay below 0.49 V/0.0103 A. Thus the zero-valued initial samples
are the unenergized network state, not evidence of an FMI coupling failure.

### Consistent loaded-start initialization

`sources/consistent_loaded_start_50us.json` exercises the loaded initialization
explicitly. It converts the three-phase node convention to phase voltages,
computes the passive R-L/load operating point, and initializes the native load,
line histories, and FMU terminal voltage before the logger records t=0. The
initial PCC KCL residual is 0 A; the maximum phase voltage drop from the grid
to the PCC is 10.10 V. In the earlier initialization, the t=0 logged branch
currents were zero despite a 1.154 A phase-current magnitude implied by the
PCC voltage and load. After initialization, the first-step explicit coupling
lag peaks at 1.169 A and falls below 0.0468 A after 1 ms. The initial KCL
consistency and the first-step explicit coupling lag are separate quantities.
The FMU continuous states remain the exported FMI start values; this is a
phasor-consistent network initialization, not a periodic steady-state solve
for the converter controller.
Full-resolution reruns of all three short cases are preserved under
`build-fmi/fmi-validation/consistent-initialization/final-reruns/`. Recreate
them and refresh the paper metrics with
`python3 dpsim-fmi/examples/inverter/run_paper_validation.py --reuse --rerun-consistent-start-cases`.
The measured windows give loaded-start PCC voltage/current RMS of 14.431 kV /
1.746 A; balanced sag 10.105 kV / 1.238 A; and unbalanced sag 9.749 kV /
1.185 A. They are 50 µs explicit-coupling results, with one network pass and
50 FMI internal steps per macro step. These are converter-case metrics, not
Newton-converged results.

To reproduce the full-resolution log after building `FMI_EMT_Ph3_Mixed`, run
it with `sources/consistent_loaded_start_50us.json` from
`build-fmi/fmi-validation/consistent-initialization/`. The result is written
to `logs/fmi_mixed.csv` for the notebook comparison.

## Bounded paper validation

Open [Paper_FMI_EMT_Results.ipynb](Paper_FMI_EMT_Results.ipynb) for the
reproducible result tables and figures, including the per-case multirate
internal-step and coupling-error logs.

Run `matlab/validation/regenerate_physical_references.m` in MATLAB to regenerate
the three 0.2 s physical traces under `build-fmi/fmi-validation/`. Then
`run_paper_validation.py --reuse` regenerates the comparison metrics and
figures from those traces and the existing DPsim logs. Without the MATLAB CSVs,
the script still generates the numerical-fixture and timestep results and
skips the independent physical-reference metrics.
Those Simulink traces cover the established long load-step and 170–190 ms
voltage-sag cases. They do not provide an independent physical reference for
the short 5–10 ms balanced and unbalanced sag cases above; those use DPsim
time-step refinement as numerical evidence only.

Running `run_paper_validation.py` without `--reuse` reruns the weak-grid sag at
1, 25, and 50 µs. It generates the comparison metrics and figures under
`build-fmi/fmi-validation/gate6-weak-grid/paper/`. The 1 µs case is the same
DPsim/FMU/network at the fine timestep, so this quantifies timestep and
single-pass sensitivity; it is not an independent physical reference. The
1 µs case executes ordinary nonlinear Newton at each step. The 25 and 50 µs
cases use the explicit event-aligned path with 25 or 50 internal 1 µs steps.

For the 5–10 ms sag window, maximum/RMS errors against the 1 µs run are 3.89
kV/283 V (PCC phase A) and 0.101 A/0.0108 A (converter phase A) at 25 µs;
they are 4.35 kV/442 V and 0.197 A/0.0230 A at 50 µs. Both maximum PCC
voltage errors occur at the 5 ms sag-onset sample. The peak coupling lag after
1 ms is 0.194 A and 0.281 A. These event transients remain visible in the
generated error plot and should not be interpreted as converged Newton results.

The importer now reports cumulative per-step `fmi_evaluation_requests` and
`fmi_evaluation_wall_time_s`. These count public importer evaluation requests
(including cache hits) and time spent in those requests; they exclude FMI event
handling and total simulator overhead. The 20 ms weak-grid runs used 71,973
requests and 4.51 s evaluation time at 1 µs, versus 50,082 requests and 4.30 /
4.33 s at 25 / 50 µs. Measured process times were 7.10, 6.30, and 6.33 s,
respectively.
The explicit cases report a large normalized `newton_residual` because that
field is the mismatch to the implicit residual while the network is intentionally
solved once; it is not comparable to a converged Newton residual.

The analytical fixture's base and 10× tighter Newton tolerance runs both
requested 4,000 FMI evaluations across 100 steps, converged in at most two
iterations, and remained within 2e-9 of the independent trapezoidal state
recurrence. Maximum relative Norton conductance/history errors were 4.72e-11 /
3.06e-11 at baseline and 5.89e-11 / 6.83e-11 when tightened. Evaluation time
was 5.11 / 5.07 ms. Tightening these already strict tolerances did not materially
change the trajectory or work count for this fixture.

### Finite-impedance Newton versus frozen-current baseline

`FMI_FiniteImpedanceCouplingTest` and `run_paper_validation.py --reuse` provide
a short, reproducible R-L comparison using the same nonlinear sampled FMI test
model, zero-current start, 50 Hz source, and grid (5 Ω and 5 mH per phase).
The 1.5 ms macro-step Newton run iterates the MNA/Norton coupling through
event-aligned FMI substeps. The explicit baseline solves MNA once per macro
step using the last accepted FMU terminal current, then advances the FMU over
the solved PCC-voltage interval. A 50 µs Newton run is the fine-step DPsim
reference; it is not an independent physical reference. This fixture isolates
the coupling policy and does not represent the switched inverter.

Metrics and plots are written to
`build-fmi/fmi-validation/gate6-weak-grid/paper/finite_impedance_metrics.{json,csv}`,
`finite_impedance_waveforms.png`, and `finite_impedance_errors.png`. The table
reports both explicit-versus-Newton error at the same macro step and each
macro-step result versus the 50 µs reference. `fmi_coupling_error_A` is the
explicit method's lag diagnostic; the explicit run's `newton_residual` is not a
Newton convergence measure. Newton converges with at most 7 iterations and a
reported terminal-current coupling lag below 2 µA. The frozen-current run has
a 0.657 A peak lag; versus Newton at the same 1.5 ms step, its post-start
maximum phase errors are 4.645 V PCC voltage (4.206 V RMS error) and 0.1723 A
converter current (0.1254 A RMS error). Against the 50 µs DPsim Newton trace,
the explicit peak errors rise to 10.604 V and 0.4118 A because that comparison
also includes macro-step discretization error. The 50 µs trace is a fine-step
DPsim reference, not an independent physical reference.

MATLAB/Simulink reference regeneration is available, although MATLAB remains
alive after writing each output and was stopped after validating the complete
CSV files. MATLAB R2025a exposes `exportToFMU`; a scoped FMI 3.0 Model Exchange
export of an averaged three-phase R-L output filter with sampled converter
voltage commands and continuous current states reached code generation, then
failed because Simulink Coder is not installed. The exporter also rejects
multi-instance FMI 3.0 Model Exchange export, so that option had to be disabled
before reaching the missing-toolbox failure. No deterministic averaged-GFL FMU
or independent second converter FMU was produced. The switched FMU's replay
failure and serialized-state abort remain the external limitation. Keep
multirate Newton out of the main paper claim/title; retain the fixture
validation and identify realistic converter operation as explicit
event-aligned coupling.
