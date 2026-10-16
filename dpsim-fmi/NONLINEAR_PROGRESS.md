# Nonlinear FMI work — larger coupled-converter stage unfinished

## Latest verified checkpoint (2026-10-16)

- All **21 CTests pass** after the final runtime changes, including sampled
  control timing, coincident command/tick, snapshot isolation, skipped optional
  completed-step callback, current-jump rejection and legacy linear imports.
- Both primary comparison notebooks execute successfully again:
  SG coupled EMT3ph MNA vs independent MATLAB trapezoid/ode15s; nonlinear MSL
  magnetic circuit vs native Modelica. Executed notebooks, standalone HTML and
  plots saved under outputs/fmi-nonlinear/sg/comparison and modelica.
- Official fully-discrete PV standalone event replay now matches all120001
  native samples to2.2e-12A and FMPy to1.2e-12A. New notebook/HTML under
  outputs/fmi-nonlinear/pv, source examples/FMI_PV_Event_Compatibility.ipynb.
- This is **not** a successful coupled PV MNA import. The full overnight goal
  remains unfinished; do not mark it complete.
- Continuous SPS network export rejects its non-inlined S-function. Replacement
  experiment confirms model update restores the S-function. Found vendor source
  sfun_spid_contc.c in installed toolbox (2631lines); it has switched descriptor
  dynamics and nSwitch+1 zero crossings, not one fixed linear network.
  S.nb={x65,y14,s1,u10}; power_analyze gives55 independent and10 dependent states.
  Any future exporter workaround must preserve switching and descriptor/state
  transformations. Do not freeze a single A/B/C/D matrix; the experiment script
  now asserts no switches before attempting such replacement.
- FMI state snapshots only wrap hybrid probes when the capability is supplied.
  Their semantics must be validated for the exporter. Community OM restore
  fails for the continuous magnetic FMU, so continuous-only models use explicit
  state/input assignment and no snapshots. Original PV restore is incomplete for
  pre-event getters; that model is not safe for Newton probing as exported.
- Next required work: select an actually exportable continuous/hybrid converter
  example or implement a faithful exporter workaround with switching preserved,
  then create real voltage/current ports and a coupled DPsim/reference comparison.
  A closed fully-discrete event replay cannot satisfy that requirement.
- Working branch remains feature/FMU_import, local, uncommitted and unpushed.

Work directory: `/Volumes/SSDMatthias/development/dpsim`, local branch
`feature/FMU_import`; no commits/pushes requested. New scripts, results and
explicit temporary files stay under `/Volumes/SSDMatthias/development`.

Latest scope (supersedes the original goal tool's MATLAB-only wording): the user
explicitly re-authorized Internet downloads/clones and Modelica examples/tools
under development. Complete in order: controlled SG, more complex bundled
Simulink example (prefer GFL/multiple loops), then Modelica portability.
No fixed time/token budget. Partial results do not complete the goal.

## Starting point

- Existing linear scalar/abc FMI import and 13 CTests previously pass.
- SG: ten-state nonlinear Model Exchange FMU plus six-state external network.
- Inputs/assets copied to `outputs/fmi-nonlinear/sg/matlab-origin`.
- Existing nonlinear references use ode15s, not trapezoidal. They include valid
  initialization, balance and energy checks; nonlinear transient import does not
  yet exist.
- Xcode 16.4 on SSD is used by MATLAB. DPsim build uses Command Line Tools via
  `env DEVELOPER_DIR=/Library/Developer/CommandLineTools` (do not change global
  xcode-select).
- MATLAB executable `/Volumes/SSDMatthias/MATLAB_R2025a.app/bin/matlab`.
- Python `dpsim-python/bin/python`; formatter paths are in session notes.

## Architecture selected

- Persistent generic FMI Model Exchange runtime (continuous states + scalar
  Float64 inputs/outputs), shared archive/XML/native library support with linear
  extraction. Strong rejection of unsupported events/clocks/discrete behavior.
- `NonlinearSSN`: implicit trapezoidal residual; finite-difference Jacobians;
  eliminate local state increments into a 3x3 Norton Jacobian and history.
- New opt-in `MNANonlinearCompInterface` allows the variable-matrix MNA solver
  to check the resulting state/current residuals and repeat restamping/solving.
  Existing SSN and legacy synchronous-generator iterations remain unchanged.
- Accept states and call FMI completed-integrator-step only after convergence;
  expose inputs, states, outputs and iteration/residual diagnostics as Attributes.
- External benchmark network is a general six-state current-input/voltage-output
  SSN branch with grid-voltage forcing and scheduled load conductance. This
  avoids artificial capacitor-history errors when a load changes at a boundary.
- Handle disturbances with exact interval splitting: endpoint before the event
  uses the preceding input, then rebase derivatives for the following interval.

## Acceptance

Independent MATLAB trapezoidal implementation with matched parameters, states,
event timing, h and converged Newton iterations. Compare against DPsim and
ode15s; verify tolerance tightening, timestep refinement and regressions.
Provisional maximum and RMS normalized errors <0.1% of declared engineering
bases, with absolute errors also reported. Produce executed notebooks/plots and
reproduction commands. No physical-parameter fitting or shifted traces.

## Next work

Implement runtime/API, opt-in MNA Newton hook, nonlinear SSN, generator network
example and MATLAB trapezoidal reference; build/test and validate before moving
to the larger exporter benchmarks.

## Checkpoint: controlled generator passed (2026-10-15 evening)

Implemented:

- ModelExchange.h/.cpp shared persistent FMI runtime; LinearModel now uses it.
- NonlinearSSN.h/.cpp, full implicit trapezoidal condensed Newton, explicit
  current sign and variable mappings, FMI/control/state/output Attributes.
- Opt-in MNANonlinearCompInterface and 55-line aggregate core changes in
  MnaSolver/MnaSolverDirect variable-matrix path. Existing SSN untouched.
- LinearNetwork3Ph helper + FMI_EMT_Ph3_Nonlinear benchmark executable.
- MATLAB create_sg_trapezoidal_reference + sg_full_rhs_reference (full 16-state
  Newton with complex-step Jacobian). 50/25 us references finish in <1 minute.
- compare_nonlinear_sg.py checks all 34 signals, KCL, convergence/refinement.
- New synthetic cubic-coupled FMI fixture and 2-instance integration test;
  explicit unconverged-step rejection test. All 15 CTests passed.
- Native EMT_Ph3_SSN_GFM_JacobianCheck and induction-motor load-switch example
  build and run successfully. Logs under outputs/fmi-nonlinear/native-regressions.

Results at outputs/fmi-nonlinear/sg:

- matlab/: matched 50/25 us references and benchmark.json with engineering bases.
- dpsim-50us, dpsim-25us, dpsim-50us-tight/: completed actual FMU/MNA runs.
- comparison.json: all 34 signals pass max/RMS normalized error <0.1%.
- 50 us matched currents max ~1.37e-6 A, bus voltage ~1.92e-6 V; two Newton
  iterations per accepted step; max scaled residual .002261; logged KCL zero.
- 25 us matched phase-A current ~5.50e-7 A, bus voltage ~1.08e-6 V.
- Tightening Newton tolerances 10x leaves six-decimal CSV identical.
- Against ode15s, 50 us worst normalized peak .0006713 (load power); 25 us
  .00016785. Phase-A bus voltage .170539 V -> .042637 V (factor ~4).
- Source notebook examples/FMI_EMT_Ph3_Nonlinear_MATLAB_Comparison.ipynb.
- Executed notebook, standalone HTML, metrics.csv, controls/waveforms/errors/
  refinement PNGs in outputs/fmi-nonlinear/sg/comparison. Visually inspected
  controls plot; all references overlay correctly.

Important runtime discovery: MATLAB FMU requests completed-step events despite
zero event indicators. Runtime now handles event iterations that leave
continuous states unchanged and refreshes nominals when requested; actual state
reinitialization/time/state events/clocks remain unsupported. inputEvent also
uses this common event path.

Last source edits after those passes (require rebuild/recheck before final):
NonlinearSSN supports n=0 in condensation and throws if post-step is called
without the opt-in solver convergence check. No numerical behavior change to SG.
Formatting of those last few lines still needed.

## Stage 2 discovery in progress

Downloaded official MathWorks example via openExample with explicit workDir:
`/Volumes/SSDMatthias/development/fmi-external-benchmarks/matlab/PVArrayGridAverageModelExample/PVArrayGridAverageModel.slx`
Official ID: simscapeelectricalsps/PVArrayGridAverageModelExample. Native model
is a 100-kW grid-connected PV average model, MPPT + boost + PLL/DC-link/dq VSC
controls; root contains 'VSC Control ' and '600V Utlity Grid', a 100kVA
260V/25kV transformer, load and powergui. It is Fixed-step/FixedStepDiscrete,
stopTime 6s. Do not silently convert away the discrete control behavior.
Inspection batch outputs/fmi-nonlinear/pv-inspection.log: diagram update passed,
but getInitialState fails because model is configured accelerator. Retry with
SimulationMode='normal' (numerical execution mode only) before getInitialState.
Batch script is outputs/fmi-nonlinear/inspect_pv_example.m; downloaded example
supports MATLAB license/products installed (39 products incl Simscape Electrical).

MATLAB model examples generally downloaded on demand; installed toolbox demo
folders don't contain most examples. `openExample(...,'workDir',SSD-path)` is
supported and keeps downloads inside development. Installed catalog JSON can
be parsed rather than grepped: help/docCatalog/topexample/sps.json and
cshapi_helptarget/sps.json. Avoid raw grep of minified JSON (huge output).

### Modelica portability stage

- User explicitly authorized downloads/clones from the Internet under SSD development.
- Cloned MSL v4.0.0 and community native arm64 OpenModelica under
  `/Volumes/SSDMatthias/development/fmi-external-benchmarks`; exact commits in
  `examples/modelica/README.md`.
- Native compiler successfully exports FMI **3.0** Model Exchange; no FMI2 adapter needed.
- Port adaptation of `ShowSaturatingInductor`, three phases, six imported states,
  successfully runs in the same `NonlinearSSN` with no magnetic equations in C++.
- Independent reference: original MSL example instances, fixed-step RK4 100 us.
- Peak current error 58.1 uA at 10 ms,14.9 uA at 5 ms; RMS ratio3.99.
  Fine runs~0.52uA hit logger quantization. KCL and convergence pass.
- New generic prescribed-voltage benchmark driver uses explicit integer number
  of Simulation steps to avoid accumulated-time rounding omitting the last sample.
- New build script explicitly loads local Complex/ModelicaServices/MSL dependencies
  with uses=false to avoid online package auto-download attempts.
- PV average source is discrete with extensive PLL/PI/MPPT states; probing export
  faithfully rather than replacing sampled control behavior.

### Sampled-control support (in progress)

- Runtime now retains default tunable parameters and supports opaque internal
  discrete variables. Boundary inputs/outputs remain continuous scalar Float64.
- `nextTimeEvent()` exposes FMI scheduled time events; strictly advancing,
  grid-aligned time events execute only after accepted endpoints.
- Coincident external reference changes and internal ticks share one event-mode
  update with right-side inputs, avoiding a one-period control delay.
- Added nonlinear sampled-plant independent recurrence and noncommensurate-step
  tests; these passed. Added sampled-LTI rejection test, mapping correction made,
  awaiting rerun.
- Actual official PV FMU initializes with zero continuous states, no inputs,
  three current outputs and next internal event50us. Standalone event-protocol
  replay underway; explicitly not claimed as MNA coupling.
- A separate continuous-powergui PV variant is being exported with all original
  sampled control algorithms retained. It is a documented numerical plant variant,
  not the original fully discrete model. Original model files unchanged.
- Original PV reference rerun with LimitDataPoints off to avoid MATLAB's default
  1000-sample output truncation.

### Continued 2026-10-16: real PV event investigation

- Original PV reference now saved over all120001samples,0..6s,50us.
- Installed FMPy0.3.32+lxml6.1.3 into separate SSD
  `/Volumes/SSDMatthias/development/fmi-external-benchmarks/python-reference`
  without changing DPsim Python packages. FMPy native ME/Euler replay exactly
  matches Simulink PV reference to~2.2e-12 A (no continuous states).
- Our old discrete diagnostic differed up to0.57A. Protocol isolation under
  `outputs/fmi-nonlinear/pv/protocol_variants.py` demonstrates that **GetFloat64
  before a time event** changes this particular export's subsequent trajectory.
  Normal/no-pre-get, repeated SetTime after event, empty SetFloat64 and snapped
  time variants match reference. Pre-get differs0.50A by0.05s.
- Even FMPy's getFMUState/setFMUState around a pre-get does NOT restore this PV
  trajectory (0.44A discrepancy by0.05s). Do not claim that snapshots make this
  exported SPS discrete model safe for provisional MNA/Newton evaluation.
- Runtime now has prepareStep (no getters), acceptedValues cache, proper
  needsCompletedIntegratorStep handling (defaultfalse), optional FMI snapshots
  for hybrid probes where capability is present. Continuous-only models skip
  snapshots because community OpenModelica advertises state snapshots but its
  restore returnederror3 during magnetic-circuit probes. This must remain
  transparent; ordinary continuous-state assignment remains the supported path.
- FMI_DiscreteProbe changed to prepareStep/completeStep/acceptedValues, which
  should match standard/FMPy zero-state event execution. Rebuild completed,
  full6s replay and primary notebooks still need rerun for newest changes.
- Constructor/commit reads accepted output/derivative once and caches them;
  old-endpoint evaluations reuse that cache, preventing repeated accepted reads.
- Continuous-powergui PV variant FAILED export: non-inlined `sfun_spid_contc`
  State-Space under powergui unsupported by MATLAB Model Exchange exporter.
- power_analyze('model','structure') gives exact55-state linear network
  A55x55,B55x10,C14x55,D14x10,x0; switch matrices also55x55. Native source utility
  is a whole feeder + transformers (not just ideal source), so do NOT replace
  it with a guessed simple R/L source.
- Current experiment `outputs/fmi-nonlinear/probe_replace_sps.m`, session2346,
  makes a new copy in external-benchmarks/matlab/PVContinuousStateSpaceProbe,
  clears only copied powergui mask/init callbacks, replaces generated S-function
  by built-in State-Space with exact power_analyze matrices, then checks whether
  update/simulation preserves that replacement and attempts export. May fail or
  be overwritten; inspect `pv-replace-sps.log`. Original files unchanged.
- Before latest prepareStep change: all19CTests passed; SG notebook passed;
  Modelica notebook failed only due unnecessary OpenModelica snapshot restore.
  Latest change restricts snapshots to hybrid models and should restore Modelica;
  must rerun both notebooks.
