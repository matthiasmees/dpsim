# IEEE 39 bus inverter comparison

This case compares one FMI 3 averaged inverter model in DPsim Model Exchange (ME), DPsim Co-Simulation (CS), and native Simulink. The same standalone FMU exposes both ME and CS from the same `IEEE39AvgInverter.cpp` source. The Simulink replay model independently implements the same current-control equations.

## Grid and event

The network uses the hardcoded three-phase IEEE 39 bus data already used by the DPsim `IEEE39_Ph3` example: 39 buses, 46 branches, 10 synchronous generators, 60 Hz. The inverter is placed at bus 30. To preserve the standard pre-event active and reactive power at bus 30, generator 30 remains at half of its published (P/Q) operating point and the FMI inverter supplies the other half: 125 MW and 80.881 Mvar. All other generators, loads, and lines remain unchanged. A balanced 100 MW load is switched onto bus 39 at 170 ms.

Runs use a 50 µs EMT network step for 0.2 s (4,000 network steps). The average converter uses an instantaneous (P/Q) current reference in the stationary αβ frame, a 5 ms first-order current response, and a 1.2 pu peak-current limit. FMI output current follows the passive sign convention. It is an averaged grid-following model, not the switched MathWorks inverter model.

## Native Simulink reference

`matlab/IEEE39_AvgInverter_Native_MatchedInit.slx` runs the same current equations with a 1 µs fixed-step RK4 solver. For each DPsim mode, the native model receives that run's inverter-terminal voltage trace and starts with the matching initial current. This is an open-loop component replay reference; only the ME-versus-CS pair is a closed-loop comparison of the full DPsim network. The notebook reports these scopes separately.

The original switched MathWorks Simulink model was also exported as FMI 3 CS during development. MATLAB R2025a produced a tool-coupling FMU that declares `needsExecutionTool=true` and contains no standalone `binaries/` model library, so DPsim cannot execute that export. The included dual-interface averaged model is the runnable ME/CS model used for this case.

## Measured run

All three executions use a 50 µs fixed step and a 0.2 s duration. The native Simulink inverter is replayed against the corresponding DPsim terminal-voltage trace. Its first call is a compile/cache warm-up; reported native runtimes are from subsequent calls. Thus native Simulink is a matched open-loop component reference, while ME and CS are closed-loop network runs.

## Timing profile and interpretation

The runner enables DPsim's sequential task profile and direct-solver phase profile for both ME and CS. `*_TaskTimes.csv` reports each component's `MnaPreStep` and `MnaPostStep` average/minimum/maximum and total time. `*_TaskTimes_SchedulerSteps.csv` gives per-step scheduler wall time and task time. `*_SolverPhaseTimes.csv` reports RHS assembly, refactorization, numeric solve, and related solver phases. `*_SimulationStepTimes.csv` contains per-step `Simulation::step()` wall time. Profile instrumentation is enabled in both modes and is included in reported runtimes.

In the profiled 0.2 s run, ME took 3.86 s total (0.28 s initialization, 3.59 s simulation loop) and CS took 3.28 s total (0.27 s initialization, 3.01 s loop). ME was about 18% slower overall. Its `Gen30_AveragedInverter.MnaPreStep` averaged 40.9 µs per step (164 ms total) versus 2.31 µs (9.24 ms) for CS. ME made 112,000 model evaluations (28 per step) taking 147 ms total. CS made 4,000 `fmi3DoStep` calls.

The full solver task averaged 334 µs per step for ME and 245 µs for CS. That task includes the FMI nonlinear convergence checks and repeat linearizations in ME. The sparse numeric solve itself was effectively the same: 6.40 µs per step for ME and 6.45 µs for CS. RHS assembly was also nearly equal at about 189 µs per step. ME performed 8,000 partial refactorizations versus 4,000 for CS, because ME repeats matrix linearization during its Newton coupling. The main ME overhead is therefore FMI model evaluations plus the extra nonlinear coupling work, not the sparse numeric solve. After one discarded Simulink compile/cache warm-up, native Simulink replay took 1.30 s for the ME trace and 1.09 s for the CS trace.

After the 20 ms startup window, the maximum phase-current difference was 0.560 A for ME versus its matching native replay, 1.617 A for CS versus its matching native replay, and 1.494 A between CS and ME. At bus 30, the maximum phase-voltage difference between the closed-loop DPsim runs was 0.331 V over the post-startup window.

## Reproduce

From the repository root:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_ieee39_comparison.py --run
```

The command builds the FMU and runner, executes ME and CS, runs both native Simulink replays, and refreshes the metrics and figures. MATLAB R2025a with Simulink is required for `--run`; to analyze existing traces only, omit the flag.

Generated CSV logs, run timing manifests, component/task timing breakdowns, JSON/CSV metrics, and PNG figures are written under `build-fmi-cosim/ieee39-inverter/` (or `build-fmi/ieee39-inverter/` when that is the configured build). The raw per-component pre/post measurements are in each run's `*_TaskTimes.csv`; aggregated tables are `component_timing_breakdown.csv`, `solver_timing_breakdown.csv`, `scheduler_step_timing.csv`, and `simulation_step_timing.csv`. The comparison notebook is `IEEE39_Inverter_Comparison.ipynb`.

## Interpretation

The ME configuration uses the explicit single-pass network policy because the purpose is a deterministic mode comparison at the same 50 µs grid step. Its logged residual is a one-pass coupling mismatch, not a Newton convergence criterion. CS makes one `fmi3DoStep` call per network step. The notebook reports current errors against the native replay under identical terminal-voltage traces, and bus-voltage differences between the two closed-loop DPsim runs.
