# IEEE 39 bus inverter comparison

This case compares one FMI 3 averaged inverter model in DPsim Model Exchange (ME), DPsim Co-Simulation (CS), and native Simulink. The same standalone FMU exposes both ME and CS from the same `IEEE39AvgInverter.cpp` source. The Simulink replay model independently implements the same current-control equations.

## Grid and event

The network uses the hardcoded three-phase IEEE 39 bus data already used by the DPsim `IEEE39_Ph3` example: 39 buses, 46 branches, 10 synchronous generators, 60 Hz. The inverter is placed at bus 30. To preserve the standard pre-event active and reactive power at bus 30, generator 30 remains at half of its published (P/Q) operating point and the FMI inverter supplies the other half: 125 MW and 80.881 Mvar. All other generators, loads, and lines remain unchanged. A balanced 100 MW load is switched onto bus 39 at 2.5 s, the midpoint of the run.

Runs use a 50 µs EMT network step for 5 s (100,000 network steps). The average converter uses an instantaneous (P/Q) current reference in the stationary αβ frame, a 5 ms first-order current response, and a 1.2 pu peak-current limit. FMI output current follows the passive sign convention. It is an averaged grid-following model, not the switched MathWorks inverter model.

## Native Simulink reference

`matlab/IEEE39_AvgInverter_Native_MatchedInit.slx` runs the same current equations as a native Simulink inverter. There is one reference run, driven by the DPsim ME run's recorded inverter-terminal voltage, initialized from that operating point, and integrated with the fixed-step ODE4 solver at 50 µs in Accelerator mode. The comparison overlays this single Simulink current trace with both DPsim ME and CS currents and plots each DPsim error against that same Simulink run. Terminal-voltage traces and their errors are included as well. This is a component-level native Simulink reference; ME and CS remain closed-loop runs of the IEEE39 network.

The original switched MathWorks Simulink model was also exported as FMI 3 CS during development. MATLAB R2025a produced a tool-coupling FMU that declares `needsExecutionTool=true` and contains no standalone `binaries/` model library, so DPsim cannot execute that export. The included dual-interface averaged model is the runnable ME/CS model used for this case.

## Measured run

All three executions use a 50 µs fixed step and a 5 s duration (100,000 intervals). The load step occurs at 2.5 s. The native Simulink inverter uses one fixed ME PCC-voltage trace, and both DPsim current results are compared to that same native output. Simulink uses a 0.05 s Accelerator compile/cache warm-up; the reported native timing is for the subsequent full 5 s run. In the current run, Simulink's measured simulation time was 3.21 s; the MATLAB batch took 107.32 s including startup, compilation, warm-up, and the measured run. Thus native Simulink is a matched component reference, while ME and CS are closed-loop network runs.

## Timing profile and interpretation

The runner enables DPsim's sequential task profile and direct-solver phase profile for both ME and CS. `*_TaskTimes.csv` reports each component's `MnaPreStep` and `MnaPostStep` average/minimum/maximum and total time. `*_TaskTimes_SchedulerSteps.csv` gives per-step scheduler wall time and task time. `*_SolverPhaseTimes.csv` reports RHS assembly, refactorization, numeric solve, and related solver phases. `*_SimulationStepTimes.csv` contains per-step `Simulation::step()` wall time. Profile instrumentation is enabled in both modes and is included in reported runtimes.

The notebook and timing CSVs report the current 5 s run, including total/init/loop timing, all component pre-step and post-step tasks, solver phases, and per-step wall time. The profile instrumentation is enabled for both DPsim modes. Simulink uses the same 50 µs fixed step and 5 s stop time; its short compile/cache warm-up is excluded from the reported run time.

The generated `ieee39_metrics.json` and `ieee39_metrics.csv` contain current and voltage errors over the post-startup period and the 20 ms window after the 2.5 s load step. These are refreshed with the latest run by the reproduction command below.

## Reproduce

For a same-network comparison of **native DPsim, FMI Model Exchange, and FMI
Co-Simulation**, run:

```sh
python3 -u dpsim-fmi/examples/ieee39_inverter/run_ieee39_dpsim_native_me_cs_sweep.py
```

This executes fresh, powerflow-initialized 10 s closed-loop cases with the 100
MW bus-39 load step at 5 s. The default sweep is 10, 25, 50, 100, 500, 1000,
and 5000 µs; the 5000 µs case is a coarse-step stress run. Native DPsim uses
the same averaged inverter equations and starts its current state from the
FMU's P/Q current reference at the solved bus-30 voltage. The full-resolution
logs are reduced to compact traces and RMS metrics before being discarded.
The notebook `IEEE39_DPsim_Comparison_10s_Powerflow.ipynb` reports matched-step
ME/CS errors relative to native DPsim, runtime, initialization checks, rolling
error traces, and a load-step voltage/current zoom. Outputs are written under
`build-fmi/ieee39-inverter/dpsim-native-me-cs-10s-powerflow/`.
The standard 10–1000 µs data completed under the previous two-iteration ME
ceiling; the 5 ms ME stress case required residual-based Newton with the
default 20-iteration ceiling after the two-iteration run failed its residual
check. The reproduction runner now uses residual-based Newton with that
default ceiling for all regenerated cases.

From the repository root:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_ieee39_comparison.py --run
```

The command builds the FMU and runner, executes ME and CS, runs one native Simulink reference, and refreshes the metrics and figures. MATLAB R2025a with Simulink is required for `--run`; to analyze existing traces only, omit the flag.

Generated CSV logs, run timing manifests, component/task timing breakdowns, JSON/CSV metrics, and PNG figures are written under `build-fmi-cosim/ieee39-inverter/` (or `build-fmi/ieee39-inverter/` when that is the configured build). The raw per-component pre/post measurements are in each run's `*_TaskTimes.csv`; aggregated tables are `component_timing_breakdown.csv`, `solver_timing_breakdown.csv`, `scheduler_step_timing.csv`, and `simulation_step_timing.csv`. The comparison notebook is `IEEE39_Inverter_Comparison.ipynb`.

For the longer CS-versus-ME timestep sweep, run:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_ieee39_timestep_sweep.py \
  --steps-us 500 625 1000 --stop-time-s 20 --load-step-s 10 \
  --output-name timestep-sweep-20s --reuse
```

This runs the 20 s pairs (or reuses completed pairs) and saves metrics and figures under `build-fmi/ieee39-inverter/timestep-sweep-20s/`. The event is placed at the midpoint. A 625 µs step is used instead of 750 µs because 625 µs divides both the 20 s stop and 10 s event times exactly.

To generate the native Simulink references for every step and duration in both sweeps, run:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_native_sweep_references.py
```

This requires MATLAB R2025a with Simulink. Each native case replays the paired ME terminal-voltage trace using the same fixed step and logged duration. The notebook plots ME and CS current errors and terminal-voltage errors as centered 0.5 s three-phase RMS percentages of the native Simulink reference. The replay is a component-level inverter reference, not a native IEEE39 network simulation.

To cross-check the FMI Co-Simulation host itself against Simulink, run:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_simulink_cosim_sweep.py
```

This imports the runnable FMI 3 CS FMU through Simulink's FMU Import block and runs it beside an independently coded native Simulink inverter. Both receive the exact logged PCC voltage samples from each matching DPsim CS case. The sweep covers 25, 50, 100, 200, and 500 µs for 5 s and 500, 625, and 1000 µs for 20 s. The load steps remain at 2.5 s and 10 s respectively. This is a component-level replay that isolates host execution; it is not a Simulink implementation of the IEEE39 network. Outputs and timing are saved under `build-fmi/ieee39-inverter/simulink-fmu-cosim/`. The comparison notebook reports both imported-FMU-vs-DPsim-CS host error and native-vs-FMU Simulink equation/integration error.

## Interpretation

For ordinary steps, the ME `NonlinearSSN` component uses the fixed-step trapezoidal state update and DPsim MNA Newton iterations. Convergence requires both scaled residuals to be at most 1, using default state tolerances of 1e-11 absolute / 1e-10 relative and current tolerances of 1e-9 A absolute / 1e-10 relative. The `newton_residual` log is the maximum normalized state/current residual. The `SinglePassExplicit` setting applies when an interval crosses an FMI internal time event; it does not replace Newton on ordinary intervals. The CS branch makes one `fmi3DoStep` call per network step and has no master-side Newton iteration. Native Simulink uses fixed-step ODE4, so it has no Newton acceptance epsilon. The notebook reports ME and CS current and terminal-voltage errors against the single native Simulink reference, plus selected bus-voltage differences between the closed-loop DPsim runs. The Simulink reference uses the ME voltage trace; CS voltage differences therefore also reflect the small ME/CS network-coupling difference.

## Closed-loop Simulink IEEE39 model

`matlab/IEEE39_ClosedLoop.slx` contains the same closed-loop IEEE39 EMT network as `IEEE39_AvgInverter.cpp`. The FMU Import block runs the same FMI 3 Co-Simulation inverter FMU at bus 30. A discrete S-function hosts the DPsim MNA network built from the same case tables, machine/load/line/transformer parameters, load-flow initialization, 100 MW bus-39 load step, and time grid. This makes the network topology and MNA solver identical by construction. It is a Simulink-hosted DPsim network, not a second, independently implemented Simscape network; use it to compare FMI block scheduling and host behavior, not as an independent network-solver benchmark.

The fixed-step model retains all 117 phase voltages, inverter currents, network step times, and initialization time. It is run at 25, 50, 100, 200, and 500 µs for 5 s, and at 500, 625, and 1000 µs for 20 s, with the 100 MW step at each run's midpoint. The matching DPsim CS traces and wall timings already exist under the same build tree. The short validation case at 50 µs for 0.1 s uses a 0.05 s load step.

Rebuild and run the closed-loop model from the repository root:

```sh
cmake --build build-fmi --target IEEE39_SimulinkNetwork
./dpsim-fmi/examples/ieee39_inverter/matlab/build_ieee39_network_sfun.sh
python3 dpsim-fmi/examples/ieee39_inverter/run_simulink_closed_loop_sweep.py
python3 dpsim-fmi/examples/ieee39_inverter/analyze_simulink_closed_loop.py
```

The runner uses MATLAB R2025a by default. Set `MATLAB_BIN` to select a different executable and `MATLAB_ROOT` when compiling the S-function. On the installed host, MATLAB's bundled C++ runtime is older than DPsim's, so the runner preloads `/usr/lib64/libstdc++.so.6`; override this with `DPsim_MATLAB_LIBSTDCXX` if needed. The FMU, S-function, and DPsim network library must be built for the same platform and ABI.

For fresh matched DPsim CS references rather than reusing the prior sweep logs, run `python3 dpsim-fmi/examples/ieee39_inverter/run_matched_dpsim_cs_sweep.py` before the analyzer. It writes a separate `dpsim-cs-matched-rerun/` tree. The focused results notebook is `IEEE39_ClosedLoop_Simulink_Comparison.ipynb`; it displays the Simulink-only load-step waveforms, all-timestep error plots, and detailed timing table.

The short 50 µs check matched 2,001 samples. Against the DPsim CS trace, the maximum bus-30 phase-voltage error was 0.402 V (0.000116% of 345 kV base), and maximum inverter-current error was 0.000211 A (0.000042% of the inverter peak-current base). Across the 5 s sweep, maximum bus-30 error ranged from 0.025% to 0.483% of base. Across 20 s runs it ranged from 7.61% at 500 µs to 14.1% at 1000 µs; long-run host trajectories therefore do not remain identical even with the same network kernel. The exact per-case comparison and component-step timings are in `build-fmi/ieee39-inverter/simulink-closed-loop-sweep-exact-order/closed_loop_comparison_metrics.json`; the 0.5 s RMS error plots overlay every tested timestep for each duration.

The measured 5 s 50 µs Simulink run took 16.45 s wall time, including 12.80 s in the DPsim network S-function; the matched DPsim CS process took 76.62 s, including its CSV logger and per-step scheduler instrumentation. At 20 s / 1000 µs, these were 3.86 s total and 2.75 s in the network S-function, versus 14.90 s in the instrumented DPsim run. These full-process times are not a clean speed benchmark because the two hosts write their traces differently. Compare network step timings and error metrics for solver work and fidelity, and treat long-run host drift as a result rather than a numerical identity claim.

### Full closed-loop native-reference sweep (20 s)

`IEEE39_20s_AllCases_NativeSimulink_Comparison.ipynb` compares independently started, closed-loop runs at 25, 50, 100, 200, 500, 625, and 1000 µs for 20 s, with a 100 MW load step at 10 s. Native Simulink and Simulink FMU CS each create their own DPsim IEEE39 network instance through the MNA S-function; DPsim ME/CS create their own network directly. No run replays a saved voltage/current trajectory. The native inverter's initial current is calculated at simulation startup from its live initialized bus-30 voltage, rather than from a stored phase-voltage vector. The DPsim ME/CS FMU uses electrical voltage/current start values matching the same solved load-flow operating point, then evolves from that state in closed loop. These are scalar initial conditions, not replayed dynamic traces.

Native Simulink remains a matched-behavior comparison, but it is not a method-independent physical truth: it advances the inverter with fixed-step ODE4 while the network uses sampled current exchange; DPsim ME uses trapezoidal state integration coupled through MNA Newton iterations, and CS runs hold the exchanged voltage over each macro step. Equal communication steps therefore do not make these discrete-time maps identical. To assess integration convergence separately, the analysis also compares every trajectory against the **25 µs DPsim ME run**, labeled as a fine-step DPsim convergence reference (not an independent physical reference). The earlier 5 s “native” comparison replayed the DPsim ME PCC-voltage trace and was not a closed-loop IEEE39 native run; it must not be used to rank closed-loop ME accuracy.

Run the three sweep commands from the repository root to regenerate the data:

```sh
python3 -u dpsim-fmi/examples/ieee39_inverter/run_ieee39_timestep_sweep.py \
  --steps-us 25 50 100 200 500 625 1000 --stop-time-s 20 --load-step-s 10 \
  --output-name timestep-sweep-20s-live-init
python3 -u dpsim-fmi/examples/ieee39_inverter/run_ieee39_native_closed_loop_sweep.py \
  --stop-time-s 20 --steps-us 25 50 100 200 500 625 1000 --load-step-s 10 \
  --output-name native-closed-loop-20s-live-init
python3 -u dpsim-fmi/examples/ieee39_inverter/run_simulink_closed_loop_sweep.py \
  --stop-time-s 20 --steps-us 25 50 100 200 500 625 1000 \
  --output-name simulink-closed-loop-20s-live-init \
  --reference-sweep-dir timestep-sweep-20s-live-init
```

The notebook's analysis cell reads these saved runs, creates native-behavior comparison and fine-step ME convergence tables/plots, reports t=0 initialization consistency and ME Newton iterations/residuals, shows load-step waveforms, compares Simulink CS directly with DPsim CS, and reports wall, network-step, DPsim component-task, and solver-phase timings. Native-relative errors are normalized by native Simulink RMS after 20 ms; fine-reference errors are normalized by the corresponding 25 µs DPsim ME RMS. Time-series error curves use 0.5 s phase-a RMS windows. Results and analysis products are under `build-fmi/ieee39-inverter/`; rerunning the notebook does not repeat the simulations.

In the corrected sweep, ME is farther from native Simulink at coarse steps (at 1000 µs, post-startup current RMS error is 28.33% of native RMS for ME and 7.37% for DPsim CS). Against the 25 µs DPsim ME refinement trajectory, however, the 1000 µs errors are 2.83% for ME and 24.66% for CS on inverter current; bus-30 voltage errors are 2.31% and 6.45%. At 625 µs, the corresponding current errors are 1.10% and 14.80%. ME therefore converges much more closely to the fine-step ME trajectory; the native-relative plot reflects the ODE4/held-sample method difference. All ME runs used two Newton iterations per logged nonzero step; the maximum normalized residual was 0.238639 at 1000 µs, below the acceptance limit of 1.

The timing difference is also expected from the coupling method. At 25 µs, DPsim ME took 729.4 s and DPsim CS took 602.0 s. The ME component `Solve` task accounted for 274.1 s versus 199.0 s in CS, while the inverter `MnaPreStep` task took 32.2 s versus 1.9 s. ME requested about 28 FMI evaluations per network step to form and iterate the nonlinear tangent; CS made one `fmi3DoStep` call per step. The Newton residuals passed throughout, so there is no evidence that the ME result degradation came from nonconvergence.

## IEEE39 DPsim comparison with subcomponents disabled

The native/ME/CS 10 s comparison is in `IEEE39_DPsim_Comparison_10s_Powerflow.ipynb`. The sweep runs 10, 25, 50, 100, 500, 1000, and 5000 µs, initializes from the solved IEEE39 powerflow, and applies the 100 MW bus-39 load step at 5 s. The EMT PiLine, RXLoad, and transformer components use their existing models with child subcomponents disabled; this reduces scheduler overhead without changing the model equations. The VBR generators and all inverter modes are unchanged. Data are saved under `build-fmi/ieee39-inverter/dpsim-native-me-cs-10s-nosubcomponents/`.

Regenerate the runs and notebook figures from the repository root with:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_ieee39_dpsim_native_me_cs_sweep.py
```

In the rerun, 10 µs solver times were 708.85 s native, 782.37 s ME, and 654.05 s CS, about 10–11% faster than the prior runs. At 25 µs they were 281.17 s, 313.27 s, and 261.96 s, about 11–12% faster. ME tracks native to numerical precision through 1 ms; CS errors grow with timestep and are 5.11% voltage / 134.51% current RMS at 5 ms. A trial with SSN passive-network components made the 10 µs CS result diverge (72.19% voltage and 71.63% current RMS error), so those results are excluded from the notebook.

## Corrected IEEE39 FMI Co-Simulation coupling

The preceding IEEE39 CS sweeps used single-pass explicit coupling. Their current
error rose with timestep because the FMU advanced using the previous accepted
bus voltage: 0.20% at 10 µs, 1.96% at 100 µs, 19.94% at 1 ms, and 134.5% at
5 ms. The bus/branch mapping was not the cause. The powerflow and EMT builder
share the same MATPOWER bus and branch tables, and the Simulink S-function
reuses this same IEEE39 builder. All three methods also recorded identical
initial bus-30 voltage and inverter current.

The IEEE39 CS FMU now implements FMI get/set/free-state functions. The `cs`
case opts into RollbackNewton and replays each network interval with four
midpoint voltage samples; `cs-explicit` remains available in the executable
to reproduce single-pass coupling. The FMU keeps its exact exponential
current-state update, while native DPsim and ME use trapezoidal integration.
This removes the previous one-step lag without forcing the CS FMU onto the ME
integrator.

A fresh 10 s, 100 µs run with a 100 MW load step at 5 s completed 100,001
samples for native, ME, and CS. ME-to-native full-run RMS errors were
1.18e-10% for voltage and 4.02e-9% for current. CS-to-native errors were
0.000143% for voltage and 0.0223% for current. CS converged in two iterations
per step, with maximum scaled current residual 0. The errors were 0.0000365%
voltage / 0.0223% current before the load step and 0.000199% / 0.0223% after
it. Solver times were 70.7 s native, 78.8 s ME, and 72.2 s CS.

A five-second sweep at 100, 500, and 1000 µs gives CS voltage/current RMS
errors of 0.0000365% / 0.0223%, 0.00091% / 0.558%, and 0.00361% / 2.26%,
respectively. Each CS run converged in two Newton iterations. The 5 ms stress
case still has 64.2% current RMS error with rollback enabled; it is too coarse
for this 60 Hz network and should be excluded from normal-step conclusions.
The fresh plots and metrics are in
`IEEE39_RollbackNewton_Validation.ipynb` and
`build-fmi/ieee39-inverter/ieee39-rollback-5s-sweep/metrics.json`. The 10 s
result is in `build-fmi/ieee39-inverter/ieee39-rollback-10s-100us/metrics.json`.
The previous `IEEE39_DPsim_Comparison_10s_Powerflow.ipynb` is marked as
historical for CS because its saved curves use single-pass explicit coupling.

Rebuild and rerun the corrected checks from the repository root:

```sh
bash dpsim-fmi/examples/ieee39_inverter/sources/build_fmu.sh build-fmi/ieee39-inverter-fmu
cmake --build build-fmi --target FMI_IEEE39_AvgInverter -j4
python3 -u dpsim-fmi/examples/ieee39_inverter/run_ieee39_dpsim_native_me_cs_sweep.py \
  --steps-us 100 500 1000 --duration-s 5 --load-step-s 2.5 \
  --output-name ieee39-rollback-5s-sweep
python3 -u dpsim-fmi/examples/ieee39_inverter/run_ieee39_dpsim_native_me_cs_sweep.py \
  --steps-us 100 --duration-s 10 --load-step-s 5 \
  --output-name ieee39-rollback-10s-100us
```
