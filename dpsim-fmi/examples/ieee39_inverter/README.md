# IEEE 39 bus inverter comparison

This case compares one FMI 3 averaged inverter model in DPsim Model Exchange (ME), DPsim Co-Simulation (CS), and native Simulink. The same standalone FMU exposes both ME and CS from the same `IEEE39AvgInverter.cpp` source. The Simulink replay model independently implements the same current-control equations.

## Grid and event

The network uses the hardcoded three-phase IEEE 39 bus data already used by the DPsim `IEEE39_Ph3` example: 39 buses, 46 branches, 10 synchronous generators, 60 Hz. The inverter is placed at bus 30. To preserve the standard pre-event active and reactive power at bus 30, generator 30 remains at half of its published (P/Q) operating point and the FMI inverter supplies the other half: 125 MW and 80.881 Mvar. All other generators, loads, and lines remain unchanged. A balanced 100 MW load is switched onto bus 39 at 170 ms.

Runs use a 50 µs EMT network step for 0.2 s (4,000 network steps). The average converter uses an instantaneous (P/Q) current reference in the stationary αβ frame, a 5 ms first-order current response, and a 1.2 pu peak-current limit. FMI output current follows the passive sign convention. It is an averaged grid-following model, not the switched MathWorks inverter model.

## Native Simulink reference

`matlab/IEEE39_AvgInverter_Native_MatchedInit.slx` runs the same current equations with a 1 µs fixed-step RK4 solver. For each DPsim mode, the native model receives that run's inverter-terminal voltage trace and starts with the matching initial current. This is an open-loop component replay reference; only the ME-versus-CS pair is a closed-loop comparison of the full DPsim network. The notebook reports these scopes separately.

The original switched MathWorks Simulink model was also exported as FMI 3 CS during development. MATLAB R2025a produced a tool-coupling FMU that declares `needsExecutionTool=true` and contains no standalone `binaries/` model library, so DPsim cannot execute that export. The included dual-interface averaged model is the runnable ME/CS model used for this case.

## Measured run

In the recorded 0.2 s run on `perf/ieee39-ieee118-optimization`, DPsim ME took 1.33 s and DPsim CS took 0.75 s for 4,000 network steps. The native 1 µs Simulink replay took 20.81 s for the ME voltage trace and 2.55 s for the CS voltage trace (the combined MATLAB batch wall time was 110.15 s, including MATLAB startup and model setup). These runtimes are specific to the host and should not be treated as portable performance claims.

After the 20 ms startup window, the maximum phase-current difference was 0.560 A (0.159% of peak current base) for ME versus its matching native replay, 1.617 A (0.459%) for CS versus its matching native replay, and 1.494 A (0.424%) between CS and ME. Around the 170 ms load step, the CS-versus-ME phase-current maximum difference was 1.483 A. At bus 30, the maximum phase-voltage difference between the closed-loop DPsim runs was 0.331 V over the post-startup window. Voltage and current percentages use their respective peak bases.

## Reproduce

From the repository root:

```sh
python3 dpsim-fmi/examples/ieee39_inverter/run_ieee39_comparison.py --run
```

The command builds the FMU and runner, executes ME and CS, runs both native Simulink replays, and refreshes the metrics and figures. MATLAB R2025a with Simulink is required for `--run`; to analyze existing traces only, omit the flag.

Generated CSV logs, run timing manifests, JSON/CSV metrics, and PNG figures are written under `build-fmi-cosim/ieee39-inverter/`. The comparison notebook is `IEEE39_Inverter_Comparison.ipynb`.

## Interpretation

The ME configuration uses the explicit single-pass network policy because the purpose is a deterministic mode comparison at the same 50 µs grid step. Its logged residual is a one-pass coupling mismatch, not a Newton convergence criterion. CS makes one `fmi3DoStep` call per network step. The notebook reports current errors against the native replay under identical terminal-voltage traces, and bus-voltage differences between the two closed-loop DPsim runs.
