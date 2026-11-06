# Small closed-loop inverter grid

This example is a three-phase voltage source, PI feeder, inverter at the PCC, and switched three-phase PCC load. The load increases by 50 MW at 2.5 s in each 5 s comparison. Simulations solve the grid and inverter together at every step; saved traces are outputs for analysis and are never replayed as inputs.

The retained notebooks are:

- `SmallGrid_DPsim_Comparison.ipynb` compares native DPsim, FMI Model Exchange (ME), and FMI Co-Simulation (CS), using native DPsim as the error reference.
- `SmallGrid_Post_Powerflow_Simulink_Comparison.ipynb` compares DPsim ME/CS, native Simulink, and Simulink CS after an independently calculated balanced AC operating point. Native Simulink is the error reference at each matched timestep.

Both reports use 25, 50, 100, 500, and 1000 µs. The previously excluded 200 µs and 625 µs cases are not part of these reports.

## Regenerate the DPsim comparison

From the repository root, run:

```bash
python dpsim-fmi/examples/small_inverter_grid/run_smallgrid_comparison.py
python dpsim-fmi/examples/small_inverter_grid/run_smallgrid_comparison.py --sweep-rest
python dpsim-fmi/examples/small_inverter_grid/run_native_inverter_sweep.py
```

The runner builds `FMI_SmallGrid_Inverter`, performs short smoke cases before the 5 s runs, validates sample grids and finite values, and writes metrics and traces below `build-fmi/small-inverter-grid/results/`.

## Regenerate the post-power-flow comparison

```bash
python dpsim-fmi/examples/small_inverter_grid/run_post_powerflow_sweep.py
python dpsim-fmi/examples/small_inverter_grid/run_post_powerflow_simulink_comparison.py
```

The first script runs DPsim ME and CS with power-flow initialization and the FMU's periodic 60 Hz current state. The second independently computes the Simulink operating point, creates a separately staged CS FMU with matching starting voltage, runs short smoke checks and fresh 5 s closed-loop Simulink native/CS cases, and regenerates the compact metrics, figures, and notebook. Outputs are under `build-fmi/small-inverter-grid/results/post_powerflow_fix/`.

To regenerate reports from completed simulation traces without rerunning the simulations, use `--analyze-only` with either `run_post_powerflow_sweep.py` or `run_post_powerflow_simulink_comparison.py`.

## Power-flow initialization support

The FMI ME and CS importers expose an opt-in `initializeFromPowerflow` option. When enabled, the wrapper reads the initialized terminal voltage from the connected DPsim nodes, calls FMI `Reset`, and re-enters FMI initialization mode with that voltage and the configured control inputs. The FMU's own initialization equations determine its internal states. The option requires the FMU to export `Reset`; it does not infer arbitrary internal FMU states from network P/Q. Existing initialization behavior is unchanged when the option is disabled.

In the post-power-flow case, the source–PI-line–PQ equivalent initializes the EMT network. The FMU current state is set to the periodic solution of its 5 ms response filter: at 60 Hz the gain is 0.46865 and phase lag is 62.05°. The corresponding filtered steady-state P/Q is used for the network initialization.

## Case parameters

| Parameter | Value |
|---|---:|
| Frequency | 60 Hz |
| PCC line-line RMS voltage | 362.2155 kV |
| Inverter active/reactive power | 125 MW / 80.881 Mvar |
| Inverter current limit | 1.20 pu |
| Inverter response time | 5 ms |
| PI-line series R/L per phase | 0.30 Ω / 20 mH |
| PI-line shunt capacitance per phase | 1 µF |
| Load step | 50 MW at 2.5 s |
| Main run | 5 s, fixed timestep |
