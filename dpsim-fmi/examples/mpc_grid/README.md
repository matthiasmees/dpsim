# MPC-controlled grid load imported in DPsim

This small FMI 3.0 Model Exchange example tests whether a sampled electrical
controller with scheduled events and rollback can be embedded in DPsim's
three-phase EMT network solve. The FMU represents an averaged, three-phase
converter and its L filter as a grid-connected **current-consuming** component.
The grid source, feeder impedance, and parallel native resistor load remain
native DPsim components.

The controller is a one-step finite-control-set MPC. At every 50 µs tick it
predicts the three phase currents for each of the eight two-level bridge
switching vectors (`u_inv ∈ {-350, +350} V` per leg), then minimizes the
predicted three-phase current-tracking error plus a switching penalty. Its
reference is a balanced 4 A peak, 50 Hz current in phase with the grid voltage.
The filter equations are `L di/dt = v_PCC - u_inv - R i`, with `L = 25 mH`
and `R = 0.4 Ω`. Positive terminal current flows from DPsim's PCC into the FMU,
so the mapped current sign is passive-positive.

## Build and run

From the repository root, with the `WITH_FMI=ON` build already configured:

```sh
bash dpsim-fmi/examples/mpc_grid/build_fmu.sh
mkdir -p build-fmi/mpc-grid-case
cd build-fmi/mpc-grid-case
../dpsim-fmi/FMI_EMT_Ph3_Mixed ../mpc-grid-fmu/mpc_grid_mpc.fmu \
  ../../dpsim-fmi/examples/mpc_grid/config.json
```

The default command-line case runs one 20 ms cycle. The executed 10-second notebook changes the reference at 5 s and includes pre-step, transient, and post-step tracking and power analysis. The sample
period matches the network step, so each MPC update is an FMI scheduled time
event at a network boundary. Rollback is enabled and implemented by copying the
small FMI instance state; DPsim can therefore retain its Newton-coupled solve.
The result logger records terminal current, current reference, selected bridge
voltage, MPC cost, controller tick count, PCC voltage, feeder/load currents,
Newton iterations, and FMI evaluation counters under `build-fmi/mpc-grid-case`.
From the repository root, `python3 dpsim-fmi/examples/mpc_grid/plot_results.py`
prints tracking/KCL/convergence metrics and saves `mpc_grid_results.png` there.

## Scope and limitations

This is a portable, hand-built FMI test FMU and not a MATLAB/Simulink export.
It tests importer compatibility and basic electrical behavior, not a validated
industrial MPC, switching-loss model, DC-link dynamics, protection logic, or
high-fidelity PWM converter. The pre-existing `microgrid_network` FMU was not
used because its Modelica electrical pins are internal and it exposes no
current-driven terminal interface for direct DPsim MNA connection.
