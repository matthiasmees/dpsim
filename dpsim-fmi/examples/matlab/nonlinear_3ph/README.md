# Nonlinear three-phase synchronous-generator benchmark

This is an explicit-equation **EMT** test grid in ordinary Simulink. The model
retains stator and field flux transients, electrical rotor position, rotor speed,
voltage control and prime-mover control. It is intentionally small enough to
inspect. The nonlinear FMI adapter and validated trapezoidal comparison are now
available; see [the nonlinear importer guide](../../../NONLINEAR.md).

It is a teaching machine with a field winding and saliency (`Ld != Lq`), without
damper windings, magnetic saturation, harmonics or manufacturer parameter fitting.
The AVR and governor are simple benchmark controllers, not IEEE exciter/turbine
standard models. The load is switched constant **impedance**, not constant power:
6/12 kW are nominal powers at 230 V phase RMS; actual load power varies with v².

The machine formulation is in the same field-winding/no-damper class described by
[MathWorks' Synchronous Machine Model 1.0](https://www.mathworks.com/help/simscape-electrical/ref/synchronousmachinemodel1.0.html).
This implementation uses a **power-invariant** Park transform and its own stated
sign conventions/parameters. It is not an implementation or validated numerical
replica of that Simscape block. MATLAB Function blocks are populated with the
[documented programmatic API](https://www.mathworks.com/help/stateflow/api/stateflow.emchart.html);
continuous states live in ordinary Integrator blocks.

## Grid

```text
mechanical reference -> droop governor -> mechanical power
                                              |
voltage reference -> PI AVR + exciter -> synchronous generator
                                              |
                                     three-phase local bus
                                        |           |
                              grounded-wye load   R-L feeder
                                        |           |
                                      ground    stiff 50 Hz grid
```

A grounded-wye capacitor is connected to the bus. It represents a small shunt
capacitance and supplies explicit bus-voltage states. It is part of the physical
benchmark, not a numerical capacitor added by an importer. The machine neutral,
load and grid neutral are grounded. Positive abc sequence is used throughout.

| Quantity | Value |
| --- | --- |
| Voltage | 230 V phase-to-neutral RMS, about 398 V line-to-line RMS |
| Frequency / power base | 50 Hz / 15 kVA |
| Initial generator active power | 8 kW, at its terminals |
| Initial local load | 6 kW at nominal voltage |
| Initial grid exchange | about 2 kW exported, including feeder losses |
| Feeder, per phase | 0.3 ohm + 5 mH |
| Bus shunt, per phase | 40 µF |
| Stator resistance | 0.2 ohm |
| Ld / Lq / L0 | 60 / 45 / 10 mH |
| Field self / mutual inductance | 1.5 / 0.27 H, power-invariant coordinates |
| Field resistance | 1.5 ohm |
| H / damping | 1.5 s / 2 pu power per pu speed deviation |
| AVR Kp / Ki / anti-windup gain | 6 / 12 s⁻¹ / 10 s⁻¹ |
| Voltage measurement / exciter time constants | 20 / 50 ms |
| Excitation command limits | 0 to 5 pu |
| Governor droop / time constant | 5% / 0.2 s |
| Prime-mover power limits | 0 to 1.2 pu |

The field voltage base is chosen so that 1 pu produces nominal no-load voltage
at nominal speed. The coupled inductance matrix is positive definite. All these
values are defined in `sg_benchmark_parameters.m`.

## Experiment

- **0–1 s:** remain on the initialized balanced periodic operating orbit.
- **1 s:** add 6 kW of nominal resistive load, giving 12 kW total nominal load.
- **2 s:** raise the AVR reference from 1.00 to 1.02 pu.
- **4 s:** stop. Log EMT and control signals on a 50 µs output grid.

The load event tests bus voltage, stator/field transients, torque, rotor speed
and governor response. The voltage-reference event makes the AVR response easy
to distinguish from the load disturbance. A stiff grid constrains terminal
frequency; the plotted frequency deviation is **rotor speed**, not an independent
grid-frequency disturbance. The generator need not supply the entire load
increment because the feeder/grid shares it. A regulated bus above grid voltage
can create substantial reactive exchange.

Limiters are present, but whether they are reached is reported rather than
assumed. This is a continuous-control benchmark with scheduled external steps;
it does not yet exercise digital sample clocks, protection relays or breaker
state machines.

## Run

In MATLAB R2025a or later:

```matlab
addpath('/Volumes/SSDMatthias/development/dpsim/dpsim-fmi/examples/matlab')
outputDir = create_sg_grid_fmu_3ph; % build and simulate
```

To also export the controlled generator as FMI 3.0 Model Exchange:

```matlab
outputDir = create_sg_grid_fmu_3ph([], true);
```

For a chosen **fresh** directory:

```matlab
outputDir = create_sg_grid_fmu_3ph('/Volumes/SSDMatthias/my_sg_3ph_test', true);
```

Close the two generated models before repeating. Existing models are preserved;
use a fresh directory. Simscape and Optimization Toolbox are not required.
Export additionally needs the standalone FMU exporter and supported C compiler;
see [exportToFMU](https://www.mathworks.com/help/fmuexport/ref/exporttofmu.html).

Open `dpsim_sg_3ph_reference.slx` for the complete grid. Its
`Controlled_generator` subsystem contains separate `Machine`, `AVR`, and
`Governor` subsystems. Each has an `Equations` MATLAB Function block, visible
continuous `States` Integrator and embedded numeric `Parameters_SI` Constant.
The saved models do not depend on function files or base-workspace variables to
simulate. Keep this source directory when rebuilding models or running the
independent MATLAB reference.

The script saves:

- `dpsim_sg_3ph_me.slx`: ten-state controlled generator with explicit ports.
- `dpsim_sg_3ph_reference.slx`: complete sixteen-state grid and step sources.
- `dpsim_sg_3ph_me.fmu` if export is enabled.
- `reference_simulink.csv`: adaptive Simulink `ode15s` output.
- `reference_uniform.csv`: Simulink traces interpolated onto a 50 µs grid,
  with exact right-hand scheduled input values restored at step times.
- `reference_matlab_ode15s.csv`: separately wired MATLAB ODE reference, solved
  more tightly, restarted exactly at both disturbance times.
- `validation.json`: integration/wiring comparison, initialization/balance
  checks, speed/control ranges and instantaneous energy conservation residual.
- `sg_3ph_reference.png` and `.fig`: EMT zooms and control/power responses.
- `sg_3ph_reference.mat`, `sg_3ph_parameters.mat`, `signal_names.json` and
  `port_mapping.json`: numeric traces, parameters and the electrical contract.

The two solvers use the same component equation kernels, so their agreement
checks solver integration and Simulink wiring, not independent machine-model
validation. Energy conservation separately checks Park scaling, electrical
power and torque signs. The reference tolerances are assertions in
`sg_validate_reference.m`, not a guarantee for other parameter sets.

On macOS, check the native exported FMU against the initial Simulink outputs:

```sh
python3 dpsim-fmi/tests/validate_sg_fmu.py /absolute/path/to/benchmark-output
```

This loads the host FMU binary in-process, supplies all five initial inputs,
checks ten continuous states, compares all 23 outputs and verifies that only
rotor electrical position has a nonzero derivative on the initial periodic
orbit. It writes `fmu_validation.json`. This is an export/initialization check,
not a nonlinear DPsim transient simulation.

## States and equations

The controlled generator has ten states:

```text
machine:   psi_d, psi_q, psi_0, psi_f, theta_e, speed_pu  (6)
AVR:       filtered_voltage_squared_pu, integral, field_supply_pu (3)
governor:  mechanical_power_pu (1)
```

The network adds `v_bus_a/b/c` and `i_line_a/b/c` (6). The MATLAB reference order
is machine, AVR, governor, network. FMI state ordering must be discovered from the
export, not assumed to match this vector.

Let `T(theta)` be orthonormal, with rows proportional to cos(theta+phase),
−sin(theta+phase), and the zero-sequence basis. Stator currents `id,iq,i0` are
positive **into the machine**:

```text
[psi_d; psi_f] = [Ld M; M Lf] [id; if]
psi_q = Lq iq;  psi_0 = L0 i0
v_dq0 = T(theta) v_abc

psi_d' = vd − Rs id + we psi_q
psi_q' = vq − Rs iq − we psi_d
psi_0' = v0 − Rs i0
psi_f' = vf − Rf if
we = wb speed_pu; theta' = we
Pgap = −we (psi_d iq − psi_q id)
speed_pu' = (Pm_pu − Pgap/Sbase − D(speed_pu−1)) / (2 H speed_pu)
i_gen_abc = −T(theta)' [id; iq; i0]
```

Nonlinearities include the state-dependent Park transform, speed×flux terms,
bilinear electromagnetic torque, speed-dependent swing equation, square-root
voltage measurement and saturated control commands. Magnetic inductances remain
constant in dq0; saturation is a separate future extension.

The AVR filters `sum(v_abc²)/(3 VphaseBase²)`, takes its square root, applies
PI control with back-calculation anti-windup and limits, then a first-order
exciter. The governor applies a droop correction and power limits followed by
a first-order prime mover.

Network equations are explicit EMT KCL/KVL:

```text
Cbus v_bus' = i_gen − i_line − Gload v_bus
Lline i_line' = v_bus − v_grid − Rline i_line
```

Load conductance is constant on each event interval. No phasor approximation
is used during the transient simulation. Phasors are used only to calculate the
initial balanced operating point with `fzero`, avoiding an artificial startup.

## Future DPsim import contract

The FMU contains the **machine and both controllers**, not the feeder, shunt,
load or event scheduler. Inputs are scalar `u_a/u_b/u_c` in V, `vref_pu`, and
`pmref_pu`. Outputs `i_a/i_b/i_c` in A are positive **out of the generator into
the bus**; this is the opposite sign of the passive RLC benchmark. Other scalar
outputs expose speed, flux, field, power and controller diagnostics.

The initial terminal voltages, reference values and port-current polarity are
saved explicitly. The FMU has fixed initial states for the calculated operating
point, not a general automatic load-flow initialization for arbitrary grids.
A changed topology/operating point requires recalculated initial conditions.
The generator terminal is a grounded abc bundle, not a floating neutral.

**The current `LinearModel` / `LinearSSN` importer cannot simulate this FMU.** It
must reject the extra control inputs/nonlinear equations. Do not fit fixed
A/B/C/D matrices and treat them as the nonlinear model.

A future adapter must retain the FMI Model Exchange instance, evaluate its
nonlinear state derivatives/current outputs inside a converged DPsim implicit
step, and couple the voltage/current Jacobian to network MNA. Controller states
are part of that same coupled state update. The network's scheduled load step
and the voltage-reference change must land exactly on solver timestamps. An
initial target would be implicit trapezoidal integration with Newton iterations
and rollback of trial FMI states; no Co-Simulation `doStep()` is needed.
The positive injection sign must be reflected in MNA stamps. That importer work
is intentionally deferred until these MATLAB traces are accepted.
