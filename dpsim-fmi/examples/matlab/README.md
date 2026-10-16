# MATLAB three-phase RLC export

Run in the MATLAB Command Window:

```matlab
addpath('/Volumes/SSDMatthias/development/dpsim/dpsim-fmi/examples/matlab')
outputDir = create_rlc_fmu_3ph;
```

The function chooses a fresh timestamped output directory beneath MATLAB's
current folder. Supply an absolute path to select another destination:

```matlab
outputDir = create_rlc_fmu_3ph('/Volumes/SSDMatthias/my_3ph_test');
```

Close the two generated benchmark models before rerunning the function. Existing
model/FMU files are preserved; use a fresh output directory for another run.

This is three independent series RLC branches connected in **grounded wye**,
with the same per-phase values as the scalar benchmark: 10 ohm, 20 mH, 100 uF.
The external source is balanced abc at 230 V phase-to-neutral RMS and 50 Hz.
There is no mutual coupling or floating-neutral constraint in this first case.

The export is an FMI 3.0 **Model Exchange** FMU with six continuous states and
separately named scalar ports:

| Role | Variables | Unit |
| --- | --- | --- |
| Voltage inputs | `u_a`, `u_b`, `u_c` | V |
| Current outputs | `i_a`, `i_b`, `i_c` | A |
| Diagnostic outputs | `v_C_a`, `v_C_b`, `v_C_c` | V |

Voltage is positive-terminal minus negative-terminal voltage. Positive current
enters the positive terminal. For the DPsim adapter the terminal order is
`{negative_abc, positive_abc}`, with the negative side grounded in this example.
The voltage source remains external to the FMU so the network supplies its input.

At time zero, `[u_a, u_b, u_c]` is approximately `[0, -281.691, 281.691]` V.
All six dynamic states start at zero. The three-phase SSN adapter seeds its
voltage history with the actual starting phase voltages; assuming three zero inputs would give an incorrect first trapezoidal step.

Outputs include:

- `dpsim_rlc_3ph_me.slx` and `dpsim_rlc_3ph_me.fmu`.
- `dpsim_rlc_3ph_reference.slx`, the balanced-source reference model.
- `reference_simulink.csv` and `reference_trapezoidal.csv`, including all three
  input voltages, currents and capacitor voltages.
- Matrix CSVs, `rlc_3ph_parameters.mat`, and `port_mapping.json` describing the
  intended electrical mapping, phase angles, initial inputs and timestep.
- `rlc_3ph_reference.fig` with the three-phase reference traces.

The MATLAB state vector is `[i_a, i_b, i_c, v_C_a, v_C_b, v_C_c]`. FMI state
names/order still need to be discovered through the model description/runtime;
they must not be assumed by the importer.

The references run for 0.1 s with a 10 microsecond trapezoidal step and tightly
toleranced Simulink `ode45`. Phase voltage/current sums are checked for balance.

For references without FMU compilation:

```matlab
create_rlc_fmu_3ph([], false)
```

The C++ importer now supports this three-channel mapping through
`LinearModel::loadMulti` and the general `LinearSSN` adapter. Run the
`FMI_EMT_Ph3_LinearSSN` target with the generated FMU, then open
`../FMI_EMT_Ph3_LinearSSN_MATLAB_Comparison.ipynb` for the comparison. See the
module README for build and run commands. The adapter preserves full matrix
coupling; this MATLAB benchmark itself uses uncoupled RLC phases.

## Nonlinear generator benchmark

For the next test, see [nonlinear_3ph/README.md](nonlinear_3ph/README.md) and
`create_sg_grid_fmu_3ph.m`. They build and simulate a sixteen-state three-phase
synchronous-generator grid with an AVR, droop governor, load step and voltage
reference step. The optional FMU contains the ten-state controlled generator.
It requires a future nonlinear DPsim FMI adapter; the existing LTI importer
cannot run it.
