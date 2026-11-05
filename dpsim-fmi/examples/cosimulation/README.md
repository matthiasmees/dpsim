# FMI Co-Simulation EMT example

## Co-Simulation counterparts

These runners cover the one-FMU and parallel two-FMU electrical coupling
patterns used by the existing FMI examples. They provide a small reusable EMT
grid for each port shape; they do not reproduce each source model's original
event schedule or validation scenario.

The `dual_interface_rl` case includes one FMI 3 test FMU with both Model
Exchange and Co-Simulation interfaces, a finite-impedance sag grid, a native
Simulink state-space reference, and a comparison notebook. Use it when the
same archive must be exercised through both DPsim import paths.

| Existing FMI example | Co-Simulation runner |
|---|---|
| `linear_scalar` | `FMI_EMT_CoSimulation` |
| `linear_three_phase`, `inverter`, `synchronous_generator` | `FMI_EMT_Ph3_CoSimulation` |
| `two_converters`, `two_generators` | `FMI_EMT_Ph3_TwoCoSimulations` |
| `microgrid_network`, `powerfactory_generator` candidates | Use a runner above if the export exposes a scalar or abc voltage/current port; these candidate FMUs have not been exported as Co-Simulation here. |

`FMI_EMT_CoSimulation` is the single-phase counterpart. It connects one FMI
voltage/current pair to an ideal 230 V RMS, 50 Hz boundary:

```sh
build-fmi/dpsim-fmi/FMI_EMT_CoSimulation /absolute/path/model.fmu [u i]
```

The optional names replace the default `u` and `i`; an optional stop time may
follow the names.

`FMI_EMT_Ph3_CoSimulation` connects a three-phase FMI 2.0 or FMI 3.0
Co-Simulation FMU to a finite R-L grid in DPsim. Its voltage inputs are in V,
current outputs in A, ordered abc. The default mapping is `u_a/u_b/u_c` to
`i_a/i_b/i_c`; the default sign is `-1`, meaning positive FMU current is
injected into the grid node. Set the sign to `+1` if positive output current
flows into the FMU.

Run with the default settings:

```sh
build-fmi/dpsim-fmi/FMI_EMT_Ph3_CoSimulation /absolute/path/model.fmu
```

Override port names, current sign, stop time and step:

```sh
build-fmi/dpsim-fmi/FMI_EMT_Ph3_CoSimulation model.fmu \
  vA vB vC iA iB iC -1 0.02 50e-6
```

`FMI_EMT_Ph3_TwoCoSimulations first.fmu second.fmu [stop_time_s [step_s]]`
places two independent Co-Simulation instances in parallel at the same PCC.
This covers the coupling structure used by the two-converter and two-generator
examples; each instance must be initialized from its own FMU archive.

The component uses explicit fixed-step coupling: the MNA solve stamps the last
completed FMU current, then advances the FMU once with the previous PCC
voltage held over the communication interval. Co-Simulation FMUs do not expose
their internal states to DPsim Newton iterations through FMI `doStep`, so the
example makes no implicit-convergence or rollback claim. The grid has 1 ohm
series resistance and 2 mH series inductance. Logs are written under `logs/`.

The checked-in inverter, synchronous-generator and RLC FMUs are Model Exchange
exports. This example is ready for matching Co-Simulation exports from those
models, but does not relabel or convert their existing archives. FMI 2.0
Co-Simulation exports require the FMI 2.0 `linux64`/`darwin64`/`win64` binary
folder convention; FMI 3.0 uses architecture-OS folder names. MATLAB R2025a
and Simulink are installed here, but the scoped export check could not build a
replacement FMU because Simulink Coder is unavailable in this installation.
Changing `ModelExchange` to `CoSimulation` in an archive is not a valid
conversion; each model must be re-exported with its source tool and Co-Simulation
interface enabled.
