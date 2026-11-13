# Three-phase EMT linear RL Co-Simulation stability benchmark

This is the ABC EMT version of the linear RL stability example. It compares
native DPsim, FMI 2.0 Model Exchange imported into DPsim's MNA solve, and two
FMI 2.0 Co-Simulation policies: single-pass explicit exchange and
rollback/Newton iteration.

| Setting | Value | Applied consistently |
|---|---:|---|
| Source | Balanced 50 Hz, 10 V RMS phase-to-neutral; 1 ms amplitude ramp | Same waveform in all variants |
| Grid resistance | 1 ohm per phase | Outside the FMU; same native/MNA resistance |
| FMU resistance | 0.1 ohm per phase | Inside the FMU or matching native resistor |
| Inductance | 20 microhenry per phase | Same branch equation in all variants |
| Initial current | `i_a=i_b=i_c=0 A` | Native inductors and both FMU states |
| Initial terminal voltage | `[0, 0, 0] V` at `t=0` | Source, ABC nodes, ME inputs, and CS input starts |
| Network/CS step | 50 microseconds nominal | One network interval per CS `doStep()` |
| Stop time | 20 ms (one 50 Hz cycle) | All variants |
| CS coupling | Single-pass explicit, held terminal voltage | Last accepted FMU current is stamped into MNA |
| Compensated CS check | Rollback/Newton at 50 microseconds | Replays the same FMU interval from a saved state as MNA iterates the interface voltage/current |

Each phase has identical independent parameters. For each phase the physical
equation is `L di/dt = v_s - (Rg + Rf)i`, with pole `-(Rg+Rf)/L = -55000
1/s`. The source waveform changes the forcing term but not the homogeneous
error recurrence of explicit CS:

```
e_(n+1) = q(h)e_n
q(h) = (1 + Rg/Rf) exp(-Rf h/L) - Rg/Rf
```

Therefore the predicted explicit-CS boundary remains 40.134139 microseconds
for each of phases a, b, and c. Balanced sinusoidal excitation makes the
example a genuine three-phase EMT waveform; its amplitude ramps from zero over
1 ms so the common initial voltage is exactly zero. This forcing does not
change the multiplier. All three phases remain independently modeled and are
logged in ABC order.

The single-pass policy is unstable above the analytical boundary. The same CS
FMU also advertises FMI state checkpoint support, which lets DPsim's
`RollbackNewton` mode restore the start of an interval and replay it while the
network solves for a self-consistent terminal voltage/current. At 50
microseconds this compensated run stays bounded (measured peak 12.856253 A
versus 12.856276 A for native and ME), while single-pass explicit coupling
diverges to about `7.61e61 A`. Rollback/Newton's max/RMS waveform differences
from native are 0.670/0.172 A. The feature requires an FMU with
`canGetAndSetFMUstate=true`; reducing the step below the boundary is the
simpler alternative when replay is unavailable.

## Build, export and run

From the repository root:

```sh
python3 dpsim-fmi/examples/linear_rl_stability_3ph/sources/build_reference_fmus.py
cmake --build build-fmi --target FMI_EMT_Ph3_LinearRLStability -j4
python3 dpsim-fmi/examples/linear_rl_stability_3ph/run_benchmark.py
```

The ME and CS archives are independently built from `sources/LinearRL3Ph.cpp`
and their corresponding extracted XML descriptions. The CS FMU integrates the
held-input RL equation exactly over each FMI communication interval. Results,
CSV traces, metrics, and plots are written under
`build-fmi/linear_rl_stability_3ph/results/`.

The installed Simulink FMU export path cannot generate archives in this
environment because Simulink Coder is not installed. These FMUs are the
reproducible C++ FMI 2.0 fallback; they were imported and executed by DPsim.

Open `LinearRL3Ph_Stability_Results.ipynb` from the repository root after
running the benchmark to see the settings, measured sweep, and updated figure. The
scalar derivation and the Simulink Coder export limitation for the earlier
single-phase version are documented in the accompanying theory PDF.
