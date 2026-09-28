# MMC without an SSN state-space model

The new `CPS::EMT::Ph3::MMC_MNA` class represents the averaged electrical
power circuit of `SSN_MMC` directly as MNA equations. Only the two MMCs are
replaced in the copied P2P example. The AC network, transformers, loads,
DC line, grounding, parameters, operating points, and power setpoint steps
remain the same. The existing DC SSN components remain part of the network.

MMC means **modular multilevel converter**, SSN means **state-space nodal**,
and MNA means **modified nodal analysis**. Both examples use DPSim's global
MNA solver; the change concerns how the MMC contributes its equations.

## Interactive notebook

Open [EMT_MMC_MNA_Comparison.ipynb](EMT_MMC_MNA_Comparison.ipynb), select the
**Python (DPsim)** kernel, and execute the cells from top to bottom:

1. Load the helper functions.
2. Choose shared settings and update the example executables.
3. Run the original SSN example and display its plots.
4. Run the direct-MNA example and display its plots.
5. Compare both runs, plot differences, and display error metrics.

The notebook uses `EMT_MMC_MNA_Notebook.py` and the comparison script in this
same directory. Each experiment writes to a new directory under
`outputs/mmc_mna/notebook_runs/`. At the default 20 µs step, earlier full
10-second runs took approximately 3 min 44 s for SSN and 1 min 35 s for MNA,
plus build and plotting time. Actual wall times are measured and displayed.

The helper removes the active Python environment's library directory from
`LD_LIBRARY_PATH` and `LIBRARY_PATH` for C++ build and simulation subprocesses.
This avoids the Conda/system library conflict found on this machine. Other
library paths and the notebook kernel's own environment are preserved.

## Model and numerical implementation

The reference model uses a limited number of harmonics in rotating
coordinates. This modeling assumption is retained. The new implementation
adds the following unknowns to the global MNA system:

| Quantity | Count | Direct equation |
|---|---:|---|
| Difference currents `iDelta_d/q` | 2 | Inductive KVL with `Larm/2 + Lreactor` |
| Circulating currents `iSigma_z/d/q` | 3 | Inductive KVL with `Larm` |
| Capacitor voltage components `vCDelta_d/q/Zd/Zq`, `vCSigma_d/q/z` | 7 | Capacitive KCL with `Csm/N` |

Modulation and coordinate transformations appear as controlled voltage or
current sources. Terms containing ω, 2ω, and 3ω represent the corresponding
rotating coordinates. The class stamps these contributions directly into
the global MNA matrix. It does not construct A/B/C/D matrices, use an SSN
base class, or reduce the power circuit to a local state-space Norton source.

Each inductive or capacitive equation uses the companion form:

```text
L/(theta*h) * i_new - electrical voltage contributions_new
  = L/(theta*h) * i_old + (1-theta)/theta * voltage balance_old

C/(theta*h) * v_new - electrical current contributions_new
  = C/(theta*h) * v_old + (1-theta)/theta * current balance_old
```

`theta=0.5` gives trapezoidal integration. Stored capacitor voltages, inductor
currents, and controller integrals remain physically necessary.

PI controllers, the PLL, optional measurement filters, and the optional Padé
delay are updated using scalar predictor/corrector steps. Integrator values
are predicted for the following implicit network step. Sources are linearized
at the last known electrical operating point: explicitly predicting the
electrical quantities could prematurely select the other side of a modulation
limit during a step. Their immediate dependence on currents and terminal
voltages is included through incremental conductances. Central differences
of the five algebraic modulation outputs are taken with respect to the
twelve electrical quantities and the three d/q/DC terminal voltage quantities.
This linearizes nonlinear controlled sources in the MNA circuit without
constructing a continuous or discrete state-space model. Purely delayed
voltage feedforward is unstable in this inductive network. The solved
electrical quantities then correct the controller integrals. Nominal grid
rotation is calculated analytically from the step count; the PLL maintains
its own separate angle deviation.

The P/Vd feedforward used in the example retains the exact scalar filter
recursion with a held input, the 40 µs sampling interval, and the current
reference held between samples. As in the reference model, the sampling
period is rounded to an integer multiple of the EMT step.

The loaded initial state is obtained from a static Newton solution of the
KVL/KCL balances and active controller conditions. Its Jacobian is used only
during initialization and is then discarded. With DC voltage control, the
specified active power fixes the otherwise undetermined loaded operating
point. No SSN object is instantiated for initialization or simulation.

This method differs numerically from the reference model's local joint
linearization of the power circuit and control system. Finite time-step
errors are therefore expected. The results include a comparison with the
EMT step halved. The implemented optional measurement filters require
`h < tau`; the optional modulation delay requires `h < delay/6`. These
options are not enabled in the P2P example.

## Build

From the repository directory, with the DPSim dependencies available:

```bash
cmake -S . -B build-mmc-mna -DCMAKE_CXX_FLAGS=-O2
cmake --build build-mmc-mna --target EMT_SSN_MMC_Matlab_P2P_Results EMT_MMC_MNA_Matlab_P2P_Results -j 4
```

The separate build directory keeps the original build intact. On this
machine, the `fmt` and `spdlog` packages matching the original build are
under `/usr/local/lib64/cmake`. They were selected explicitly while using
the Conda environment:

```bash
cmake -S . -B build-mmc-mna -DCMAKE_CXX_FLAGS=-O2 \
  -Dfmt_DIR=/usr/local/lib64/cmake/fmt \
  -Dspdlog_DIR=/usr/local/lib64/cmake/spdlog \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DFETCHCONTENT_BASE_DIR="$PWD/_deps" \
  -DFETCH_GRID_DATA=OFF -DDPSIM_BUILD_DOC=OFF
```

All dependencies in `_deps` must already be present for this local offline
configuration.

## Simulation and plots

The four arguments of both examples are the EMT step, final time, theta,
and the number of steps between logged samples. Separate working directories
keep the results apart. For example, from the repository directory:

```bash
mkdir -p outputs/mmc_mna/reference_20us outputs/mmc_mna/mna_20us
(cd outputs/mmc_mna/reference_20us && ../../../build-mmc-mna/dpsim/examples/cxx/EMT_SSN_MMC_Matlab_P2P_Results 20e-6 10 0.5 10 > run.log 2>&1)
(cd outputs/mmc_mna/mna_20us && ../../../build-mmc-mna/dpsim/examples/cxx/EMT_MMC_MNA_Matlab_P2P_Results 20e-6 10 0.5 10 > run.log 2>&1)

python examples/Notebooks/Components/EMT_MMC_MNA_Comparison.py \
  outputs/mmc_mna/reference_20us/logs/EMT_SSN_MMC_Matlab_P2P_Results/EMT_SSN_MMC_Matlab_P2P_Results.csv \
  outputs/mmc_mna/mna_20us/logs/EMT_MMC_MNA_Matlab_P2P_Results/EMT_MMC_MNA_Matlab_P2P_Results.csv \
  --output outputs/mmc_mna/comparison_20us
```

These command-line examples reuse the named output directories. Use new
directory names to preserve earlier results, or use the notebook, which
creates a fresh directory for every run.

Python requires NumPy with `trapezoid` and Matplotlib. On this machine,
`/home/ssc/miniforge3/envs/dpsim-python/bin/python` provides these packages.

For the refinement comparison, repeat both simulations with
`10e-6 10 0.5 20` and pass `--reference-dt 10e-6 --mna-dt 10e-6` to the
comparison script. This retains the regular 200 µs logging interval. For
very fast transients, the last example argument can be set to `1`.

The analysis creates PNG and PDF figures for both stations, difference
plots, close-ups of the steps at 4/6/8 s, AC waveforms, and `metrics.csv`,
`metrics.json`, and `report.md`. Errors are reported both in physical units
and relative to fixed engineering bases. This avoids misleading relative
errors when reactive power or currents are near zero.
`--max-error-percent` can enforce a project-specific acceptance limit.

Both examples log the already incremented simulation time after
`sim.step()`. For every row after the initial row, the analysis subtracts
the corresponding EMT step from the timestamp. Unequal time grids are
linearly interpolated only within the shared time interval.

## Verified results (22 September 2026)

Both variants were simulated for 10 s, including all setpoint steps at
4/6/8 s, using EMT steps of 20 µs and 10 µs. All four runs completed normally
and produced finite values. Initial powers and DC currents agree up to
numerical rounding. Both comparisons meet the acceptance limit chosen
for this evaluation: 0.1% of each engineering base for all 16 main signals.

Maximum observed difference across both stations:

| Quantity | 20 µs | 10 µs |
|---|---:|---:|
| Active power | 0.618977 MW | 0.149966 MW |
| Reactive power | 0.736587 Mvar | 0.322217 Mvar |
| DC voltage | 25.4703 V | 4.24987 V |
| DC current | 0.027154 A | 0.006812 A |
| Capacitor energy | 103.226 J | 30.8074 J |

The largest normalized error decreases from **0.073659%** to **0.032222%**.
Remaining differences occur mainly during fast transients after the
setpoint steps. These figures refer to the 200 µs logging interval and
fixed engineering bases; they are not continuous-time error bounds or
instantaneous relative errors near zero.

The build, a `-Wall -Werror` check, the identity of the copied network
scenario, and the physical energy balance of the electrical couplings
were checked. Optional control modes, measurement filters, and Padé
delays were not independently validated by these P2P runs.

The following links point to the earlier archived comparison results;
the notebook generates new plots and reports for each experiment:

- [20 µs result report](../../../outputs/mmc_mna/comparison_20us/report.md)
- [10 µs result report](../../../outputs/mmc_mna/comparison_10us/report.md)
- [Time-step comparison](../../../outputs/mmc_mna/timestep_comparison.png)
