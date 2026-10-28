# SSN GFM load-step and timing benchmark

Build from the repository root with:

```sh
cmake --build build --target EMT_SSN_GFM_Test -j2
```

Run `build/dpsim/examples/cxx/EMT_SSN_GFM_Test` from the repository root. The
default buffered real-time logger writes CSV and log files below
`logs/EMT_SSN_GFM_LoadStep_EMT/`. Use
`build/dpsim/examples/cxx/EMT_SSN_GFM_Test --csv-logger` for a conventional
per-step CSV logger comparison; it writes to the sibling `_CSVLogger`
directory. The scenario is 50 Hz, 50 µs, 3 s total, with a 3 kW resistive load
switched on at 1 s and off at 2 s.

The connection between the converter and load bus is an EMT SSN PiLine with
0.20 Ω series resistance, 0.1 mH series inductance, and 0.1 µF shunt
capacitance per phase. The power-flow PiLine uses the same parameters. The
GFM reactive-power reference includes the PiLine's calculated initial reactive
demand (about 26.6 var), so the voltage excitation integrator does not see a
persistent Q error at the initial operating point.
Both loads use the balanced fixed-impedance
`EMT::Ph3::SSN::ImpedanceLoad` model, parameterized from their P/Q targets at
nominal voltage.

The configured horizon is 3 s. This run records 60,001 measurements from 0
through 2.99995 s and executes 60,000 scheduler steps; the final-time boundary
is not an additional integration step. The electrical CSV starts with `time`,
then all GFM quantities, followed by load and remaining network quantities.

To plot the run's waveforms and cumulative stage timings, execute
`examples/Notebooks/Grids/EMT_SSN_GFM_LoadStep_Benchmark_Results.ipynb` after
running the benchmark.

## Output files

- `EMT_SSN_GFM_LoadStep_EMT.csv`: 60,001 sampled rows, including the initial
  point. Columns include all 17 local GFM states, the derived dq and abc
  interface/filter quantities, references, powers, PiLine voltage/current and
  R/L/C parameters, load impedance, and load terminal voltage and current. All
  GFM columns come first after time.
  Electrical units are in each column name.
- `per_step_task_timing.csv`: one row for every executed time step and one
  wall-time column for every scheduled task. It separates each component's
  pre-step and post-step work, the complete solver task, and data logging.
- `per_step_component_timing.csv`: one row per time step with separate GFM,
  SSN PiLine, base-load, and switched-load pre-step, post-step, and total
  component times, plus combined load, solver-task, and data-logger task time.
- `performance_summary.csv`: simulation wall-clock duration, modeled duration,
  time step, and executed-step count.
- `EMT_SSN_GFM_LoadStep_EMT_Solver.log`: DPsim's LU/refactorization and linear
  solve aggregate timing summary. The factorization timing covers the solver
  backend's combined L/U factor build; KLU does not expose separate L and U
  build timers. Those operations are inside the per-step solver task duration.

The load impedance columns are the per-phase equivalent R and X at the nominal
line-to-line RMS voltage, calculated from each load's configured P+jQ. Load P/Q
  columns report configured targets; terminal voltage and current are measured
  each step. The GFM model has no current or voltage limiter state in this
  implementation; virtual impedance
is a fixed configured value (zero in this islanded baseline), rather than a
dynamic state. Its controller states are the power filters, VSG omega/theta,
voltage magnitude/excitation state, voltage/current error integrators, delay
voltage states, capacitor phase voltages, and filter phase currents.

## Timing interpretation

The per-step task CSV is the step-resolved profile for component pre-step,
solver/integration, component post-step, and logging. DPsim's solver log adds
cumulative, average, and maximum matrix-refactorization and linear-solve times
for the same run. The active GFM state-space model requires matrix
recomputation at each step. The selected KLU backend reports each numerical
factorization as one operation, so it cannot separate internal L and U
construction. For an apples-to-apples comparison, keep the model, time step,
build configuration, solver backend, and logging columns unchanged.

The component profile groups each component's pre-step and post-step tasks so
GFM, PiLine, base-load, and switched-load time can be compared directly. It also
reports solver/integration and logger task time. The buffered logger reduces
per-step formatting and I/O work, then writes the CSV at shutdown; the notebook
compares logger-task, calculation, flush, and end-to-end times for both modes.
