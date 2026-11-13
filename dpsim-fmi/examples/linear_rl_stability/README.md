# Linear RL FMU stability benchmark

This benchmark compares one DC series circuit through three DPsim paths:

| Variant | DPsim representation | Coupling / integration |
|---|---|---|
| Native | EMT source, (R_g), (R_f), and (L) components | MNA with trapezoidal inductor integration |
| FMI ME | FMI 2.0 model exchange port for (R_f+L) | LinearSSN state included in DPsim MNA |
| FMI CS | FMI 2.0 co-simulation port for the same (R_f+L) | `doStep()` once per network step; voltage held between calls and last returned current stamped into MNA |

Parameters are (V_s=10\,\mathrm{V}), (R_g=1\,\Omega), (R_f=0.1\,\Omega), and (L=20\,\mu\mathrm{H}). All variants start from zero inductor current and use the same DC step at (t=0). The nominal network and CS communication step is 50 µs.

## Export the FMUs

The Simulink model and export script are in `matlab/`. With MATLAB R2025a,
Simulink, and the Simulink FMU Exporter installed, run from the repository root:

```sh
matlab -batch "addpath('dpsim-fmi/examples/linear_rl_stability/matlab'); export_linear_rl_fmus('build-fmi/linear_rl_stability/fmus')"
```

This creates the Simulink source model and attempts to create FMI 2.0 ME and CS
archives. The script copies successful exports and their model descriptions into
`sources/`; it never changes an ME archive into a CS archive by editing XML.

On the development host used for this benchmark, MATLAB R2025a and Simulink's
`exportToFMU` are installed, but the export stopped at code generation because
**Simulink Coder is not installed**. The Simulink source model is retained in
`sources/LinearRLStability.slx`. Until Simulink Coder is available, build the
standards-based fallback FMUs from the checked-in FMI 2.0 C++ source and XML:

```sh
python3 dpsim-fmi/examples/linear_rl_stability/sources/build_reference_fmus.py
```

This compiles and packages separate FMI 2.0 Model Exchange and Co-Simulation
archives under `sources/`. The ME archive exposes the differential state and
derivative. The CS archive advances that same state with the exact solution for
constant input over each `doStep` interval. The package contains the required
FMI 2.0 `modelDescription.xml` and Linux `binaries/linux64` implementation.

## Build and execute

```sh
cmake -S . -B build-fmi -DWITH_FMI=ON -DDPSIM_BUILD_EXAMPLES=ON -DDPSIM_FMI_BUILD_TESTS=ON
cmake --build build-fmi --target FMI_LinearRLStability -j4
python3 dpsim-fmi/examples/linear_rl_stability/run_benchmark.py
```

The runner executes all three variants at 50 µs and at the listed sweep points
around the analytical boundary. It writes raw traces, `summary.csv`,
`metrics.json`, and `linear_rl_stability.png` under
`build-fmi/linear_rl_stability/results/`. Generated logs and plots stay under
the configured build tree. The regression checks for the analytical boundary
are:

```sh
python3 dpsim-fmi/tests/test_linear_rl_stability.py
```

## Stability analysis

The physical circuit obeys

\[
L\,\frac{di}{dt}=V_s-(R_g+R_f)i,
\qquad
i(t)=\frac{V_s}{R_g+R_f}\left(1-e^{-(R_g+R_f)t/L}\right).
\]

For explicit CS, DPsim solves the network using the last completed FMU current.
The resulting PCC voltage is held during the FMU's next communication interval.
With an exact internal solution of the FMU's linear RL equation, this gives

\[
i_{k+1}=e^{-R_f h/L}i_k +
\left(1-e^{-R_f h/L}\right)\frac{V_s-R_g i_k}{R_f},
\qquad
q=\left(1+\frac{R_g}{R_f}\right)e^{-R_fh/L}-\frac{R_g}{R_f}.
\]

The fixed point is stable when \(|q|<1\). Since \(R_g>R_f\), the first
loss of stability is at \(q=-1\):

\[
h_{\mathrm{crit}}=\frac{L}{R_f}\ln\left(\frac{R_g+R_f}{R_g-R_f}\right)
=40.1341\,\mu\mathrm{s}.
\]

Thus the 50 µs explicit CS case is unstable even though the underlying RL
circuit is passive and linear. At 40.1341 µs the ideal recurrence is marginal;
the sweep includes points just below and above it. Native and ME remain stable
at 50 µs because their circuit state is integrated inside the DPsim MNA solve.
This result describes the selected non-iterative explicit coupling policy, not
every possible FMI CS master algorithm.

The checked-in fallback CS FMU advances the held-input RL equation with its
exact exponential solution on every `doStep`, which matches the stability
recurrence. The Simulink source model uses `ode15s`; the export attempt on this
host stopped before code generation. The runner compares every output trace
against the continuous analytical response and reports the measured results
alongside the discrete multiplier.

## Measured results on the development host

The fallback FMI 2.0 ME and CS FMUs were imported and simulated with the native
DPsim benchmark executable. At 50 µs, native and ME produced identical traces
and settled to 9.090909 A. Their peak absolute error against the continuous
solution was 3.246589 A at the first communication step; the final error after
5 ms was below (10^{-7}\,\mathrm{A}). The startup peak reflects trapezoidal
integration at a step longer than the 18.18 µs electrical time constant.

| CS step | Analytical multiplier \(q\) | Measured peak current (A) | Result over 5 ms |
|---:|---:|---:|---|
| 40.00 µs | −0.993962 | 18.13 A | Bounded, slow alternating convergence |
| 40.13 µs | −0.999814 | 18.18 A | Near marginal stability |
| 40.25 µs | −1.005212 | 26.32 A | Growing oscillation |
| 41.00 µs | −1.038880 | 874.86 A | Diverging |
| 50.00 µs | −1.433191 | (2.71\times10^{16}\) A | Rapid divergence |

The runner writes the full per-step native, ME, and CS measurements to
`build-fmi/linear_rl_stability/results/summary.csv`; the plot separates native
and ME from the diverging CS trace so both remain readable.
