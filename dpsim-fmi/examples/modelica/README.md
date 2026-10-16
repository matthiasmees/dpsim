# Nonlinear MSL circuit portability benchmark

`NonlinearThreePhase.mo` adapts **Modelica Standard Library 4.0.0**
`Modelica.Electrical.Analog.Examples.ShowSaturatingInductor`. Each phase retains
its parallel saturating/linear inductors (2/1/0.5 H, nominal current 1 A), zero
initial currents and original 1.25 V, 1 rad/s waveform. Three original example
instances with phase offsets 0 and ±120 degrees form the independent native
reference. The port adaptation replaces only the source with a SignalVoltage
connected to an FMI voltage input and exposes the sum of the two branch currents.

The exported FMI 3.0 Model Exchange has six continuous states. `NonlinearSSN`
evaluates these native equations and condenses the implicit trapezoidal residual
into DPsim MNA. `FMI_EMT_Ph3_Prescribed` supplies an ideal voltage constraint and
logs source current for KCL checks. It contains no magnetic equations or fitted
matrices. The low source frequency follows the original standard example.
Nonlinear phase currents need not sum to zero.

## Reproduce on the SSD

Dependencies already cloned below `/Volumes/SSDMatthias/development/fmi-external-benchmarks`:

- [Modelica Standard Library](https://github.com/modelica/ModelicaStandardLibrary/tree/v4.0.0),
  tag v4.0.0, commit `6626538a236a34a390655b55d5e8c906bfe31772`.
- [Community native macOS arm64 OpenModelica distribution](https://github.com/yannrichet/openmodelica-macos-arm64),
  commit `c1ffa03cdf0c6240f048ebd5525e094a5cde0449`, reporting
  `v1.28.0-dev-213-g821cdd7bd3-cmake`. This is a community binary distribution,
  not an official signed OpenModelica macOS installer.

From the DPsim repository root:

```sh
env DEVELOPER_DIR=/Library/Developer/CommandLineTools \
  dpsim-python/bin/python dpsim-fmi/examples/modelica/build_benchmark.py \
  --omc /Volumes/SSDMatthias/development/fmi-external-benchmarks/openmodelica-macos-arm64/bin/omc \
  --msl /Volumes/SSDMatthias/development/fmi-external-benchmarks/ModelicaStandardLibrary \
  --output /Volumes/SSDMatthias/development/dpsim/outputs/fmi-nonlinear/modelica

env DEVELOPER_DIR=/Library/Developer/CommandLineTools \
  cmake --build build-macos-check --parallel 3 --target FMI_EMT_Ph3_Prescribed
```

Open `dpsim-fmi/examples/FMI_EMT_Ph3_Modelica_Comparison.ipynb` and run all cells.
It reuses the saved native RK4 reference (100 µs), reruns imported cases at
10 ms, 5 ms, 200 µs and 100 µs, checks all three currents and source-current KCL,
and writes plots and metrics under `outputs/fmi-nonlinear/modelica`.
The native simulation stops at the original example's 6.2832 s; comparison
runs stop at 6.28 s, an integer number of each chosen DPsim timestep.

The initial verified runs give peak current errors of 58.1 µA at 10 ms,
14.9 µA at 5 ms and about 0.52 µA at 200 µs. RMS error falls by 3.99 when the
coarse step is halved. Fine errors are limited by DPsim's six-decimal CSV format.
This establishes compatibility for this smooth nonlinear circuit, not universal
support for Modelica exports. See `../../NONLINEAR.md` for supported FMI features.

MSL sources remain in their separate checkout with their original license.
The wrapper imports those library components; it does not redistribute their
implementation. OpenModelica's warning about unpatched standard ModelicaServices
is retained in the build log; this electrical benchmark uses its numerical
constants and no service-dependent animation or external file functions.
