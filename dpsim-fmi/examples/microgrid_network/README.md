# OpenModelica microgrid network candidate

This directory retains the upstream `grid.network.fmu` as a genuine external
FMI 2.0 Model Exchange candidate. It contains two inverter subsystems, three
filter networks and 48 continuous states. The FMU and its extracted
`modelDescription.xml` are stored in `sources/` together with the explicit
interface description in `sources/config.json`.

`FMI_MicrogridNetwork_Candidate.cpp` is an import and interface probe. It is
intentionally not presented as a working EMT network example: the upstream
FMU only ships `linux64` and `win64` binaries, so it cannot execute in the
current macOS/ARM64 build. In addition, its six terminal currents are FMI
inputs while the useful terminal voltages are local Modelica variables rather
than `causality="output"` variables. The current DPsim `NonlinearSSN` port
expects the opposite voltage-input/current-output direction, so no unsafe
mapping or XML modification is made.

The probe is still useful on Linux or Windows to verify that the archive is
valid and that DPsim can parse and instantiate the FMI 2.0 Model Exchange
instance. A future current-driven FMU adapter can promote the local terminal
voltages explicitly and then connect this candidate to a native MNA network.

Source: <https://github.com/upb-lea/openmodelica-microgrid-gym>
