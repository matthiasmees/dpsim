# FMI examples

The runnable cases are `linear_scalar`, `linear_three_phase`, `synchronous_generator`, `inverter`, and `two_generators`. `microgrid_network` and `powerfactory_generator` are additional external FMI 2.0 Model Exchange candidates with capability probes and interface notebooks; their upstream archives are retained unchanged, but they require unavailable host binaries and explicit electrical-interface adapters. Each case contains `sources/` (final FMU/XML/source model), `matlab/` (saved reference results where applicable), and a self-contained notebook. `common/` contains shared C++ headers.

Build with CMake as usual. Targets also include `FMI_Capabilities`, `FMI_EMT_LinearSSN`, `FMI_EMT_Ph3_LinearSSN`, `FMI_EMT_Ph3_Nonlinear`, `FMI_EMT_Ph3_Mixed`, `FMI_EMT_Ph3_TwoConverters`, `FMI_EMT_Ph3_TwoGenerators`, `FMI_MicrogridNetwork_Candidate`, and `FMI_PowerFactoryGenerator_Candidate`. Executables remain in the configured build directory; notebooks write their generated outputs under `build-macos-check/fmi-examples/`. Regression tests live in `dpsim-fmi/tests`.

The inverter retains the final successful V3 model. Its comparison notebook documents accepted errors and the limits of the refinement evidence.

The `two_generators` case also validates model-supplied nonzero operating-point initialization, independent control inputs, genuine shared-network coupling, singleton rejection and instance-order isolation against full-state MATLAB references.
