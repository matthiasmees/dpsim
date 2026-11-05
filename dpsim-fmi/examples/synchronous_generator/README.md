# Synchronous Generator

Build the C++ example with CMake in the repository build directory. `sources/` contains the final FMU, extracted `modelDescription.xml`, and Simulink model. `matlab/` contains saved reference results and model definitions where needed. The comparison notebook contains its comparison logic and writes generated runs and plots under `build-macos-check/fmi-examples/synchronous_generator/`.

This is the complex EMT validation case: the MATLAB-exported FMI 3.0 Model
Exchange FMU contains ten continuous generator/AVR/governor states and exposes
three-phase voltage inputs and current outputs. It has a native ARM64 macOS
binary and was rerun through DPsim at 50 microseconds after the FMI 2.0
backend changes. Open
`FMI_EMT_Ph3_Nonlinear_MATLAB_Comparison.ipynb` with the DPsim Python kernel
for the full MATLAB/DPsim waveform and error plots.
