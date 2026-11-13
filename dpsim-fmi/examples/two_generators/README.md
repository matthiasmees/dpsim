# Coupled generator FMUs

`FMI_EMT_Ph3_TwoGenerators.cpp` connects two separately instantiated nonlinear
three-phase FMUs to one shared EMT bus. Build the `FMI_EMT_Ph3_TwoGenerators`
target with `WITH_FMI=ON` and `DPSIM_BUILD_EXAMPLES=ON`. The archives and model
descriptions are in `sources/`; simulation logs are written under the configured
build directory.
