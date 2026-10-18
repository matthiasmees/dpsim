# FMI importer regression tests

Build with `DPSIM_FMI_BUILD_TESTS=ON`; run `ctest --test-dir build-macos-check -R '^fmi_' --output-on-failure`. The C++ fixtures cover linear/nonlinear import, sampled controls, events, isolation and intentional error conditions. Build products are generated under the configured build directory.
