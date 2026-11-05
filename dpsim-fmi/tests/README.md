# FMI importer regression tests

Build with `DPSIM_FMI_BUILD_TESTS=ON`; run `ctest --test-dir build-macos-check -R '^fmi_' --output-on-failure`. The C++ fixtures cover linear/nonlinear import, sampled controls, events, isolation and intentional error conditions. Build products are generated under the configured build directory.

The `fmi_multirate_distinct_fmu_clocks` regression couples two separate
sampled FMUs to one EMT bus with a 2 ms network step. Their internal event
periods are 1 ms and 0.5 ms; the test checks that each executes the expected
number of event-aligned substeps and matches an independent single-rate
trapezoidal recurrence. The maximum state-vector norm error must stay below
`2e-9` for both instances.
