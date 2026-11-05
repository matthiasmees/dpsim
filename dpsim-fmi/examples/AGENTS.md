# FMI example layout

For future development, keep examples small and retain only final working versions:

- One C++ example and one self-contained comparison notebook per case.
- Store the final FMU, extracted modelDescription.xml and source model in each case's sources/ directory.
- Keep meaningful MATLAB reference results/model definitions in matlab/.
- Generate binaries, simulation logs, plots and derived metrics under the configured build directory.
- Keep regression tests in dpsim-fmi/tests, not examples.
- Do not accumulate analysis utilities, intermediate exports, caches or failed development steps in examples.

Write example C++ in the native DPsim example style: clear parameter names and visible sections for parameters, nodes, components, topology, initialization, logging, events and simulation as applicable. Keep electrical connections easy to follow; explain FMU port units/signs and special initialization choices.
