# FMI examples

The examples focus on reusable FMI import and co-simulation workflows:

- `linear_scalar` and `linear_three_phase` demonstrate Model Exchange FMUs in
  the EMT linear SSN path.
- `synchronous_generator` and `two_generators` demonstrate nonlinear
  three-phase Model Exchange, operating-point initialization and multiple
  coupled FMU instances.
- `cosimulation` contains scalar and three-phase Co-Simulation examples.
- `python` contains import quickstarts and explicit port-mapping guidance.

The `sources/` directories include the FMU archives and model descriptions
needed to run the examples. General examples use portable FMI inputs; generated
simulation logs and plots belong in the configured build directory. `common/`
contains shared C++ helpers. Regression fixtures and tests live in
`dpsim-fmi/tests`.
