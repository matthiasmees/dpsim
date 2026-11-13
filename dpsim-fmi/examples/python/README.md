# Python FMU import examples

These notebooks use DPsim's `dpsimpy.fmi` bindings and the portable three-phase
RLC FMU in `../linear_three_phase/sources/`. Build DPsim with `WITH_FMI=ON` and
`WITH_PYBIND=ON`, then open the notebooks from the repository checkout. They
look for the `build-fmi` Python extension and keep simulation output there.

- [Import_FMU_Quickstart.ipynb](Import_FMU_Quickstart.ipynb) shows automatic
  three-phase port discovery and the `import_fmu(path)` / `connect(node)` flow.
- [Import_FMU_Mapping_and_Settings.ipynb](Import_FMU_Mapping_and_Settings.ipynb)
  shows explicit port mapping, current polarity, settings, and coupling logs.

The examples select `SinglePassExplicit` for reproducible fixed-step coupling;
they demonstrate import and network connection, not implicit macro-step
convergence.

For FMI Co-Simulation exports, `dpsimpy.fmi.import_cosimulation_fmu(path, ... )`
creates an explicit one- or three-channel electrical adapter. Configure the
voltage input names, current output names and current sign to match the FMU.
The importer exposes `dpsimpy.fmi.CoSimulation` directly for applications that
want to control communication steps themselves. The checked-in three-phase
archive is Model Exchange and cannot be passed to this Co-Simulation API.
