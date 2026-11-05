# Python FMU import examples

These notebooks use DPsim's `dpsimpy.fmi` bindings and the retained inverter
FMU in `../inverter/sources/`. Build DPsim with `WITH_FMI=ON` and
`WITH_PYBIND=ON`, then open the notebooks from the repository checkout. They
look for the `build-fmi` Python extension and keep simulation output there.

- [Import_FMU_Quickstart.ipynb](Import_FMU_Quickstart.ipynb) shows automatic
  three-phase port discovery and the `import_fmu(path)` / `connect(node)` flow.
- [Import_FMU_Mapping_and_Settings.ipynb](Import_FMU_Mapping_and_Settings.ipynb)
  shows explicit port mapping, current polarity, settings, and coupling logs.

The inverter's exported runtime does not reproduce a rollback/replay interval
deterministically. Both short examples therefore select
`SinglePassExplicit`; their results demonstrate import and network connection,
not implicit macro-step convergence.

For FMI Co-Simulation exports, `dpsimpy.fmi.import_cosimulation_fmu(path, ... )`
creates an explicit one- or three-channel electrical adapter. Configure the
voltage input names, current output names and current sign to match the FMU.
The importer exposes `dpsimpy.fmi.CoSimulation` directly for applications that
want to control communication steps themselves. The checked-in inverter
archive is Model Exchange and cannot be passed to this Co-Simulation API.
