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
