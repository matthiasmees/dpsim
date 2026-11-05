# Available FMU capability characterization

The repository currently contains one switched MathWorks inverter export and
separate linear/generator exports. No second independent MathWorks switched
converter export is available in the workspace; the existing two-generator
example uses a separate generator export instead of cloning an FMU.

The table below was obtained with the `FMI_Capabilities` executable, including
the runtime rollback replay probe.

| FMU | Event indicators | Get/set state | Serialize state | Rollback replay | Recommended mode |
|---|---:|---:|---:|---|---|
| `linear_scalar/sources/dpsim_rlc_me.fmu` | no | yes | yes | deterministic | `rollback_newton` |
| `linear_three_phase/sources/dpsim_rlc_3ph_me.fmu` | no | yes | yes | deterministic | `rollback_newton` |
| `synchronous_generator/sources/dpsim_sg_3ph_me.fmu` | no | yes | yes | deterministic | `rollback_newton` |
| `two_generators/sources/dpsim_sg_3ph_second.fmu` | no | yes | yes | deterministic | `rollback_newton` |
| `inverter/sources/GridTiedInverterPortedV3.fmu` | yes | yes | yes | nondeterministic | `single_pass_explicit` |

All listed exports advertise the FMI singleton-per-process restriction. The
two-generator test therefore relies on two independent exported model
identities, while the switched inverter is not claimed to support two genuine
instances in one process.
