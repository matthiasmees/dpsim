# PowerFactory/Modelica synchronous-generator candidate

This directory retains `genInterfacePQ.fmu` from the ALSETLab
PowerFactory/Modelica FMI tutorial. It is a useful candidate for the next
operating-point initialization step: the model contains a GENROU machine,
SCRX excitation and P/Q interface filters, with ten continuous states.

The C++ probe verifies the FMI 2.0 Model Exchange metadata and reports its
interface. It is not advertised as a runnable macOS EMT example because the
upstream archive ships only a Windows binary. Its `vi`/`vr` inputs and `p`/`q`
outputs also require an explicit complex-power or dq adapter; they must not be
silently treated as abc voltage/current channels. The FMU does not advertise
FMI state rollback, so it is not yet suitable for implicit multirate Newton
coupling.

Source: <https://github.com/ALSETLab/Modelica-FMI-PowerFactory-ModelImport-Tutorial>
