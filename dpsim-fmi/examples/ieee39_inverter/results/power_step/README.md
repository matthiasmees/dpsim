# IEEE39 active-power command-step comparison

This result uses the three-phase EMT IEEE39 network and the averaged bus-30
inverter. At 0.1 s, the active-power command changes from 125 MW to 250 MW;
reactive power remains 80.881 Mvar. The three cases use the same IEEE39
power-flow initialization, inverter voltage, and current state. Their logged
initial three-phase current RMS is 237.313923 A in every case.

| Setting | Value |
|---|---:|
| Event | Active-power reference step, 125 → 250 MW |
| Event time | 0.1 s |
| Reactive-power reference | 80.881 Mvar, constant |
| Simulation stop | 0.3 s |
| Communication/network steps | 50 µs, 1 ms, 5 ms |
| Variants | Native DPsim, FMI Model Exchange, FMI Co-Simulation (explicit ZOH) |
| Converter current response | 5 ms first-order lag, shared model |
| Current limit | 1.20 pu, shared model |
| Initialization | Common IEEE39 power flow and equal inverter initial voltage/current |

The plot shows the one-cycle current RMS envelope around the event. The power
command creates an obvious response: at 50 µs, ME current rises from roughly
111 A immediately before the command to 196 A after settling. The corresponding
5 ms explicit-CS trace rises from roughly 128 A to 227 A.

This case does **not** demonstrate an event-triggered divergence. At 5 ms, the
explicit-CS current error relative to ME is already 130.94% before the event,
130.97% during the response, and 130.95% afterward. The step makes the physical
response visible, but does not materially increase the existing explicit
coupling error. At 1 ms, the corresponding errors are about 19.2% in all three
windows. The traces remain finite and settle; call this a large bounded
coupling error, not numerical divergence.

## Reproduce

From the DPsim repository root, build the FMU and runner, then run the three
modes at any one of the listed steps. For the 5 ms example:

```sh
cmake --build build-fmi --target FMI_IEEE39_AvgInverter -j4
for mode in native me cs-explicit; do
  mkdir -p build-fmi/ieee39-inverter/power-step-cosim-demo/dt_5000us/$mode
  (cd build-fmi/ieee39-inverter/power-step-cosim-demo/dt_5000us/$mode && \
    /path/to/build-fmi/dpsim-fmi/FMI_IEEE39_AvgInverter \
      "$mode" /path/to/build-fmi/ieee39-inverter-fmu/IEEE39AvgInverter.fmu \
      0.3 0.005 0.1 30 power-step 250 0)
done
```

Replace `/path/to` with the repository's absolute path. `error_summary.csv`
contains current RMS errors relative to ME in pre-event, response, and settled
windows. `initialization_and_response.csv` records matched starting current
and representative response samples. The detailed raw solver traces are in
the corresponding `build-fmi/ieee39-inverter/power-step-cosim-demo` folders.
