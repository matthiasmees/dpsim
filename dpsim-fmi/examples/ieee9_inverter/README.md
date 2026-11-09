# IEEE9 inverter: native DPsim, FMI ME, and FMI CS

This case replaces the IEEE9 bus-3 synchronous generator with the same
averaged inverter in native DPsim, FMI Model Exchange, and FMI Co-Simulation.
Every case initializes its EMT network from the solved IEEE9 powerflow. The
inverter starts at the bus-3 powerflow voltage and its configured 85 MW /
-10.9 Mvar operating point. All three methods have identical initial PCC
voltage and inverter current. A 25 MW resistive load is switched onto bus 8 at
the requested event time.

## Topology and model checks

The mapping agrees with the standard IEEE9 case: lines 5-4, 6-4, 7-5, 9-6,
7-8, and 8-9; transformers 1-4, 2-7, and 3-9. The powerflow and EMT models use
the same bus numbering and connections. The bus-3 generator is represented in
the powerflow by its equivalent P/Q injection and replaced in EMT by the
inverter. Lines, transformers, loads, and passive-current signs match the
standard case. No topology remapping error was found.

The long-run failure in the earlier case came from omitting the IEEE9
excitation systems and turbine governors on the remaining synchronous
generators. That left a sustained generation deficit after the load step. The
EMT generators now include the reference IEEE9 DC1-style exciters and TGOV1
governors. This is needed for a meaningful 100-second test.

The three inverter realizations share the same 5 ms current-response model,
current limit, and initial P/Q current. Native DPsim and ME use trapezoidal
integration. The CS FMU retains its own exact exponential current-state update.
For the opt-in rollback/Newton CS case, each outer communication interval is
replayed using four FMU substeps with midpoint voltage samples linearly
interpolated from the previous PCC voltage to the candidate endpoint. This
keeps the CS method's own integrator while avoiding a full-step jump to the
candidate endpoint. Single-pass explicit CS remains available as
`cs-explicit`, but it is not the stable closed-loop comparison mode here.

## Validated results

All methods completed a fresh 100-second run at 100 µs, with the load step at
50 s and 1,000,000 samples per run. They started from exactly the same recorded
bus-3 voltage and inverter current. ME agreed with native DPsim to 0.000077%
PCC-voltage RMS and 0.0000014% current RMS over the full run. CS remained
finite and physically bounded through 100 s; its full-run waveform errors
were 7.46% in PCC voltage and 1.85% in inverter current. Its 100 ms rolling
voltage error rose gradually to 12.84% near the end, while current error
reached 3.01%. This is accumulated trajectory/phase drift, not numerical
runaway: the CS inverter current magnitude stayed close to native (about
0.02% mean difference in the saved trace), and the PCC voltage space-vector
magnitude difference was about 0.4% near the end. The dominant discrepancy is
the accumulated difference between the CS exact-exponential state update and
the native/ME trapezoidal update in a long closed-loop network simulation.

A separate fresh 5-second 100 µs run with the load step at 2.5 s gives a useful
short-horizon check: CS versus native is 0.313% PCC-voltage RMS and 0.0768%
current RMS; ME versus native is 0.0000039% and 0.000000087%, respectively.
The longer test confirms that CS is stable for 100 µs, but its waveform does
not remain identical to native DPsim over 100 seconds. Interpret these as
same-setup implementation comparisons, not independent physical-accuracy
measurements.

The 100-second run summary and compact traces are in
`build-fmi/ieee9-inverter/ieee9-corrected-100s-100us/metrics.json` and its
per-method `*_trace.csv` and `*_rolling_error.csv` files. The 5-second short
check is in `build-fmi/ieee9-inverter/ieee9-cs-substep4-5s/metrics.json`.

## Reproduce

From the repository root:

```bash
cmake --build build-fmi --target FMI_IEEE9_AvgInverter -j4
bash dpsim-fmi/examples/ieee9_inverter/sources/build_fmu.sh build-fmi/ieee9-inverter-fmu
python3 dpsim-fmi/examples/ieee9_inverter/run_ieee9_dpsim_native_me_cs_sweep.py \
  --steps-us 100 --duration-s 100 --load-step-s 50 \
  --output-name ieee9-corrected-100s-100us
```

The focused notebook `IEEE9_CoSimulation_RollbackNewton_Validation.ipynb`
visualizes the corrected 5-second, 100 µs check. The earlier full-component
10-second notebook contains data from the old single-pass CS case and should
not be used as the corrected comparison.
