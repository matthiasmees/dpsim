# Converter import extension selected after exporter probes

The controlled-generator benchmark and21 regressions remain the baseline.
The next work must preserve the additional converter model's plant and controls.

## Verified exporter evidence

- Native MATLAB exporter: GridTiedInverter with the original physical network
  and controllers, changing only Simscape UseLocalSolver to off, exports FMI3 ME
  with64 continuous states and40 event indicators. Original local-solver-on model
  exported zero continuous states and cannot establish DPsim plant integration.
- DistributionSTATCOMAverageModel with continuous powergui retains its sampled
  controls. Its network has zero physical switches,42 S-function states; the
  generated complete model has48 continuous states and5 event indicators.
- FMI Kit S-function target preserves non-inlined network source compilation.
  It exports FMI2 ME, not FMI3 ME. Its simple continuous integrator builds and
  simulates correctly on this host after local ARM/R2025a compatibility fixes.
  The full STATCOM library builds without MATLAB-runtime dependencies, but its
  instantiation currently crashes; that must be fixed and verified before use.

## Runtime changes, explicit opt-in

Files:

- include/dpsim-fmi/ModelExchange.h: optional settings for endpoint state events
  and continuous-state reset handling; expose committed continuous state.
- src/ModelExchange.cpp: load/evaluate FMI event indicators, detect crossings
  only after a converged full MNA step, and combine state/time/input events in
  one event-mode update. Refresh actual continuous states after event iteration.
- src/NonlinearSSN.cpp: publish the refreshed post-event state and history.
  Retain terminal-current consistency checks: a reset that changes electrical
  current needs an instantaneous network solve and must still fail explicitly.
- examples/FMI_EMT_Ph3_Prescribed.cpp: JSON flags enabling the above policies for
  the converter benchmark, while previous configs keep their existing behavior.
- tests: a small independent phase-wrap fixture verifies retained physical
  current, reset-state publication, no duplicate event, and rejection by default.

Fixed-step endpoint event handling quantizes event location to at most one DPsim
step. It is not adaptive root localization. The converter comparison must
therefore demonstrate the requested physical-signal accuracy under timestep
refinement; merely running the model does not establish success. All existing
scalar/abc linear, nonlinear and native SSN behavior must remain valid.

## Exporter fallback if the native64-state converter is impractical

The switchless48-state STATCOM is a suitable original averaged converter example.
If FMI Kit is used, preserve its source model and generated equations. A clearly
identified FMI3 ABI adapter may delegate ME functions to the generated FMI2 code;
it must not call doStep or embed a plant integrator. Such an adapter must preserve
state/event metadata, types and value references, and be independently compared
with the original FMI2 and MATLAB simulations. Renaming an FMI2 XML version or
reporting an event-only replay as a coupled simulation is not acceptable.

No general native equation/MNA reconstruction, VILLAS coupling or dynamic-phasor
conversion is part of these changes. All new files/results stay on the SSD.
