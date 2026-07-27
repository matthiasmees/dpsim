# SSN MMC Phase 3 Handover

## Status

Phase 3 is partially implemented and intentionally uncommitted. Stages 3B
through 3E and the scalar-capacitor part of 3F pass against the real EMT MMC.
The scalar PiLine 3F case is the first failing layer, so Phase 3 completion is
not claimed.

## MMC inspection and external command

`EMT::Ph3::SSN_MMC` has terminals `[ABC, dc+, dc-]` and SSN vectors
`u=[va,vb,vc,vdcp,vdcn]`, `y=[ia,ib,ic,idcp,idcn]`. AC current is positive out
of the converter. `idc=3*iSigma_z` is positive from dc+ through the converter
to dc-; nodal injections are `idcp=-idc`, `idcn=+idc`.

The internal output-current controller produces `vMDelta_d/q` in peak phase
volts in the nominal positive-sequence dq frame. The plant consumes it in
`evaluateStateDerivative()` through:

```text
mDelta_d = -2 vMDelta_d / (vdcp-vdcn)
mDelta_q = -2 vMDelta_q / (vdcp-vdcn)
```

The new `ControlSource` defaults to `InternalControllers`. Only the explicit
`ExternalDifferentialVoltage` selection substitutes the external 2x1 command.
Internal circulating-current, zero-sequence, capacitor, and electrical states
are unchanged. External mode holds the now-unused internal output-current
integrators. Non-finite commands throw.

Existing power definitions are:

```text
Pac = 1.5(vd*id+vq*iq)
Qac = 1.5(-vd*iq+vq*id)
Pdc = (vdcp-vdcn)*idc
power_balance_error = Pac-Pdc
```

This instantaneous diagnostic contains losses and stored-energy rate; it is
not asserted to be zero during transients. No validated plant blocking exists.
Station `Blocked` therefore means controller-command inhibition and restoration
of the plant's default internal control source.

## Station architecture

`EMT::Ph3::SSN_MMCStation` is a `SimSignalComp`. It references exactly one
existing `SSN_MMC`; it owns no terminal, node, MNA stamp, or electrical state.
It coordinates:

- `Signal::ExternallyAngledDQAdapter`;
- three `DQSymSecondOrderFilter` instances for P, Q, and Vdc;
- active-power, DC-voltage, and reactive-power `DQSymOuterController` blocks;
- `DQSymCurrentController`;
- `DQSymModulation`.

Modes are `ActivePowerReactivePower` and `DCVoltageReactivePower`. States are
`Blocked`, `Ready`, and `Enabled`. Rejection diagnostics distinguish not-ready,
non-finite measurements/references, excessive enable error, command mismatch,
and invalid plant mode.

## Transform and bases

The angle and angular frequency are externally supplied in radians and rad/s.
No DQsym-equivalent PLL is claimed. `ExternallyAngledDQAdapter` uses:

- phase sequence a-b-c;
- positive increasing angle;
- phase-a cosine maximum at `theta=0`;
- amplitude-invariant `2/3` abc-to-dq scaling;
- `va=d*cos(theta)-q*sin(theta)` inverse transform;
- no zero sequence;
- current positive out of the converter.

Controller bases are:

```text
Vbase_ac = sqrt(2/3)*Vac_line-line_RMS
Ibase_ac = (2/3)*Pnom/Vbase_ac
Vbase_dc = Vdc_pole-to-pole_nominal
Pbase = Pnom
```

The DQsym modulation output is dimensionless. The plant command conversion is:

```text
vMDelta_dq [V peak phase] = vdq_ref_pu * Vbase_ac
mDelta_dq = -2*vMDelta_dq/Vdc
```

At nominal Vdc the DQsym modulation magnitude equals
`2*|vMDelta_dq|/Vdc`; the plant applies the required negative insertion-index
sign internally.

## Signal flow and timing

Declared tasks are:

1. measurements and abc/dq/P/Q calculation;
1. filters and selected outer loops;
1. dq current controller;
1. modulation and plant-command conversion.

Each downstream task declares dependencies on upstream modified attributes.
The command task modifies the plant's dynamic external command attribute.
The electrical solution of step `k` feeds these tasks; the held command is
consumed by the MMC SSN model rebuild for step `k+1`. This avoids an
algebraic loop through the current-step MNA solution.

Active-power mode follows the pinned cascade:

```text
Pmeas-Pref -> active PI -> Vdc_ref
(Vdc_meas-Vdc_ref)/Vdc_nom -> DC-voltage PI -> id_ref
Qref-Qmeas -> reactive PI -> iq_ref
```

DC-voltage mode bypasses the active PI and uses the external Vdc reference.
At the station/plant boundary, positive generated Q requires negative `iq`, so
the reactive PI scalar output is negated. In direct DC-voltage mode the
validated plant requires negative `id` for DC overvoltage, so the DQsym
`Vdc_meas-Vdc_ref` PI output is also negated. The active-power cascade retains
its independently validated direct orientation.

## Initialization and enable

Initialization sets externally supplied angle/frequency, dq measurements,
P/Q/Vdc filter states (`x1=measurement`, `x2=0`), and references to the
operating point. Outer-loop integrator states reproduce initial Vdc, id, and
iq references. Current-loop states reproduce the initial converter command
after subtracting PCC feedforward and R/omega-L terms. `Ready` holds that
command.

Enable requires finite measurements/references, active and reactive errors
within the requested per-unit tolerance, and held-command agreement within the
requested voltage tolerance. No reset, tracking PI, ramp, precharge, or
undocumented bumpless-transfer equation is implemented.

## Parameters and limits

Pinned DQsym gains and limits remain:

- active PI `Kp=0.5/3`, `Ki=1`, Vdc reference `[0.8,1.2] pu`;
- Vdc PI `Kp=4`, `Ki=100`, id reference `[-2,2] pu`;
- reactive PI `Kp=0.5/3`, `Ki=1`, iq reference `[-0.25,0.25] pu`;
- current PI `Kp=0.6`, `Ki=6`, axis voltage `[-2,2] pu`;
- measurement filter `f0=1000 Hz`, `zeta=1`;
- modulation axes `[-2,2]`;
- Vdc normalization clamp `[0.75,1.25]*Vdc_nom`.

Axis saturation and Phase 2 conditional-integration anti-windup are preserved.
No vector limiter was added.

Every PI now exposes unsaturated output plus independent upper- and
lower-saturation flags. Both current-controller axes and both modulation axes
expose the same directional diagnostics; the station republishes the
directional outer-loop, current-loop, and modulation flags.

## Validation

Passing targets:

```text
DQSym_Controller_Validation
MMC_ExternalCommand_Validation
MMCStation_OpenLoop_Validation
MMCStation_CurrentControl_Validation
MMCStation_OuterControl_Validation
EMT_SSN_MMCStation_ClosedLoop
EMT_SSN_MMCStation_DCCapacitor
```

The external-command validation covers default-source invariance, explicit
selection, non-finite rejection, a physically initialized derivative response,
transform round trip and P/Q signs, Ready-to-Enabled behavior, and task
dependency declarations. The station targets cover signal-level command
scaling, current-loop connection/saturation, and both outer-loop hierarchies.

The full stiff-source EMT target covers the initialized internal/external
transition, positive and negative id/iq steps, return to the operating point,
saturation entry/recovery, reactive-power closure, and the complete
`ActivePowerReactivePower` cascade. The capacitor target closes
`DCVoltageReactivePower` around a scalar 1.5 mF DC capacitor with positive and
negative 10 A disturbances and a reactive-power step.

Default MMC regressions built and ran after these changes:

```text
EMT_SSN_MMC_ActivePowerControl
EMT_SSN_MMC_Cap_VoltageControl
EMT_SSN_MMC_PiLine
EMT_SSN_MMC_Resistor_NegRefStep
EMT_SSN_MMC_Trafo
```

Selected numerical results (2026-07-27):

- initialized command jump `1.53503 V`; first enabled current-vector norm
  `0.026149 A`;
- stiff-case AC KCL `1.71e-13 A`; dc+ and dc- KCL `1.14e-13 A`;
- stiff-case stored energy `14.386969..14.528057 MJ`;
- id peaks `[-0.008648,+0.008502] pu`, iq peaks
  `[-0.009298,+0.009044] pu`, with final id/iq below `4.1e-5 pu`;
- reactive closure Q peaks `[-0.001520,+0.001683] pu`, maximum active-power
  disturbance `0.001299 pu`;
- active closure P peaks `[-0.005642,+0.005265] pu`, final P
  `-4.74e-5 pu`;
- capacitor Vdc `439895.824..440115.497 V`, final `440000.261 V`, KCL
  `1.08e-8 A`, stored energy `14.514657..14.523677 MJ`;
- legacy PiLine `Vdc=440000.001258 V`, dc+ KCL `1.87e-10 A`, dc- KCL
  `2.53e-10 A`.

The full-EMT tests currently bound KCL and energy and use the plant's
instantaneous `Pac-Pdc` diagnostic. A derived transient balance that separates
arm losses and `dE/dt` has not yet met the requested acceptance evidence.

## Outstanding acceptance criteria

`EMT_SSN_MMCStation_DCPiLine` is implemented but fails. With the validated
Phase-1 stiff remote-pole boundary and symmetric `+/-100 V` pole-to-pole
disturbances, its first failed invariant is local DC KCL after the state has
already diverged. At theta `0.55` the observed extrema were:

```text
local Vdc = [-9.440e18, +5.115e18] V
remote Vdc = [439900, 440100] V
line current = 9.441e16 A
stored energy max = 3.328e33 J
local KCL residuals = [14432,14848] A
```

The earlier theta `0.5` run diverged on the same order (`Vdc` around
`1e17 V`, line current around `1e16 A`). Therefore slight theta damping does
not identify this as the known trapezoidal alternating mode. The earliest
unresolved layer is the closed Vdc-controller/MMC/PiLine interaction. Gains
were not retuned and extra damping was not added.

Generalized-theta support was added to EMT SSN discretization with separate
old/new input weights:

```text
Ad = (I-theta*dt*A)^-1 (I+(1-theta)*dt*A)
Bd,new = (I-theta*dt*A)^-1 theta*dt*B
Bd,old = (I-theta*dt*A)^-1 (1-theta)*dt*B
```

The default remains theta `0.5`; setters reject values outside `[0.5,1]`.
Only the failing station PiLine study explicitly selects `0.55`.

Still required before completion:

- diagnose the closed Vdc/MMC/PiLine instability without controller retuning;
- validate directional saturation flags in full EMT for every limited axis;
- derive and assert transient AC/DC power balance including losses and
  stored-energy derivative;
- quantify settling time, overshoot, and high-frequency oscillation with
  explicit windowed assertions.

Do not advance to a point-to-point link until these Phase 3 items pass.

## Upstream boundary

Exact MathWorks SPS PLL, `abc to dq0`, and `dq0 to abc` internals remain
unavailable. No class is called or claimed to be an exact DQsym PLL or SPS
transform. No undocumented PI reset or tracking semantics are implemented.
