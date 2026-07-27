# SSN MMC Phase 2 Controller Design

## Authoritative sources

The controller source is TU Delft's DQsym repository:

```text
https://github.com/control-protection-grids-tudelft/DP
immutable revision: 36f47adb7c10c8b109e25a2180f20e08c645cfc4
```

Immutable source links use:

```text
https://github.com/control-protection-grids-tudelft/DP/tree/36f47adb7c10c8b109e25a2180f20e08c645cfc4
```

The extracted equations below come from:

- `DP_v1.0/P2PHVDCMMc/P2PHVDCMMCcode.m`;
- `P2PHVDCMMC.slx`, expanded XML systems `27380` (active-power loop),
  `27486` (DC-voltage loop), `27749` (reactive-power loop), `27402`
  (current loop), `27454/27455` (PI with anti-windup), and `27562`
  (PLL, transforms, filtering and P/Q measurement).

The Simulink PLL and second-order filter are references to MathWorks
Specialized Power Systems library blocks. DQsym stores their complete mask
parameters but not their internal implementation.

## Verified controller cascade

```text
P_ref - P_meas
  DQsym uses e_P = P_meas - P_ref
       -> active-power PI -> Vdc_ref [pu * Vnom_dc]

Vdc_ref, Vdc_meas
  e_Vdc = (Vdc_meas - Vdc_ref) / Vnom_dc
       -> DC-voltage PI -> id_ref [pu]

Qref, Qmeas
  e_Q = Qref - Qmeas
       -> reactive-power PI -> iq_ref [pu]

id_ref-id, iq_ref-iq
       -> axis PI controllers
       + PCC-voltage feedforward
       + R and omega-L decoupling
       -> vd_conv, vq_conv [pu]
       -> DC-voltage normalization and axis saturation
       -> inverse dq transform -> three-phase voltage command
```

Thus DQsym does **not** map active-power error directly to `id_ref`.
`P_ref` produces `Vdc_ref`, which is consumed by the DC-voltage loop.

## Units, bases and signs

DQsym defines:

```text
Vb_ac = sqrt(2/3) Vnom_sec
Ib_ac = (2/3) Pnom / Vb_ac
Vb_dc = Vnom_dc
Ib_dc = Pnom / Vb_dc
```

The controller uses amplitude-normalized dq quantities. Current is positive
out of the converter. DQsym annotations state:

```text
P = vd id + vq iq
Q = vq id - vd iq
positive id -> positive generated active power
negative iq -> positive generated reactive power
```

With the PLL aligned to the PCC voltage, `vq=0`; hence
`id=P/vd` and `iq=-Q/vd`.

## PI and enable semantics

For the current-loop PI subsystem:

```text
u_raw[k] = Kp e[k] + x[k]
u[k] = clamp(u_raw[k], u_min, u_max)
```

The integrator uses the Simulink trapezoidal discrete integrator. Conditional
integration is applied:

```text
hold = (u_raw >= u_max and e > 0)
    or (u_raw <= u_min and e < 0)

g[k] = hold ? 0 : Ki e[k]
x[k] = clamp(x[k-1] + Ts/2 (g[k-1] + g[k]), u_min, u_max)
```

DQsym describes this explicitly as suspending integration when saturation
and controller error have the same sign.

The outer-loop enable switch substitutes the reference for the measurement
when disabled. Therefore the error becomes exactly zero and the integrator
retains its previous state:

```text
measurement_selected = enable ? measurement : reference
```

This is a state-hold disable, not a reset. It is bumpless only when the held
integrator state already equals the desired enabled output. DQsym provides no
integrator-reset input, tracking input, reset equation, or general bumpless
transfer algorithm. Phase 2 therefore does not invent those features.

Feedforward is not an input to the DQsym PI block. It is summed externally in
the current controller.

## Measurement filter

DQsym explicitly instantiates a Specialized Power Systems
`Second-Order Filter`, type `Lowpass`, for both dq voltage and dq current:

```text
f0 = Fn_filter = 1000 Hz
zeta = Zeta_filter = 1
Ts = 40 us
Initialize = off
```

The continuous low-pass state model represented in DPsim is:

```text
omega0 = 2 pi f0
xdot_1 = x_2
xdot_2 = omega0^2 (u-x_1) - 2 zeta omega0 x_2
y = x_1
```

It is discretized with DPsim's trapezoidal state-space rule. Explicit
steady-state initialization uses `x_1=u_initial`, `x_2=0`, preventing a
first-step transient. This is not replaced by a first-order approximation.

## Outer loops

Active-power loop:

```text
e_P = enable ? (Pmeas-Pref) : 0
vdc_ref_pu = clamp(Kp_P e_P + integral(Ki_P e_P), 0.8, 1.2)
Vdc_ref = Vnom_dc vdc_ref_pu
Kp_P=0.5/3, Ki_P=1, initial integrator=1 pu
```

DC-voltage loop:

```text
e_Vdc = enable ? (Vdc_meas-Vdc_ref)/Vnom_dc : 0
id_ref = clamp(Kp_Vdc e_Vdc + integral(Ki_Vdc e_Vdc), -2, 2)
Kp_Vdc=4, Ki_Vdc=100
```

Reactive-power loop:

```text
e_Q = enable ? (Qref-Qmeas) : 0
iq_ref = clamp(Kp_Q e_Q + integral(Ki_Q e_Q), -0.25, 0.25)
Kp_Q=0.5/3, Ki_Q=1
```

## Current controller and modulation

The verified current errors and voltage equations are:

```text
e_d = id_ref-id
e_q = iq_ref-iq

vd_conv = vd + Rff id - omega_pu Lff iq + PI_d(e_d)
vq_conv = vq + Rff iq + omega_pu Lff id + PI_q(e_q)

Rff = Rarm_pu/2
Lff = Larm_pu/2
omega_pu = frequency_hz/Fnom
```

DQsym uses `Kp_I=0.6`, `Ki_I=6` and axis limits `[-2,2] pu`.

Before inverse transformation:

```text
Vdc_limited = clamp(mean(Vdc), 0.75 Vnom_dc, 1.25 Vnom_dc)
kdc = Vnom_dc / Vdc_limited
kscale = 1 / ((Vnom_dc/2) sqrt(3/2) / Vnom_sec)
u_dq = axis_clamp(kscale kdc [vd_conv,vq_conv], -2, 2)
m = hypot(u_d, u_q)
u_abc = inverse_dq(u_dq, theta)
```

The model uses axis saturation, not vector-magnitude saturation.

## PLL and transform boundary

DQsym instantiates the MathWorks SPS `PLL (3ph)` with:

```text
initial [angle_deg, frequency_hz] = [0, Fnom]
PI parameters = [180, 3200, 1]
minimum frequency = 0.9 Fnom
maximum frequency rate = 12
input-filter cutoff = 25 Hz
sample time = 40 us
AGC = on
```

It then uses the SPS `abc to dq0` and `dq0 to abc` blocks with the PLL `wt`
output. The upstream repository does not contain the internal PLL, AGC, or
transform equations. Consequently an exact replacement PLL and an
independent Park-orientation test cannot be implemented from this revision
alone. Existing DPsim `Signal::PLL` is not claimed to be equivalent.

## Implementable Phase 2 scope

Independently verifiable from the pinned source:

- conditional-integration PI;
- DQsym second-order measurement filter;
- active-power, reactive-power and DC-voltage outer-loop equations;
- dq current PI, feedforward and decoupling;
- DC normalization, axis saturation and modulation scaling.

Unresolved and intentionally not implemented:

- PI reset semantics;
- arbitrary tracking/bumpless-transfer logic beyond DQsym state-hold disable;
- exact SPS PLL/AGC equations;
- a replacement Park transform claimed to be bit-equivalent to SPS.
