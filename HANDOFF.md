# Scalar-DC SSN MMC Development Handoff

## Phase 2 controller implementation

Authoritative source:

```text
https://github.com/control-protection-grids-tudelft/DP
commit 36f47adb7c10c8b109e25a2180f20e08c645cfc4
```

Detailed source locations and derivations are in
`docs/development/SSN_MMC_PHASE2_CONTROL_DESIGN.md`.

### Classes and files

- `Signal::DQSymPIController`
- `Signal::DQSymSecondOrderFilter`
- `Signal::DQSymOuterController`
- `Signal::DQSymCurrentController`
- `Signal::DQSymModulation`

Public declarations:
`dpsim-models/include/dpsim-models/Signal/DQSymControllerBlocks.h`.

Implementation:
`dpsim-models/src/Signal/DQSymControllerBlocks.cpp`.

Focused validation:
`dpsim/examples/cxx/Circuits/DQSym_Controller_Validation.cpp`.

The public header is included by `dpsim-models/Components.h`; the source and
validation target are registered in the respective CMake source lists. These
classes are signal components, own no electrical nodes and provide no MNA
stamps.

### Equations and hierarchy

DQsym uses the cascade:

```text
eP = Pmeas-Pref
eP -> active-power PI -> Vdc_ref

eVdc = (Vdc_meas-Vdc_ref)/Vnom_dc
eVdc -> DC-voltage PI -> id_ref

eQ = Qref-Qmeas
eQ -> reactive-power PI -> iq_ref
```

PI and conditional-integration anti-windup:

```text
u_raw = Kp e + x + feedforward
u = clamp(u_raw, lower, upper)

hold = (u_raw >= upper and e > 0)
    or (u_raw <= lower and e < 0)

g = hold ? 0 : Ki e
x[k] = clamp(x[k-1] + Ts/2 (g[k-1]+g[k]), lower, upper)
```

The same internal PI update routine is used by the reusable scalar PI and
both current-controller axes.

Current controller:

```text
ed = id_ref-id
eq = iq_ref-iq

vd_ref = vd + Rff id - (f/Fnom) Lff iq + PI_d(ed)
vq_ref = vq + Rff iq + (f/Fnom) Lff id + PI_q(eq)
```

DQsym parameters are `Kp_I=0.6`, `Ki_I=6`, `Rff=Rarm_pu/2`,
`Lff=Larm_pu/2`, and per-axis limits `[-2,2] pu`.

The measurement block is the DQsym second-order low-pass filter:

```text
omega0 = 2 pi f0
xdot1 = x2
xdot2 = omega0^2 (u-x1) - 2 zeta omega0 x2
y = x1
```

DQsym specifies `f0=1000 Hz`, `zeta=1`, and controller period `40 us`.
DPsim discretizes this state model with its trapezoidal state-space routine.

Modulation/reference normalization:

```text
Vdc_limited = clamp(Vdc, 0.75 Vnom_dc, 1.25 Vnom_dc)
kdc = Vnom_dc/Vdc_limited
kscale = 1/((Vnom_dc/2) sqrt(3/2)/Vnom_sec)
udq = axis_clamp(kscale kdc vdq_ref, -2, 2)
m = hypot(ud,uq)
```

### Units and signs

- `Vdc`: scalar pole-to-pole volts.
- AC voltage base: line-line RMS volts; controller dq voltages are
  amplitude-based per unit.
- Currents and current references: amplitude-based per unit.
- Frequency input: hertz; the decoupling multiplier is `f/Fnom`.
- PI limits and states use the unit of the corresponding controller output.
- Modulation dq and abc outputs are dimensionless reference quantities.
- Current is positive out of the converter.
- `P=vd*id+vq*iq`; positive `id` means positive generated active power.
- `Q=vq*id-vd*iq`; negative `iq` means positive generated reactive power.
- The implemented inverse-dq convention is
  `va=d cos(theta)-q sin(theta)` with phases shifted by `-2pi/3`.
  It matches the traced zero-angle DQsym connection but is not claimed to be
  an exact replacement for the internal SPS transform.

### Enable and initialization

DQsym outer-loop disable substitutes the reference for the measurement,
making the error zero. The implementation holds the integrator state while
disabled. Re-enable is bumpless when error is zero and the retained state is
already the required output.

PI initialization explicitly supplies the integrator state. Filter
steady-state initialization uses `x1=input_initial`, `x2=0`, producing no
first-step transient. Current-controller d/q integrator states are set
explicitly.

DQsym provides no reset or tracking-reset input. No reset equation or general
tracking-based bumpless transfer was invented.

### Validation

Commands:

```text
git diff --check
cmake --build build --target DQSym_Controller_Validation --parallel 2
./build/dpsim/examples/cxx/DQSym_Controller_Validation
cmake --build build --target dpsim-models --parallel 2
cmake --build build --target EMT_SSN_MMC_PiLine --parallel 2
```

Validated tolerances/results:

- PI proportional/trapezoidal response: `1e-12`, observed error `0`.
- Saturation and anti-windup state hold: exact.
- Disabled-state hold: exact.
- Zero-error state-hold re-enable: `1e-12`, observed error `0`.
- Critically damped filter analytical comparison: tolerance `3e-4`,
  observed error `4.26545e-5`.
- Filter steady-state initialization: `1e-12`, observed error `0`.
- dq feedforward and decoupling: `1e-12`, observed error `0`.
- Modulation scaling: `1e-12`, observed error `0`.
- Zero-angle inverse-dq phases: `1e-12`, maximum observed error
  `4.44e-16`.
- Axis saturation: exact.
- All checked states and outputs remain finite.

### Unresolved authoritative boundaries

DQsym references MathWorks Specialized Power Systems `PLL (3ph)`,
`abc to dq0`, and `dq0 to abc` blocks. The pinned upstream repository
contains mask parameters and connections but not their internal equations.
No exact DQsym-equivalent PLL or transform is claimed.

DQsym also supplies no PI integrator-reset, tracking-reset, or general
bumpless-transfer equation. Only the verified state-hold disable behavior is
implemented.

Do not start station composition by silently substituting the existing DPsim
PLL or by claiming exact SPS transform equivalence.
