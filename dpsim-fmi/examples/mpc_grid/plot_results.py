#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""Summarize and plot the short MPC-FMU/DPsim experiment."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "build-fmi/mpc-grid-case"
LOG = RESULTS / "logs/fmi_mixed.csv"
if not LOG.is_file():
    raise FileNotFoundError(f"Run the DPsim case first; missing {LOG}")
d = pd.read_csv(LOG, skipinitialspace=True)
d.columns = d.columns.str.strip()


def phases(prefix):
    return d[[f"{prefix}_{p}" for p in range(3)]].to_numpy(dtype=float)


t = d.time.to_numpy()
mask = t >= 1e-3
current = np.column_stack([d[f"i_{p}"] for p in "abc"])
reference = np.column_stack([d[f"i_ref_{p}"] for p in "abc"])
error = current - reference
kcl = phases("line_current") + phases("load_current") + phases("converter_current")
metrics = {
    "samples": int(len(d)),
    "duration_s": float(t[-1] - t[0]),
    "controller_ticks_final": int(d.controller_ticks.iloc[-1]),
    "current_tracking_rms_error_after_1ms_A": float(np.sqrt(np.mean(error[mask] ** 2))),
    "current_tracking_max_error_after_1ms_A": float(np.max(np.abs(error[mask]))),
    "max_pcc_kcl_residual_A": float(np.max(np.abs(kcl))),
    "max_newton_iterations": int(d.newton_iterations.max()),
    "max_reported_newton_residual": float(d.newton_residual.max()),
    "max_fmi_coupling_error_A": float(d.fmi_coupling_error_A.max()),
    "mean_converter_active_power_after_1ms_W": float(
        np.mean(np.sum(phases("pcc_voltage")[mask] * current[mask], axis=1))
    ),
}
print(json.dumps(metrics, indent=2))
RESULTS.mkdir(parents=True, exist_ok=True)
(RESULTS / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
axes[0].plot(t * 1e3, d.grid_voltage_0, label="grid phase a")
axes[0].plot(t * 1e3, d.pcc_voltage_0, label="PCC phase a")
axes[0].set_ylabel("Voltage [V]")
axes[0].legend()
for phase in "abc":
    axes[1].plot(t * 1e3, d[f"i_{phase}"], label=f"current {phase}")
    axes[1].plot(
        t * 1e3, d[f"i_ref_{phase}"], ls="--", alpha=0.7, label=f"reference {phase}"
    )
axes[1].set_ylabel("Current [A]")
axes[1].legend(ncol=3, fontsize=8)
axes[2].step(t * 1e3, d.controller_ticks, where="post", label="MPC ticks")
axes[2].set_ylabel("Tick count")
axes[2].set_xlabel("Time [ms]")
axes[2].legend()
for ax in axes:
    ax.grid(alpha=0.25)
fig.suptitle("FCS-MPC-controlled three-phase component imported into DPsim")
fig.tight_layout()
fig.savefig(RESULTS / "mpc_grid_results.png", dpi=150)
