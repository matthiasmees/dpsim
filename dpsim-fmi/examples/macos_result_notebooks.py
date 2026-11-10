#!/usr/bin/env python3
"""Generate the checked-in macOS comparison notebooks from saved sweep data."""
from pathlib import Path
import nbformat as nbf

HERE = Path(__file__).resolve().parent


def notebook(title: str, intro: str, setup: str, analysis: str, output: Path) -> None:
    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    }
    nb.cells = [
        nbf.v4.new_markdown_cell(intro),
        nbf.v4.new_code_cell(setup),
        nbf.v4.new_markdown_cell("## Error and timing summary"),
        nbf.v4.new_code_cell(analysis),
    ]
    output.write_text(nbf.writes(nb), encoding="utf-8")


def dpsim_notebook(case: str, label: str, bus: str, output: Path) -> None:
    setup = f"""from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

HERE = Path.cwd().resolve()
ROOT = next((p for p in (HERE, *HERE.parents) if (p / "dpsim-fmi").is_dir()), None)
if ROOT is None:
    raise FileNotFoundError("Open this notebook from a DPsim checkout")
BUILD = ROOT / "build-fmi"
if not BUILD.exists():
    BUILD = ROOT / "build-macos-fmi-latest"
RESULTS = BUILD / "{case}-inverter" / "macos-{case}-10s"
metrics = json.loads((RESULTS / "metrics.json").read_text())
steps = [int(x) for x in metrics["configuration"]["steps_us"]]
print("macOS results:", RESULTS)
print(metrics["configuration"])

def saved_path(value):
    path = Path(value)
    if path.exists():
        return path
    marker = "/build-fmi/"
    if marker in str(value):
        return BUILD / str(value).split(marker, 1)[1]
    return path
"""
    analysis = f"""rows = []
for step in steps:
    case = metrics["cases"][str(step)]
    for method in ("me", "cs"):
        run = case["runtimes"][method]
        full = case["comparisons"][method]["full_run"]
        rows.append({{
            "timestep [µs]": step,
            "method": method.upper(),
            "PCC voltage RMS error [%]": full["voltage_rms_error_percent"],
            "inverter current RMS error [%]": full["current_rms_error_percent"],
            "solver time [s]": run["solver_wall_s"],
            "process time [s]": run["process_wall_s"],
            "max Newton iterations": run.get("newton_metrics", {{}}).get("max_iterations"),
            "max scaled residual": run.get("newton_metrics", {{}}).get("max_scaled_residual"),
        }})
summary = pd.DataFrame(rows).sort_values(["timestep [µs]", "method"])
display(summary.style.format({{c: "{{:.6g}}" for c in summary.select_dtypes("number").columns}}))

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
for method, color in (("ME", "tab:blue"), ("CS", "tab:orange")):
    part = summary[summary.method == method]
    axes[0].plot(part["timestep [µs]"], part["PCC voltage RMS error [%]"], "o-", color=color, label=method)
    axes[1].plot(part["timestep [µs]"], part["inverter current RMS error [%]"], "o-", color=color, label=method)
    axes[2].plot(part["timestep [µs]"], part["process time [s]"], "o-", color=color, label=method)
for ax, title, ylabel in zip(axes,
    ("PCC voltage error", "Inverter current error", "Process runtime"),
    ("RMS error [% of native DPsim RMS]", "RMS error [% of native DPsim RMS]", "seconds")):
    ax.set_xscale("log"); ax.set_xlabel("Fixed timestep [µs]"); ax.set_ylabel(ylabel)
    ax.set_title(title); ax.grid(True, which="both", alpha=0.3); ax.legend()
plt.show()

step = 100 if 100 in steps else steps[len(steps)//2]
case = metrics["cases"][str(step)]
event = metrics["configuration"]["load_step_s"]
me = pd.read_csv(saved_path(case["comparisons"]["me"]["load_step_zoom_trace"]))
cs = pd.read_csv(saved_path(case["comparisons"]["cs"]["load_step_zoom_trace"]))
t = (me.time_s - event) * 1000
fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True, constrained_layout=True)
axes[0].plot(t, me.native_v_a_V / 1000, "k--", label="Native DPsim")
axes[0].plot(t, me.me_v_a_V / 1000, color="tab:blue", label="DPsim ME")
axes[0].plot((cs.time_s - event) * 1000, cs.cs_v_a_V / 1000, color="tab:orange", label="DPsim CS")
axes[1].plot(t, me.native_i_a_A, "k--", label="Native DPsim")
axes[1].plot(t, me.me_i_a_A, color="tab:blue", label="DPsim ME")
axes[1].plot((cs.time_s - event) * 1000, cs.cs_i_a_A, color="tab:orange", label="DPsim CS")
axes[0].set_ylabel("PCC phase-a voltage [kV]"); axes[1].set_ylabel("Inverter phase-a current [A]")
axes[1].set_xlabel("Time from load step [ms]")
for ax in axes: ax.axvline(0, color="0.3", ls=":"); ax.grid(True, alpha=0.3); ax.legend()
fig.suptitle(f"{label}: {{step}} µs load-step response")
plt.show()
"""
    intro = f"""# {label} macOS comparison

Fresh Apple Silicon results for the {label} native DPsim, FMI Model Exchange (ME), and FMI Co-Simulation (CS) cases. The network, load event, timesteps, initialization, and metric definitions match the Linux comparison. The saved traces and metrics are under `build-fmi/{case}-inverter/macos-{case}-10s/`."""
    notebook(label, intro, setup, analysis, output)


small_setup = """from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

HERE = Path.cwd().resolve()
ROOT = next((p for p in (HERE, *HERE.parents) if (p / "dpsim-fmi").is_dir()), None)
if ROOT is None:
    raise FileNotFoundError("Open this notebook from a DPsim checkout")
BUILD = ROOT / "build-fmi"
if not BUILD.exists():
    BUILD = ROOT / "build-macos-fmi-latest"
DATA = BUILD / "small-inverter-grid" / "results"
sweep = json.loads((DATA / "metrics_sweep.json").read_text())
steps = [50, 100, 500, 1000]
records = {}
for us in steps:
    key = f"dt_{us}us"
    traces = {
        "DPsim native": pd.read_csv(DATA / "dpsim" / "native" / key / "dpsim_native.csv"),
        "DPsim ME": pd.read_csv(DATA / "dpsim" / "me" / key / "dpsim_me.csv"),
        "DPsim CS": pd.read_csv(DATA / "dpsim" / "cs" / key / "dpsim_cs.csv"),
        "Simulink native": pd.read_csv(DATA / "simulink" / "native" / key / "simulink_native.csv"),
        "Simulink CS": pd.read_csv(DATA / "simulink" / "cs" / key / "simulink_cs.csv"),
    }
    for name, frame in traces.items():
        if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
            raise ValueError(f"non-finite values in {name} at {us} µs")
    records[us] = traces
print("macOS results:", DATA)
"""
small_analysis = """rows = []
for us in steps:
    tr = records[us]
    for method, reference, label in (
        ("DPsim ME", "DPsim native", "DPsim ME vs native"),
        ("DPsim CS", "DPsim native", "DPsim CS vs native"),
        ("Simulink CS", "Simulink native", "Simulink CS vs native"),
    ):
        for prefix, unit, quantity in (("v", "V", "PCC voltage"), ("i", "A", "Inverter current")):
            cols = [f"{prefix}_{p}_{unit}" for p in "abc"]
            a = tr[method][cols].to_numpy(); b = tr[reference][cols].to_numpy()
            rms = float(np.sqrt(np.mean((a - b) ** 2)))
            denom = float(np.sqrt(np.mean(b ** 2)))
            rows.append({"timestep [µs]": us, "comparison": label, "quantity": quantity,
                         "RMS error [unit]": rms, "RMS error [% reference RMS]": 100 * rms / denom})
errors = pd.DataFrame(rows)
display(errors.pivot_table(index=["timestep [µs]", "comparison"], columns="quantity", values="RMS error [% reference RMS]").style.format("{:.6g}"))
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
for comparison, group in errors.groupby("comparison"):
    for ax, quantity in zip(axes, ("PCC voltage", "Inverter current")):
        part = group[group.quantity == quantity]
        ax.plot(part["timestep [µs]"], part["RMS error [% reference RMS]"], "o-", label=comparison)
for ax, title in zip(axes, ("PCC voltage error", "Inverter current error")):
    ax.set_xscale("log"); ax.set_xlabel("Fixed timestep [µs]"); ax.set_ylabel("RMS error [%]")
    ax.set_title(title); ax.grid(True, which="both", alpha=0.3); ax.legend()
plt.show()

us = 100
tr = records[us]
fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True, constrained_layout=True)
for name, style in (("DPsim native", "k--"), ("DPsim ME", "b-"), ("DPsim CS", "orange")):
    axes[0].plot(tr[name].time_s, tr[name].v_a_V / 1000, style, label=name)
    axes[1].plot(tr[name].time_s, tr[name].i_a_A, style, label=name)
axes[0].set_ylabel("PCC phase-a voltage [kV]"); axes[1].set_ylabel("Inverter phase-a current [A]")
axes[1].set_xlabel("Time [s]")
for ax in axes: ax.axvline(2.5, color="0.3", ls=":"); ax.grid(True, alpha=0.3); ax.legend()
fig.suptitle("Small grid macOS response at 100 µs")
plt.show()
"""
small_intro = """# Small-grid macOS comparison

Fresh Apple Silicon results for the small inverter grid. The notebook compares DPsim native, FMI ME, and FMI CS against the native DPsim baseline, and also compares the Simulink native and imported-FMU CS replay at the same timesteps. The saved traces are under `build-fmi/small-inverter-grid/results/`."""
notebook(
    "Small-grid macOS comparison",
    small_intro,
    small_setup,
    small_analysis,
    HERE / "small_inverter_grid" / "SmallGrid_macOS_Comparison.ipynb",
)

dpsim_notebook(
    "ieee39",
    "IEEE39 inverter macOS comparison",
    "bus 30",
    HERE / "ieee39_inverter" / "IEEE39_macOS_Comparison.ipynb",
)
dpsim_notebook(
    "ieee9",
    "IEEE9 inverter macOS comparison",
    "bus 3",
    HERE / "ieee9_inverter" / "IEEE9_macOS_Comparison.ipynb",
)
