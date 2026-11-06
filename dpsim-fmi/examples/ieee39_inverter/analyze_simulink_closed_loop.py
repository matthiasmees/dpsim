#!/usr/bin/env python3
"""Compare closed-loop Simulink runs against matched DPsim IEEE39 CS runs."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BUILD = next(
    (
        p
        for p in (ROOT / "build-fmi-cosim", ROOT / "build-fmi")
        if (p / "CMakeCache.txt").is_file()
    ),
    ROOT / "build-fmi-cosim",
)
DATA = BUILD / "ieee39-inverter"
SIMULINK = DATA / "simulink-closed-loop-sweep-exact-order"
VOLTAGE_BASE_V = 345e3
CURRENT_BASE_A = math.sqrt(2) * math.hypot(125e6, 80.881e6) / (math.sqrt(3) * 345e3)


def metric(values: np.ndarray, base: float) -> dict[str, float]:
    return {
        "max_abs": float(np.max(np.abs(values))),
        "rms_all_phases": float(np.sqrt(np.mean(values * values))),
        "max_abs_percent_base": float(100 * np.max(np.abs(values)) / base),
    }


def runner_timing(path: Path) -> dict[str, float]:
    text = path.read_text(errors="replace")
    match = re.search(
        r"runtime_s=([0-9.eE+-]+)\s+initialization_s=([0-9.eE+-]+)\s+simulation_s=([0-9.eE+-]+)",
        text,
    )
    if not match:
        raise ValueError(f"Cannot parse DPsim timing from {path}")
    return {
        "total_runtime_s": float(match.group(1)),
        "initialization_s": float(match.group(2)),
        "simulation_s": float(match.group(3)),
    }


def main() -> None:
    manifest = json.loads((SIMULINK / "manifest.json").read_text())
    summaries = []
    traces: dict[tuple[float, int], pd.DataFrame] = {}
    for entry in manifest:
        folder = SIMULINK / entry["sweep_dir"] / entry["case_dir"]
        sim = pd.read_csv(folder / "simulink_closed_loop_pcc.csv")
        dpsim_case = (
            DATA
            / "dpsim-cs-matched-rerun"
            / entry["sweep_dir"]
            / entry["case_dir"]
            / "cs"
        )
        dpsim_file = (
            dpsim_case / "logs" / "IEEE39_Inverter_cs" / "IEEE39_Inverter_cs.csv"
        )
        dpsim = pd.read_csv(dpsim_file, skipinitialspace=True)
        dpsim.columns = dpsim.columns.str.strip()
        n = min(len(sim), len(dpsim))
        if abs(len(sim) - len(dpsim)) > 1:
            raise ValueError(
                f"Sample count mismatch for {entry['label']}: Simulink={len(sim)}, DPsim={len(dpsim)}"
            )
        dt = float(entry["step_s"])
        if abs(sim.time.iloc[-1] - entry["stop_time_s"]) > dt * 1e-6:
            raise ValueError(f"Simulink endpoint mismatch for {entry['label']}")
        if not np.allclose(
            sim.time.iloc[:n], dpsim.time.iloc[:n], rtol=0, atol=max(5.1e-6, dt * 1e-6)
        ):
            raise ValueError(f"Time grids differ for {entry['label']}")
        t = sim.time.to_numpy()[:n]
        currents_s = sim[["i_a", "i_b", "i_c"]].to_numpy()[:n]
        currents_d = dpsim[[f"inverter.i_{p}" for p in range(3)]].to_numpy()[:n]
        v30_s = sim[["bus30_v_a", "bus30_v_b", "bus30_v_c"]].to_numpy()[:n]
        v30_d = dpsim[[f"N30.V_{p}" for p in range(3)]].to_numpy()[:n]
        v39_s = sim[["bus39_v_a", "bus39_v_b", "bus39_v_c"]].to_numpy()[:n]
        v39_d = dpsim[[f"N39.V_{p}" for p in range(3)]].to_numpy()[:n]
        errors = {
            "inverter_current_A": currents_s - currents_d,
            "bus30_voltage_V": v30_s - v30_d,
            "bus39_voltage_V": v39_s - v39_d,
        }
        event = float(entry["load_step_s"])
        windows = {
            "full_run": np.ones(n, dtype=bool),
            "pre_event": t < event,
            "first_20ms_after_event": (t >= event) & (t <= event + 0.02),
            "post_event": t >= event,
        }
        base = {
            "inverter_current_A": CURRENT_BASE_A,
            "bus30_voltage_V": VOLTAGE_BASE_V,
            "bus39_voltage_V": VOLTAGE_BASE_V,
        }
        case_metrics = json.loads(
            (folder / "simulink_closed_loop_timing.json").read_text()
        )
        dpsim_time = runner_timing(dpsim_case / "runner.log")
        row: dict = {
            "sweep_dir": entry["sweep_dir"],
            "case_dir": entry["case_dir"],
            "stop_time_s": float(entry["stop_time_s"]),
            "load_step_s": event,
            "step_us": dt * 1e6,
            "samples": n,
            "simulink_samples": len(sim),
            "dpsim_samples": len(dpsim),
            "simulink_only_terminal_samples": max(0, len(sim) - len(dpsim)),
            "simulink_wall_s": case_metrics["simulink_wall_time_s"],
            "simulink_network_init_s": case_metrics["network_initialization_s"],
            "simulink_network_step_mean_s": case_metrics["network_step_mean_s"],
            "simulink_network_step_max_s": case_metrics["network_step_max_s"],
            "simulink_network_steps_total_s": case_metrics["network_step_total_s"],
            "dpsim_total_wall_s": dpsim_time["total_runtime_s"],
            "dpsim_init_s": dpsim_time["initialization_s"],
            "dpsim_simulation_s": dpsim_time["simulation_s"],
        }
        for signal, values in errors.items():
            for window, mask in windows.items():
                row[f"{signal}.{window}"] = metric(values[mask], base[signal])
        summaries.append(row)
        trace = pd.DataFrame(
            {
                "time_s": t,
                "step_us": dt * 1e6,
                "current_phase_a_error_percent": 100
                * errors["inverter_current_A"][:, 0]
                / CURRENT_BASE_A,
                "bus30_phase_a_voltage_error_percent": 100
                * errors["bus30_voltage_V"][:, 0]
                / VOLTAGE_BASE_V,
                "bus39_phase_a_voltage_error_percent": 100
                * errors["bus39_voltage_V"][:, 0]
                / VOLTAGE_BASE_V,
            }
        )
        window = max(1, round(0.5 / dt))
        for key in (
            "current_phase_a_error_percent",
            "bus30_phase_a_voltage_error_percent",
            "bus39_phase_a_voltage_error_percent",
        ):
            trace[key.replace("error_percent", "0p5s_rms_error_percent")] = np.sqrt(
                trace[key].pow(2).rolling(window, center=True, min_periods=1).mean()
            )
        traces[(float(entry["stop_time_s"]), round(dt * 1e6))] = trace
        trace.to_csv(folder / "closed_loop_error_trace.csv", index=False)

    (SIMULINK / "closed_loop_comparison_metrics.json").write_text(
        json.dumps(
            {
                "reference": "DPsim FMI Co-Simulation on the same IEEE39 MNA network, same FMU and parameters",
                "voltage_percent_base": VOLTAGE_BASE_V,
                "current_percent_base_A": CURRENT_BASE_A,
                "cases": summaries,
            },
            indent=2,
        )
        + "\n"
    )
    flat = []
    for case in summaries:
        for name, value in case.items():
            if isinstance(value, dict):
                flat.append(
                    {
                        "sweep_dir": case["sweep_dir"],
                        "step_us": case["step_us"],
                        "stop_time_s": case["stop_time_s"],
                        "signal_window": name,
                        **value,
                    }
                )
    pd.DataFrame(flat).to_csv(
        SIMULINK / "closed_loop_comparison_metrics.csv", index=False
    )

    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(summaries)))
    for stop in sorted({case["stop_time_s"] for case in summaries}):
        selected = [case for case in summaries if case["stop_time_s"] == stop]
        fig, axes = plt.subplots(
            3, 1, figsize=(12, 9), sharex=True, constrained_layout=True
        )
        keys = [
            "current_phase_a_0p5s_rms_error_percent",
            "bus30_phase_a_voltage_0p5s_rms_error_percent",
            "bus39_phase_a_voltage_0p5s_rms_error_percent",
        ]
        titles = [
            "Inverter phase-a current error",
            "Bus 30 phase-a voltage error",
            "Bus 39 phase-a voltage error",
        ]
        for case, color in zip(selected, colors):
            trace = traces[(case["stop_time_s"], round(case["step_us"]))]
            label = f"{case['step_us']:g} µs"
            for axis, key in zip(axes, keys):
                axis.plot(
                    trace.time_s, trace[key], color=color, linewidth=1.1, label=label
                )
        for axis, title in zip(axes, titles):
            axis.set_title(title)
            axis.set_ylabel("0.5 s RMS error [% of base]")
            axis.grid(True, alpha=0.35)
            axis.axvline(
                stop / 2,
                color="0.35",
                linestyle=":",
                linewidth=1.1,
                label="100 MW load step",
            )
        axes[0].legend(ncol=4, fontsize=9)
        axes[-1].set_xlabel("Time [s]")
        fig.suptitle(f"Closed-loop Simulink IEEE39 vs DPsim CS — {stop:g} s run")
        fig.savefig(SIMULINK / f"closed_loop_error_{stop:g}s.png", dpi=160)
        plt.close(fig)
    print(f"Analyzed {len(summaries)} matched closed-loop cases: {SIMULINK}")


if __name__ == "__main__":
    main()
