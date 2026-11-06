#!/usr/bin/env python3
"""Analyze matched 20 s IEEE39 closed-loop sweeps and ME refinement error.

Native Simulink is retained as a behavior comparison. The 25 us DPsim ME run
is separately used as a fine-step DPsim convergence reference; it is not an
independent physical reference.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STEPS_US = (25, 50, 100, 200, 500, 625, 1000)
METHODS = ("Simulink FMU CS", "DPsim ME", "DPsim CS")
DP_SWEEP = "timestep-sweep-20s-live-init"
NATIVE_SWEEP = "native-closed-loop-20s-live-init"
SIMULINK_CS_SWEEP = "simulink-closed-loop-20s-live-init"
COLORS = plt.cm.viridis(np.linspace(0.05, 0.95, len(STEPS_US)))
STYLES = {"Simulink FMU CS": "-", "DPsim ME": "--", "DPsim CS": ":"}


def find_root(start: Path) -> Path:
    for path in (start, *start.parents):
        if (path / "build-fmi" / "ieee39-inverter").is_dir():
            return path
    raise FileNotFoundError(
        "Could not locate build-fmi/ieee39-inverter from the current directory"
    )


def read_trace(path: Path, method: str, step_us: int) -> dict:
    frame = pd.read_csv(path, skipinitialspace=True)
    frame.columns = frame.columns.str.strip()
    if method.startswith("DPsim"):
        phases = ("a", "b", "c")
        trace = {
            "time": frame["time"].to_numpy(dtype=float),
            "current": np.column_stack([frame[f"inverter.i_{i}"] for i in range(3)]),
            "bus30": np.column_stack([frame[f"N30.V_{i}"] for i in range(3)]),
            "bus39": np.column_stack([frame[f"N39.V_{i}"] for i in range(3)]),
        }
    else:
        trace = {
            "time": frame["time"].to_numpy(dtype=float),
            "current": frame[[f"i_{p}" for p in "abc"]].to_numpy(dtype=float),
            "bus30": frame[[f"bus30_v_{p}" for p in "abc"]].to_numpy(dtype=float),
            "bus39": frame[[f"bus39_v_{p}" for p in "abc"]].to_numpy(dtype=float),
        }
    for key in ("current", "bus30", "bus39"):
        if not np.isfinite(trace[key]).all():
            raise ValueError(f"Non-finite values in {path} ({key})")
    if len(trace["time"]) < 2 or np.any(np.diff(trace["time"]) <= 0):
        raise ValueError(f"Invalid time samples in {path}")
    trace["source"] = str(path)
    trace["step_us"] = step_us
    return trace


def load_case(data: Path, step_us: int) -> dict:
    key = f"dt-{step_us:03d}us"
    native = read_trace(
        data / NATIVE_SWEEP / "timestep-sweep-20s" / key / "native_closed_loop_pcc.csv",
        "Native Simulink",
        step_us,
    )
    simcs = read_trace(
        data / SIMULINK_CS_SWEEP / DP_SWEEP / key / "simulink_closed_loop_pcc.csv",
        "Simulink FMU CS",
        step_us,
    )
    root = data / DP_SWEEP / key
    me = read_trace(
        root / "me" / "logs" / "IEEE39_Inverter_me" / "IEEE39_Inverter_me.csv",
        "DPsim ME",
        step_us,
    )
    cs = read_trace(
        root / "cs" / "logs" / "IEEE39_Inverter_cs" / "IEEE39_Inverter_cs.csv",
        "DPsim CS",
        step_us,
    )
    traces = {
        "Native Simulink": native,
        "Simulink FMU CS": simcs,
        "DPsim ME": me,
        "DPsim CS": cs,
    }
    dt = step_us * 1e-6
    for name, trace in traces.items():
        if not np.allclose(
            np.diff(trace["time"]), dt, atol=max(5.1e-6, dt * 1e-6), rtol=0
        ):
            raise ValueError(f"{name} time increments do not match {step_us} us")
    # DPsim loggers can omit the final endpoint. Compare identical sample indices,
    # with one common time vector, instead of extrapolating a trace.
    count = min(len(t["time"]) for t in traces.values())
    common_time = native["time"][:count]
    for name, trace in traces.items():
        if not np.allclose(
            trace["time"][:count], common_time, atol=max(5.1e-6, dt * 1e-6), rtol=0
        ):
            raise ValueError(
                f"{name} sample grid differs from native Simulink at {step_us} us"
            )
        for field in ("current", "bus30", "bus39"):
            trace[field] = trace[field][:count]
        trace["time"] = common_time
    return {"traces": traces, "time": common_time, "step_us": step_us, "dt": dt}


def run_timing(data: Path, step_us: int, method: str) -> dict:
    key = f"dt-{step_us:03d}us"
    if method == "Native Simulink":
        path = (
            data
            / NATIVE_SWEEP
            / "timestep-sweep-20s"
            / key
            / "native_closed_loop_timing.json"
        )
        raw = json.loads(path.read_text())
        return {
            "method": method,
            "step_us": step_us,
            "wall_s": raw["simulink_wall_time_s"],
            "initialization_s": raw["network_initialization_s"],
            "network_step_mean_us": raw["network_step_mean_s"] * 1e6,
            "network_step_max_us": raw["network_step_max_s"] * 1e6,
            "network_step_total_s": raw["network_step_total_s"],
            "steps": raw["samples"] - 1,
        }
    if method == "Simulink FMU CS":
        path = (
            data
            / SIMULINK_CS_SWEEP
            / DP_SWEEP
            / key
            / "simulink_closed_loop_timing.json"
        )
        raw = json.loads(path.read_text())
        return {
            "method": method,
            "step_us": step_us,
            "wall_s": raw["simulink_wall_time_s"],
            "initialization_s": raw["network_initialization_s"],
            "network_step_mean_us": raw["network_step_mean_s"] * 1e6,
            "network_step_max_us": raw["network_step_max_s"] * 1e6,
            "network_step_total_s": raw["network_step_total_s"],
            "steps": raw["samples"] - 1,
        }
    mode = "me" if method == "DPsim ME" else "cs"
    folder = data / DP_SWEEP / key / mode
    match = re.search(
        rf"IEEE39_Inverter_{mode} runtime_s=([0-9.eE+-]+) initialization_s=([0-9.eE+-]+) simulation_s=([0-9.eE+-]+) steps=(\d+)",
        (folder / "runner.log").read_text(),
    )
    if not match:
        raise ValueError(f"Could not parse DPsim timing from {folder/'runner.log'}")
    runtime, init, sim, steps = match.groups()
    return {
        "method": method,
        "step_us": step_us,
        "wall_s": float(runtime),
        "initialization_s": float(init),
        "network_step_mean_us": np.nan,
        "network_step_max_us": np.nan,
        "network_step_total_s": float(sim),
        "steps": int(steps),
    }


def binned_rms_percent(
    time: np.ndarray, error: np.ndarray, reference: np.ndarray, width_s=0.5
) -> tuple[np.ndarray, np.ndarray]:
    bins = np.floor((time - time[0]) / width_s).astype(int)
    out_t, out_e = [], []
    for bucket in np.unique(bins):
        mask = bins == bucket
        ref_rms = np.sqrt(np.mean(reference[mask] ** 2))
        err_rms = np.sqrt(np.mean(error[mask] ** 2))
        out_t.append(float(np.mean(time[mask])))
        out_e.append(100 * err_rms / max(ref_rms, 1e-12))
    return np.asarray(out_t), np.asarray(out_e)


def build(start: Path | None = None) -> dict:
    root = find_root((start or Path.cwd()).resolve())
    data = root / "build-fmi" / "ieee39-inverter"
    output = data / "native-simulink-20s-allcases-analysis"
    output.mkdir(parents=True, exist_ok=True)
    cases = {step: load_case(data, step) for step in STEPS_US}
    metric_rows, time_rows, direct_rows, fine_rows, initial_rows = [], [], [], [], []
    signal_fields = {
        "Inverter current": "current",
        "Bus 30 PCC voltage": "bus30",
        "Bus 39 load-bus voltage": "bus39",
    }
    fine_case = cases[25]
    fine_trace = fine_case["traces"]["DPsim ME"]
    for step, case in cases.items():
        ref = case["traces"]["Native Simulink"]
        time = case["time"]
        for method, trace in case["traces"].items():
            for signal, field in (
                ("Inverter current", "current"),
                ("Bus 30 PCC voltage", "bus30"),
            ):
                difference = trace[field][0] - ref[field][0]
                initial_rows.append(
                    {
                        "step_us": step,
                        "method": method,
                        "signal": signal,
                        "initial_time_s": float(trace["time"][0]),
                        "native_initial_rms": float(
                            np.sqrt(np.mean(ref[field][0] ** 2))
                        ),
                        "candidate_initial_rms": float(
                            np.sqrt(np.mean(trace[field][0] ** 2))
                        ),
                        "max_abs_phase_difference": float(np.max(np.abs(difference))),
                        "max_abs_phase_difference_percent_of_native_rms": 100
                        * float(np.max(np.abs(difference)))
                        / max(float(np.sqrt(np.mean(ref[field][0] ** 2))), 1e-12),
                    }
                )
        for method in METHODS:
            candidate = case["traces"][method]
            for signal, field in signal_fields.items():
                err = candidate[field] - ref[field]
                ref_sig = ref[field]
                keep = time >= 0.02
                rms_ref = float(np.sqrt(np.mean(ref_sig[keep] ** 2)))
                rms_err = float(np.sqrt(np.mean(err[keep] ** 2)))
                peak_percent = float(
                    100 * np.max(np.abs(err[keep])) / max(rms_ref, 1e-12)
                )
                max_t = int(np.argmax(np.max(np.abs(err[keep]), axis=1)))
                metric_rows.append(
                    {
                        "method": method,
                        "step_us": step,
                        "signal": signal,
                        "compared_samples": int(keep.sum()),
                        "comparison_stop_s": float(time[-1]),
                        "rms_error_percent_of_native_rms": 100
                        * rms_err
                        / max(rms_ref, 1e-12),
                        "max_abs_phase_percent_of_native_rms": peak_percent,
                        "max_phase_error_abs": float(np.max(np.abs(err[keep]))),
                        "native_rms": rms_ref,
                        "max_error_time_s": float(time[keep][max_t]),
                    }
                )
                # Phase-a time-local RMS percentage, smoothed only for display.
                bt, bp = binned_rms_percent(
                    time[keep], err[keep, 0], ref_sig[keep, 0], width_s=0.5
                )
                time_rows.extend(
                    {
                        "method": method,
                        "step_us": step,
                        "signal": signal,
                        "time_s": float(t),
                        "phase_a_rms_error_percent": float(p),
                    }
                    for t, p in zip(bt, bp)
                )
        for signal, field in signal_fields.items():
            err = (
                case["traces"]["Simulink FMU CS"][field]
                - case["traces"]["DPsim CS"][field]
            )
            ref_sig = ref[field]
            keep = time >= 0.02
            ref_rms = float(np.sqrt(np.mean(ref_sig[keep] ** 2)))
            err_rms = float(np.sqrt(np.mean(err[keep] ** 2)))
            direct_rows.append(
                {
                    "step_us": step,
                    "signal": signal,
                    "rms_error_percent_of_native_rms": 100
                    * err_rms
                    / max(ref_rms, 1e-12),
                    "max_abs_phase_percent_of_native_rms": 100
                    * float(np.max(np.abs(err[keep])))
                    / max(ref_rms, 1e-12),
                    "max_phase_error_abs": float(np.max(np.abs(err[keep]))),
                }
            )
        # Compare every independently initialized closed-loop result against
        # the finest available DPsim ME trajectory on exactly matching times.
        stride = step // 25
        if step % 25:
            raise ValueError(
                f"{step} us is not an integer multiple of the 25 us reference grid"
            )
        fine_indices = np.arange(0, len(fine_trace["time"]), stride)
        for method, candidate in case["traces"].items():
            count = min(len(candidate["time"]), len(fine_indices))
            candidate_time = candidate["time"][:count]
            reference_time = fine_trace["time"][fine_indices[:count]]
            if not np.allclose(
                candidate_time, reference_time, atol=max(5.1e-6, step * 1e-12), rtol=0
            ):
                raise ValueError(
                    f"{method} time grid cannot be aligned to 25 us DPsim ME at {step} us"
                )
            keep = candidate_time >= 0.02
            for signal, field in signal_fields.items():
                ref_sig = fine_trace[field][fine_indices[:count]][keep]
                err = candidate[field][:count][keep] - ref_sig
                ref_rms = float(np.sqrt(np.mean(ref_sig**2)))
                err_rms = float(np.sqrt(np.mean(err**2)))
                fine_rows.append(
                    {
                        "method": method,
                        "step_us": step,
                        "signal": signal,
                        "reference": "DPsim ME 25 us fine-step trajectory",
                        "reference_is_independent_physical_model": False,
                        "compared_samples": int(keep.sum()),
                        "rms_error_percent_of_fine_me_rms": 100
                        * err_rms
                        / max(ref_rms, 1e-12),
                        "max_abs_phase_percent_of_fine_me_rms": 100
                        * float(np.max(np.abs(err)))
                        / max(ref_rms, 1e-12),
                        "fine_me_rms": ref_rms,
                    }
                )
    metrics = pd.DataFrame(metric_rows)
    time_errors = pd.DataFrame(time_rows)
    direct = pd.DataFrame(direct_rows)
    fine_metrics = pd.DataFrame(fine_rows)
    initial_metrics = pd.DataFrame(initial_rows)
    timings = pd.DataFrame(
        [
            run_timing(data, step, method)
            for step in STEPS_US
            for method in ("Native Simulink", *METHODS)
        ]
    )
    metrics.to_csv(output / "native_reference_error_metrics.csv", index=False)
    time_errors.to_csv(
        output / "native_reference_error_time_series_0p5s.csv", index=False
    )
    direct.to_csv(output / "simulink_cs_vs_dpsim_cs_metrics.csv", index=False)
    fine_metrics.to_csv(output / "fine_me25_reference_metrics.csv", index=False)
    initial_metrics.to_csv(output / "initialization_consistency.csv", index=False)
    timings.to_csv(output / "simulation_timing_summary.csv", index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axs = plt.subplots(
        2, 1, figsize=(15, 10), sharex=True, constrained_layout=True
    )
    plotted = [
        ("Inverter current", 0, "Phase-a current RMS error [% of native phase-a RMS]"),
        (
            "Bus 30 PCC voltage",
            1,
            "Bus 30 phase-a voltage RMS error [% of native phase-a RMS]",
        ),
    ]
    for title, axis_index, ylabel in plotted:
        ax = axs[axis_index]
        for method in METHODS:
            for color, step in zip(COLORS, STEPS_US):
                group = time_errors[
                    (time_errors.method == method)
                    & (time_errors.step_us == step)
                    & (time_errors.signal == title)
                ]
                ax.plot(
                    group.time_s,
                    group.phase_a_rms_error_percent,
                    color=color,
                    linestyle=STYLES[method],
                    linewidth=1.5,
                    label=f"{method}, {step} µs",
                )
        ax.axhline(0, color="black", linewidth=0.7)
        ax.axvline(
            10,
            color="0.25",
            linestyle="--",
            linewidth=1,
            label="100 MW load step" if axis_index == 0 else None,
        )
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 20)
        ax.legend(ncol=4, fontsize=7, loc="upper left")
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle(
        "20 s closed-loop IEEE39: every timestep compared with native Simulink"
    )
    fig.savefig(output / "native_reference_error_vs_time.png", dpi=170)
    plt.close(fig)

    # The imported FMI block is much closer to the native Simulink build than
    # the DPsim-hosted runs, so show its seven traces on a dedicated scale too.
    fig, axs = plt.subplots(2, 1, figsize=(14, 7), sharex=True, constrained_layout=True)
    for ax, signal, ylabel in (
        (
            axs[0],
            "Inverter current",
            "Phase-a current RMS error [% of native phase-a RMS]",
        ),
        (
            axs[1],
            "Bus 30 PCC voltage",
            "Bus 30 phase-a voltage RMS error [% of native phase-a RMS]",
        ),
    ):
        for color, step in zip(COLORS, STEPS_US):
            group = time_errors[
                (time_errors.method == "Simulink FMU CS")
                & (time_errors.step_us == step)
                & (time_errors.signal == signal)
            ]
            ax.plot(
                group.time_s,
                group.phase_a_rms_error_percent,
                color=color,
                linewidth=1.6,
                label=f"{step} µs",
            )
        ax.axhline(0, color="black", linewidth=0.7)
        ax.axvline(10, color="0.25", linestyle="--", linewidth=1)
        ax.set_ylabel(ylabel)
        ax.legend(ncol=4, fontsize=8)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle(
        "Simulink FMU CS error relative to native Simulink (zoomed, all timesteps)"
    )
    fig.savefig(output / "native_reference_error_zoom_simulink_cs.png", dpi=170)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for ax, (signal, _) in zip(axs, signal_fields.items()):
        part = metrics[metrics.signal == signal]
        for method in METHODS:
            m = part[part.method == method]
            ax.plot(
                m.step_us,
                m.max_abs_phase_percent_of_native_rms,
                marker="o",
                label=method,
            )
        ax.set_xscale("log", base=2)
        ax.set_xticks(STEPS_US, labels=[str(x) for x in STEPS_US])
        ax.set_xlabel("Timestep [µs]")
        ax.set_ylabel("Maximum absolute phase error [% of native RMS]")
        ax.set_title(signal)
        ax.legend(fontsize=8)
    fig.suptitle(
        "Maximum 20 s error against native Simulink (startup included after 20 ms)"
    )
    fig.savefig(output / "native_reference_max_error_vs_timestep.png", dpi=170)
    plt.close(fig)

    fig, axs = plt.subplots(2, 1, figsize=(14, 8), sharex=True, constrained_layout=True)
    for ax, signal in zip(axs, ("Inverter current", "Bus 30 PCC voltage")):
        for color, step in zip(COLORS, STEPS_US):
            part = time_errors[
                (time_errors.method == "DPsim CS")
                & (time_errors.step_us == step)
                & (time_errors.signal == signal)
            ]
            other = time_errors[
                (time_errors.method == "Simulink FMU CS")
                & (time_errors.step_us == step)
                & (time_errors.signal == signal)
            ]
            merged = part.merge(
                other,
                on=["step_us", "signal", "time_s"],
                suffixes=("_dpsim", "_simulink"),
            )
            ax.plot(
                merged.time_s,
                merged.phase_a_rms_error_percent_dpsim
                - merged.phase_a_rms_error_percent_simulink,
                color=color,
                label=f"{step} µs",
            )
        ax.axhline(0, color="black", linewidth=0.7)
        ax.axvline(10, color="0.25", linestyle="--")
        ax.set_ylabel(f"Difference in 0.5 s error envelope\n{signal} [% of native RMS]")
        ax.legend(ncol=4, fontsize=8)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle(
        "DPsim CS versus Simulink FMU CS error envelopes, both relative to native Simulink"
    )
    fig.savefig(output / "simulink_cs_vs_dpsim_cs_error_vs_time.png", dpi=170)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.7), constrained_layout=True)
    for ax, signal in zip(axs, ("Inverter current", "Bus 30 PCC voltage")):
        part = direct[direct.signal == signal]
        ax.plot(
            part.step_us,
            part.rms_error_percent_of_native_rms,
            marker="o",
            label="RMS trace difference",
        )
        ax.plot(
            part.step_us,
            part.max_abs_phase_percent_of_native_rms,
            marker="x",
            linestyle="--",
            label="Maximum phase difference",
        )
        ax.set_xscale("log", base=2)
        ax.set_xticks(STEPS_US, labels=[str(x) for x in STEPS_US])
        ax.set_xlabel("Timestep [µs]")
        ax.set_ylabel("Simulink CS − DPsim CS [% of native RMS]")
        ax.set_title(signal)
        ax.legend(fontsize=8)
    fig.suptitle("Direct trace difference: Simulink FMU CS versus DPsim CS")
    fig.savefig(output / "simulink_cs_vs_dpsim_cs_error_vs_timestep.png", dpi=170)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for ax, signal in zip(axs, signal_fields):
        part = fine_metrics[fine_metrics.signal == signal]
        for method in ("Native Simulink", "Simulink FMU CS", "DPsim CS", "DPsim ME"):
            group = part[part.method == method].sort_values("step_us")
            ax.plot(
                group.step_us,
                group.rms_error_percent_of_fine_me_rms,
                marker="o",
                label=method,
            )
        ax.set_xscale("log", base=2)
        ax.set_xticks(STEPS_US, labels=[str(x) for x in STEPS_US])
        ax.set_xlabel("Timestep [µs]")
        ax.set_ylabel("RMS error [% of 25 µs DPsim ME RMS]")
        ax.set_title(signal)
        ax.legend(fontsize=8)
    fig.suptitle(
        "Closed-loop trajectory convergence against 25 µs DPsim ME (fine-step reference)"
    )
    fig.savefig(output / "fine_me25_reference_error_vs_timestep.png", dpi=170)
    plt.close(fig)

    # Independent, close-up waveforms around the 10 s event for the middle 100 µs case.
    case = cases[100]
    t = case["time"]
    mask = (t >= 9.9) & (t <= 10.1)
    fig, axs = plt.subplots(2, 1, figsize=(13, 7), sharex=True, constrained_layout=True)
    for name, trace in case["traces"].items():
        axs[0].plot(t[mask], trace["current"][mask, 0], label=name, linewidth=1)
        axs[1].plot(t[mask], trace["bus30"][mask, 0], label=name, linewidth=1)
    axs[0].set_ylabel("Phase-a inverter current [A]")
    axs[1].set_ylabel("Bus 30 phase-a voltage [V]")
    axs[1].set_xlabel("Time [s]")
    axs[0].legend(ncol=2, fontsize=8)
    fig.suptitle("10 s, 100 MW load step: closed-loop phase-a waveforms at 100 µs")
    fig.savefig(output / "load_step_waveforms_100us.png", dpi=170)
    plt.close(fig)

    # DPsim component and solver details are retained as compact, traceable tables.
    profiles = []
    for step in STEPS_US:
        for method in ("me", "cs"):
            prefix = f"IEEE39_Inverter_{method}"
            folder = data / DP_SWEEP / f"dt-{step:03d}us" / method / "logs" / prefix
            task = pd.read_csv(folder / f"{prefix}_TaskTimes.csv")
            task["step_us"] = step
            task["method"] = method.upper()
            profiles.append(task)
    task_profiles = pd.concat(profiles, ignore_index=True)
    task_profiles.to_csv(
        output / "dpsim_component_task_timings_all_steps.csv", index=False
    )
    phase_profiles = []
    for step in STEPS_US:
        for method in ("me", "cs"):
            prefix = f"IEEE39_Inverter_{method}"
            folder = data / DP_SWEEP / f"dt-{step:03d}us" / method / "logs" / prefix
            phase = pd.read_csv(folder / f"{prefix}_SolverPhaseTimes.csv")
            phase["step_us"] = step
            phase["method"] = method.upper()
            phase_profiles.append(phase)
    solver_profiles = pd.concat(phase_profiles, ignore_index=True)
    solver_profiles.to_csv(
        output / "dpsim_solver_phase_timings_all_steps.csv", index=False
    )
    newton_rows = []
    for step in STEPS_US:
        path = (
            data
            / DP_SWEEP
            / f"dt-{step:03d}us"
            / "me"
            / "logs"
            / "IEEE39_Inverter_me"
            / "IEEE39_Inverter_me.csv"
        )
        frame = pd.read_csv(path, skipinitialspace=True)
        frame.columns = frame.columns.str.strip()
        active = frame[frame["time"] > 0]
        iterations = active["newton_iterations"].to_numpy(dtype=float)
        residual = active["newton_residual"].to_numpy(dtype=float)
        newton_rows.append(
            {
                "step_us": step,
                "network_steps": int(round(20.0 / (step * 1e-6))),
                "solver_samples_observed": len(active),
                "mean_newton_iterations": float(np.mean(iterations)),
                "max_newton_iterations": int(np.max(iterations)),
                "steps_with_2_iterations": int(np.count_nonzero(iterations == 2)),
                "steps_over_normalized_residual_tolerance_1": int(
                    np.count_nonzero(residual > 1)
                ),
                "max_normalized_newton_residual": float(np.max(residual)),
                "acceptance_threshold": 1.0,
                "fmi_evaluation_requests_total": int(
                    frame["fmi_evaluation_requests"].max()
                ),
                "fmi_evaluation_requests_per_network_step": float(
                    frame["fmi_evaluation_requests"].max() / round(20.0 / (step * 1e-6))
                ),
                "fmi_evaluation_wall_time_total_s": float(
                    frame["fmi_evaluation_wall_time_s"].max()
                ),
                "max_fmi_internal_steps": int(frame["fmi_internal_steps"].max()),
            }
        )
    newton_summary = pd.DataFrame(newton_rows)
    newton_summary.to_csv(output / "dpsim_me_newton_iteration_summary.csv", index=False)
    return {
        "root": root,
        "output": output,
        "metrics": metrics,
        "time_errors": time_errors,
        "direct": direct,
        "timings": timings,
        "fine_metrics": fine_metrics,
        "initial_metrics": initial_metrics,
        "task_profiles": task_profiles,
        "solver_profiles": solver_profiles,
        "newton_summary": newton_summary,
        "cases": cases,
    }


if __name__ == "__main__":
    result = build()
    print(f"Wrote native-reference analysis to {result['output']}")
    print(
        result["metrics"]
        .groupby(["method", "signal"])["rms_error_percent_of_native_rms"]
        .max()
        .to_string()
    )
