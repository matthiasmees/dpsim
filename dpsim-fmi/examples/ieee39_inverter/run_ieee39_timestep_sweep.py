#!/usr/bin/env python3
"""Run matched DPsim ME/CS pairs over a communication-step sweep."""
from __future__ import annotations

import json
import argparse
import re
import subprocess
import time
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
OUT_ROOT = BUILD / "ieee39-inverter"
RUNNER = BUILD / "dpsim-fmi" / "FMI_IEEE39_AvgInverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
STEP_US = (25, 50, 100, 200, 500)
STOP_S = 5.0
LOAD_STEP_S = 2.5
OUT = OUT_ROOT / "timestep-sweep"
CURRENT_BASE_PEAK_A = np.sqrt(2) * np.hypot(125e6, 80.881e6) / (np.sqrt(3) * 345e3)


def load_run(folder: Path, mode: str, dt_s: float) -> pd.DataFrame:
    csv = folder / "logs" / f"IEEE39_Inverter_{mode}" / f"IEEE39_Inverter_{mode}.csv"
    frame = pd.read_csv(csv, skipinitialspace=True)
    frame.columns = frame.columns.str.strip()
    expected = int(round(STOP_S / dt_s)) + 1
    # Simulation logs may omit the final endpoint row while still executing
    # exactly STOP_S / dt_s solver steps. Accept [0, STOP-dt] or [0, STOP].
    if (
        len(frame) not in (expected, expected - 1)
        or frame.time.iloc[-1] < STOP_S - 1.01 * dt_s
        or frame.time.iloc[-1] > STOP_S + dt_s * 1e-6
    ):
        raise ValueError(
            f"{csv}: expected {expected} or {expected-1} samples over a {STOP_S} s run"
        )
    if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ValueError(f"{csv}: non-finite simulation values")
    # DPsim writes timestamps with six decimal places. For sub-millisecond
    # steps, formatting alone can perturb adjacent logged intervals by 5 us.
    time_tolerance = max(dt_s * 1e-6, 5.1e-6)
    if not np.allclose(np.diff(frame.time), dt_s, rtol=0, atol=time_tolerance):
        raise ValueError(f"{csv}: time grid differs from requested {dt_s:g} s")
    return frame


def phase_vector(frame: pd.DataFrame, name: str) -> np.ndarray:
    return np.column_stack([frame[f"{name}{phase}"] for phase in range(3)])


def summarize(
    error: np.ndarray, mask: np.ndarray, base: float
) -> dict[str, float | int]:
    selected = error[mask]
    return {
        "samples": int(mask.sum()),
        "max_abs_phase": float(np.max(np.abs(selected))),
        "rms_all_phases": float(np.sqrt(np.mean(selected * selected))),
        "max_three_phase_norm": float(np.max(np.linalg.norm(selected, axis=1))),
        "max_abs_percent_base": float(100 * np.max(np.abs(selected)) / base),
    }


def run_pair(step_us: int) -> tuple[dict, pd.DataFrame]:
    dt_s = step_us * 1e-6
    pair_dir = OUT / f"dt-{step_us:03d}us"
    pair_dir.mkdir(parents=True, exist_ok=True)
    frames = {}
    run_timing = {}
    for mode in ("me", "cs"):
        folder = pair_dir / mode
        folder.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        result = subprocess.run(
            [str(RUNNER), mode, str(FMU), str(STOP_S), str(dt_s), str(LOAD_STEP_S)],
            cwd=folder,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        elapsed = time.perf_counter() - started
        (folder / "runner.log").write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f"{step_us} us {mode} failed:\n{result.stdout}")
        parsed = {
            key: re.search(rf"{key}=([0-9.eE+-]+)", result.stdout)
            for key in ("runtime_s", "initialization_s", "simulation_s", "steps")
        }
        if any(v is None for v in parsed.values()):
            raise RuntimeError(
                f"Could not parse {mode} timing from {folder}/runner.log"
            )
        run_timing[mode] = {
            "total_runtime_s": float(parsed["runtime_s"].group(1)),
            "initialization_runtime_s": float(parsed["initialization_s"].group(1)),
            "simulation_runtime_s": float(parsed["simulation_s"].group(1)),
            "process_wall_time_s": elapsed,
        }
        frames[mode] = load_run(folder, mode, dt_s)

    if not np.allclose(
        frames["me"].time, frames["cs"].time, rtol=0, atol=max(dt_s * 1e-6, 5.1e-6)
    ):
        raise ValueError(f"ME and CS time grids differ at {step_us} us")
    t = frames["me"].time.to_numpy()
    current_error = phase_vector(frames["cs"], "inverter.i_") - phase_vector(
        frames["me"], "inverter.i_"
    )
    terminal_voltage_error = phase_vector(frames["cs"], "inverter.v_") - phase_vector(
        frames["me"], "inverter.v_"
    )
    bus39_voltage_error = phase_vector(frames["cs"], "N39.V_") - phase_vector(
        frames["me"], "N39.V_"
    )
    windows = {
        "post_startup": t >= 0.02,
        "pre_event": (t >= 0.02) & (t < LOAD_STEP_S),
        "load_step_0_to_20ms": (t >= LOAD_STEP_S) & (t <= LOAD_STEP_S + 0.02),
        "post_event": t > LOAD_STEP_S + 0.02,
        "final_1s": t >= STOP_S - 1.0,
    }
    summary = {
        "step_us": step_us,
        "step_s": dt_s,
        "stop_time_s": STOP_S,
        "load_step_time_s": LOAD_STEP_S,
        "network_steps": int(round(STOP_S / dt_s)),
        "me_last_logged_time_s": float(frames["me"].time.iloc[-1]),
        "cs_last_logged_time_s": float(frames["cs"].time.iloc[-1]),
        "me_timing": run_timing["me"],
        "cs_timing": run_timing["cs"],
        "current_error_A": {
            k: summarize(current_error, mask, CURRENT_BASE_PEAK_A)
            for k, mask in windows.items()
        },
        "terminal_voltage_error_V": {
            k: summarize(terminal_voltage_error, mask, 345e3)
            for k, mask in windows.items()
        },
        "bus39_voltage_error_V": {
            k: summarize(bus39_voltage_error, mask, 345e3)
            for k, mask in windows.items()
        },
        "me_fmi_evaluation_requests": int(
            frames["me"].fmi_evaluation_requests.iloc[-1]
        ),
        "me_fmi_evaluation_time_s": float(
            frames["me"].fmi_evaluation_wall_time_s.iloc[-1]
        ),
        "cs_do_step_calls": int(frames["cs"].fmi_do_steps.max()),
    }

    # Keep a compact 10 ms RMS error envelope for comparing divergence over 5 s.
    bin_id = np.minimum(np.floor(t / 0.01).astype(int), int(STOP_S / 0.01) - 1)
    envelope = (
        pd.DataFrame(
            {
                "time_s": t,
                "current_error_sq": np.sum(current_error * current_error, axis=1),
                "terminal_voltage_error_sq": np.sum(
                    terminal_voltage_error * terminal_voltage_error, axis=1
                ),
                "bus39_voltage_error_sq": np.sum(
                    bus39_voltage_error * bus39_voltage_error, axis=1
                ),
                "bin": bin_id,
            }
        )
        .groupby("bin", as_index=False)
        .mean()
    )
    envelope["step_us"] = step_us
    envelope["time_s"] = (envelope.bin + 0.5) * 0.01
    envelope["current_error_rms_A"] = np.sqrt(envelope.current_error_sq)
    envelope["terminal_voltage_error_rms_V"] = np.sqrt(
        envelope.terminal_voltage_error_sq
    )
    envelope["bus39_voltage_error_rms_V"] = np.sqrt(envelope.bus39_voltage_error_sq)
    envelope = envelope[
        [
            "step_us",
            "time_s",
            "current_error_rms_A",
            "terminal_voltage_error_rms_V",
            "bus39_voltage_error_rms_V",
        ]
    ]

    event = (t >= LOAD_STEP_S - 0.01) & (t <= LOAD_STEP_S + 0.03)
    event_trace = pd.DataFrame(
        {
            "step_us": step_us,
            "time_s": t[event],
            "phase_a_current_error_A": current_error[event, 0],
            "current_error_norm_A": np.linalg.norm(current_error[event], axis=1),
            "phase_a_terminal_voltage_error_V": terminal_voltage_error[event, 0],
        }
    )
    trace = pd.concat([envelope, event_trace], ignore_index=True, sort=False)
    trace.to_csv(pair_dir / "error_trace.csv", index=False)
    return summary, trace


def flatten_summaries(results: list[dict]) -> pd.DataFrame:
    rows = []
    for result in results:
        if "current_error_A" not in result:
            continue
        for signal, windows in (
            ("current_error_A", result["current_error_A"]),
            ("terminal_voltage_error_V", result["terminal_voltage_error_V"]),
            ("bus39_voltage_error_V", result["bus39_voltage_error_V"]),
        ):
            for window, values in windows.items():
                rows.append(
                    {
                        "step_us": result["step_us"],
                        "step_s": result["step_s"],
                        "network_steps": result["network_steps"],
                        "signal": signal,
                        "window": window,
                        **values,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    global OUT, STOP_S, LOAD_STEP_S
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--steps-us",
        type=int,
        nargs="+",
        default=list(STEP_US),
        help="communication steps to run (default: 25 50 100 200 500)",
    )
    parser.add_argument(
        "--stop-time-s",
        type=float,
        default=STOP_S,
        help="duration of each matched ME/CS run",
    )
    parser.add_argument(
        "--load-step-s",
        type=float,
        default=None,
        help="load step time (default: midpoint of the run)",
    )
    parser.add_argument(
        "--output-name",
        default=None,
        help="output folder name under the build's ieee39-inverter directory",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed timestep pairs and aggregate the full sweep",
    )
    args = parser.parse_args()
    STOP_S = args.stop_time_s
    LOAD_STEP_S = STOP_S / 2 if args.load_step_s is None else args.load_step_s
    if STOP_S <= 0 or not 0 < LOAD_STEP_S < STOP_S:
        parser.error(
            "stop time must be positive and load step must fall inside the run"
        )
    if any(step <= 0 for step in args.steps_us):
        parser.error("all communication steps must be positive")
    OUT = OUT_ROOT / (
        args.output_name
        or ("timestep-sweep" if STOP_S == 5.0 else f"timestep-sweep-{STOP_S:g}s")
    )
    if not RUNNER.is_file() or not FMU.is_file():
        raise FileNotFoundError(
            "Build the FMI runner/FMU first with run_ieee39_comparison.py --run"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    old_trace_path = OUT / "error_traces.csv"
    old_traces = (
        pd.read_csv(old_trace_path) if old_trace_path.is_file() else pd.DataFrame()
    )
    for step_us in args.steps_us:
        pair_dir = OUT / f"dt-{step_us:03d}us"
        metrics_path = pair_dir / "metrics.json"
        trace_path = pair_dir / "error_trace.csv"
        if args.reuse and metrics_path.is_file() and trace_path.is_file():
            print(f"Reusing completed {step_us} us pair", flush=True)
            continue
        print(
            f"Running matched ME/CS pair at {step_us} us for {STOP_S:g} s; event at {LOAD_STEP_S:g} s",
            flush=True,
        )
        try:
            summary, data = run_pair(step_us)
            summary["status"] = "complete"
            metrics_path.write_text(json.dumps(summary, indent=2) + "\n")
            (pair_dir / "failure.json").unlink(missing_ok=True)
        except Exception as error:
            failure = {
                "step_us": step_us,
                "step_s": step_us * 1e-6,
                "stop_time_s": STOP_S,
                "load_step_time_s": LOAD_STEP_S,
                "status": "failed",
                "error": str(error),
            }
            (pair_dir / "failure.json").write_text(json.dumps(failure, indent=2) + "\n")
            print(f"Step {step_us} us failed: {error}", flush=True)

    results = []
    envelopes = []
    for step_us in args.steps_us:
        pair_dir = OUT / f"dt-{step_us:03d}us"
        metrics_path = pair_dir / "metrics.json"
        trace_path = pair_dir / "error_trace.csv"
        if metrics_path.is_file():
            summary = json.loads(metrics_path.read_text())
            results.append(summary)
            if trace_path.is_file():
                envelopes.append(pd.read_csv(trace_path))
            elif not old_traces.empty:
                envelopes.append(old_traces[old_traces.step_us == step_us])
        elif (pair_dir / "failure.json").is_file():
            results.append(json.loads((pair_dir / "failure.json").read_text()))
    if not any("current_error_A" in result for result in results):
        raise RuntimeError("No completed timestep pairs are available for analysis")
    (OUT / "sweep_metrics.json").write_text(
        json.dumps(
            {
                "description": "Matched closed-loop IEEE39 ME and CS simulations; single-pass ME coupling; compare each CS trace with ME at identical dt.",
                "step_us": args.steps_us,
                "stop_time_s": STOP_S,
                "load_step_time_s": LOAD_STEP_S,
                "load_step_MW": 100,
                "results": results,
            },
            indent=2,
        )
        + "\n"
    )
    pd.concat(envelopes, ignore_index=True).to_csv(
        OUT / "error_traces.csv", index=False
    )
    flat = flatten_summaries(results)
    flat.to_csv(OUT / "timestep_sweep_metrics.csv", index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    post = flat[(flat.window == "post_startup")]
    for signal, label in (
        ("current_error_A", "Current"),
        ("terminal_voltage_error_V", "Inverter terminal voltage"),
        ("bus39_voltage_error_V", "Bus 39 voltage"),
    ):
        part = post[post.signal == signal]
        axis = (
            ax[0, 0]
            if signal == "current_error_A"
            else ax[0, 1] if signal == "terminal_voltage_error_V" else ax[1, 0]
        )
        axis.plot(part.step_us, part.rms_all_phases, marker="o", label=f"{label} RMS")
        axis.plot(
            part.step_us,
            part.max_abs_phase,
            marker="x",
            linestyle="--",
            label=f"{label} max phase",
        )
        axis.set_xscale("log", base=2)
        axis.set_xticks(args.steps_us, labels=[str(x) for x in args.steps_us])
        axis.set_xlabel("Communication step [µs]")
        axis.set_ylabel("Error [A]" if signal == "current_error_A" else "Error [V]")
        axis.set_title(f"CS minus ME: {label}, post-startup")
        axis.legend(fontsize=8)
    ax[1, 1].set_axis_off()
    ax[1, 1].text(
        0.02,
        0.95,
        f"All pairs: {STOP_S:g} s, {min(args.steps_us)} to {max(args.steps_us)} µs sweep.\n"
        f"100 MW bus 39 load step at {LOAD_STEP_S:g} s.\n"
        "Ordinary ME steps use DPsim MNA Newton; CS takes one fmi3DoStep\n"
        "per network step. Errors use each pair's own identical time grid.",
        va="top",
        fontsize=11,
    )
    fig.suptitle("IEEE39 timestep sweep: CS trajectory difference from ME", fontsize=14)
    fig.savefig(OUT / "timestep_sweep_summary.png", dpi=180)
    plt.close(fig)

    traces = pd.concat(envelopes, ignore_index=True)
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), constrained_layout=True)
    env = traces.dropna(subset=["current_error_rms_A"])
    for step_us, group in env.groupby("step_us"):
        ax[0].plot(group.time_s, group.current_error_rms_A, label=f"{step_us} µs")
        ax[1].plot(
            group.time_s, group.terminal_voltage_error_rms_V, label=f"{step_us} µs"
        )
    ax[0].set(
        title="CS–ME three-phase current error, 10 ms RMS bins",
        xlabel="Time [s]",
        ylabel="Current error [A RMS]",
    )
    ax[1].set(
        title="CS–ME inverter terminal-voltage error, 10 ms RMS bins",
        xlabel="Time [s]",
        ylabel="Voltage error [V RMS]",
    )
    for axis in ax:
        axis.axvline(LOAD_STEP_S, color="black", linestyle=":", alpha=0.6)
        axis.legend(ncol=5, fontsize=8)
    fig.savefig(OUT / "timestep_sweep_error_traces.png", dpi=180)
    plt.close(fig)
    print(
        json.dumps(
            {
                "metrics": str(OUT / "timestep_sweep_metrics.csv"),
                "summary_plot": str(OUT / "timestep_sweep_summary.png"),
                "error_plot": str(OUT / "timestep_sweep_error_traces.png"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
