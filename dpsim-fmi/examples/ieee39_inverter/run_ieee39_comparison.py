#!/usr/bin/env python3
"""Run and analyze the IEEE39 averaged-inverter ME/CS/Simulink comparison."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = next(
    (
        candidate
        for candidate in (ROOT / "build-fmi-cosim", ROOT / "build-fmi")
        if (candidate / "CMakeCache.txt").is_file()
    ),
    ROOT / "build-fmi-cosim",
)
DATA = BUILD / "ieee39-inverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
RUNNER = BUILD / "dpsim-fmi" / "FMI_IEEE39_AvgInverter"
DT = 50e-6
STOP = 5.0
LOAD_STEP = STOP / 2
P_SET = 125e6
Q_SET = 80.881e6
CURRENT_BASE_PEAK = np.sqrt(2) * np.hypot(P_SET, Q_SET) / (np.sqrt(3) * 345e3)
VOLTAGE_BASE_PEAK = np.sqrt(2 / 3) * 345e3


def load_dpsim(mode: str) -> pd.DataFrame:
    folder = DATA / f"run-{mode}"
    path = folder / "logs" / f"IEEE39_Inverter_{mode}" / f"IEEE39_Inverter_{mode}.csv"
    frame = pd.read_csv(path, skipinitialspace=True)
    frame.columns = frame.columns.str.strip()
    if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ValueError(f"Non-finite values found in {path}")
    if len(frame) != int(round(STOP / DT)) + 1 or not np.isclose(
        frame.time.iloc[-1], STOP, atol=DT * 1e-6
    ):
        raise ValueError(f"Unexpected final time in {path}: {frame.time.iloc[-1]}")
    return frame


def load_native() -> pd.DataFrame:
    path = DATA / "native_simulink_reference.csv"
    frame = pd.read_csv(path)
    if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ValueError(f"Non-finite values found in {path}")
    return frame


def vector(frame: pd.DataFrame, prefix: str, phases=(0, 1, 2)) -> np.ndarray:
    return np.column_stack([frame[f"{prefix}{phase}"] for phase in phases])


def metrics(
    error: np.ndarray, reference: np.ndarray, mask: np.ndarray, base_peak: float
) -> dict:
    e = error[mask]
    r = reference[mask]
    return {
        "samples": int(mask.sum()),
        "max_abs": float(np.max(np.abs(e))),
        "rms": float(np.sqrt(np.mean(e * e))),
        "reference_rms": float(np.sqrt(np.mean(r * r))),
        "max_abs_percent_base": float(100 * np.max(np.abs(e)) / base_peak),
    }


def timing_artifacts(mode: str) -> dict[str, Path]:
    folder = DATA / f"run-{mode}"
    profile_dir = folder / "logs" / f"IEEE39_Inverter_{mode}"
    task_csv = profile_dir / f"IEEE39_Inverter_{mode}_TaskTimes.csv"
    solver_matches = list(profile_dir.glob("*_SolverPhaseTimes.csv"))
    paths = {
        "component_task_profile_csv": task_csv,
        "solver_phase_profile_csv": (
            solver_matches[0] if len(solver_matches) == 1 else None
        ),
        "scheduler_step_profile_csv": profile_dir
        / f"IEEE39_Inverter_{mode}_TaskTimes_SchedulerSteps.csv",
        "simulation_step_profile_csv": profile_dir
        / f"IEEE39_Inverter_{mode}_SimulationStepTimes.csv",
    }
    missing = [
        str(path) for path in paths.values() if path is None or not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Missing detailed timing artifacts: " + ", ".join(missing)
        )
    return paths


def run_case(mode: str) -> dict:
    folder = DATA / f"run-{mode}"
    folder.mkdir(parents=True, exist_ok=True)
    log = folder / "runner.log"
    started = time.perf_counter()
    result = subprocess.run(
        [str(RUNNER), mode, str(FMU), str(STOP), str(DT), str(LOAD_STEP)],
        cwd=folder,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    wall = time.perf_counter() - started
    log.write_text(result.stdout)
    fields = {
        name: re.search(rf"{name}=([0-9.eE+-]+)", result.stdout)
        for name in ("runtime_s", "initialization_s", "simulation_s")
    }
    if any(match is None for match in fields.values()):
        raise RuntimeError(f"Could not parse simulation time from {log}")
    profile = timing_artifacts(mode)
    timing = {
        "mode": mode,
        "total_runtime_s": float(fields["runtime_s"].group(1)),
        "initialization_runtime_s": float(fields["initialization_s"].group(1)),
        "simulation_runtime_s": float(fields["simulation_s"].group(1)),
        "process_wall_time_s": wall,
        "communication_step_s": DT,
        "stop_time_s": STOP,
        "load_step_time_s": LOAD_STEP,
        "network_steps": int(round(STOP / DT)),
        "log": str(log),
        **{name: str(path) for name, path in profile.items()},
    }
    (folder / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
    return timing


def run_all() -> None:
    subprocess.run(
        ["bash", str(CASE / "sources" / "build_fmu.sh"), str(FMU.parent)],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        ["cmake", "--build", str(BUILD), "--target", "FMI_IEEE39_AvgInverter", "-j4"],
        cwd=ROOT,
        check=True,
    )
    timings = {mode: run_case(mode) for mode in ("me", "cs")}
    native_out = DATA / "native_simulink_reference.csv"
    reference_input = (
        DATA / "run-me" / "logs" / "IEEE39_Inverter_me" / "IEEE39_Inverter_me.csv"
    )
    matlab_dir = CASE / "matlab"
    matlab_expr = (
        f"addpath('{matlab_dir}'); "
        # One reference simulation is used for both DPsim overlays and errors.
        # Discard the first run as MATLAB/Simulink compile/cache warm-up.
        f"run_native_replay('{reference_input}','{native_out}',0.05,false); "
        f"run_native_replay('{reference_input}','{native_out}');"
    )
    timing_path = DATA / "native_timing.json"
    timing_path.touch()
    start_ns = time.time_ns()
    start = time.perf_counter()
    process = subprocess.Popen(
        ["matlab", "-batch", matlab_expr],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    deadline = time.monotonic() + 1800
    while time.monotonic() < deadline:
        status = process.poll()
        if status is not None:
            output = process.communicate()[0]
            if status != 0:
                raise subprocess.CalledProcessError(status, process.args, output=output)
            break
        if timing_path.exists() and timing_path.stat().st_mtime_ns > start_ns:
            # MATLAB writes the manifest after Simulink has completed. Its
            # headless UI cleanup can remain alive, so stop it after results are safe.
            process.terminate()
            try:
                output = process.communicate(timeout=10)[0]
            except subprocess.TimeoutExpired:
                process.kill()
                output = process.communicate()[0]
            break
        time.sleep(1)
    else:
        process.kill()
        output = process.communicate()[0]
        raise TimeoutError(
            "MATLAB did not finish the native Simulink reference within 30 minutes"
        )
    wall = time.perf_counter() - start
    (DATA / "native_batch.log").write_text(output)
    timings["native_replay_batch_wall_time_s"] = wall
    timings["native_simulink_warmup_runs"] = 1
    timings["timing_profile_overhead"] = (
        "Both DPsim modes use the same SequentialScheduler task profiler and "
        "solver phase profiler; reported runtime includes this instrumentation."
    )
    timings["native_simulink"] = json.loads(timing_path.read_text())
    native = timings["native_simulink"]
    if not np.isclose(native["fixed_step_s"], DT, rtol=0, atol=1e-12):
        raise ValueError("Native Simulink timestep does not match DPsim")
    if not np.isclose(native["stop_time_s"], STOP, rtol=0, atol=DT * 1e-6):
        raise ValueError("Native Simulink duration does not match DPsim")
    if native["sample_count"] != int(round(STOP / DT)) + 1:
        raise ValueError(
            "Native Simulink sample count does not match the DPsim time grid"
        )
    (DATA / "timings.json").write_text(json.dumps(timings, indent=2) + "\n")


def write_timing_tables(timings: dict) -> None:
    component_rows, solver_rows, step_frames, simulation_step_frames = [], [], [], []
    for mode in ("me", "cs"):
        profile = timings[mode]
        tasks = pd.read_csv(profile["component_task_profile_csv"])
        tasks["mode"] = mode.upper()
        tasks["category"] = np.select(
            [
                tasks.task.str.endswith(".MnaPreStep"),
                tasks.task.str.endswith(".MnaPostStep"),
                tasks.task.str.endswith(".Solve"),
            ],
            ["component_pre_step", "component_post_step", "solver_task"],
            default="other_task",
        )
        tasks["avg_us"] = tasks.avg_ns / 1e3
        tasks["total_ms"] = tasks.total_ns / 1e6
        component_rows.append(tasks.rename(columns={"task": "component_task"}))

        phases = pd.read_csv(profile["solver_phase_profile_csv"])
        phases["mode"] = mode.upper()
        solver_rows.append(phases)

        per_step = pd.read_csv(profile["scheduler_step_profile_csv"])
        per_step.insert(0, "mode", mode.upper())
        step_frames.append(per_step)
        simulation_steps = pd.read_csv(profile["simulation_step_profile_csv"])
        simulation_steps.insert(0, "mode", mode.upper())
        simulation_step_frames.append(simulation_steps)

    pd.concat(component_rows, ignore_index=True).to_csv(
        DATA / "component_timing_breakdown.csv", index=False
    )
    pd.concat(solver_rows, ignore_index=True).to_csv(
        DATA / "solver_timing_breakdown.csv", index=False
    )
    pd.concat(step_frames, ignore_index=True).to_csv(
        DATA / "scheduler_step_timing.csv", index=False
    )
    pd.concat(simulation_step_frames, ignore_index=True).to_csv(
        DATA / "simulation_step_timing.csv", index=False
    )


def analyze() -> dict:
    me, cs = load_dpsim("me"), load_dpsim("cs")
    native = load_native()
    time_me = me.time.to_numpy()
    time_cs = cs.time.to_numpy()
    if not np.allclose(time_me, time_cs, atol=1e-12, rtol=0):
        raise ValueError("ME and CS outputs are not on the same time grid")
    t = time_me
    masks = {
        "post_startup": t >= 0.02,
        "pre_event": (t >= 0.02) & (t < LOAD_STEP),
        "load_step_window": (t >= LOAD_STEP) & (t <= LOAD_STEP + 0.02),
        "recovery": (t > LOAD_STEP + 0.02) & (t <= STOP),
    }
    current_me = vector(me, "inverter.i_")
    current_cs = vector(cs, "inverter.i_")
    native_i = np.column_stack(
        [np.interp(t, native.time, native[f"i_{phase}"]) for phase in "abc"]
    )
    native_v = np.column_stack(
        [np.interp(t, native.time, native[f"v_{phase}"]) for phase in "abc"]
    )
    current_diffs = {
        "me_vs_native_simulink": current_me - native_i,
        "cs_vs_native_simulink": current_cs - native_i,
        "cs_vs_me": current_cs - current_me,
    }
    voltage_me = vector(me, "inverter.v_")
    voltage_cs = vector(cs, "inverter.v_")
    voltage_diffs = {
        "me_vs_native_simulink": voltage_me - native_v,
        "cs_vs_native_simulink": voltage_cs - native_v,
    }
    metrics_out = {
        "system": {
            "name": "IEEE 39 bus, three-phase EMT",
            "buses": 39,
            "branches": 46,
            "frequency_Hz": 60,
            "network_step_s": DT,
            "stop_time_s": STOP,
            "network_steps": int(round(STOP / DT)),
            "domain_time_grid": f"DPsim ME, DPsim CS, and native Simulink all use 50 us over {STOP:g} s ({int(round(STOP / DT)):,} intervals).",
            "disturbance": f"100 MW balanced load step at bus 39 at {LOAD_STEP:.3f} s",
            "device": "average PQ-controlled current source at bus 30; half of original generator-30 P/Q is retained as synchronous generation",
            "device_setpoint": {
                "P_W": P_SET,
                "Q_var": Q_SET,
                "current_response_time_s": 0.005,
                "current_limit_pu": 1.2,
            },
            "current_base_peak_A": CURRENT_BASE_PEAK,
            "native_reference_scope": "one native Simulink inverter run driven by the DPsim ME PCC-voltage trace; both DPsim current traces and both terminal-voltage traces are compared against this same reference",
        },
        "timings": (
            json.loads((DATA / "timings.json").read_text())
            if (DATA / "timings.json").exists()
            else {}
        ),
        "current_error": {},
        "terminal_voltage_error": {},
        "bus_voltage_me_vs_cs": {},
        "coupling_diagnostics": {},
    }
    if metrics_out["timings"].get("me", {}).get("component_task_profile_csv"):
        write_timing_tables(metrics_out["timings"])
    for name, err in current_diffs.items():
        metrics_out["current_error"][name] = {
            window: metrics(
                err,
                native_i if "native_simulink" in name else current_me,
                mask,
                CURRENT_BASE_PEAK,
            )
            for window, mask in masks.items()
            if mask.any()
        }
    for name, err in voltage_diffs.items():
        metrics_out["terminal_voltage_error"][name] = {
            window: metrics(err, native_v, mask, VOLTAGE_BASE_PEAK)
            for window, mask in masks.items()
            if mask.any()
        }
    for bus in (1, 10, 20, 30, 31, 39):
        name = f"N{bus}.V_"
        v_me = vector(me, name)
        v_cs = vector(cs, name)
        err = v_cs - v_me
        metrics_out["bus_voltage_me_vs_cs"][str(bus)] = {
            window: metrics(err, v_me, mask, VOLTAGE_BASE_PEAK)
            for window, mask in masks.items()
            if mask.any()
        }
    metrics_out["coupling_diagnostics"] = {
        "me_max_newton_iterations": int(me.newton_iterations.max()),
        "me_mean_newton_iterations": float(me.newton_iterations[t >= 0.02].mean()),
        "me_max_reported_residual": float(me.newton_residual.max()),
        "me_max_internal_steps_per_macro_step": int(me.fmi_internal_steps.max()),
        "me_fmi_evaluation_requests_total": int(me.fmi_evaluation_requests.iloc[-1]),
        "me_fmi_evaluation_time_total_s": float(me.fmi_evaluation_wall_time_s.iloc[-1]),
        "me_fmi_evaluations_per_network_step": float(
            me.fmi_evaluation_requests.iloc[-1] / round(STOP / DT)
        ),
        "me_fmi_evaluation_time_per_network_step_us": float(
            1e6 * me.fmi_evaluation_wall_time_s.iloc[-1] / round(STOP / DT)
        ),
        "cs_do_step_calls": int(cs.fmi_do_steps.max()),
        "interpretation": "Ordinary ME network steps use DPsim MNA Newton convergence with scaled state/current residuals <= 1; SinglePassExplicit only handles macrosteps crossing FMI internal time events. CS performs one fmi3DoStep per 50 us network step and no master-side Newton iterations.",
    }
    json_path = DATA / "ieee39_metrics.json"
    csv_path = DATA / "ieee39_metrics.csv"
    json_path.write_text(json.dumps(metrics_out, indent=2) + "\n")
    flat = []
    for comparison, by_window in metrics_out["current_error"].items():
        for window, values in by_window.items():
            flat.append({"comparison": comparison, "window": window, **values})
    for bus, by_window in metrics_out["bus_voltage_me_vs_cs"].items():
        for window, values in by_window.items():
            flat.append(
                {"comparison": f"CS_vs_ME_bus{bus}_voltage", "window": window, **values}
            )
    for comparison, by_window in metrics_out["terminal_voltage_error"].items():
        for window, values in by_window.items():
            flat.append({"comparison": comparison, "window": window, **values})
    pd.DataFrame(flat).to_csv(csv_path, index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    rms_window_s = 0.5  # exactly 30 cycles at 60 Hz, suppresses carrier ripple
    rms_window_samples = int(round(rms_window_s / DT))

    def moving_rms(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if values.ndim == 1:
            mean_square = (
                pd.Series(values * values)
                .rolling(
                    rms_window_samples, center=True, min_periods=rms_window_samples
                )
                .mean()
                .to_numpy()
            )
        else:
            mean_square = (
                pd.Series(np.sum(values * values, axis=1))
                .rolling(
                    rms_window_samples, center=True, min_periods=rms_window_samples
                )
                .mean()
                .to_numpy()
            )
        return np.sqrt(mean_square)

    native_v_on_grid = np.column_stack(
        [np.interp(t, native.time, native[f"v_{phase}"]) for phase in "abc"]
    )
    native_i_on_grid = np.column_stack(
        [np.interp(t, native.time, native[f"i_{phase}"]) for phase in "abc"]
    )
    rms_terminal_me = moving_rms(voltage_me[:, 0]) / 1e3
    rms_terminal_cs = moving_rms(voltage_cs[:, 0]) / 1e3
    rms_terminal_native = moving_rms(native_v_on_grid[:, 0]) / 1e3
    rms_bus39_me = moving_rms(vector(me, "N39.V_")[:, 0]) / 1e3
    rms_bus39_cs = moving_rms(vector(cs, "N39.V_")[:, 0]) / 1e3
    rms_current_me = moving_rms(current_me)
    rms_current_cs = moving_rms(current_cs)
    rms_current_native = moving_rms(native_i_on_grid)
    event_zoom = (t >= LOAD_STEP - 0.01) & (t <= LOAD_STEP + 0.03)

    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    ax[0, 0].plot(t, rms_terminal_me, label="DPsim ME")
    ax[0, 0].plot(t, rms_terminal_cs, "--", label="DPsim CS")
    ax[0, 0].plot(t, rms_terminal_native, ":", label="Native Simulink input")
    ax[0, 0].set(
        title="Inverter terminal phase-A voltage (0.5 s moving RMS)",
        xlabel="Time [s]",
        ylabel="Voltage [kV RMS]",
    )
    ax[0, 0].ticklabel_format(axis="y", style="plain", useOffset=False)
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot(t, rms_bus39_me, label="ME")
    ax[0, 1].plot(t, rms_bus39_cs, "--", label="CS")
    ax[0, 1].set(
        title="Bus 39 phase-A voltage (0.5 s moving RMS)",
        xlabel="Time [s]",
        ylabel="Voltage [kV RMS]",
    )
    ax[0, 1].legend()
    ax[1, 0].plot(t, rms_current_me, label="ME")
    ax[1, 0].plot(t, rms_current_cs, "--", label="CS")
    ax[1, 0].plot(t, rms_current_native, ":", label="Native Simulink")
    ax[1, 0].set(
        title="Converter current (three-phase 0.5 s moving RMS)",
        xlabel="Time [s]",
        ylabel="Current [A RMS]",
    )
    ax[1, 0].legend()
    ax[1, 1].plot(t[event_zoom], current_me[event_zoom, 0], label="ME")
    ax[1, 1].plot(t[event_zoom], current_cs[event_zoom, 0], "--", label="CS")
    ax[1, 1].plot(
        t[event_zoom], native_i_on_grid[event_zoom, 0], ":", label="Native Simulink"
    )
    ax[1, 1].set(
        title="Phase-A current around 2.5 s load step",
        xlabel="Time [s]",
        ylabel="Current [A]",
    )
    ax[1, 1].legend(fontsize=8)
    for panel in ax.flat:
        panel.axvline(LOAD_STEP, color="black", alpha=0.5, linestyle=":")
    fig.suptitle("IEEE 39 bus averaged inverter: ME, CS, native Simulink", fontsize=14)
    fig.savefig(DATA / "ieee39_waveforms.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    ax[0, 0].plot(
        t, moving_rms(voltage_diffs["me_vs_native_simulink"]), label="ME vs Simulink"
    )
    ax[0, 0].plot(
        t, moving_rms(voltage_diffs["cs_vs_native_simulink"]), label="CS vs Simulink"
    )
    ax[0, 0].set(
        title="Inverter terminal-voltage error (0.5 s moving RMS)",
        xlabel="Time [s]",
        ylabel="Three-phase error [V RMS]",
    )
    ax[0, 0].legend()
    ax[0, 1].plot(
        t, moving_rms(current_diffs["me_vs_native_simulink"]), label="ME vs Simulink"
    )
    ax[0, 1].plot(
        t, moving_rms(current_diffs["cs_vs_native_simulink"]), label="CS vs Simulink"
    )
    ax[0, 1].set(
        title="Converter-current error (0.5 s moving RMS)",
        xlabel="Time [s]",
        ylabel="Three-phase error [A RMS]",
    )
    ax[0, 1].legend()
    ax[1, 0].plot(
        t[event_zoom],
        current_diffs["me_vs_native_simulink"][event_zoom, 0],
        label="ME minus Simulink",
    )
    ax[1, 0].plot(
        t[event_zoom],
        current_diffs["cs_vs_native_simulink"][event_zoom, 0],
        label="CS minus Simulink",
    )
    ax[1, 0].set(
        title="Phase-A current error around load step",
        xlabel="Time [s]",
        ylabel="Error [A]",
    )
    ax[1, 0].legend(fontsize=8)
    ax[1, 1].plot(t[event_zoom], current_me[event_zoom, 0], label="ME")
    ax[1, 1].plot(t[event_zoom], current_cs[event_zoom, 0], "--", label="CS")
    ax[1, 1].plot(
        t[event_zoom], native_i_on_grid[event_zoom, 0], ":", label="Native Simulink"
    )
    ax[1, 1].set(
        title="Current overlay around load step",
        xlabel="Time [s]",
        ylabel="Phase-A current [A]",
    )
    ax[1, 1].legend(fontsize=8)
    for panel in ax.flat:
        panel.axvline(LOAD_STEP, color="black", alpha=0.5, linestyle=":")
    fig.suptitle("IEEE 39 bus comparison errors", fontsize=14)
    fig.savefig(DATA / "ieee39_errors.png", dpi=180)
    plt.close(fig)
    print(
        json.dumps(
            {
                "metrics": str(json_path),
                "csv": str(csv_path),
                "waveforms": str(DATA / "ieee39_waveforms.png"),
                "errors": str(DATA / "ieee39_errors.png"),
            },
            indent=2,
        )
    )
    return metrics_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        action="store_true",
        help="rerun DPsim ME/CS and native Simulink references",
    )
    args = parser.parse_args()
    if args.run:
        run_all()
    analyze()


if __name__ == "__main__":
    main()
