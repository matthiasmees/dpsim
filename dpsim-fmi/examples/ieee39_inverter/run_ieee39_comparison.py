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
BUILD = ROOT / "build-fmi-cosim"
DATA = BUILD / "ieee39-inverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
RUNNER = BUILD / "dpsim-fmi" / "FMI_IEEE39_AvgInverter"
DT = 50e-6
STOP = 0.2
LOAD_STEP = 0.17
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
    if not np.isclose(frame.time.iloc[-1], STOP, atol=DT * 1e-6):
        raise ValueError(f"Unexpected final time in {path}: {frame.time.iloc[-1]}")
    return frame


def load_native(mode: str) -> pd.DataFrame:
    path = DATA / f"run-{mode}" / "native_replay.csv"
    frame = pd.read_csv(path)
    if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ValueError(f"Non-finite values found in {path}")
    return frame


def vector(frame: pd.DataFrame, prefix: str, phases=(0, 1, 2)) -> np.ndarray:
    return np.column_stack([frame[f"{prefix}{phase}"] for phase in phases])


def metrics(error: np.ndarray, reference: np.ndarray, mask: np.ndarray,
            base_peak: float) -> dict:
    e = error[mask]
    r = reference[mask]
    return {
        "samples": int(mask.sum()),
        "max_abs": float(np.max(np.abs(e))),
        "rms": float(np.sqrt(np.mean(e * e))),
        "reference_rms": float(np.sqrt(np.mean(r * r))),
        "max_abs_percent_base": float(100 * np.max(np.abs(e)) / base_peak),
    }


def run_case(mode: str) -> dict:
    folder = DATA / f"run-{mode}"
    folder.mkdir(parents=True, exist_ok=True)
    log = folder / "runner.log"
    started = time.perf_counter()
    result = subprocess.run(
        [str(RUNNER), mode, str(FMU), str(STOP), str(DT)],
        cwd=folder, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=True,
    )
    wall = time.perf_counter() - started
    log.write_text(result.stdout)
    match = re.search(r"runtime_s=([0-9.eE+-]+)", result.stdout)
    if not match:
        raise RuntimeError(f"Could not parse simulation time from {log}")
    timing = {
        "mode": mode,
        "simulation_runtime_s": float(match.group(1)),
        "process_wall_time_s": wall,
        "communication_step_s": DT,
        "stop_time_s": STOP,
        "network_steps": int(round(STOP / DT)),
        "log": str(log),
    }
    (folder / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
    return timing


def run_all() -> None:
    subprocess.run(["bash", str(CASE / "sources" / "build_fmu.sh")], cwd=ROOT, check=True)
    subprocess.run(["cmake", "--build", str(BUILD), "--target", "FMI_IEEE39_AvgInverter", "-j4"], cwd=ROOT, check=True)
    timings = {mode: run_case(mode) for mode in ("me", "cs")}
    native_out = [DATA / f"run-{mode}" / "native_replay.csv" for mode in ("me", "cs")]
    inputs = [DATA / f"run-{mode}" / "logs" / f"IEEE39_Inverter_{mode}" / f"IEEE39_Inverter_{mode}.csv" for mode in ("me", "cs")]
    matlab_dir = CASE / "matlab"
    matlab_expr = (
        f"addpath('{matlab_dir}'); "
        f"run_native_replay('{inputs[0]}','{native_out[0]}'); "
        f"run_native_replay('{inputs[1]}','{native_out[1]}');"
    )
    timing_paths = [DATA / f"run-{mode}" / "native_timing.json" for mode in ("me", "cs")]
    for path in timing_paths:
        path.touch()
    start_ns = time.time_ns()
    start = time.perf_counter()
    process = subprocess.Popen(["matlab", "-batch", matlab_expr], cwd=ROOT,
                               text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    deadline = time.monotonic() + 600
    while time.monotonic() < deadline:
        status = process.poll()
        if status is not None:
            output = process.communicate()[0]
            if status != 0:
                raise subprocess.CalledProcessError(status, process.args, output=output)
            break
        if all(path.exists() and path.stat().st_mtime_ns > start_ns for path in timing_paths):
            # MATLAB writes both manifests after Simulink has completed. Its
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
        raise TimeoutError("MATLAB did not finish both native Simulink replays within 10 minutes")
    wall = time.perf_counter() - start
    (DATA / "native_batch.log").write_text(output)
    timings["native_replay_batch_wall_time_s"] = wall
    timings["native_replay"] = {
        mode: json.loads((DATA / f"run-{mode}" / "native_timing.json").read_text())
        for mode in ("me", "cs")
    }
    (DATA / "timings.json").write_text(json.dumps(timings, indent=2) + "\n")


def analyze() -> dict:
    me, cs = load_dpsim("me"), load_dpsim("cs")
    native_me, native_cs = load_native("me"), load_native("cs")
    time_me = me.time.to_numpy()
    time_cs = cs.time.to_numpy()
    if not np.allclose(time_me, time_cs, atol=1e-12, rtol=0):
        raise ValueError("ME and CS outputs are not on the same time grid")
    t = time_me
    masks = {
        "post_startup": t >= 0.02,
        "pre_event": (t >= 0.02) & (t < LOAD_STEP),
        "load_step_window": (t >= LOAD_STEP) & (t <= 0.19),
        "recovery": (t > 0.19) & (t <= STOP),
    }
    current_me = vector(me, "inverter.i_")
    current_cs = vector(cs, "inverter.i_")
    native_me_i = np.column_stack([
        np.interp(t, native_me.time, native_me[f"i_{phase}"], left=np.nan, right=np.nan)
        for phase in "abc"
    ])
    native_cs_i = np.column_stack([
        np.interp(t, native_cs.time, native_cs[f"i_{phase}"], left=np.nan, right=np.nan)
        for phase in "abc"
    ])
    current_diffs = {
        "me_vs_native_replay": current_me - native_me_i,
        "cs_vs_native_replay": current_cs - native_cs_i,
        "cs_vs_me": current_cs - current_me,
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
            "disturbance": f"100 MW balanced load step at bus 39 at {LOAD_STEP:.3f} s",
            "device": "average PQ-controlled current source at bus 30; half of original generator-30 P/Q is retained as synchronous generation",
            "device_setpoint": {"P_W": P_SET, "Q_var": Q_SET, "current_response_time_s": 0.005, "current_limit_pu": 1.2},
            "current_base_peak_A": CURRENT_BASE_PEAK,
            "native_reference_scope": "open-loop Simulink replay using each DPsim mode's own inverter-terminal voltage trace",
        },
        "timings": json.loads((DATA / "timings.json").read_text()) if (DATA / "timings.json").exists() else {},
        "current_error": {},
        "bus_voltage_me_vs_cs": {},
        "coupling_diagnostics": {},
    }
    for name, err in current_diffs.items():
        metrics_out["current_error"][name] = {
            window: metrics(err, current_me if name != "cs_vs_native_replay" else native_cs_i,
                            mask, CURRENT_BASE_PEAK)
            for window, mask in masks.items() if mask.any()
        }
    for bus in (1, 10, 20, 30, 31, 39):
        name = f"N{bus}.V_"
        v_me = vector(me, name)
        v_cs = vector(cs, name)
        err = v_cs - v_me
        metrics_out["bus_voltage_me_vs_cs"][str(bus)] = {
            window: metrics(err, v_me, mask, VOLTAGE_BASE_PEAK)
            for window, mask in masks.items() if mask.any()
        }
    metrics_out["coupling_diagnostics"] = {
        "me_max_newton_iterations": int(me.newton_iterations.max()),
        "me_mean_newton_iterations": float(me.newton_iterations[t >= 0.02].mean()),
        "me_max_reported_residual": float(me.newton_residual.max()),
        "me_max_internal_steps_per_macro_step": int(me.fmi_internal_steps.max()),
        "cs_do_step_calls": int(cs.fmi_do_steps.max()),
        "interpretation": "ME uses SinglePassExplicit as the network coupling policy; its residual is the one-pass mismatch, not a converged Newton residual. CS performs one fmi3DoStep per 50 us network step.",
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
            flat.append({"comparison": f"CS_vs_ME_bus{bus}_voltage", "window": window, **values})
    pd.DataFrame(flat).to_csv(csv_path, index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    ax[0,0].plot(t, me["N30.V_0"]/1e3, label="ME")
    ax[0,0].plot(t, cs["N30.V_0"]/1e3, "--", label="CS")
    ax[0,0].set(title="Bus 30 phase-A voltage", xlabel="Time [s]", ylabel="Voltage [kV]")
    ax[0,0].legend()
    ax[0,1].plot(t, me["N39.V_0"]/1e3, label="ME")
    ax[0,1].plot(t, cs["N39.V_0"]/1e3, "--", label="CS")
    ax[0,1].set(title="Bus 39 phase-A voltage", xlabel="Time [s]", ylabel="Voltage [kV]")
    ax[0,1].legend()
    ax[1,0].plot(t, current_me[:,0], label="ME")
    ax[1,0].plot(t, current_cs[:,0], "--", label="CS")
    ax[1,0].plot(native_me.time, native_me.i_a, ":", label="Simulink replay (ME trace)")
    ax[1,0].set(title="Converter phase-A current", xlabel="Time [s]", ylabel="Current [A]")
    ax[1,0].legend(fontsize=8)
    ax[1,1].plot(t, current_diffs["me_vs_native_replay"][:,0], label="ME minus Simulink replay")
    ax[1,1].plot(t, current_diffs["cs_vs_native_replay"][:,0], label="CS minus Simulink replay")
    ax[1,1].set(title="Phase-A current error", xlabel="Time [s]", ylabel="Error [A]")
    ax[1,1].legend(fontsize=8)
    for panel in ax.flat:
        panel.axvline(LOAD_STEP, color="black", alpha=0.5, linestyle=":")
    fig.suptitle("IEEE 39 bus averaged inverter: ME, CS, native Simulink", fontsize=14)
    fig.savefig(DATA / "ieee39_waveforms.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for bus in (30, 39):
        name = f"N{bus}.V_"
        difference = vector(cs, name) - vector(me, name)
        ax[0,0].plot(t, difference[:,0], label=f"Bus {bus} phase A")
    ax[0,0].set(title="CS minus ME bus voltage", xlabel="Time [s]", ylabel="Error [V]")
    ax[0,0].legend()
    ax[0,1].plot(t, np.linalg.norm(current_diffs["cs_vs_me"], axis=1))
    ax[0,1].set(title="CS minus ME converter current", xlabel="Time [s]", ylabel="Vector error [A]")
    post_startup = t >= 0.02
    ax[1,0].plot(t[post_startup], np.linalg.norm(current_diffs["me_vs_native_replay"][post_startup], axis=1), label="ME")
    ax[1,0].plot(t[post_startup], np.linalg.norm(current_diffs["cs_vs_native_replay"][post_startup], axis=1), label="CS")
    ax[1,0].set(title="Post-startup current error vs same-trace native Simulink", xlabel="Time [s]", ylabel="Vector error [A]")
    ax[1,0].legend()
    zoom = (t >= 0.165) & (t <= 0.19)
    ax[1,1].plot(t[zoom], current_me[zoom,0], label="ME")
    ax[1,1].plot(t[zoom], current_cs[zoom,0], label="CS")
    ax[1,1].plot(t[zoom], native_cs_i[zoom,0], ":", label="Simulink replay (CS trace)")
    ax[1,1].set(title="Load-step detail", xlabel="Time [s]", ylabel="Phase-A current [A]")
    ax[1,1].legend(fontsize=8)
    for panel in ax.flat:
        panel.axvline(LOAD_STEP, color="black", alpha=0.5, linestyle=":")
    fig.suptitle("IEEE 39 bus comparison errors", fontsize=14)
    fig.savefig(DATA / "ieee39_errors.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"metrics": str(json_path), "csv": str(csv_path),
                      "waveforms": str(DATA / "ieee39_waveforms.png"),
                      "errors": str(DATA / "ieee39_errors.png")}, indent=2))
    return metrics_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="rerun DPsim ME/CS and native Simulink references")
    args = parser.parse_args()
    if args.run:
        run_all()
    analyze()


if __name__ == "__main__":
    main()
