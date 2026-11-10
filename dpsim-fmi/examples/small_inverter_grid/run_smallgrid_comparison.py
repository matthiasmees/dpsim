#!/usr/bin/env python3
"""Run the focused source--PI-line--inverter comparison in four domains."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
DATA = BUILD / "small-inverter-grid" / "results"
FMU = CASE / "sources" / "IEEE39AvgInverter.fmu"
RUNNER = BUILD / "dpsim-fmi" / "FMI_SmallGrid_Inverter"
MATLAB_DIR = CASE / "matlab"
DT_SMOKE, STOP_SMOKE, EVENT_SMOKE = 25e-6, 0.002, 0.001


def cases(
    stop_time: float, dt: float, load_step: float, sweep_rest: bool = False
) -> list[dict]:
    if sweep_rest:
        return [
            {
                "label": f"{step * 1e6:g}us, {stop_time:g}s load step",
                "sweep_dir": f"dt_{int(round(step * 1e6))}us",
                "step_s": step,
                "stop_time_s": stop_time,
                "load_step_s": load_step,
            }
            for step in (50e-6, 100e-6, 500e-6, 1000e-6)
        ]
    result = [
        {
            "label": "25us startup smoke",
            "sweep_dir": "smoke_25us",
            "step_s": DT_SMOKE,
            "stop_time_s": STOP_SMOKE,
            "load_step_s": EVENT_SMOKE,
        }
    ]
    if stop_time > STOP_SMOKE:
        result.append(
            {
                "label": f"{dt * 1e6:g}us, {stop_time:g}s load step",
                "sweep_dir": "main",
                "step_s": dt,
                "stop_time_s": stop_time,
                "load_step_s": load_step,
            }
        )
    return result


def run_dpsim(mode: str, entry: dict) -> dict:
    out = DATA / "dpsim" / mode / entry["sweep_dir"]
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    proc = subprocess.run(
        [
            str(RUNNER),
            mode,
            str(FMU),
            str(entry["stop_time_s"]),
            str(entry["step_s"]),
            str(entry["load_step_s"]),
        ],
        cwd=out,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    process_wall = time.perf_counter() - started
    (out / "runner.log").write_text(proc.stdout)
    line = next(
        (x for x in proc.stdout.splitlines() if x.startswith("small_grid ")), None
    )
    if line is None:
        raise RuntimeError(f"DPsim timing summary missing from {out / 'runner.log'}")
    values = dict(part.split("=", 1) for part in line.split() if "=" in part)
    log_name = f"SmallGrid_Inverter_{mode}"
    trace_path = out / "logs" / log_name / f"{log_name}.csv"
    frame = pd.read_csv(trace_path, skipinitialspace=True)
    frame.columns = frame.columns.str.strip()
    endpoint_tol = max(1e-8, entry["step_s"] * 1e-5)
    frame = frame[frame.time <= entry["stop_time_s"] + endpoint_tol].reset_index(
        drop=True
    )
    if len(frame) != round(entry["stop_time_s"] / entry["step_s"]) + 1:
        raise ValueError(f"Unexpected DPsim samples in {trace_path}: {len(frame)}")
    if not np.isfinite(frame.select_dtypes(include=[np.number]).to_numpy()).all():
        raise ValueError(f"Non-finite DPsim output in {trace_path}")
    trace = pd.DataFrame({"time_s": frame.time})
    for phase, suffix in enumerate("abc"):
        trace[f"v_{suffix}_V"] = frame[f"pcc.v_{phase}"]
        trace[f"i_{suffix}_A"] = frame[f"inverter.i_{phase}"]
    trace.to_csv(out / f"dpsim_{mode}.csv", index=False)
    timing = {
        "mode": mode,
        "label": entry["label"],
        "step_s": entry["step_s"],
        "stop_time_s": entry["stop_time_s"],
        "load_step_s": entry["load_step_s"],
        "samples": len(frame),
        "process_wall_time_s": process_wall,
        **{key: float(val) for key, val in values.items() if key.endswith("_s")},
    }
    timing["trace"] = str(out / f"dpsim_{mode}.csv")
    (out / f"dpsim_{mode}_timing.json").write_text(json.dumps(timing, indent=2) + "\n")
    return timing


def run_matlab(entry_list: list[dict], mode: str, run_tag: str = "") -> None:
    out = DATA / "simulink" / mode
    out.mkdir(parents=True, exist_ok=True)
    tag_suffix = f"_{run_tag}" if run_tag else ""
    manifest = DATA / f"manifest_{mode}{tag_suffix}.json"
    manifest.write_text(json.dumps(entry_list, indent=2) + "\n")
    model = BUILD / "small-inverter-grid" / "models" / f"SmallGrid_{mode}.slx"
    model.parent.mkdir(parents=True, exist_ok=True)
    # The MATLAB runner writes ``<mode>_complete.marker`` in its output
    # directory. Keep the marker name independent of the manifest's run tag;
    # the tag only isolates the manifest/log files and the case directories.
    complete = out / f"{mode}_complete.marker"
    complete.unlink(missing_ok=True)
    for entry in entry_list:
        case_dir = out / entry["sweep_dir"]
        for suffix in (
            f"simulink_{mode}.mat",
            f"simulink_{mode}.csv",
            f"simulink_{mode}_timing.json",
        ):
            (case_dir / suffix).unlink(missing_ok=True)
    expression = (
        f"addpath('{MATLAB_DIR}'); run_smallgrid_sweep('{mode}',"
        f"'{manifest}','{out}','{FMU}','unused','{model}');"
    )
    proc = subprocess.Popen(
        [os.environ.get("MATLAB_BIN", "matlab"), "-batch", expression],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    deadline = time.monotonic() + 3600
    while time.monotonic() < deadline:
        if complete.is_file():
            proc.terminate()
            break
        if proc.poll() is not None:
            break
        time.sleep(1)
    if proc.poll() is None:
        proc.terminate()
    try:
        output = proc.communicate(timeout=8)[0]
    except subprocess.TimeoutExpired:
        proc.kill()
        output = proc.communicate()[0]
    (DATA / f"matlab_{mode}{tag_suffix}.log").write_text(output)
    expected = [
        out / entry["sweep_dir"] / f"simulink_{mode}_timing.json"
        for entry in entry_list
    ]
    missing = [str(path) for path in expected if not path.is_file()]
    if missing or not complete.is_file():
        raise RuntimeError(
            f"MATLAB {mode} run did not complete; missing {missing}\n{output[-4000:]}"
        )


def validate_and_summarize(
    entry_list: list[dict], output_name: str = "metrics.json"
) -> dict:
    summary = {}
    for entry in entry_list:
        key = entry["sweep_dir"]
        traces = {}
        for mode in ("me", "cs"):
            path = DATA / "dpsim" / mode / key / f"dpsim_{mode}.csv"
            traces[f"DPsim {mode.upper()}"] = pd.read_csv(path)
        for mode in ("native", "cs"):
            path = DATA / "simulink" / mode / key / f"simulink_{mode}.csv"
            traces[f"Simulink {mode.upper()}"] = pd.read_csv(path)
        reference = traces["Simulink NATIVE"]
        t = reference.time_s.to_numpy()
        metrics = {}
        signal_health = {}
        dt = entry["step_s"]
        for name, trace in traces.items():
            if len(trace) != len(reference) or not np.allclose(
                trace.time_s, t, rtol=0, atol=max(1e-12, dt * 1e-8)
            ):
                raise ValueError(f"Time grid mismatch for {key}: {name}")
            numeric = trace.select_dtypes(include=[np.number]).to_numpy()
            if not np.isfinite(numeric).all():
                raise ValueError(f"Non-finite {name} output in {key}")
            signal_health[name] = {
                "finite_all_signals": True,
                "max_abs_pcc_voltage_V": float(
                    trace[[f"v_{phase}_V" for phase in "abc"]].abs().to_numpy().max()
                ),
                "max_abs_inverter_current_A": float(
                    trace[[f"i_{phase}_A" for phase in "abc"]].abs().to_numpy().max()
                ),
            }
            if name == "Simulink NATIVE":
                continue
            fields = {}
            for phase in "abc":
                for prefix in ("v", "i"):
                    signal = f"{prefix}_{phase}_{'V' if prefix == 'v' else 'A'}"
                    scale = np.sqrt(np.mean(reference[signal].to_numpy() ** 2))
                    err = trace[signal].to_numpy() - reference[signal].to_numpy()
                    fields[signal] = {
                        "rms_error": float(np.sqrt(np.mean(err**2))),
                        "rms_error_percent_reference_rms": float(
                            100 * np.sqrt(np.mean(err**2)) / scale
                        ),
                        "max_abs_error": float(np.max(np.abs(err))),
                    }
            metrics[name] = fields
        summary[key] = {
            "case": entry,
            "signal_health": signal_health,
            "metrics_vs_simulink_native": metrics,
            "traces": {
                k: str(
                    DATA
                    / ("dpsim" if k.startswith("DPsim") else "simulink")
                    / (
                        "me"
                        if k == "DPsim ME"
                        else "cs" if k in ("DPsim CS", "Simulink CS") else "native"
                    )
                    / key
                    / (
                        "dpsim_" + ("me" if k == "DPsim ME" else "cs") + ".csv"
                        if k.startswith("DPsim")
                        else "simulink_"
                        + ("native" if k == "Simulink NATIVE" else "cs")
                        + ".csv"
                    )
                )
                for k in traces
            },
        }
    out = DATA / output_name
    out.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stop-time", type=float, default=5.0)
    parser.add_argument("--dt", type=float, default=25e-6)
    parser.add_argument("--load-step", type=float, default=2.5)
    parser.add_argument("--skip-matlab", action="store_true")
    parser.add_argument("--skip-dpsim", action="store_true")
    parser.add_argument(
        "--sweep-rest",
        action="store_true",
        help="run the remaining 5 s timestep cases: 50 to 1000 microseconds",
    )
    args = parser.parse_args()
    if args.sweep_rest:
        args.stop_time = 5.0
        args.load_step = 2.5
    if args.dt <= 0 or args.stop_time <= 0 or not 0 < args.load_step < args.stop_time:
        parser.error(
            "dt and stop time must be positive and the load step must be in range"
        )
    DATA.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["cmake", "--build", str(BUILD), "--target", "FMI_SmallGrid_Inverter", "-j2"],
        check=True,
    )
    entry_list = cases(args.stop_time, args.dt, args.load_step, args.sweep_rest)
    if not args.skip_dpsim:
        for mode in ("me", "cs"):
            for entry in entry_list:
                print(f"Running DPsim {mode.upper()} {entry['label']}", flush=True)
                run_dpsim(mode, entry)
    if not args.skip_matlab:
        for mode in ("native", "cs"):
            print(f"Running Simulink {mode.upper()} cases", flush=True)
            run_matlab(entry_list, mode, "sweep" if args.sweep_rest else "")
    metrics = validate_and_summarize(
        entry_list, "metrics_sweep.json" if args.sweep_rest else "metrics.json"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
