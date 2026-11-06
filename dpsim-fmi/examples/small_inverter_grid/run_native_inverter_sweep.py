#!/usr/bin/env python3
"""Run or summarize the native DPsim inverter on the matched small-grid sweep."""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[2]
BUILD = ROOT / "build-fmi"
DATA = BUILD / "small-inverter-grid" / "results"
METHOD_RUNNER = CASE / "run_smallgrid_comparison.py"
STEPS_US = (25, 50, 100, 500, 1000)


def comparison_cases() -> list[dict]:
    return [
        {
            "label": f"{dt}us, 5s load step",
            "sweep_dir": "main" if dt == 25 else f"dt_{dt}us",
            "native_dir": f"dt_{dt}us",
            "step_s": dt * 1e-6,
            "stop_time_s": 5.0,
            "load_step_s": 2.5,
        }
        for dt in STEPS_US
    ]


def load_runner_module():
    spec = importlib.util.spec_from_file_location("smallgrid_runner", METHOD_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {METHOD_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def method_traces(key: str) -> tuple[dict[str, pd.DataFrame], dict]:
    if key == "main":
        item = json.loads((DATA / "metrics.json").read_text())[key]
    else:
        item = json.loads((DATA / "metrics_sweep.json").read_text())[key]
    traces = {name: pd.read_csv(path) for name, path in item["traces"].items()}
    return traces, item


def run_native(entry: dict) -> None:
    runner = load_runner_module()
    native_entry = dict(entry)
    native_entry["sweep_dir"] = entry["native_dir"]
    print(f"Running native DPsim {entry['label']}", flush=True)
    result = runner.run_dpsim("native", native_entry)
    print(json.dumps(result), flush=True)


def summarize() -> dict:
    report = {}
    for entry in comparison_cases():
        key = entry["sweep_dir"]
        traces, existing_metrics = method_traces(key)
        native_path = (
            DATA / "dpsim" / "native" / entry["native_dir"] / "dpsim_native.csv"
        )
        native = pd.read_csv(native_path)
        reference = traces["Simulink NATIVE"]
        expected_time = reference["time_s"].to_numpy()
        dt = entry["step_s"]
        if len(native) != len(reference) or not np.allclose(
            native.time_s.to_numpy(), expected_time, rtol=0, atol=max(1e-12, dt * 1e-8)
        ):
            raise ValueError(f"Native DPsim time grid does not match {key}")
        traces["DPsim native"] = native
        signals = {}
        for method, frame in traces.items():
            if not np.isfinite(
                frame.select_dtypes(include=[np.number]).to_numpy()
            ).all():
                raise ValueError(f"Non-finite values in {method} at {key}")
            if method == "Simulink NATIVE":
                continue
            metrics = {}
            for prefix, unit in (("v", "V"), ("i", "A")):
                columns = [f"{prefix}_{phase}_{unit}" for phase in "abc"]
                delta = frame[columns].to_numpy() - reference[columns].to_numpy()
                ref_values = reference[columns].to_numpy()
                rms_error = float(np.sqrt(np.mean(delta**2)))
                ref_rms = float(np.sqrt(np.mean(ref_values**2)))
                metrics["PCC voltage" if prefix == "v" else "Inverter current"] = {
                    f"rms_error_{unit}": rms_error,
                    "rms_error_percent_reference_rms": 100 * rms_error / ref_rms,
                    f"max_abs_error_{unit}": float(np.max(np.abs(delta))),
                }
            signals[method] = metrics

        raw_path = (
            DATA
            / "dpsim"
            / "native"
            / entry["native_dir"]
            / "logs"
            / "SmallGrid_Inverter_native"
            / "SmallGrid_Inverter_native.csv"
        )
        raw = pd.read_csv(raw_path, skipinitialspace=True)
        steps = raw[raw.time > 0]
        iterations = steps.newton_iterations.astype(int)
        timing = json.loads(
            (
                DATA
                / "dpsim"
                / "native"
                / entry["native_dir"]
                / "dpsim_native_timing.json"
            ).read_text()
        )
        native_relative_to_me = {}
        for prefix, unit in (("v", "V"), ("i", "A")):
            cols = [f"{prefix}_{phase}_{unit}" for phase in "abc"]
            delta = native[cols].to_numpy() - traces["DPsim ME"][cols].to_numpy()
            ref = traces["DPsim ME"][cols].to_numpy()
            err = float(np.sqrt(np.mean(delta**2)))
            native_relative_to_me[prefix] = {
                f"rms_error_{unit}": err,
                "rms_error_percent_dpsim_me_rms": 100
                * err
                / float(np.sqrt(np.mean(ref**2))),
            }
        report[key] = {
            "case": {
                k: entry[k] for k in ("label", "step_s", "stop_time_s", "load_step_s")
            },
            "native_trace": str(native_path),
            "metrics_vs_simulink_native": signals,
            "native_vs_dpsim_me": native_relative_to_me,
            "newton": {
                "accepted_steps": int(len(steps)),
                "min_iterations": int(iterations.min()),
                "median_iterations": float(iterations.median()),
                "p95_iterations": float(iterations.quantile(0.95)),
                "max_iterations": int(iterations.max()),
                "iteration_histogram": {
                    str(int(k)): int(v)
                    for k, v in iterations.value_counts().sort_index().items()
                },
                "max_scaled_residual": float(steps.newton_residual.max()),
            },
            "timing": timing,
            "other_methods": existing_metrics.get("method_status", {}),
        }
    output = DATA / "metrics_native_dpsim.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="reuse native DPsim runs already present on disk",
    )
    args = parser.parse_args()
    if not args.analyze_only:
        subprocess.run(
            [
                "cmake",
                "--build",
                str(BUILD),
                "--target",
                "FMI_SmallGrid_Inverter",
                "-j2",
            ],
            check=True,
        )
        for entry in comparison_cases():
            run_native(entry)
    print(json.dumps(summarize(), indent=2))


if __name__ == "__main__":
    main()
