#!/usr/bin/env python3
"""Rerun DPsim ME and CS at every retained step with periodic PF starts."""
from __future__ import annotations

import csv
import argparse
import json
import math
import subprocess
import time
from pathlib import Path

from run_powerflow_initialization_comparison import rebuild_fmu

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
FMU = CASE / "sources" / "IEEE39AvgInverter.fmu"
EXECUTABLE = BUILD / "dpsim-fmi" / "FMI_SmallGrid_Inverter"
OUT = BUILD / "small-inverter-grid" / "results" / "post_powerflow_fix"
STEP_US = (25, 50, 100, 500, 1000)
STOP = 5.0
LOAD_STEP = 2.5
BASE_DT = 25e-6
VISUAL_POINTS = 8000


def case_dir(stage: str, mode: str, step_us: int) -> Path:
    return OUT / stage / mode / f"dt_{step_us}us"


def log_name(mode: str) -> str:
    return f"SmallGrid_Inverter_{mode}_powerflow"


def timing_profiles(raw_trace: Path, mode: str) -> dict:
    directory = raw_trace.parent
    stem = log_name(mode)
    task_path = directory / f"{stem}_TaskTimes.csv"
    phase_path = directory / f"{stem}_SolverPhaseTimes.csv"
    tasks = []
    with task_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if int(row["calls"]) == 0:
                continue
            tasks.append(
                {
                    "task": row["task"],
                    "avg_us": float(row["avg_ns"]) / 1000,
                    "total_s": float(row["total_ns"]) / 1e9,
                    "calls": int(row["calls"]),
                }
            )
    phases = []
    with phase_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if int(row["calls"]) == 0:
                continue
            phases.append(
                {
                    "phase": row["phase"],
                    "avg_us": float(row["avg_us"]),
                    "total_s": float(row["total_s"]),
                    "calls": int(row["calls"]),
                }
            )
    return {"component_tasks": tasks, "solver_phases": phases}


def run_one(stage: str, mode: str, step_us: int, duration: float) -> dict:
    dt = step_us * 1e-6
    load_step = LOAD_STEP if duration > LOAD_STEP else duration / 2
    directory = case_dir(stage, mode, step_us)
    directory.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = subprocess.run(
        [
            str(EXECUTABLE),
            mode,
            str(FMU),
            str(duration),
            str(dt),
            str(load_step),
            "powerflow",
        ],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    process_wall = time.perf_counter() - started
    (directory / "runner.log").write_text(completed.stdout)
    summary = next(
        (
            line
            for line in completed.stdout.splitlines()
            if line.startswith("small_grid ")
        ),
        None,
    )
    if summary is None:
        raise RuntimeError(f"Missing DPsim timing summary for {mode} {step_us} µs")
    timing = dict(field.split("=", 1) for field in summary.split() if "=" in field)
    raw = directory / "logs" / log_name(mode) / f"{log_name(mode)}.csv"
    if not raw.is_file():
        raise FileNotFoundError(raw)
    expected = round(duration / dt) + 1
    compact_stride = max(1, math.ceil(expected / VISUAL_POINTS))
    compact_path = directory / f"{mode}_powerflow.csv"
    count = 0
    iteration_values: list[int] = []
    maximum_residual = 0.0
    initial_voltage: list[float] | None = None
    initial_current: list[float] | None = None
    v_transient: list[float] = []
    i_transient: list[float] = []
    v_steady: list[float] = []
    i_steady: list[float] = []
    with raw.open(newline="") as source, compact_path.open("w", newline="") as target:
        reader = csv.DictReader(source, skipinitialspace=True)
        reader.fieldnames = [name.strip() for name in (reader.fieldnames or [])]
        fields = set(reader.fieldnames)
        if not {
            "time",
            *(f"pcc.v_{p}" for p in range(3)),
            *(f"inverter.i_{p}" for p in range(3)),
        }.issubset(fields):
            raise RuntimeError(f"Missing voltage/current columns in {raw}")
        is_me = mode == "me"
        writer_fields = [
            "time_s",
            *(f"v_{p}_V" for p in "abc"),
            *(f"i_{p}_A" for p in "abc"),
        ]
        if is_me:
            writer_fields += ["newton_iterations", "newton_residual"]
        writer = csv.DictWriter(target, fieldnames=writer_fields)
        writer.writeheader()
        for index, row in enumerate(reader):
            t = float(row["time"])
            voltages = [float(row[f"pcc.v_{p}"]) for p in range(3)]
            currents = [float(row[f"inverter.i_{p}"]) for p in range(3)]
            if not all(math.isfinite(value) for value in [t, *voltages, *currents]):
                raise RuntimeError(f"Non-finite output in {raw} at row {index}")
            if initial_voltage is None:
                initial_voltage, initial_current = voltages, currents
            vmag = math.sqrt(2.0 / 3.0 * sum(value * value for value in voltages))
            imag = math.sqrt(2.0 / 3.0 * sum(value * value for value in currents))
            if t <= min(duration, 0.1):
                v_transient.append(vmag)
                i_transient.append(imag)
            if 0.5 <= t < min(duration, load_step) - 0.1:
                v_steady.append(vmag)
                i_steady.append(imag)
            output = {
                "time_s": t,
                **{f"v_{p}_V": voltages[k] for k, p in enumerate("abc")},
                **{f"i_{p}_A": currents[k] for k, p in enumerate("abc")},
            }
            if is_me:
                iterations = int(float(row["newton_iterations"]))
                residual = float(row["newton_residual"])
                if t > 0:
                    iteration_values.append(iterations)
                maximum_residual = max(maximum_residual, residual)
                output["newton_iterations"] = iterations
                output["newton_residual"] = residual
            if index % compact_stride == 0 or index == expected - 1:
                writer.writerow(output)
            count += 1
    if count != expected:
        raise RuntimeError(f"Expected {expected} samples but read {count} from {raw}")
    deviation = {}
    for signal, early, steady in (
        ("voltage", v_transient, v_steady),
        ("current", i_transient, i_steady),
    ):
        reference = sum(steady) / len(steady) if steady else math.nan
        deviation[signal] = (
            100
            * math.sqrt(sum((value - reference) ** 2 for value in early) / len(early))
            / reference
            if early and math.isfinite(reference) and reference
            else None
        )
    result = {
        "mode": mode.upper(),
        "initialization": "powerflow_periodic_state",
        "step_us": step_us,
        "step_s": dt,
        "stop_s": duration,
        "load_step_s": load_step,
        "samples": count,
        "process_wall_s": process_wall,
        "simulation_init_s": float(timing["init_s"]),
        "simulation_wall_s": float(timing["wall_s"]),
        "simulation_run_s": float(timing["simulation_s"]),
        "initial_pcc_abc_V": initial_voltage,
        "initial_current_abc_A": initial_current,
        "initial_100ms_magnitude_deviation_percent": deviation,
        "raw_trace": str(raw),
        "compact_trace": str(compact_path),
        "newton": None,
        "timing_profiles": timing_profiles(raw, mode),
    }
    if iteration_values:
        ordered = sorted(iteration_values)
        result["newton"] = {
            "steps": len(ordered),
            "min": ordered[0],
            "median": ordered[len(ordered) // 2],
            "p95": ordered[min(len(ordered) - 1, math.ceil(0.95 * len(ordered)) - 1)],
            "max": ordered[-1],
            "max_scaled_residual": maximum_residual,
        }
    (directory / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    startup_text = (
        f"{deviation['voltage']:.6g}%/{deviation['current']:.6g}%"
        if deviation["voltage"] is not None
        else "n/a"
    )
    print(
        f"{stage:5} {mode.upper():2} {step_us:4} µs: "
        f"{count:7} samples, {result['simulation_wall_s']:.3f}s, "
        f"startup dev V/I={startup_text}",
        flush=True,
    )
    return result


def full_reference_metrics(records: list[dict]) -> dict:
    by_case = {
        (entry["mode"].lower(), entry["step_us"]): entry
        for entry in records
        if entry["stop_s"] == STOP
    }
    reference_path = Path(by_case[("me", 25)]["raw_trace"])
    report = {}
    for step_us in STEP_US:
        stride = step_us // 25
        method_results = {}
        for mode in ("me", "cs"):
            case = by_case[(mode, step_us)]
            source_path = Path(case["raw_trace"])
            error_sq = {"voltage": 0.0, "current": 0.0}
            reference_sq = {"voltage": 0.0, "current": 0.0}
            max_error = {"voltage": 0.0, "current": 0.0}
            sample_count = 0
            with (
                reference_path.open(newline="") as ref_stream,
                source_path.open(newline="") as src_stream,
            ):
                ref_rows = csv.DictReader(ref_stream, skipinitialspace=True)
                src_rows = csv.DictReader(src_stream, skipinitialspace=True)
                ref_rows.fieldnames = [
                    name.strip() for name in (ref_rows.fieldnames or [])
                ]
                src_rows.fieldnames = [
                    name.strip() for name in (src_rows.fieldnames or [])
                ]
                source_index = -1
                for ref_index, ref in enumerate(ref_rows):
                    if ref_index % stride != 0:
                        continue
                    try:
                        value = next(src_rows)
                    except StopIteration as error:
                        raise RuntimeError(
                            f"Trace shorter than reference: {source_path}"
                        ) from error
                    source_index += 1
                    sample_count += 1
                    if abs(float(ref["time"]) - float(value["time"])) > 1e-12:
                        raise RuntimeError(
                            f"Time grids differ at {step_us} µs for {mode}"
                        )
                    for kind, suffix in (("voltage", "V"), ("current", "A")):
                        for phase in "abc":
                            phase_index = "abc".index(phase)
                            column = (
                                f"pcc.v_{phase_index}"
                                if kind == "voltage"
                                else f"inverter.i_{phase_index}"
                            )
                            reference_value = float(ref[column])
                            difference = float(value[column]) - reference_value
                            error_sq[kind] += difference * difference
                            reference_sq[kind] += reference_value * reference_value
                            max_error[kind] = max(max_error[kind], abs(difference))
            metrics = {}
            for kind in ("voltage", "current"):
                denominator = max(1, sample_count * 3)
                rms_error = math.sqrt(error_sq[kind] / denominator)
                rms_ref = math.sqrt(reference_sq[kind] / denominator)
                metrics[kind] = {
                    "rms_error": rms_error,
                    "rms_error_percent_25us_me_rms": (
                        100 * rms_error / rms_ref if rms_ref else None
                    ),
                    "max_abs_error": max_error[kind],
                }
            method_results[mode.upper()] = metrics
        pair = {}
        me_path = Path(by_case[("me", step_us)]["raw_trace"])
        cs_path = Path(by_case[("cs", step_us)]["raw_trace"])
        with me_path.open(newline="") as a, cs_path.open(newline="") as b:
            for kind, suffix in (("voltage", "V"), ("current", "A")):
                err_sq = ref_sq = 0.0
                sample_count = 0
                a.seek(0)
                b.seek(0)
                ra = csv.DictReader(a, skipinitialspace=True)
                rb = csv.DictReader(b, skipinitialspace=True)
                ra.fieldnames = [name.strip() for name in (ra.fieldnames or [])]
                rb.fieldnames = [name.strip() for name in (rb.fieldnames or [])]
                for row_a, row_b in zip(ra, rb):
                    sample_count += 1
                    for phase in "abc":
                        phase_index = "abc".index(phase)
                        col = (
                            f"pcc.v_{phase_index}"
                            if kind == "voltage"
                            else f"inverter.i_{phase_index}"
                        )
                        difference = float(row_b[col]) - float(row_a[col])
                        err_sq += difference * difference
                        ref_sq += float(row_a[col]) ** 2
                denominator = max(1, sample_count * 3)
                rms_error = math.sqrt(err_sq / denominator)
                rms_ref = math.sqrt(ref_sq / denominator)
                pair[kind] = {
                    "rms_error": rms_error,
                    "rms_error_percent_me_rms": (
                        100 * rms_error / rms_ref if rms_ref else None
                    ),
                }
        report[str(step_us)] = {
            "ME_vs_25us_ME": method_results["ME"],
            "CS_vs_25us_ME": method_results["CS"],
            "CS_vs_ME_same_step": pair,
        }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="reuse the fresh sweep traces already on disk",
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    if args.analyze_only:
        for stage, duration in (("smoke", 0.05), ("full", STOP)):
            for step_us in STEP_US:
                for mode in ("me", "cs"):
                    path = case_dir(stage, mode, step_us) / "metrics.json"
                    item = json.loads(path.read_text())
                    item["timing_profiles"] = timing_profiles(
                        Path(item["raw_trace"]), mode
                    )
                    if mode == "me":
                        iterations: list[int] = []
                        max_residual = 0.0
                        with Path(item["raw_trace"]).open(newline="") as stream:
                            rows = csv.DictReader(stream, skipinitialspace=True)
                            rows.fieldnames = [
                                name.strip() for name in (rows.fieldnames or [])
                            ]
                            for row in rows:
                                if float(row["time"]) <= 0:
                                    continue
                                iterations.append(int(float(row["newton_iterations"])))
                                max_residual = max(
                                    max_residual, float(row["newton_residual"])
                                )
                        ordered = sorted(iterations)
                        item["newton"] = {
                            "steps": len(ordered),
                            "min": ordered[0],
                            "median": ordered[len(ordered) // 2],
                            "p95": ordered[
                                min(
                                    len(ordered) - 1, math.ceil(0.95 * len(ordered)) - 1
                                )
                            ],
                            "max": ordered[-1],
                            "max_scaled_residual": max_residual,
                        }
                    records.append(item)
    else:
        subprocess.run(
            [
                "cmake",
                "--build",
                str(BUILD),
                "--target",
                "FMI_SmallGrid_Inverter",
                "-j1",
            ],
            check=True,
        )
        rebuild_fmu()
        # Fresh short checks precede the 5 s sweep for every mode and step.
        for step_us in STEP_US:
            for mode in ("me", "cs"):
                records.append(run_one("smoke", mode, step_us, 0.05))
        for step_us in STEP_US:
            for mode in ("me", "cs"):
                records.append(run_one("full", mode, step_us, STOP))
    report = {
        "configuration": {
            "stop_s": STOP,
            "load_step_s": LOAD_STEP,
            "steps_us": list(STEP_US),
            "modes": ["ME", "CS"],
            "initialization": "powerflow with periodic current state",
            "samples": {str(us): round(STOP / (us * 1e-6)) + 1 for us in STEP_US},
        },
        "runs": records,
        "metrics_vs_fine_step": full_reference_metrics(records),
    }
    report_path = OUT / "metrics.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
