#!/usr/bin/env python3
"""Run native, FMI ME and FMI CS small-grid cases for the 10 s timestep sweep.

Raw solver logs are streamed and removed after their metrics are computed. The
script retains compact traces and exact full-resolution RMS error summaries
under build-fmi/small-inverter-grid/results/dpsim_10s_1to1000us_native_pf_init/.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import time
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
EXE = BUILD / "dpsim-fmi" / "FMI_SmallGrid_Inverter"
FMU = CASE / "sources" / "IEEE39AvgInverter.fmu"
OUT = BUILD / "small-inverter-grid" / "results" / "dpsim_10s_1to1000us_native_pf_init"
STEPS_US = (10, 25, 50, 100, 500, 1000)
MODES = ("native", "me", "cs")
STOP_S = 10.0
LOAD_STEP_S = 5.0
ERROR_WINDOW_S = 0.1
MAX_PLOT_POINTS = 10000


def stem(mode: str) -> str:
    return f"SmallGrid_Inverter_{mode}_powerflow"


def raw_path(directory: Path, mode: str) -> Path:
    name = stem(mode)
    return directory / "logs" / name / f"{name}.csv"


def run_case(mode: str, step_us: int) -> dict:
    dt = step_us * 1e-6
    directory = OUT / f"dt_{step_us}us" / mode
    directory.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = subprocess.run(
        [str(EXE), mode, str(FMU), str(STOP_S), str(dt), str(LOAD_STEP_S), "powerflow"],
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
        raise RuntimeError(f"Missing solver summary for {mode} at {step_us} us")
    timing = dict(item.split("=", 1) for item in summary.split() if "=" in item)
    raw = raw_path(directory, mode)
    if not raw.is_file():
        raise FileNotFoundError(raw)
    return {
        "mode": mode,
        "step_us": step_us,
        "step_s": dt,
        "duration_s": STOP_S,
        "load_step_s": LOAD_STEP_S,
        "initialization": "powerflow (network powerflow plus FMU state reinitialization)",
        "samples_expected": round(STOP_S / dt) + 1,
        "process_wall_s": process_wall,
        "init_s": float(timing["init_s"]),
        "solver_wall_s": float(timing["wall_s"]),
        "simulation_s": float(timing["simulation_s"]),
        "raw_path": str(raw),
        "directory": str(directory),
    }


def names(row: dict) -> list[float]:
    return [float(row[f"pcc.v_{i}"]) for i in range(3)] + [
        float(row[f"inverter.i_{i}"]) for i in range(3)
    ]


def write_compact_native(record: dict) -> tuple[Path, dict]:
    raw = Path(record["raw_path"])
    compact = Path(record["directory"]) / "native_trace.csv"
    expected = record["samples_expected"]
    stride = max(1, math.ceil(expected / MAX_PLOT_POINTS))
    count = 0
    first_v = first_i = None
    with raw.open(newline="") as source, compact.open("w", newline="") as target:
        reader = csv.DictReader(source, skipinitialspace=True)
        writer = csv.writer(target)
        writer.writerow(
            ["time_s", "v_a_V", "v_b_V", "v_c_V", "i_a_A", "i_b_A", "i_c_A"]
        )
        for index, row in enumerate(reader):
            t = float(row["time"])
            values = names(row)
            if not math.isfinite(t) or not all(math.isfinite(x) for x in values):
                raise ValueError(f"Non-finite native output at sample {index}")
            if first_v is None:
                first_v, first_i = values[:3], values[3:]
            if index % stride == 0 or index == expected - 1:
                writer.writerow([t, *values])
            count += 1
    if count != expected:
        raise RuntimeError(f"Native trace has {count} samples, expected {expected}")
    return compact, {
        "samples": count,
        "initial_pcc_abc_V": first_v,
        "initial_current_abc_A": first_i,
    }


def compare_candidate(reference: dict, candidate: dict, mode: str) -> dict:
    ref_path = Path(reference["raw_path"])
    src_path = Path(candidate["raw_path"])
    directory = Path(candidate["directory"])
    compact = directory / f"{mode}_trace.csv"
    error_trace = directory / f"{mode}_rolling_error.csv"
    dt = candidate["step_s"]
    window = max(1, round(ERROR_WINDOW_S / dt))
    stride = max(1, math.ceil(candidate["samples_expected"] / MAX_PLOT_POINTS))
    ring: deque[tuple[float, float]] = deque()
    sum_err_v = sum_ref_v = sum_err_i = sum_ref_i = 0.0
    full = {"v_err_sq": 0.0, "v_ref_sq": 0.0, "i_err_sq": 0.0, "i_ref_sq": 0.0}
    peak = {"v": 0.0, "i": 0.0}
    before = {
        "v_err_sq": 0.0,
        "v_ref_sq": 0.0,
        "i_err_sq": 0.0,
        "i_ref_sq": 0.0,
        "n": 0,
    }
    after = {"v_err_sq": 0.0, "v_ref_sq": 0.0, "i_err_sq": 0.0, "i_ref_sq": 0.0, "n": 0}
    count = 0
    first = None
    with (
        ref_path.open(newline="") as ref_stream,
        src_path.open(newline="") as src_stream,
        compact.open("w", newline="") as compact_stream,
        error_trace.open("w", newline="") as error_stream,
    ):
        rr = csv.DictReader(ref_stream, skipinitialspace=True)
        cr = csv.DictReader(src_stream, skipinitialspace=True)
        compact_writer = csv.writer(compact_stream)
        error_writer = csv.writer(error_stream)
        compact_writer.writerow(
            ["time_s", "v_a_V", "v_b_V", "v_c_V", "i_a_A", "i_b_A", "i_c_A"]
        )
        error_writer.writerow(
            ["time_s", "voltage_error_percent_100ms", "current_error_percent_100ms"]
        )
        for index, (r, c) in enumerate(zip(rr, cr)):
            t = float(c["time"])
            tr = float(r["time"])
            rv, cv = names(r), names(c)
            if abs(t - tr) > max(1e-12, dt * 1e-8):
                raise RuntimeError(
                    f"Time grids differ for {mode} at sample {index}: {t} vs {tr}"
                )
            if not math.isfinite(t) or not all(math.isfinite(x) for x in (*rv, *cv)):
                raise ValueError(f"Non-finite output for {mode} at sample {index}")
            if first is None:
                first = {"pcc_abc_V": cv[:3], "current_abc_A": cv[3:]}
            if index % stride == 0 or index == candidate["samples_expected"] - 1:
                compact_writer.writerow([t, *cv])
            v_err = sum((cv[j] - rv[j]) ** 2 for j in range(3))
            v_ref = sum(rv[j] ** 2 for j in range(3))
            i_err = sum((cv[j + 3] - rv[j + 3]) ** 2 for j in range(3))
            i_ref = sum(rv[j + 3] ** 2 for j in range(3))
            full["v_err_sq"] += v_err
            full["v_ref_sq"] += v_ref
            full["i_err_sq"] += i_err
            full["i_ref_sq"] += i_ref
            peak["v"] = max(peak["v"], math.sqrt(v_err / 3))
            peak["i"] = max(peak["i"], math.sqrt(i_err / 3))
            window_item = (v_err, v_ref, i_err, i_ref)
            ring.append(window_item)
            sum_err_v += v_err
            sum_ref_v += v_ref
            sum_err_i += i_err
            sum_ref_i += i_ref
            if len(ring) > window:
                old = ring.popleft()
                sum_err_v -= old[0]
                sum_ref_v -= old[1]
                sum_err_i -= old[2]
                sum_ref_i -= old[3]
            if index % stride == 0 or index == candidate["samples_expected"] - 1:
                # Repeated add/remove updates can leave an O(machine-epsilon)
                # negative sum after the true squared error reaches zero.
                ev = (
                    100 * math.sqrt(max(0.0, sum_err_v) / sum_ref_v)
                    if sum_ref_v > 0
                    else 0.0
                )
                ei = (
                    100 * math.sqrt(max(0.0, sum_err_i) / sum_ref_i)
                    if sum_ref_i > 0
                    else 0.0
                )
                error_writer.writerow([t, ev, ei])
            section = before if t < LOAD_STEP_S else after
            section["v_err_sq"] += v_err
            section["v_ref_sq"] += v_ref
            section["i_err_sq"] += i_err
            section["i_ref_sq"] += i_ref
            section["n"] += 1
            count += 1
    if count != candidate["samples_expected"]:
        raise RuntimeError(
            f"{mode} trace has {count} samples, expected {candidate['samples_expected']}"
        )

    def score(data: dict, kind: str) -> dict:
        err = data[f"{kind}_err_sq"]
        ref = data[f"{kind}_ref_sq"]
        return {
            "rms_error": math.sqrt(err / (3 * max(1, data.get("n", count)))),
            "rms_error_percent": 100 * math.sqrt(err / ref) if ref else None,
        }

    report = {
        "mode": mode.upper(),
        "samples": count,
        "initial": first,
        "full_run": {"voltage": score(full, "v"), "current": score(full, "i")},
        "before_load_step": {
            "voltage": score(before, "v"),
            "current": score(before, "i"),
        },
        "after_load_step": {"voltage": score(after, "v"), "current": score(after, "i")},
        "peak_instantaneous_phase_rms_error": peak,
        "compact_trace": str(compact),
        "rolling_error_trace": str(error_trace),
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--steps-us", nargs="+", type=int, choices=STEPS_US, default=list(STEPS_US)
    )
    parser.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    parser.add_argument(
        "--resume", action="store_true", help="reuse completed run logs and metrics"
    )
    args = parser.parse_args()
    if "native" not in args.modes:
        raise ValueError("Each timestep needs native DPsim as the error reference")
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "configuration": {
            "duration_s": STOP_S,
            "load_step_s": LOAD_STEP_S,
            "steps_us": list(STEPS_US),
            "methods": list(MODES),
            "initialization": "network powerflow; FMU state reinitialized from powerflow",
            "error_window_s": ERROR_WINDOW_S,
        },
        "cases": {},
    }
    manifest_path = OUT / "metrics.json"
    if manifest_path.is_file():
        manifest["cases"].update(json.loads(manifest_path.read_text()).get("cases", {}))
    for step_us in args.steps_us:
        print(f"\n=== {step_us} us ===", flush=True)
        case: dict = {}
        for mode in args.modes:
            record_path = OUT / f"dt_{step_us}us" / mode / "run.json"
            if args.resume and record_path.is_file():
                record = json.loads(record_path.read_text())
                if not Path(record["raw_path"]).is_file():
                    # Reuse compact run only if already fully analyzed below.
                    if (
                        mode == "native"
                        and Path(record["directory"], "native_trace.csv").is_file()
                    ):
                        case[mode] = record
                        continue
                else:
                    case[mode] = record
                    continue
            record = run_case(mode, step_us)
            record_path.write_text(json.dumps(record, indent=2) + "\n")
            case[mode] = record
            print(
                f"{mode:6} {record['samples_expected']:>9,} samples; solver {record['solver_wall_s']:.3f}s; process {record['process_wall_s']:.3f}s",
                flush=True,
            )
        if "native" not in case:
            raise RuntimeError(f"Native reference missing for {step_us} us")
        reference = case["native"]
        native_compact = Path(reference["directory"]) / "native_trace.csv"
        if Path(reference["raw_path"]).is_file():
            _, native_meta = write_compact_native(reference)
            reference["compact_trace"] = str(native_compact)
            reference["native_metadata"] = native_meta
        results = {}
        for mode in args.modes:
            if mode == "native":
                continue
            candidate = case[mode]
            if (
                Path(reference["raw_path"]).is_file()
                and Path(candidate["raw_path"]).is_file()
            ):
                results[mode] = compare_candidate(reference, candidate, mode)
                (Path(candidate["directory"]) / "comparison.json").write_text(
                    json.dumps(results[mode], indent=2) + "\n"
                )
                print(
                    f"{mode.upper()} vs native: V {results[mode]['full_run']['voltage']['rms_error_percent']:.6g}% RMS, I {results[mode]['full_run']['current']['rms_error_percent']:.6g}% RMS",
                    flush=True,
                )
                # Full-resolution CSVs are reproducible; retain compact traces only.
                Path(candidate["raw_path"]).unlink()
        if Path(reference["raw_path"]).is_file():
            Path(reference["raw_path"]).unlink()
        for mode in args.modes:
            directory = Path(case[mode]["directory"])
            for path in directory.glob("logs/**/*TaskTimes*.csv"):
                path.unlink()
            for path in directory.glob("logs/**/*SimulationStepTimes.csv"):
                path.unlink()
            for path in directory.glob("logs/**/*SolverPhaseTimes.csv"):
                path.unlink()
        manifest["cases"][str(step_us)] = {
            "native": {k: v for k, v in reference.items() if k != "raw_path"},
            "comparisons": results,
        }
        (OUT / "metrics.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nResults: {OUT / 'metrics.json'}", flush=True)


if __name__ == "__main__":
    main()
