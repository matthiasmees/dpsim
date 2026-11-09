#!/usr/bin/env python3
"""Smoke-check then run the 5 ms, 10 s small-grid coarse-step stress case."""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

import run_dpsim_10s_sweep as sweep

ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "build-fmi"
CASE = Path(__file__).resolve().parent
EXE = BUILD / "dpsim-fmi" / "FMI_SmallGrid_Inverter"
FMU = CASE / "sources" / "IEEE39AvgInverter.fmu"
OUT = BUILD / "small-inverter-grid" / "results" / "dpsim_5ms_stress"
STEP_US = 5000


def smoke_check() -> None:
    dt = STEP_US * 1e-6
    smoke = OUT / "smoke_100ms"
    for mode in ("native", "me", "cs"):
        directory = smoke / mode
        directory.mkdir(parents=True, exist_ok=True)
        completed = subprocess.run(
            [str(EXE), mode, str(FMU), "0.1", str(dt), "0.05", "powerflow"],
            cwd=directory,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=True,
        )
        (directory / "runner.log").write_text(completed.stdout)
        name = f"SmallGrid_Inverter_{mode}_powerflow"
        raw = directory / "logs" / name / f"{name}.csv"
        with raw.open(newline="") as stream:
            reader = csv.DictReader(stream, skipinitialspace=True)
            count = 0
            for row in reader:
                values = (
                    [float(row["time"])]
                    + [float(row[f"pcc.v_{phase}"]) for phase in range(3)]
                    + [float(row[f"inverter.i_{phase}"]) for phase in range(3)]
                )
                if not all(math.isfinite(value) for value in values):
                    raise RuntimeError(f"Non-finite smoke output in {mode}")
                count += 1
        if count != 21:
            raise RuntimeError(f"Expected 21 smoke samples for {mode}, got {count}")
        print(f"smoke {mode.upper()}: {count} finite samples, 100 ms", flush=True)


def full_run() -> None:
    sweep.OUT = OUT
    OUT.mkdir(parents=True, exist_ok=True)
    records = {mode: sweep.run_case(mode, STEP_US) for mode in sweep.MODES}
    native = records["native"]
    _, native_metadata = sweep.write_compact_native(native)
    native["compact_trace"] = str(Path(native["directory"]) / "native_trace.csv")
    native["native_metadata"] = native_metadata
    comparisons = {}
    for mode in ("me", "cs"):
        candidate = records[mode]
        comparisons[mode] = sweep.compare_candidate(native, candidate, mode)
        (Path(candidate["directory"]) / "comparison.json").write_text(
            json.dumps(comparisons[mode], indent=2) + "\n"
        )
        print(
            f"{mode.upper()} vs native: V "
            f"{comparisons[mode]['full_run']['voltage']['rms_error_percent']:.6g}% RMS, I "
            f"{comparisons[mode]['full_run']['current']['rms_error_percent']:.6g}% RMS",
            flush=True,
        )

    for mode in sweep.MODES:
        record = records[mode]
        Path(record["raw_path"]).unlink()
        directory = Path(record["directory"])
        for pattern in (
            "logs/**/*TaskTimes*.csv",
            "logs/**/*SimulationStepTimes.csv",
            "logs/**/*SolverPhaseTimes.csv",
        ):
            for path in directory.glob(pattern):
                path.unlink()

    manifest = {
        "configuration": {
            "duration_s": sweep.STOP_S,
            "load_step_s": sweep.LOAD_STEP_S,
            "steps_us": [STEP_US],
            "methods": list(sweep.MODES),
            "initialization": "network powerflow; native and FMU filter states initialized to the same 60 Hz periodic state",
            "error_window_s": sweep.ERROR_WINDOW_S,
            "case_type": "coarse-step EMT stress test; not a waveform-accuracy reference",
        },
        "cases": {
            str(STEP_US): {
                "native": {
                    key: value for key, value in native.items() if key != "raw_path"
                },
                "comparisons": comparisons,
            }
        },
    }
    path = OUT / "metrics.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Results: {path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke-only", action="store_true", help="run only the 100 ms checks"
    )
    args = parser.parse_args()
    smoke_check()
    if not args.smoke_only:
        full_run()


if __name__ == "__main__":
    main()
