#!/usr/bin/env python3
"""Generate matched native Simulink references for every stored timestep pair."""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = next(
    (
        path
        for path in (ROOT / "build-fmi-cosim", ROOT / "build-fmi")
        if (path / "CMakeCache.txt").is_file()
    ),
    ROOT / "build-fmi-cosim",
)
DATA = BUILD / "ieee39-inverter"


def matlab_quote(value: str) -> str:
    return value.replace("'", "''")


def main() -> None:
    runs = (
        ("timestep-sweep", 5.0, (25, 50, 100, 200, 500)),
        ("timestep-sweep-20s", 20.0, (500, 625, 1000)),
    )
    cases = []
    for sweep_name, stop_time, step_sizes in runs:
        sweep_dir = DATA / sweep_name
        for step_us in step_sizes:
            pair_dir = sweep_dir / f"dt-{step_us:03d}us"
            input_csv = (
                pair_dir
                / "me"
                / "logs"
                / "IEEE39_Inverter_me"
                / "IEEE39_Inverter_me.csv"
            )
            if not input_csv.is_file():
                raise FileNotFoundError(input_csv)
            output_csv = pair_dir / "native_simulink_reference.csv"
            cases.append(
                {
                    "label": f"{sweep_name}, {step_us} us",
                    "input_csv": str(input_csv.resolve()),
                    "output_csv": str(output_csv.resolve()),
                    "stop_time_s": stop_time,
                }
            )

    manifest = DATA / "native_reference_sweep_manifest.json"
    manifest.write_text(json.dumps(cases, indent=2) + "\n")
    batch_log = DATA / "native_reference_sweep_batch.log"
    timing_files = [
        Path(entry["output_csv"]).with_name("native_timing.json") for entry in cases
    ]
    start = time.time()
    expression = (
        f"addpath('{matlab_quote(str((CASE / 'matlab').resolve()))}'); "
        f"run_native_sweep_references('{matlab_quote(str(manifest.resolve()))}');"
    )
    with batch_log.open("w") as log:
        process = subprocess.Popen(
            ["matlab", "-batch", expression],
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        deadline = time.monotonic() + 1800
        while time.monotonic() < deadline:
            status = process.poll()
            if status is not None:
                if status:
                    raise subprocess.CalledProcessError(status, process.args)
                break
            complete = True
            for timing_file in timing_files:
                if not timing_file.is_file() or timing_file.stat().st_mtime < start:
                    complete = False
                    break
                timing = json.loads(timing_file.read_text())
                if timing.get("reference_type") != "open_loop_voltage_replay":
                    complete = False
                    break
            if complete:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                break
            time.sleep(1)
        else:
            process.kill()
            process.wait()
            raise TimeoutError(
                "MATLAB did not finish the native sweep references in 30 minutes"
            )
    print(f"Native references written for {len(cases)} matched cases.")
    print(f"Manifest: {manifest}")
    print(f"MATLAB log: {batch_log}")


if __name__ == "__main__":
    main()
