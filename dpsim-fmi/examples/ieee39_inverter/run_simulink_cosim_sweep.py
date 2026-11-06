#!/usr/bin/env python3
"""Run the native-in-Simulink versus imported FMI CS sweep."""
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
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"


def matlab_quote(value: str) -> str:
    return value.replace("'", "''")


def main() -> None:
    runs = (
        ("timestep-sweep", 5.0, (25, 50, 100, 200, 500)),
        ("timestep-sweep-20s", 20.0, (500, 625, 1000)),
    )
    cases = []
    for sweep_dir, stop_time, step_sizes in runs:
        for step_us in step_sizes:
            case_dir = f"dt-{step_us:03d}us"
            input_csv = (
                DATA
                / sweep_dir
                / case_dir
                / "cs"
                / "logs"
                / "IEEE39_Inverter_cs"
                / "IEEE39_Inverter_cs.csv"
            )
            if not input_csv.is_file():
                raise FileNotFoundError(input_csv)
            cases.append(
                {
                    "label": f"{sweep_dir}, {step_us} us",
                    "sweep_dir": sweep_dir,
                    "case_dir": case_dir,
                    "input_csv": str(input_csv.resolve()),
                    "fmu_path": str(FMU.resolve()),
                    "step_s": step_us * 1e-6,
                    "stop_time_s": stop_time,
                }
            )

    output_root = DATA / "simulink-fmu-cosim"
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = output_root / "manifest.json"
    manifest.write_text(json.dumps(cases, indent=2) + "\n")
    model_path = CASE / "matlab" / "IEEE39_Simulink_FMU_Cosim_Validation.slx"
    batch_log = output_root / "matlab_batch.log"
    result_files = [
        output_root / c["sweep_dir"] / c["case_dir"] / "simulink_cosim_timing.json"
        for c in cases
    ]
    start = time.time()
    expression = (
        f"addpath('{matlab_quote(str((CASE / 'matlab').resolve()))}'); "
        f"run_fmu_cosim_sweep('{matlab_quote(str(manifest.resolve()))}',"
        f"'{matlab_quote(str(output_root.resolve()))}',"
        f"'{matlab_quote(str(model_path.resolve()))}');"
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
            complete = all(
                path.is_file() and path.stat().st_mtime >= start
                for path in result_files
            )
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
            raise TimeoutError("Simulink FMU CS sweep exceeded 30 minutes")
    print(f"Completed {len(cases)} side-by-side Simulink native/FMU-CS runs.")
    print(f"Results: {output_root}")
    print(f"Harness model: {model_path}")


if __name__ == "__main__":
    main()
