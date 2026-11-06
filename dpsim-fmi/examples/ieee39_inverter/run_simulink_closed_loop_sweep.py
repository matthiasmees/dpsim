#!/usr/bin/env python3
"""Run the exact IEEE39 DPsim network in Simulink with the FMI CS inverter."""
from __future__ import annotations

import json
import argparse
import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = next(
    (
        p
        for p in (ROOT / "build-fmi-cosim", ROOT / "build-fmi")
        if (p / "CMakeCache.txt").is_file()
    ),
    ROOT / "build-fmi-cosim",
)
DATA = BUILD / "ieee39-inverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
NETWORK_LIB = BUILD / "dpsim-fmi" / "libieee39_simulink_network.so"
DEFAULT_OUTPUT = DATA / "simulink-closed-loop-sweep-exact-order"
MODEL = CASE / "matlab" / "IEEE39_ClosedLoop.slx"
MATLAB_RUNNER = CASE / "matlab" / "run_ieee39_closed_loop.m"
MATLAB = os.environ.get("MATLAB_BIN", "/usr/local/MATLAB/R2025a/bin/matlab")


def quote(value: Path | str) -> str:
    return str(value).replace("'", "''")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stop-time-s",
        type=float,
        default=None,
        help="run only this duration instead of the default 5s/20s suite",
    )
    parser.add_argument(
        "--steps-us",
        type=int,
        nargs="+",
        default=None,
        help="timestep sweep for --stop-time-s",
    )
    parser.add_argument(
        "--output-name",
        default=None,
        help="fresh output directory name under the build's ieee39-inverter folder",
    )
    parser.add_argument(
        "--reference-sweep-dir",
        default=None,
        help="directory containing matching DPsim CS traces for the selected duration",
    )
    parser.add_argument(
        "--model-path",
        type=Path,
        default=MODEL,
        help="Simulink model path (use a temporary path to preserve an existing model)",
    )
    args = parser.parse_args()
    if not FMU.is_file() or not NETWORK_LIB.is_file():
        raise FileNotFoundError(
            "Build the IEEE39 FMU and Simulink network library first."
        )
    cases = []
    if args.stop_time_s is None:
        runs = (
            ("timestep-sweep", 5.0, (25, 50, 100, 200, 500)),
            ("timestep-sweep-20s", 20.0, (500, 625, 1000)),
        )
        output = DEFAULT_OUTPUT if args.output_name is None else DATA / args.output_name
    else:
        stop = args.stop_time_s
        steps = tuple(args.steps_us or (25, 50, 100, 200, 500, 625, 1000))
        sweep = args.reference_sweep_dir or f"timestep-sweep-{stop:g}s"
        runs = ((sweep, stop, steps),)
        output = DATA / (args.output_name or f"simulink-closed-loop-{stop:g}s")
        if stop <= 0 or any(step <= 0 for step in steps):
            parser.error("duration and timesteps must be positive")
    for sweep, stop, steps in runs:
        for step_us in steps:
            case_dir = f"dt-{step_us:03d}us"
            cpp_csv = (
                DATA
                / sweep
                / case_dir
                / "cs"
                / "logs"
                / "IEEE39_Inverter_cs"
                / "IEEE39_Inverter_cs.csv"
            )
            if not cpp_csv.is_file():
                raise FileNotFoundError(
                    f"Missing matched DPsim CS reference: {cpp_csv}"
                )
            cases.append(
                {
                    "label": f"{stop:g} s, {step_us} us",
                    "sweep_dir": sweep,
                    "case_dir": case_dir,
                    "step_s": step_us * 1e-6,
                    "stop_time_s": stop,
                    "load_step_s": stop / 2,
                }
            )

    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "manifest.json"
    manifest.write_text(json.dumps(cases, indent=2) + "\n")
    log_path = output / "matlab_batch.log"
    expected = []
    for case in cases:
        case_output = output / case["sweep_dir"] / case["case_dir"]
        expected.extend(
            (
                case_output / "simulink_closed_loop.mat",
                case_output / "simulink_closed_loop_pcc.csv",
                case_output / "simulink_closed_loop_timing.json",
            )
        )
    expression = (
        f"addpath('{quote(CASE / 'matlab')}'); "
        f"run_ieee39_closed_loop('{quote(manifest)}','{quote(output)}',"
        f"'{quote(FMU)}','{quote(NETWORK_LIB)}','{quote(args.model_path)}');"
    )
    env = os.environ.copy()
    preload = os.environ.get("DPsim_MATLAB_LIBSTDCXX", "/usr/lib64/libstdc++.so.6")
    env["LD_PRELOAD"] = preload + (
        ":" + env["LD_PRELOAD"] if env.get("LD_PRELOAD") else ""
    )
    started = time.time()
    with log_path.open("w") as log:
        proc = subprocess.Popen(
            [MATLAB, "-batch", expression],
            cwd=ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        deadline = time.monotonic() + 4 * 60 * 60
        while time.monotonic() < deadline:
            code = proc.poll()
            if code is not None:
                if code:
                    raise subprocess.CalledProcessError(code, proc.args)
                break
            if all(
                path.is_file()
                and path.stat().st_size > 0
                and path.stat().st_mtime >= started
                for path in expected
            ):
                # The result writer emits timing before its compact CSV, so wait
                # for every nonempty artifact before stopping MATLAB. MATLAB on
                # this host may then hang during its own exit cleanup.
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                break
            time.sleep(5)
        else:
            proc.kill()
            proc.wait()
            raise TimeoutError("Closed-loop MATLAB timestep sweep exceeded four hours")
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        tail = log_path.read_text(errors="replace")[-8000:]
        raise RuntimeError(f"MATLAB sweep is incomplete: {missing}\n{tail}")
    print(
        f"Completed {len(cases)} closed-loop Simulink IEEE39 cases in {time.time()-started:.1f} s."
    )
    print(f"Outputs: {output}")
    print(f"Model: {args.model_path}")


if __name__ == "__main__":
    main()
