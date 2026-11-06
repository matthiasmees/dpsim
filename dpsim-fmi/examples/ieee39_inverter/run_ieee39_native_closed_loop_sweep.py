#!/usr/bin/env python3
"""Run the native Simulink inverter in the full closed-loop IEEE39 network."""
from __future__ import annotations

import argparse
import json
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
MATLAB = os.environ.get("MATLAB_BIN", "/usr/local/MATLAB/R2025a/bin/matlab")
MODEL = CASE / "matlab" / "IEEE39_ClosedLoop_Native.slx"


def q(value: str | Path) -> str:
    return str(value).replace("'", "''")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stop-time-s", type=float, default=20.0)
    parser.add_argument(
        "--steps-us", type=int, nargs="+", default=[25, 50, 100, 200, 500, 625, 1000]
    )
    parser.add_argument("--load-step-s", type=float, default=None)
    parser.add_argument("--output-name", default="native-closed-loop-20s-sweep")
    parser.add_argument(
        "--model-path",
        type=Path,
        default=MODEL,
        help="Simulink model path (use a temporary path to preserve an existing model)",
    )
    args = parser.parse_args()
    stop = args.stop_time_s
    event = stop / 2 if args.load_step_s is None else args.load_step_s
    if stop <= 0 or not 0 < event < stop:
        parser.error("load step must be inside a positive duration")
    if any(step <= 0 for step in args.steps_us):
        parser.error("timesteps must be positive")
    sweep_dir = "timestep-sweep-20s" if stop == 20.0 else f"timestep-sweep-{stop:g}s"
    output = DATA / args.output_name
    output.mkdir(parents=True, exist_ok=True)
    cases = []
    for step_us in args.steps_us:
        case_dir = f"dt-{step_us:03d}us"
        cases.append(
            {
                "label": f"{stop:g} s, {step_us} us",
                "sweep_dir": sweep_dir,
                "case_dir": case_dir,
                "step_s": step_us * 1e-6,
                "stop_time_s": stop,
                "load_step_s": event,
            }
        )
    manifest = output / "manifest.json"
    manifest.write_text(json.dumps(cases, indent=2) + "\n")
    expected = [
        output / x["sweep_dir"] / x["case_dir"] / "native_closed_loop_timing.json"
        for x in cases
    ]
    expression = (
        f"addpath('{q(CASE / 'matlab')}'); "
        f"run_ieee39_native_closed_loop('{q(manifest)}','{q(output)}','{q(args.model_path)}');"
    )
    env = os.environ.copy()
    preload = os.environ.get("DPsim_MATLAB_LIBSTDCXX", "/usr/lib64/libstdc++.so.6")
    env["LD_PRELOAD"] = preload + (
        ":" + env["LD_PRELOAD"] if env.get("LD_PRELOAD") else ""
    )
    log_path = output / "matlab_batch.log"
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
            if all(p.is_file() and p.stat().st_mtime >= started for p in expected):
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
            raise TimeoutError("Native Simulink IEEE39 sweep exceeded four hours")
    missing = [str(p) for p in expected if not p.is_file()]
    if missing:
        raise RuntimeError(
            f"Missing native closed-loop outputs: {missing}\n{log_path.read_text(errors='replace')[-8000:]}"
        )
    print(
        f"Completed {len(cases)} native closed-loop cases in {time.time()-started:.1f} s."
    )
    print(f"Outputs: {output}")
    print(f"Model: {args.model_path}")


if __name__ == "__main__":
    main()
