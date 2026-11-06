#!/usr/bin/env python3
"""Rerun DPsim CS cases in a fresh folder for the closed-loop Simulink check."""
from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BUILD = next(
    (
        p
        for p in (ROOT / "build-fmi-cosim", ROOT / "build-fmi")
        if (p / "CMakeCache.txt").is_file()
    ),
    ROOT / "build-fmi-cosim",
)
DATA = BUILD / "ieee39-inverter"
RUNNER = BUILD / "dpsim-fmi" / "FMI_IEEE39_AvgInverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
OUTPUT = DATA / "dpsim-cs-matched-rerun"


def main() -> None:
    if not RUNNER.is_file() or not FMU.is_file():
        raise FileNotFoundError("Build the IEEE39 DPsim runner and FMU first.")
    runs = (
        ("timestep-sweep", 5.0, 2.5, (25, 50, 100, 200, 500)),
        ("timestep-sweep-20s", 20.0, 10.0, (500, 625, 1000)),
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    for sweep, stop, event, steps in runs:
        for step_us in steps:
            case = f"dt-{step_us:03d}us"
            folder = OUTPUT / sweep / case / "cs"
            folder.mkdir(parents=True, exist_ok=True)
            start = time.perf_counter()
            result = subprocess.run(
                [
                    str(RUNNER),
                    "cs",
                    str(FMU),
                    str(stop),
                    str(step_us * 1e-6),
                    str(event),
                ],
                cwd=folder,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            wall = time.perf_counter() - start
            (folder / "runner.log").write_text(result.stdout)
            if result.returncode:
                raise RuntimeError(
                    f"DPsim CS {stop}s/{step_us}us failed:\n{result.stdout[-6000:]}"
                )
            match = re.search(
                r"runtime_s=([0-9.eE+-]+)\s+initialization_s=([0-9.eE+-]+)\s+simulation_s=([0-9.eE+-]+)\s+steps=([0-9]+)",
                result.stdout,
            )
            if not match:
                raise RuntimeError(
                    f"Could not parse run timing for {stop}s/{step_us}us"
                )
            log = folder / "logs" / "IEEE39_Inverter_cs" / "IEEE39_Inverter_cs.csv"
            if not log.is_file():
                raise FileNotFoundError(log)
            summary = {
                "sweep_dir": sweep,
                "case_dir": case,
                "stop_time_s": stop,
                "load_step_s": event,
                "step_us": step_us,
                "total_runtime_s": float(match.group(1)),
                "initialization_s": float(match.group(2)),
                "simulation_s": float(match.group(3)),
                "steps": int(match.group(4)),
                "process_wall_s": wall,
                "log_path": str(log.relative_to(DATA)),
            }
            summaries.append(summary)
            (folder.parent / "timing.json").write_text(
                json.dumps(summary, indent=2) + "\n"
            )
            print(
                f"Done {stop:g}s / {step_us}us: {summary['total_runtime_s']:.2f}s total, {summary['steps']} steps",
                flush=True,
            )
    (OUTPUT / "sweep_timing.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(f"Fresh DPsim CS reference logs: {OUTPUT}")


if __name__ == "__main__":
    main()
