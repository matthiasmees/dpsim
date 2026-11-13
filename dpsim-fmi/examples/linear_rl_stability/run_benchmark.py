#!/usr/bin/env python3
"""Run native, FMI ME and FMI CS variants and analyze their CSV traces."""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

VS = 10.0
RG = 1.0
RF = 0.1
INDUCTANCE = 20e-6
DEFAULT_STEPS = [20, 30, 35, 38, 39, 40, 40.13, 40.25, 41, 45, 50]
FINAL_TIME = 5e-3


def stability_limit() -> float:
    """Exact CS boundary for an exact internal RL step and explicit grid ZOH."""
    return INDUCTANCE / RF * math.log((RG + RF) / (RG - RF))


def cs_multiplier(step: float) -> float:
    a = math.exp(-RF * step / INDUCTANCE)
    return a - (RG / RF) * (1.0 - a)


def analytic_current(time: float) -> float:
    total_r = RG + RF
    return VS / total_r * (1.0 - math.exp(-total_r * time / INDUCTANCE))


def _load_trace(path: Path) -> tuple[list[float], list[float]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"No CSV header in {path}")
        reader.fieldnames = [key.strip() for key in reader.fieldnames]
        time_key = next(
            (key for key in reader.fieldnames if key.lower() in {"time", "t"}), None
        )
        current_key = next(
            (key for key in reader.fieldnames if key in {"i_A", "i_a", "current_A"}),
            None,
        )
        if time_key is None or current_key is None:
            raise ValueError(
                f"Expected time and i_A columns in {path}; found {reader.fieldnames}"
            )
        rows = [(float(row[time_key]), float(row[current_key])) for row in reader]
    return [r[0] for r in rows], [r[1] for r in rows]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runner", type=Path, default=Path("build-fmi/dpsim-fmi/FMI_LinearRLStability")
    )
    parser.add_argument(
        "--me-fmu",
        type=Path,
        default=Path("dpsim-fmi/examples/linear_rl_stability/sources/LinearRL_ME.fmu"),
    )
    parser.add_argument(
        "--cs-fmu",
        type=Path,
        default=Path("dpsim-fmi/examples/linear_rl_stability/sources/LinearRL_CS.fmu"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("build-fmi/linear_rl_stability/results")
    )
    parser.add_argument("--steps-us", type=float, nargs="+", default=DEFAULT_STEPS)
    parser.add_argument("--stop-time", type=float, default=FINAL_TIME)
    args = parser.parse_args()
    import matplotlib.pyplot as plt

    repo = Path.cwd()
    runner = args.runner if args.runner.is_absolute() else repo / args.runner
    me_fmu = args.me_fmu if args.me_fmu.is_absolute() else repo / args.me_fmu
    cs_fmu = args.cs_fmu if args.cs_fmu.is_absolute() else repo / args.cs_fmu
    output = args.output if args.output.is_absolute() else repo / args.output
    for path in (runner, me_fmu, cs_fmu):
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing {path}; build the target and export both FMUs first"
            )
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    traces: dict[tuple[str, float], tuple[list[float], list[float]]] = {}
    summary_rows: list[dict[str, float | str | bool]] = []
    step_values = sorted(set(args.steps_us) | {50.0})
    for step_us in step_values:
        step = step_us * 1e-6
        for mode, fmu in (("native", "-"), ("me", str(me_fmu)), ("cs", str(cs_fmu))):
            work = output / "runs" / f"{mode}_{step_us:g}us"
            if work.exists():
                shutil.rmtree(work)
            work.mkdir(parents=True)
            subprocess.run(
                [str(runner), mode, fmu, f"{step:.17g}", f"{args.stop_time:.17g}"],
                cwd=work,
                check=True,
            )
            produced = work / "logs" / f"linear_rl_{mode}.csv"
            if not produced.is_file():
                raise FileNotFoundError(f"Runner did not produce {produced}")
            destination = raw / f"{mode}_{step_us:g}us.csv"
            shutil.copy2(produced, destination)
            times, currents = _load_trace(destination)
            traces[(mode, step_us)] = (times, currents)
            peak = max((abs(v) for v in currents), default=0.0)
            errors = [v - analytic_current(t) for t, v in zip(times, currents)]
            end_error = abs(errors[-1]) if errors else math.inf
            rms_error = (
                math.sqrt(sum(e * e for e in errors) / len(errors))
                if errors
                else math.inf
            )
            summary_rows.append(
                {
                    "mode": mode,
                    "step_us": step_us,
                    "cs_multiplier": cs_multiplier(step),
                    "cs_stable_by_analysis": abs(cs_multiplier(step)) < 1.0,
                    "peak_current_A": peak,
                    "final_current_A": currents[-1] if currents else math.nan,
                    "final_error_A": end_error,
                    "max_abs_error_A": max((abs(e) for e in errors), default=math.inf),
                    "rms_error_A": rms_error,
                    "trace_csv": str(destination.relative_to(output)),
                }
            )

    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    result = {
        "parameters": {
            "Vs_V": VS,
            "Rg_ohm": RG,
            "Rf_ohm": RF,
            "L_H": INDUCTANCE,
            "nominal_step_s": 50e-6,
            "stop_time_s": args.stop_time,
        },
        "analytical_cs_stability_limit_s": stability_limit(),
        "analytical_cs_stability_limit_us": stability_limit() * 1e6,
        "step_multiplier": "q=(1+Rg/Rf)*exp(-Rf*h/L)-Rg/Rf",
        "summary_csv": "summary.csv",
    }
    (output / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")

    fig = plt.figure(figsize=(12, 8))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
    native_axis = fig.add_subplot(grid[0, 0])
    cs_axis = fig.add_subplot(grid[0, 1])
    sweep_axis = fig.add_subplot(grid[1, :])
    nominal = 50.0
    for mode, style in (("native", "-"), ("me", "--")):
        ts, ys = traces[(mode, nominal)]
        native_axis.plot([t * 1e3 for t in ts], ys, style, label=mode.upper())
    ts = [
        idx * nominal * 1e-6
        for idx in range(int(args.stop_time / (nominal * 1e-6)) + 1)
    ]
    native_axis.plot(
        [t * 1e3 for t in ts], [analytic_current(t) for t in ts], ":", label="analytic"
    )
    native_axis.set(
        title="Native and ME at 50 µs", xlabel="Time (ms)", ylabel="Current (A)"
    )
    native_axis.grid(True, alpha=0.3)
    native_axis.legend()

    cs_times, cs_currents = traces[("cs", nominal)]
    cs_axis.plot([t * 1e3 for t in cs_times], cs_currents, color="tab:red")
    cs_axis.set_yscale("symlog", linthresh=10)
    cs_axis.set(
        title="Explicit CS at 50 µs (diverging)",
        xlabel="Time (ms)",
        ylabel="Current (A, symlog)",
    )
    cs_axis.grid(True, which="both", alpha=0.3)

    sweep = step_values
    cs_peak = [max(abs(y) for y in traces[("cs", s)][1]) for s in sweep]
    sweep_axis.semilogy(sweep, [max(v, 1e-16) for v in cs_peak], "o-")
    sweep_axis.axvline(
        stability_limit() * 1e6,
        color="tab:red",
        linestyle="--",
        label=f"analytic {stability_limit()*1e6:.3f} µs",
    )
    sweep_axis.set(
        title="Explicit CS peak current across step sweep",
        xlabel="Step (µs)",
        ylabel="Peak |i| (A)",
    )
    sweep_axis.grid(True, which="both", alpha=0.3)
    sweep_axis.legend()
    fig.tight_layout()
    fig.savefig(output / "linear_rl_stability.png", dpi=180)
    print(json.dumps(result, indent=2))
    print(f"Wrote CSV and plot results to {output}")


if __name__ == "__main__":
    main()
