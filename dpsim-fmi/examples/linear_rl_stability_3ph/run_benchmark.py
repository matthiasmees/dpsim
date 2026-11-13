#!/usr/bin/env python3
"""Run native, FMI ME and FMI CS balanced three-phase EMT stability cases."""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

RG, RF, L = 1.0, 0.1, 20e-6
DEFAULT_STEPS_US = [20, 30, 35, 38, 39, 40, 40.13, 40.25, 41, 45, 50]
STOP_TIME = 20e-3


def multiplier(step: float) -> float:
    return (1 + RG / RF) * math.exp(-RF * step / L) - RG / RF


def critical_step() -> float:
    return L / RF * math.log((RG + RF) / (RG - RF))


def trace(path: Path) -> tuple[list[float], list[list[float]], list[list[float]]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        rows = list(reader)
    phases = [[float(row[f"i_abc_A_{phase}"]) for phase in range(3)] for row in rows]
    voltages = [[float(row[f"u_abc_V_{phase}"]) for phase in range(3)] for row in rows]
    times = [float(row["time"]) for row in rows]
    return times, phases, voltages


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runner",
        type=Path,
        default=Path("build-fmi/dpsim-fmi/FMI_EMT_Ph3_LinearRLStability"),
    )
    parser.add_argument(
        "--me-fmu",
        type=Path,
        default=Path(
            "dpsim-fmi/examples/linear_rl_stability_3ph/sources/LinearRL3Ph_ME.fmu"
        ),
    )
    parser.add_argument(
        "--cs-fmu",
        type=Path,
        default=Path(
            "dpsim-fmi/examples/linear_rl_stability_3ph/sources/LinearRL3Ph_CS.fmu"
        ),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("build-fmi/linear_rl_stability_3ph/results")
    )
    parser.add_argument("--steps-us", type=float, nargs="+", default=DEFAULT_STEPS_US)
    parser.add_argument("--stop-time", type=float, default=STOP_TIME)
    args = parser.parse_args()
    import matplotlib.pyplot as plt

    root = Path.cwd()
    resolve = lambda value: value if value.is_absolute() else root / value
    runner, me_fmu, cs_fmu, output = map(
        resolve, (args.runner, args.me_fmu, args.cs_fmu, args.output)
    )
    for path in (runner, me_fmu, cs_fmu):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; build the example and FMUs first")
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    steps = sorted(set(args.steps_us) | {50.0})
    traces: dict[tuple[str, float], tuple[list[float], list[list[float]]]] = {}
    summary: list[dict[str, object]] = []

    for step_us in steps:
        step = step_us * 1e-6
        variants = [("native", "-"), ("me", str(me_fmu)), ("cs", str(cs_fmu))]
        if step_us == 50.0:
            variants.append(("cs-rollback", str(cs_fmu)))
        for mode, fmu in variants:
            work = output / "runs" / f"{mode}_{step_us:g}us"
            if work.exists():
                shutil.rmtree(work)
            work.mkdir(parents=True)
            subprocess.run(
                [str(runner), mode, fmu, f"{step:.17g}", f"{args.stop_time:.17g}"],
                cwd=work,
                check=True,
            )
            produced = work / "logs" / f"linear_rl_3ph_{mode}.csv"
            destination = raw / f"{mode}_{step_us:g}us.csv"
            shutil.copy2(produced, destination)
            times, currents, voltages = trace(destination)
            if not currents or any(abs(value) > 1e-12 for value in currents[0]):
                raise AssertionError(
                    f"{mode} at {step_us} us did not start with zero ABC current"
                )
            if any(abs(measured) > 1e-12 for measured in voltages[0]):
                raise AssertionError(
                    f"{mode} at {step_us} us has inconsistent initial ABC terminal voltage"
                )
            peak = max(abs(value) for row in currents for value in row)
            final = currents[-1]
            balance = max(abs(sum(row)) for row in currents)
            traces[(mode, step_us)] = (times, currents)
            summary.append(
                {
                    "mode": mode,
                    "step_us": step_us,
                    "cs_multiplier_per_phase": multiplier(step),
                    "cs_stable_by_analysis": abs(multiplier(step)) < 1,
                    "peak_abs_phase_current_A": peak,
                    "final_i_a_A": final[0],
                    "final_i_b_A": final[1],
                    "final_i_c_A": final[2],
                    "max_abs_abc_sum_A": balance,
                    "trace_csv": str(destination.relative_to(output)),
                }
            )

    for row in summary:
        step_us = float(row["step_us"])
        reference = traces[("native", step_us)][1]
        mode_values = traces[(str(row["mode"]), step_us)][1]
        differences = [
            value - ref
            for values, refs in zip(mode_values, reference)
            for value, ref in zip(values, refs)
        ]
        row["max_abs_delta_to_native_A"] = max(map(abs, differences), default=0.0)
        row["rms_delta_to_native_A"] = (
            math.sqrt(sum(error * error for error in differences) / len(differences))
            if differences
            else 0.0
        )

    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    metadata = {
        "model": "balanced three-phase EMT, grounded-wye independent RL branches",
        "phase_voltage_rms_V": 10.0,
        "frequency_Hz": 50.0,
        "Rg_ohm_per_phase": RG,
        "Rf_ohm_per_phase": RF,
        "L_H_per_phase": L,
        "initial_current_abc_A": [0, 0, 0],
        "nominal_step_s": 50e-6,
        "stop_time_s": args.stop_time,
        "analytical_cs_stability_limit_s": critical_step(),
        "analytical_cs_stability_limit_us": critical_step() * 1e6,
        "multiplier": "q=(1+Rg/Rf)*exp(-Rf*h/L)-Rg/Rf; same for a,b,c",
        "summary_csv": "summary.csv",
    }
    (output / "metrics.json").write_text(json.dumps(metadata, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    colors = ("tab:blue", "tab:orange", "tab:green")
    for mode, style in (("native", "-"), ("me", "--"), ("cs-rollback", ":")):
        times, values = traces[(mode, 50.0)]
        for phase, label in enumerate("abc"):
            axes[0, 0].plot(
                [t * 1e3 for t in times],
                [row[phase] for row in values],
                style,
                color=colors[phase],
                label=f"{mode.upper()} {label}",
            )
    axes[0, 0].set(
        title="Native, ME and rollback/Newton, h = 50 µs",
        xlabel="Time (ms)",
        ylabel="Phase current (A)",
    )
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].legend(ncol=3, fontsize=7)
    times, values = traces[("cs", 50.0)]
    for phase, label in enumerate("abc"):
        axes[0, 1].plot(
            [t * 1e3 for t in times],
            [row[phase] for row in values],
            color=colors[phase],
            label=label,
        )
    axes[0, 1].set_yscale("symlog", linthresh=10)
    axes[0, 1].set_yticks([-1e60, -1e40, -1e20, 0, 1e20, 1e40, 1e60])
    axes[0, 1].set_yticklabels(
        [
            r"$-10^{60}$",
            r"$-10^{40}$",
            r"$-10^{20}$",
            "0",
            r"$10^{20}$",
            r"$10^{40}$",
            r"$10^{60}$",
        ]
    )
    axes[0, 1].set(
        title="Explicit CS, h = 50 µs",
        xlabel="Time (ms)",
        ylabel="Phase current (A, symlog)",
    )
    axes[0, 1].grid(which="both", alpha=0.3)
    axes[0, 1].legend()
    for mode, marker in (("native", "o"), ("me", "s"), ("cs", "^")):
        peaks = [
            next(
                row["peak_abs_phase_current_A"]
                for row in summary
                if row["mode"] == mode and row["step_us"] == step
            )
            for step in steps
        ]
        axes[1, 0].semilogy(steps, peaks, marker + "-", label=mode.upper())
    rollback_peak = next(
        row["peak_abs_phase_current_A"]
        for row in summary
        if row["mode"] == "cs-rollback" and row["step_us"] == 50.0
    )
    axes[1, 0].semilogy(
        [50.0],
        [rollback_peak],
        "D",
        color="tab:purple",
        label="CS Rollback/Newton (50 µs)",
    )
    axes[1, 0].axvline(
        critical_step() * 1e6,
        color="red",
        linestyle="--",
        label=f"predicted {critical_step()*1e6:.3f} µs",
    )
    axes[1, 0].set(
        title="Peak current over the step sweep",
        xlabel="Step (µs)",
        ylabel="Peak |phase current| (A)",
    )
    axes[1, 0].grid(which="both", alpha=0.3)
    axes[1, 0].legend()
    q_values = [multiplier(step * 1e-6) for step in steps]
    axes[1, 1].plot(steps, q_values, "o-")
    axes[1, 1].axhline(-1, color="red", linestyle="--")
    axes[1, 1].axhline(1, color="gray", linestyle=":")
    axes[1, 1].axvline(critical_step() * 1e6, color="red", linestyle="--")
    axes[1, 1].set(
        title="Shared per-phase CS multiplier q(h)", xlabel="Step (µs)", ylabel="q"
    )
    axes[1, 1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output / "linear_rl_3ph_stability.png", dpi=180)
    print(json.dumps(metadata, indent=2))
    print(f"Wrote three-phase results to {output}")


if __name__ == "__main__":
    main()
