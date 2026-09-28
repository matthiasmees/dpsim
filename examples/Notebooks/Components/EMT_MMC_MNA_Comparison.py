#!/usr/bin/env python3
"""Compare the unchanged SSN P2P example with its direct-MNA MMC copy.

Requires numpy and matplotlib. Inputs are the two DPSim CSV files. Pass the
actual EMT steps: both examples log sim.time() AFTER step(), one step later
than the physical solution time (except the initial row). This script corrects
that shared timestamp convention before matching/interpolating the samples.
"""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STATIONS = ("rectifier", "inverter")
# Signal, label, plotting unit, unit divisor, engineering normalization.
SIGNALS = (
    ("p_ac", "Active power", "MW", 1e6, 1e9),
    ("q_ac", "Reactive power", "Mvar", 1e6, 1e9),
    ("vdc", "DC voltage", "kV", 1e3, 640e3),
    ("idc", "DC current", "kA", 1e3, 1e9 / 640e3),
    ("energy", "Capacitor energy", "MJ", 1e6, 60e6),
    ("i_sigma_z", "Circulating current (zero component)", "A", 1, 1e9 / 640e3 / 3),
    ("i_delta_d", "AC current d", "A", 1, np.sqrt(2 / 3) * 1e9 / 333e3),
    ("i_delta_q", "AC current q", "A", 1, np.sqrt(2 / 3) * 1e9 / 333e3),
)


def read_csv(path, step):
    if step <= 0:
        raise ValueError("EMT steps must be positive")
    with path.open() as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        names = [name.strip() for name in next(reader)]
        values = np.loadtxt(stream, delimiter=",", ndmin=2)
    if values.shape[0] < 2 or values.shape[1] != len(names):
        raise ValueError(f"Invalid or incomplete CSV: {path}")
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite simulation output: {path}")
    data = dict(zip(names, values.T))
    if "time" not in data:
        raise ValueError(f"Missing time column: {path}")
    data["time"] = data["time"].copy()
    data["time"][1:] -= step
    if np.any(np.diff(data["time"]) <= 0):
        raise ValueError(f"Non-monotonic physical time: {path}")
    return data


def aligned(reference, native):
    rt, nt = reference["time"], native["time"]
    tolerance = 1e-10
    if abs(rt[0] - nt[0]) > tolerance or abs(rt[-1] - nt[-1]) > max(np.diff(rt).max(), np.diff(nt).max()) + tolerance:
        raise ValueError("Runs cover different time intervals; use equal final times")
    mask = (rt >= nt[0] - tolerance) & (rt <= nt[-1] + tolerance)
    time = rt[mask]
    keys = set(reference) & set(native) - {"time"}
    same_grid = len(time) == len(nt) and np.allclose(time, nt, rtol=0, atol=tolerance)
    a = {key: reference[key][mask] for key in keys}
    b = {key: native[key] if same_grid else np.interp(time, nt, native[key]) for key in keys}
    if "reference.p_inverter" not in keys:
        raise ValueError("Missing power-reference trace")
    # Compare held references at recorded times; unequal event positions may
    # indicate an incompatible scenario, beyond the one-sample interpolation edge.
    if same_grid and not np.allclose(a["reference.p_inverter"], b["reference.p_inverter"], rtol=0, atol=1):
        raise ValueError("Power-reference sequences differ")
    return time, a, b, not same_grid


def save_figure(fig, output, name):
    fig.savefig(output / f"{name}.png", dpi=160)
    fig.savefig(output / f"{name}.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("native", type=Path)
    parser.add_argument("--reference-dt", type=float, default=20e-6)
    parser.add_argument("--mna-dt", type=float, default=20e-6)
    parser.add_argument("--output", type=Path, default=Path("outputs/mmc_mna/comparison"))
    parser.add_argument("--max-error-percent", type=float, help="Optional acceptance limit relative to engineering bases")
    args = parser.parse_args()
    reference = read_csv(args.reference, args.reference_dt)
    native = read_csv(args.native, args.mna_dt)
    time, a, b, interpolated = aligned(reference, native)
    for station in STATIONS:
        for signal, *_ in SIGNALS:
            if f"{station}.{signal}" not in a:
                raise ValueError(f"Missing required comparison signal: {station}.{signal}")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    colors = ("#1864ab", "#d9480f")

    for station in STATIONS:
        fig, axes = plt.subplots(4, 2, figsize=(13, 12), sharex=True, layout="constrained")
        for ax, (signal, label, unit, divisor, _) in zip(axes.flat, SIGNALS):
            key = f"{station}.{signal}"
            ax.plot(time, a[key] / divisor, color=colors[0], lw=1.3, label="SSN MMC")
            ax.plot(time, b[key] / divisor, color=colors[1], lw=1, ls="--", label="MMC MNA")
            if signal == "p_ac" and station == "inverter":
                ax.plot(time, a["reference.p_inverter"] / divisor, color="0.4", ls=":", lw=1, label="Setpoint")
            ax.set(title=label, ylabel=unit)
        for ax in axes[-1]: ax.set_xlabel("Time [s]")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"{station.capitalize()}: same P2P case, SSN and direct MNA")
        save_figure(fig, output, f"{station}_overview")

    fig, axes = plt.subplots(4, 2, figsize=(13, 12), sharex=True, layout="constrained")
    for ax, (signal, label, unit, divisor, _) in zip(axes.flat, SIGNALS):
        for station, color in zip(STATIONS, colors):
            key = f"{station}.{signal}"
            ax.plot(time, (b[key] - a[key]) / divisor, lw=1, color=color, label=station)
        ax.set(title=label, ylabel=f"Difference [{unit}]")
    for ax in axes[-1]: ax.set_xlabel("Time [s]")
    axes[0, 0].legend()
    fig.suptitle("Difference: MMC MNA minus SSN MMC")
    save_figure(fig, output, "differences")

    for event in (4, 6, 8):
        mask = (time >= event - 0.025) & (time <= event + 0.25)
        if mask.sum() < 2: continue
        fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True, layout="constrained")
        for col, station in enumerate(STATIONS):
            for row, index in enumerate((0, 2, 3)):
                signal, label, unit, divisor, _ = SIGNALS[index]
                key = f"{station}.{signal}"
                ax = axes[row, col]
                ax.plot(time[mask], a[key][mask] / divisor, color=colors[0], label="SSN MMC")
                ax.plot(time[mask], b[key][mask] / divisor, color=colors[1], ls="--", label="MMC MNA")
                ax.set(title=f"{station}: {label}", ylabel=unit)
                ax.axvline(event, color="0.5", lw=0.7, ls=":")
        for ax in axes[-1]: ax.set_xlabel("Time [s]")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Setpoint step at {event} s")
        save_figure(fig, output, f"step_{event}s")

    # Final 40 ms near the last operating point for actual abc traces.
    mask = time >= max(0, time[-1] - 0.04)
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True, layout="constrained")
    for ax, prefix, unit, divisor in zip(axes, ("ac.rectifier_terminal.v_abc", "ac.inverter_terminal.v_abc", "ac.right_grid.i_abc"), ("kV", "kV", "kA"), (1e3, 1e3, 1e3)):
        for phase, color in zip(range(3), ("#1864ab", "#2b8a3e", "#c92a2a")):
            key = f"{prefix}_{phase}"
            ax.plot(time[mask], a[key][mask] / divisor, color=color, label=f"SSN Phase {'ABC'[phase]}")
            ax.plot(time[mask], b[key][mask] / divisor, color=color, ls="--", label=f"MNA Phase {'ABC'[phase]}")
        ax.set(title=prefix, ylabel=unit)
    axes[0].legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("Time [s]")
    save_figure(fig, output, "ac_waveforms")

    metrics = []
    duration = time[-1] - time[0]
    if duration <= 0: raise ValueError("No simulation interval")
    for station in STATIONS:
        for signal, label, unit, divisor, base in SIGNALS:
            key = f"{station}.{signal}"
            error = b[key] - a[key]
            maximum = float(np.max(np.abs(error)))
            rms = float(np.sqrt(np.trapezoid(error * error, time) / duration))
            metrics.append(dict(signal=key, unit=unit, max_error=maximum / divisor,
                                rms_error=rms / divisor, base_si=base,
                                max_error_percent=100 * maximum / base,
                                rms_error_percent=100 * rms / base,
                                initial_error=float(error[0]) / divisor))
    worst = max(metrics, key=lambda m: m["max_error_percent"])
    passed = (None if args.max_error_percent is None else
              worst["max_error_percent"] <= args.max_error_percent)
    summary = dict(reference=str(args.reference.resolve()), native=str(args.native.resolve()),
                   reference_dt=args.reference_dt, mna_dt=args.mna_dt,
                   physical_end_time=float(time[-1]), samples=len(time), interpolated=interpolated,
                   error_definition="MNA minus SSN; percent of documented engineering base",
                   acceptance_limit_percent=args.max_error_percent, passed=passed,
                   metrics=metrics)
    (output / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=metrics[0].keys())
        writer.writeheader(); writer.writerows(metrics)
    lines = ["# MMC comparison", "", f"Physical time interval: 0–{time[-1]:.6f} s. "
             f"Time steps: SSN {args.reference_dt*1e6:g} µs, MNA {args.mna_dt*1e6:g} µs.", "",
             "Percentages use fixed engineering bases (1 GW, 1 Gvar, 640 kV, "
             "1562.5 A DC, 60 MJ, and derived current bases), rather than instantaneous values near zero.", "",
             "| Signal | Maximum error | RMS error | Maximum / base |",
             "|---|---:|---:|---:|"]
    for m in metrics:
        lines.append(f"| {m['signal']} | {m['max_error']:.6g} {m['unit']} | {m['rms_error']:.6g} {m['unit']} | {m['max_error_percent']:.5g} % |")
    if passed is not None:
        lines += ["", f"Acceptance limit: {args.max_error_percent:g} % of the corresponding engineering base. "
                  f"Result: {'passed' if passed else 'exceeded'}."]
    lines += ["", "![Inverter](inverter_overview.png)", "", "![Rectifier](rectifier_overview.png)",
              "", "![Differences](differences.png)", "", "![AC waveforms](ac_waveforms.png)", ""]
    (output / "report.md").write_text("\n".join(lines))
    print(f"Plots and metrics: {output}")
    print(f"Largest normalized error: {worst['signal']} = {worst['max_error_percent']:.6g}%")
    if args.max_error_percent is not None and worst["max_error_percent"] > args.max_error_percent:
        raise SystemExit(f"Comparison exceeds requested {args.max_error_percent:g}% limit")


if __name__ == "__main__":
    main()
