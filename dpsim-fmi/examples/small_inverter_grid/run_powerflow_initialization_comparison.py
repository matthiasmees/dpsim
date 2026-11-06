#!/usr/bin/env python3
"""Compare fresh DPsim FMI runs with XML and network powerflow FMU starts."""
from __future__ import annotations

import csv
import json
import math
import zipfile
import subprocess
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
FMU = CASE / "sources" / "IEEE39AvgInverter.fmu"
RUNNER = BUILD / "dpsim-fmi" / "FMI_SmallGrid_Inverter"
OUT = BUILD / "small-inverter-grid" / "results" / "powerflow_initialization"
DT = 25e-6
LOAD_STEP = 2.5


def rebuild_fmu() -> None:
    """Compile the example FMU from its checked-in C++ and XML sources."""
    package = BUILD / "small-inverter-grid" / "fmu-stage"
    binary = package / "binaries" / "x86_64-linux" / "IEEE39AvgInverter.so"
    package.mkdir(parents=True, exist_ok=True)
    binary.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "c++",
            "-std=c++17",
            "-O3",
            "-fPIC",
            "-shared",
            "-I",
            str(ROOT / "dpsim-fmi" / "third_party" / "fmi3"),
            str(CASE / "sources" / "IEEE39AvgInverter.cpp"),
            "-o",
            str(binary),
        ],
        check=True,
    )
    temporary_fmu = package / "IEEE39AvgInverter.fmu"
    with zipfile.ZipFile(
        temporary_fmu, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        archive.write(CASE / "sources" / "modelDescription.xml", "modelDescription.xml")
        archive.write(binary, "binaries/x86_64-linux/IEEE39AvgInverter.so")
    temporary_fmu.replace(FMU)
    print(f"Rebuilt {FMU} from the example C++ source", flush=True)


def run(mode: str, init_mode: str, stop: float, label: str) -> dict:
    case_dir = OUT / label / mode / init_mode
    case_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    load_step = LOAD_STEP if stop > LOAD_STEP else stop / 2
    proc = subprocess.run(
        [str(RUNNER), mode, str(FMU), str(stop), str(DT), str(load_step), init_mode],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    process_wall = time.perf_counter() - t0
    (case_dir / "runner.log").write_text(proc.stdout)
    timing_line = next(
        (line for line in proc.stdout.splitlines() if line.startswith("small_grid ")),
        None,
    )
    if timing_line is None:
        raise RuntimeError(f"No simulation summary in {case_dir / 'runner.log'}")
    timing = dict(item.split("=", 1) for item in timing_line.split() if "=" in item)
    log_name = f"SmallGrid_Inverter_{mode}"
    if init_mode != "manual":
        log_name += f"_{init_mode}"
    source = case_dir / "logs" / log_name / f"{log_name}.csv"
    target = case_dir / f"{mode}_{init_mode}.csv"
    fields = None
    first = None
    last = None
    count = 0
    min_v = [math.inf] * 3
    max_v = [-math.inf] * 3
    min_i = [math.inf] * 3
    max_i = [-math.inf] * 3
    peak_v = 0.0
    peak_i = 0.0
    transient_magnitude = {"v": [], "i": []}
    steady_magnitude = {"v": [], "i": []}
    with source.open(newline="") as fin, target.open("w", newline="") as fout:
        reader = csv.DictReader(fin, skipinitialspace=True)
        fields = [field.strip() for field in reader.fieldnames or []]
        reader.fieldnames = fields
        needed = (
            ["time"]
            + [f"pcc.v_{p}" for p in range(3)]
            + [f"inverter.i_{p}" for p in range(3)]
        )
        if not set(needed).issubset(fields):
            raise RuntimeError(f"Missing expected signals in {source}: {fields}")
        writer = csv.DictWriter(
            fout,
            fieldnames=["time_s"]
            + [f"v_{p}_V" for p in "abc"]
            + [f"i_{p}_A" for p in "abc"],
        )
        writer.writeheader()
        for row in reader:
            values = {key: float(row[key]) for key in needed}
            if not all(math.isfinite(value) for value in values.values()):
                raise RuntimeError(f"Non-finite logged values in {source}")
            if first is None:
                first = values
            last = values
            count += 1
            if values["time"] <= min(stop, 0.1):
                for phase in range(3):
                    v = values[f"pcc.v_{phase}"]
                    i = values[f"inverter.i_{phase}"]
                    min_v[phase] = min(min_v[phase], v)
                    max_v[phase] = max(max_v[phase], v)
                    min_i[phase] = min(min_i[phase], i)
                    max_i[phase] = max(max_i[phase], i)
                    peak_v = max(peak_v, abs(v))
                    peak_i = max(peak_i, abs(i))
                vmag = math.sqrt(
                    2.0 / 3.0 * sum(values[f"pcc.v_{phase}"] ** 2 for phase in range(3))
                )
                imag = math.sqrt(
                    2.0
                    / 3.0
                    * sum(values[f"inverter.i_{phase}"] ** 2 for phase in range(3))
                )
                transient_magnitude["v"].append(vmag)
                transient_magnitude["i"].append(imag)
            if 0.5 <= values["time"] < min(stop, LOAD_STEP) - 0.1:
                steady_magnitude["v"].append(
                    math.sqrt(
                        2.0
                        / 3.0
                        * sum(values[f"pcc.v_{phase}"] ** 2 for phase in range(3))
                    )
                )
                steady_magnitude["i"].append(
                    math.sqrt(
                        2.0
                        / 3.0
                        * sum(values[f"inverter.i_{phase}"] ** 2 for phase in range(3))
                    )
                )
            writer.writerow(
                {
                    "time_s": values["time"],
                    **{f"v_{p}_V": values[f"pcc.v_{k}"] for k, p in enumerate("abc")},
                    **{
                        f"i_{p}_A": values[f"inverter.i_{k}"]
                        for k, p in enumerate("abc")
                    },
                }
            )
    expected = round(stop / DT) + 1
    if count != expected:
        raise RuntimeError(f"{source} has {count} samples; expected {expected}")
    result = {
        "mode": mode,
        "initialization": init_mode,
        "stop_s": stop,
        "dt_s": DT,
        "load_step_s": load_step,
        "samples": count,
        "process_wall_s": process_wall,
        "initial_pcc_abc_V": [first[f"pcc.v_{p}"] for p in range(3)],
        "initial_current_abc_A": [first[f"inverter.i_{p}"] for p in range(3)],
        "first_100ms_pcc_phase_range_V": [max_v[k] - min_v[k] for k in range(3)],
        "first_100ms_current_phase_range_A": [max_i[k] - min_i[k] for k in range(3)],
        "first_100ms_peak_abs_voltage_V": peak_v,
        "first_100ms_peak_abs_current_A": peak_i,
        "initial_100ms_magnitude_rms_deviation_from_pre_event_mean_percent": {},
        "trace": str(target),
        **{key: float(value) for key, value in timing.items() if key.endswith("_s")},
    }
    for signal, values in transient_magnitude.items():
        steady = steady_magnitude[signal]
        if steady:
            reference = sum(steady) / len(steady)
            result["initial_100ms_magnitude_rms_deviation_from_pre_event_mean_percent"][
                signal
            ] = (
                100.0
                * math.sqrt(
                    sum((value - reference) ** 2 for value in values) / len(values)
                )
                / reference
            )
        else:
            result["initial_100ms_magnitude_rms_deviation_from_pre_event_mean_percent"][
                signal
            ] = None
    (case_dir / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"{mode.upper():2} {init_mode:23} {stop:g}s: "
        f"wall={result['wall_s']:.3f}s init={result['init_s']:.3f}s "
        f"samples={count} V0={result['initial_pcc_abc_V']} "
        f"I0={result['initial_current_abc_A']} "
        f"startup-deviation={result['initial_100ms_magnitude_rms_deviation_from_pre_event_mean_percent']}",
        flush=True,
    )
    return result


def trace_path(mode: str, init_mode: str) -> Path:
    return OUT / "five_second" / mode / init_mode / f"{mode}_{init_mode}.csv"


def compare_pair(mode: str) -> dict:
    unadjusted = trace_path(mode, "powerflow-unadjusted")
    initialized = trace_path(mode, "powerflow")
    windows = {
        "startup_0_to_100ms": (0.0, 0.1),
        "pre_event_0p5_to_2p4s": (0.5, 2.4),
        "post_event_2p5_to_5s": (2.5, 5.0),
    }
    accum = {
        name: {kind: [0.0, 0.0] for kind in ("voltage", "current")} for name in windows
    }
    with unadjusted.open(newline="") as a, initialized.open(newline="") as b:
        ra, rb = csv.DictReader(a), csv.DictReader(b)
        for row_a, row_b in zip(ra, rb):
            t = float(row_a["time_s"])
            if abs(t - float(row_b["time_s"])) > 1e-12:
                raise RuntimeError(f"Time grids differ for {mode} at {t}")
            window = next(
                (name for name, (lo, hi) in windows.items() if lo <= t < hi), None
            )
            if window is None:
                continue
            for kind, suffix in (("voltage", "V"), ("current", "A")):
                for phase in "abc":
                    column = f"{'v' if kind == 'voltage' else 'i'}_{phase}_{suffix}"
                    reference = float(row_a[column])
                    error = float(row_b[column]) - reference
                    accum[window][kind][0] += error * error
                    accum[window][kind][1] += reference * reference
    output = {}
    for window, values in accum.items():
        output[window] = {}
        for kind, (error_sum, ref_sum) in values.items():
            error_rms = math.sqrt(
                error_sum / max(1, 3 * (windows[window][1] - windows[window][0]) / DT)
            )
            ref_rms = math.sqrt(
                ref_sum / max(1, 3 * (windows[window][1] - windows[window][0]) / DT)
            )
            output[window][kind] = {
                "rms_difference": error_rms,
                "rms_difference_percent_of_unadjusted_rms": (
                    100.0 * error_rms / ref_rms if ref_rms else None
                ),
            }
    return output


def plot_comparison() -> Path:
    colors = {"me": "#2455a4", "cs": "#d55e00"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex="col")
    axes[0, 0].set_title("PCC three-phase voltage magnitude")
    axes[0, 1].set_title("Inverter phase-a current")
    for mode in ("me", "cs"):
        for init_mode, style, label in (
            ("powerflow-unadjusted", "--", "FMU XML state"),
            ("powerflow", "-", "FMU reinitialized from PF voltage"),
        ):
            t_values, v_values, i_values = [], [], []
            with trace_path(mode, init_mode).open(newline="") as stream:
                for index, row in enumerate(csv.DictReader(stream)):
                    if index % 20:
                        continue
                    t = float(row["time_s"])
                    v = math.sqrt(
                        2 / 3 * sum(float(row[f"v_{p}_V"]) ** 2 for p in "abc")
                    )
                    i = float(row["i_a_A"])
                    t_values.append(t)
                    v_values.append(v)
                    i_values.append(i)
            suffix = "ME" if mode == "me" else "CS"
            name = f"DPsim {suffix}: {label}"
            axes[0, 0].plot(
                t_values,
                v_values,
                linestyle=style,
                color=colors[mode],
                alpha=0.95 if init_mode == "powerflow" else 0.7,
                label=name,
            )
            axes[0, 1].plot(
                t_values,
                i_values,
                linestyle=style,
                color=colors[mode],
                alpha=0.95 if init_mode == "powerflow" else 0.7,
                label=name,
            )
            keep = [index for index, time_s in enumerate(t_values) if time_s <= 0.1]
            axes[1, 0].plot(
                [t_values[k] * 1000 for k in keep],
                [v_values[k] for k in keep],
                linestyle=style,
                color=colors[mode],
                alpha=0.95 if init_mode == "powerflow" else 0.7,
                label=name,
            )
            axes[1, 1].plot(
                [t_values[k] * 1000 for k in keep],
                [i_values[k] for k in keep],
                linestyle=style,
                color=colors[mode],
                alpha=0.95 if init_mode == "powerflow" else 0.7,
                label=name,
            )
    axes[1, 0].set_title("Startup zoom (0–100 ms)")
    axes[1, 1].set_title("Startup zoom (0–100 ms)")
    for ax in axes[0]:
        ax.axvline(
            LOAD_STEP, color="0.35", linestyle=":", linewidth=1, label="50 MW load step"
        )
        ax.set_xlim(0, 5)
    for ax in axes[1]:
        ax.set_xlim(0, 100)
        ax.set_xlabel("Time [ms]")
    axes[0, 0].set_ylabel("Voltage magnitude [V]")
    axes[1, 0].set_ylabel("Voltage magnitude [V]")
    axes[0, 1].set_ylabel("Phase-a current [A]")
    axes[1, 1].set_ylabel("Phase-a current [A]")
    for row in axes:
        for ax in row:
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=7, loc="best")
    fig.suptitle("Power-flow-based FMU initialization: DPsim ME and CS")
    fig.tight_layout()
    path = OUT / "powerflow_initialization_comparison.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    subprocess.run(
        ["cmake", "--build", str(BUILD), "--target", "FMI_SmallGrid_Inverter", "-j1"],
        check=True,
    )
    rebuild_fmu()
    OUT.mkdir(parents=True, exist_ok=True)
    # Always validate the four configurations on a brief run before the 5 s case.
    results = {"short": [], "five_second": []}
    for mode in ("me", "cs"):
        for init_mode in ("powerflow-unadjusted", "powerflow"):
            results["short"].append(run(mode, init_mode, 0.05, "short"))
    for mode in ("me", "cs"):
        for init_mode in ("powerflow-unadjusted", "powerflow"):
            results["five_second"].append(run(mode, init_mode, 5.0, "five_second"))
    results["feature_vs_unadjusted"] = {
        mode: compare_pair(mode) for mode in ("me", "cs")
    }
    results["plot"] = str(plot_comparison())
    report = OUT / "comparison_metrics.json"
    report.write_text(json.dumps(results, indent=2) + "\n")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
