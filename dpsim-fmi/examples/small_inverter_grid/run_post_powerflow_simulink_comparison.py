#!/usr/bin/env python3
"""Compare post-PF DPsim ME/CS against fresh Simulink native/CS sweeps."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import zipfile
from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
POST = BUILD / "small-inverter-grid" / "results" / "post_powerflow_fix"
OUT = POST / "simulink_comparison"
DP_METRICS = POST / "metrics.json"
STEPS_US = (25, 50, 100, 500, 1000)
STOP_S = 5.0
LOAD_STEP_S = 2.5
VISUAL_LIMIT = 6000
METHODS = ("DPsim ME", "DPsim CS", "Simulink CS")
COLORS = {
    "Native Simulink": "#202020",
    "DPsim ME": "#2463a6",
    "DPsim CS": "#d55e00",
    "Simulink CS": "#009e73",
}


def one_bus_powerflow() -> dict:
    """Solve the balanced fundamental-frequency PF without DPsim outputs."""
    vll = 362215.5
    frequency = 60.0
    omega = 2 * math.pi * frequency
    resistance, inductance, capacitance = 0.30, 0.02, 1e-6
    tau, active, reactive = 5e-3, 125e6, 80.881e6
    omega_tau = omega * tau
    denominator = 1 + omega_tau**2
    p_periodic = (active - omega_tau * reactive) / denominator
    q_periodic = (reactive + omega_tau * active) / denominator
    z_line = complex(resistance, omega * inductance)
    v_phase_nominal = vll / math.sqrt(3)

    def electrical(magnitude: float, angle: float):
        vpcc = magnitude * complex(math.cos(angle), math.sin(angle))
        i_inverter = -((complex(p_periodic, q_periodic) / (3 * vpcc)).conjugate())
        i_series = i_inverter + 1j * omega * (capacitance / 2) * vpcc
        v_source = vpcc + z_line * i_series
        return vpcc, i_inverter, i_series, v_source

    x = [v_phase_nominal, 0.0]
    for _ in range(12):
        _, _, _, source = electrical(*x)
        residual = np.array([source.real - v_phase_nominal, source.imag], dtype=float)
        if np.linalg.norm(residual, ord=np.inf) < 1e-9:
            break
        jacobian = np.empty((2, 2))
        for column, delta in ((0, 0.01), (1, 1e-7)):
            xp = x.copy()
            xp[column] += delta
            _, _, _, source_p = electrical(*xp)
            rp = np.array([source_p.real - v_phase_nominal, source_p.imag])
            jacobian[:, column] = (rp - residual) / delta
        correction = np.linalg.solve(jacobian, residual)
        x = [x[k] - correction[k] for k in range(2)]
    vpcc, i_inverter, i_series, v_source = electrical(*x)
    if abs(v_source - complex(v_phase_nominal, 0)) > 1e-5:
        raise RuntimeError(
            f"Independent one-bus PF failed to converge: residual={v_source}"
        )

    def abc_peak(phasor: complex) -> list[float]:
        return [
            math.sqrt(2)
            * (
                phasor
                * complex(
                    math.cos(-k * 2 * math.pi / 3), math.sin(-k * 2 * math.pi / 3)
                )
            ).real
            for k in range(3)
        ]

    source_phase = [
        math.atan2(
            (
                v_source
                * complex(
                    math.cos(-k * 2 * math.pi / 3), math.sin(-k * 2 * math.pi / 3)
                )
            ).imag,
            (
                v_source
                * complex(
                    math.cos(-k * 2 * math.pi / 3), math.sin(-k * 2 * math.pi / 3)
                )
            ).real,
        )
        for k in range(3)
    ]
    return {
        "method": "independent balanced one-bus AC power flow solved from source, PI-line, and periodic inverter equations",
        "uses_dpsim_initial_values": False,
        "frequency_hz": frequency,
        "vll_source_rms_V": vll,
        "periodic_inverter_P_W": p_periodic,
        "periodic_inverter_Q_var": q_periodic,
        "pcc_phase_a_phasor_rms_V": [vpcc.real, vpcc.imag],
        "source_phase_a_phasor_rms_V": [v_source.real, v_source.imag],
        "v_pcc_abc_V": abc_peak(vpcc),
        "i_inverter_abc_A": abc_peak(i_inverter),
        "i_line_abc_A": abc_peak(i_series),
        "source_phase_peak_V": math.sqrt(2) * abs(v_source),
        "source_phase_rad": source_phase,
        "source_voltage_residual_V": abs(v_source - complex(v_phase_nominal, 0)),
        "line_R_ohm": resistance,
        "line_L_H": inductance,
        "line_C_total_F": capacitance,
    }


def make_initial_condition_fmu(op: dict) -> Path:
    """Build a separate CS FMU whose starts match the independent PF point."""
    stage = BUILD / "small-inverter-grid" / "post_powerflow_simulink_fmu"
    binary = stage / "binaries" / "x86_64-linux" / "IEEE39AvgInverter.so"
    stage.mkdir(parents=True, exist_ok=True)
    binary.parent.mkdir(parents=True, exist_ok=True)
    source = (CASE / "sources" / "IEEE39AvgInverter.cpp").read_text()
    default = "model->voltage = {288466.352188, -200716.863664, -87749.488524};"
    replacement = (
        "model->voltage = {" + ", ".join(f"{x:.15g}" for x in op["v_pcc_abc_V"]) + "};"
    )
    if default not in source:
        raise RuntimeError(
            "FMU Reset voltage start changed; review specialized Simulink start generation."
        )
    source = source.replace(default, replacement, 1)
    source_file = stage / "IEEE39AvgInverter.cpp"
    source_file.write_text(source)
    subprocess.run(
        [
            "c++",
            "-std=c++17",
            "-O3",
            "-fPIC",
            "-shared",
            "-I",
            str(ROOT / "dpsim-fmi" / "third_party" / "fmi3"),
            str(source_file),
            "-o",
            str(binary),
        ],
        check=True,
    )

    xml = (CASE / "sources" / "modelDescription.xml").read_text()
    for phase, voltage, current in zip(
        "abc", op["v_pcc_abc_V"], op["i_inverter_abc_A"]
    ):
        xml, n_v = re.subn(
            rf'(<Float64 name="v_{phase}"[^>]*start=")[^"]+("[^>]*/>)',
            rf"\g<1>{voltage:.15g}\g<2>",
            xml,
            count=1,
        )
        xml, n_i = re.subn(
            rf'(<Float64 name="x_{phase}"[^>]*start=")[^"]+("[^>]*/>)',
            rf"\g<1>{current:.15g}\g<2>",
            xml,
            count=1,
        )
        if n_v != 1 or n_i != 1:
            raise RuntimeError(f"Could not update FMU XML start for phase {phase}.")
    fmu = stage / "IEEE39AvgInverter.fmu"
    with zipfile.ZipFile(fmu, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("modelDescription.xml", xml)
        archive.write(binary, "binaries/x86_64-linux/IEEE39AvgInverter.so")
    return fmu


def run_matlab_sweep(op: dict, fmu: Path, full_sweep: bool = True) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    op_path = OUT / "independent_operating_point.json"
    op_path.write_text(json.dumps(op, indent=2) + "\n")
    matlab_script = CASE / "matlab" / "run_post_powerflow_simulink_sweep.m"
    args = [str(ROOT), str(OUT / "simulink"), str(fmu), str(op_path)]
    # This workspace has no apostrophes in its paths; MATLAB string quoting is safe.
    call = ",".join("'" + item.replace("\\", "/") + "'" for item in args)
    call += ",true" if full_sweep else ",false"
    cmd = [
        "matlab",
        "-batch",
        f"addpath('{matlab_script.parent.as_posix()}'); run_post_powerflow_simulink_sweep({call});",
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def read_rows(path: Path):
    stream = path.open(newline="")
    reader = csv.DictReader(stream, skipinitialspace=True)
    reader.fieldnames = [field.strip() for field in (reader.fieldnames or [])]
    return stream, reader


def paths_for(step_us: int) -> dict[str, Path]:
    dp = POST / "full"
    sl = OUT / "simulink"
    return {
        "DPsim ME": dp
        / "me"
        / f"dt_{step_us}us"
        / "logs"
        / "SmallGrid_Inverter_me_powerflow"
        / "SmallGrid_Inverter_me_powerflow.csv",
        "DPsim CS": dp
        / "cs"
        / f"dt_{step_us}us"
        / "logs"
        / "SmallGrid_Inverter_cs_powerflow"
        / "SmallGrid_Inverter_cs_powerflow.csv",
        "Native Simulink": sl / f"dt_{step_us}us" / "native" / "simulink_native.csv",
        "Simulink CS": sl / f"dt_{step_us}us" / "cs" / "simulink_cs.csv",
    }


def get_signal(row: dict, method: str, signal: str) -> list[float]:
    if method.startswith("DPsim"):
        if signal == "v":
            return [float(row[f"pcc.v_{k}"]) for k in range(3)]
        return [float(row[f"inverter.i_{k}"]) for k in range(3)]
    prefix = "v" if signal == "v" else "i"
    return [
        float(row[f"{prefix}_{phase}_{'V' if signal == 'v' else 'A'}"])
        for phase in "abc"
    ]


def analyze_step(step_us: int, output_dir: Path) -> dict:
    paths = paths_for(step_us)
    for method, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {method} run at {step_us} us: {path}")
    streams = {}
    readers = {}
    for method, path in paths.items():
        streams[method], readers[method] = read_rows(path)
    dt = step_us * 1e-6
    window_n = max(1, round(0.1 / dt))
    queues = {method: {s: deque() for s in ("v", "i")} for method in METHODS}
    rolling_sums = {method: {s: [0.0, 0.0] for s in ("v", "i")} for method in METHODS}
    accum = {method: {s: [0.0, 0.0, 0.0] for s in ("v", "i")} for method in METHODS}
    count = 0
    compact_stride = max(1, math.ceil((round(STOP_S / dt) + 1) / VISUAL_LIMIT))
    compact_path = output_dir / f"comparison_{step_us}us.csv"
    fields = ["time_s", "native_vmag_kV", "native_i_a_A"]
    for method in METHODS:
        key = method.lower().replace(" ", "_")
        fields += [
            f"{key}_vmag_kV",
            f"{key}_i_a_A",
            f"{key}_v_error_100ms_percent",
            f"{key}_i_error_100ms_percent",
        ]
    initial = {}
    with compact_path.open("w", newline="") as compact_stream:
        writer = csv.DictWriter(compact_stream, fieldnames=fields)
        writer.writeheader()
        for index in range(round(STOP_S / dt) + 1):
            rows = {}
            for method in paths:
                try:
                    rows[method] = next(readers[method])
                except StopIteration as error:
                    raise RuntimeError(f"Short trace: {paths[method]}") from error
            times = {
                method: float(row.get("time", row.get("time_s")))
                for method, row in rows.items()
            }
            t = times["Native Simulink"]
            if any(abs(value - t) > 1e-12 for value in times.values()):
                raise RuntimeError(
                    f"Time-grid mismatch at {step_us} us, sample {index}: {times}"
                )
            signals = {
                method: {s: get_signal(rows[method], method, s) for s in ("v", "i")}
                for method in paths
            }
            if not all(
                math.isfinite(x)
                for item in signals.values()
                for arr in item.values()
                for x in arr
            ):
                raise RuntimeError(f"Non-finite output at {step_us} us, sample {index}")
            if index == 0:
                initial = {method: signals[method] for method in paths}
            native = signals["Native Simulink"]
            error_percent_now = {}
            for method in METHODS:
                error_percent_now[method] = {}
                for signal in ("v", "i"):
                    err = [
                        signals[method][signal][k] - native[signal][k] for k in range(3)
                    ]
                    q = queues[method][signal]
                    sums = rolling_sums[method][signal]
                    ref2 = sum(x * x for x in native[signal])
                    err2 = sum(x * x for x in err)
                    q.append((ref2, err2))
                    sums[0] += ref2
                    sums[1] += err2
                    if len(q) > window_n:
                        old = q.popleft()
                        sums[0] -= old[0]
                        sums[1] -= old[1]
                    if len(q) == window_n and sums[0] > 0:
                        error_percent_now[method][signal] = 100 * math.sqrt(
                            sums[1] / sums[0]
                        )
                    else:
                        error_percent_now[method][signal] = math.nan
                    accum[method][signal][0] += err2
                    accum[method][signal][1] += ref2
                    accum[method][signal][2] = max(
                        accum[method][signal][2], math.sqrt(err2 / 3)
                    )
            if index % compact_stride == 0 or index == round(STOP_S / dt):
                row = {
                    "time_s": t,
                    "native_vmag_kV": math.sqrt(2 / 3 * sum(x * x for x in native["v"]))
                    / 1000,
                    "native_i_a_A": native["i"][0],
                }
                for method in METHODS:
                    key = method.lower().replace(" ", "_")
                    row[f"{key}_vmag_kV"] = (
                        math.sqrt(2 / 3 * sum(x * x for x in signals[method]["v"]))
                        / 1000
                    )
                    row[f"{key}_i_a_A"] = signals[method]["i"][0]
                    row[f"{key}_v_error_100ms_percent"] = error_percent_now[method]["v"]
                    row[f"{key}_i_error_100ms_percent"] = error_percent_now[method]["i"]
                writer.writerow(row)
            count += 1
    for stream in streams.values():
        stream.close()
    expected = round(STOP_S / dt) + 1
    if count != expected:
        raise RuntimeError(f"Expected {expected} aligned samples, got {count}")
    metrics = {}
    for method in METHODS:
        method_metrics = {}
        for signal in ("v", "i"):
            err2, ref2, max_error = accum[method][signal]
            denominator = max(1, count * 3)
            err_rms = math.sqrt(err2 / denominator)
            ref_rms = math.sqrt(ref2 / denominator)
            method_metrics["voltage" if signal == "v" else "current"] = {
                "rms_error": err_rms,
                "rms_error_percent_native_rms": (
                    100 * err_rms / ref_rms if ref_rms else None
                ),
                "max_instantaneous_phase_error": max_error,
            }
        metrics[method] = method_metrics
    metrics["Native Simulink"] = {
        s: {
            "rms_error": 0,
            "rms_error_percent_native_rms": 0,
            "max_instantaneous_phase_error": 0,
        }
        for s in ("voltage", "current")
    }
    initial_metrics = {}
    for method, signals in initial.items():
        initial_metrics[method] = {"v_abc_V": signals["v"], "i_abc_A": signals["i"]}
    return {
        "step_us": step_us,
        "samples": count,
        "compact_trace": str(compact_path),
        "rms_error_vs_native_simulink": metrics,
        "initial_values": initial_metrics,
    }


def add_runtimes(summary: dict) -> None:
    dp_metrics = json.loads(DP_METRICS.read_text())
    run_by_key = {
        (entry["mode"], entry["step_us"]): entry
        for entry in dp_metrics["runs"]
        if entry["stop_s"] == STOP_S
    }
    for entry in summary["steps"]:
        step = entry["step_us"]
        entry["wall_time_s"] = {
            "DPsim ME": run_by_key[("ME", step)]["simulation_wall_s"],
            "DPsim CS": run_by_key[("CS", step)]["simulation_wall_s"],
        }
        for mode in ("native", "cs"):
            timing_path = (
                OUT
                / "simulink"
                / f"dt_{step}us"
                / mode
                / f"simulink_{mode}_timing.json"
            )
            timing = json.loads(timing_path.read_text())
            label = "Native Simulink" if mode == "native" else "Simulink CS"
            entry["wall_time_s"][label] = timing["simulink_wall_time_s"]


def create_figures(summary: dict) -> list[Path]:
    figure_dir = OUT / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    produced = []
    for item in summary["steps"]:
        step = item["step_us"]
        rows = list(csv.DictReader(Path(item["compact_trace"]).open(newline="")))
        t = np.array([float(r["time_s"]) for r in rows])
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex="col")
        axes[0, 0].set_title("PCC three-phase voltage magnitude")
        axes[0, 1].set_title("Inverter phase-a current")
        axes[1, 0].set_title("100 ms rolling voltage RMS error vs native Simulink")
        axes[1, 1].set_title("100 ms rolling current RMS error vs native Simulink")
        axes[0, 0].plot(
            t,
            [float(r["native_vmag_kV"]) for r in rows],
            color=COLORS["Native Simulink"],
            label="Native Simulink",
            lw=1.5,
        )
        axes[0, 1].plot(
            t,
            [float(r["native_i_a_A"]) for r in rows],
            color=COLORS["Native Simulink"],
            label="Native Simulink",
            lw=1.5,
        )
        for method in METHODS:
            key = method.lower().replace(" ", "_")
            axes[0, 0].plot(
                t,
                [float(r[f"{key}_vmag_kV"]) for r in rows],
                color=COLORS[method],
                label=method,
                lw=1,
            )
            axes[0, 1].plot(
                t,
                [float(r[f"{key}_i_a_A"]) for r in rows],
                color=COLORS[method],
                label=method,
                lw=1,
            )
            axes[1, 0].plot(
                t,
                [float(r[f"{key}_v_error_100ms_percent"]) for r in rows],
                color=COLORS[method],
                label=method,
                lw=1,
            )
            axes[1, 1].plot(
                t,
                [float(r[f"{key}_i_error_100ms_percent"]) for r in rows],
                color=COLORS[method],
                label=method,
                lw=1,
            )
        for ax in axes.flat:
            ax.axvline(
                LOAD_STEP_S, color="#555555", ls=":", lw=1, label="50 MW load step"
            )
            ax.grid(True, alpha=0.3)
        axes[0, 0].set_ylabel("Voltage [kV]")
        axes[0, 1].set_ylabel("Current [A]")
        axes[1, 0].set_ylabel("RMS error [%]")
        axes[1, 1].set_ylabel("RMS error [%]")
        axes[1, 0].set_xlabel("Time [s]")
        axes[1, 1].set_xlabel("Time [s]")
        axes[0, 0].legend(fontsize=8, ncol=2)
        axes[0, 1].legend(fontsize=8, ncol=2)
        for ax in axes[1, :]:
            ax.legend(fontsize=8, ncol=2)
        fig.suptitle(f"Post-PF small grid: four closed-loop methods, {step} µs")
        fig.tight_layout()
        path = figure_dir / f"comparison_{step}us.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        produced.append(path)

    x = np.arange(len(summary["steps"]))
    width = 0.19
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for j, method in enumerate(METHODS):
        color = COLORS[method]
        for axis, signal, offset in ((axes[0], "voltage", j), (axes[1], "current", j)):
            vals = [
                s["rms_error_vs_native_simulink"][method][signal][
                    "rms_error_percent_native_rms"
                ]
                for s in summary["steps"]
            ]
            axis.bar(x + (offset - 1) * width, vals, width, label=method, color=color)
    axes[0].set_title("Three-phase PCC voltage RMS error")
    axes[1].set_title("Three-phase inverter current RMS error")
    for ax in axes:
        ax.set_xticks(x, [f"{s['step_us']} µs" for s in summary["steps"]])
        ax.set_ylabel("Error vs native Simulink [%]")
        ax.grid(axis="y", alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = figure_dir / "rms_error_by_timestep.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    produced.append(path)

    fig, ax = plt.subplots(figsize=(10, 5))
    for method in ("DPsim ME", "DPsim CS", "Native Simulink", "Simulink CS"):
        ax.plot(
            [s["step_us"] for s in summary["steps"]],
            [s["wall_time_s"][method] for s in summary["steps"]],
            marker="o",
            label=method,
            color=COLORS[method],
        )
    ax.set_xscale("log")
    ax.set_xticks(
        [s["step_us"] for s in summary["steps"]],
        [str(s["step_us"]) for s in summary["steps"]],
    )
    ax.set_xlabel("Fixed step [µs]")
    ax.set_ylabel("Whole-run wall time [s]")
    ax.set_title("5 s runtime by timestep")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    path = figure_dir / "runtime_by_timestep.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    produced.append(path)
    return produced


def markdown_table(summary: dict) -> str:
    lines = [
        "| Step | Samples | DPsim ME V / I [%] | DPsim CS V / I [%] | Simulink CS V / I [%] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for step in summary["steps"]:
        values = []
        for method in METHODS:
            m = step["rms_error_vs_native_simulink"][method]
            values.append(
                f"{m['voltage']['rms_error_percent_native_rms']:.5g} / {m['current']['rms_error_percent_native_rms']:.5g}"
            )
        lines.append(
            f"| {step['step_us']} µs | {step['samples']:,} | {values[0]} | {values[1]} | {values[2]} |"
        )
    return "\n".join(lines)


def create_notebook(summary: dict, figures: list[Path]) -> Path:
    relfig = {p.name: p for p in figures}
    cells = []

    def md(text):
        cells.append(
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": text.splitlines(keepends=True),
            }
        )

    def image_cell(path: Path):
        import base64

        payload = base64.b64encode(path.read_bytes()).decode("ascii")
        cells.append(
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [
                    {
                        "data": {
                            "image/png": payload,
                            "text/plain": [f"<Figure: {path.name}>\n"],
                        },
                        "metadata": {},
                        "output_type": "display_data",
                    }
                ],
                "source": [
                    "from IPython.display import Image, display\n",
                    f"display(Image(filename=str(figures / '{path.name}')))\n",
                ],
            }
        )

    md(
        "# Small grid: post-PF DPsim and Simulink comparison\n\nThis notebook compares fresh closed-loop DPsim ME, DPsim CS, native Simulink, and Simulink CS runs after an independently calculated balanced AC power-flow initialization. Every case uses a 5 s run, the same timestep within each comparison, and a 50 MW PCC load step at 2.5 s.\n\nSimulink initialization was solved independently from the one-bus AC equations for the source, PI feeder, shunt capacitance, and the inverter's 60 Hz periodic current. No DPsim initial values or recorded traces were used to initialize Simulink. The 25 µs ME run is not used as the reference here: native Simulink is the error reference for each matched timestep.\n"
    )
    cells.append(
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "from pathlib import Path\n",
                "roots = [Path.cwd(), *Path.cwd().parents]\n",
                "repo = next(root for root in roots if (root / 'build-fmi/small-inverter-grid/results/post_powerflow_fix/simulink_comparison/figures').is_dir())\n",
                "figures = repo / 'build-fmi/small-inverter-grid/results/post_powerflow_fix/simulink_comparison/figures'\n",
                "print(f'Using compact saved plots in {figures}')\n",
            ],
        }
    )
    init_lines = [
        "| Step | Method | Initial voltage error [%] | Initial current error [%] |",
        "|---:|---|---:|---:|",
    ]
    for step in summary["steps"]:
        for method, checks in summary["initial_point_check"][
            str(step["step_us"])
        ].items():
            init_lines.append(
                f"| {step['step_us']} µs | {method} | {checks['voltage_initial_error_percent']:.6g} | {checks['current_initial_error_percent']:.6g} |"
            )
    md(
        "## Independent startup operating point\n\nThe balanced AC power flow was solved separately from the DPsim outputs. The periodic inverter injection accounts for the 5 ms current response at 60 Hz. Native Simulink starts its line, PCC, and inverter current states at that operating point; the Simulink CS copy of the FMU carries matching voltage and periodic-current start values. The DPsim runs use their own power-flow initialization. The table checks initial abc values against the independently computed AC point; small numerical differences are expected from the separate solvers.\n\n"
        + "\n".join(init_lines)
        + "\n"
    )
    md(
        "## RMS errors against native Simulink\n\nEach error uses all three phases over the complete 5 s trace and is normalized by the corresponding native-Simulink RMS. `Simulink CS` is included as a method alongside DPsim ME and CS.\n\n"
        + markdown_table(summary)
        + "\n"
    )
    md(
        "## Matched traces and rolling errors\n\nEach timestep has its own plot. The lower panels show 100 ms rolling RMS error percentages against native Simulink; the dotted line marks the 50 MW load step.\n"
    )
    for step in summary["steps"]:
        path = OUT / "figures" / f"comparison_{step['step_us']}us.png"
        md(f"### {step['step_us']} µs\n")
        image_cell(path)
    md(
        "## RMS error and runtime across the sweep\n\nRuntimes are whole-run wall times for the fresh 5 s simulations. Each method uses the same fixed timestep within a row; different methods use their own solver implementations.\n"
    )
    image_cell(relfig["rms_error_by_timestep.png"])
    image_cell(relfig["runtime_by_timestep.png"])
    md(
        "## Reproduce\n\nThe script independently computes the AC startup point, builds a Simulink-CS-specific copy of the FMU with matching start values, runs 50 ms smoke cases first, then the 5 s native and CS sweeps, and regenerates compact metrics and figures:\n\n```bash\npython dpsim-fmi/examples/small_inverter_grid/run_post_powerflow_simulink_comparison.py\n```\n\nThe full-resolution traces, timing records, metrics, and figures are under `build-fmi/small-inverter-grid/results/post_powerflow_fix/simulink_comparison/`. The customized FMU is staged separately under `build-fmi/`; the checked-in FMU sources are not altered by this comparison runner.\n"
    )
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3.x"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = CASE / "SmallGrid_Post_Powerflow_Simulink_Comparison.ipynb"
    path.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n")
    return path


def analyze_only() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    steps = []
    for step in STEPS_US:
        output_dir = OUT / "simulink" / f"dt_{step}us" / "analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        steps.append(analyze_step(step, output_dir))
    op = json.loads((OUT / "independent_operating_point.json").read_text())
    reference_initial = {"v": op["v_pcc_abc_V"], "i": op["i_inverter_abc_A"]}
    startup_check = {}
    for step in steps:
        for method, values in step["initial_values"].items():
            check = {}
            for signal, label in (("v", "voltage"), ("i", "current")):
                actual = values["v_abc_V" if signal == "v" else "i_abc_A"]
                reference = reference_initial[signal]
                rms_ref = math.sqrt(sum(x * x for x in reference) / 3)
                rms_error = math.sqrt(
                    sum((x - y) ** 2 for x, y in zip(actual, reference)) / 3
                )
                check[label + "_initial_error_percent"] = (
                    100 * rms_error / rms_ref if rms_ref else None
                )
            startup_check.setdefault(str(step["step_us"]), {})[method] = check
    summary = {
        "configuration": {
            "stop_s": STOP_S,
            "load_step_s": LOAD_STEP_S,
            "steps_us": list(STEPS_US),
            "methods": ["Native Simulink", *METHODS],
            "error_reference": "Native Simulink at matched timestep",
            "initialization": "Independent one-bus AC power flow for Simulink; DPsim's own EMT PF for DPsim",
        },
        "independent_operating_point": op,
        "initial_point_check": startup_check,
        "steps": steps,
    }
    add_runtimes(summary)
    summary_path = OUT / "comparison_metrics.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "simulink_full_sweep_complete.json").write_text(
        json.dumps(
            {
                "complete": True,
                "steps_us": list(STEPS_US),
                "methods": ["native", "cs"],
                "stop_s": STOP_S,
                "load_step_s": LOAD_STEP_S,
                "samples_per_step": {
                    str(step): round(STOP_S / (step * 1e-6)) + 1 for step in STEPS_US
                },
                "smoke_cases": 2,
                "smoke_stop_s": 0.05,
            },
            indent=2,
        )
        + "\n"
    )
    figures = create_figures(summary)
    notebook = create_notebook(summary, figures)
    print(f"Wrote {summary_path}\nWrote {notebook}")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="Regenerate comparison artifacts from completed Simulink runs",
    )
    parser.add_argument(
        "--smoke-only",
        action="store_true",
        help="Run only the two 50 ms Simulink checks",
    )
    args = parser.parse_args()
    if args.analyze_only:
        analyze_only()
        return
    op = one_bus_powerflow()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "independent_operating_point.json").write_text(
        json.dumps(op, indent=2) + "\n"
    )
    fmu = make_initial_condition_fmu(op)
    run_matlab_sweep(op, fmu, full_sweep=not args.smoke_only)
    if args.smoke_only:
        print("Both short Simulink checks completed; full 5 s sweep was not started.")
        return
    summary = analyze_only()
    print("RMS error summary vs matched native Simulink:")
    print(markdown_table(summary))


if __name__ == "__main__":
    main()
