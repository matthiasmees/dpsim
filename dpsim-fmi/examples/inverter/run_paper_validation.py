#!/usr/bin/env python3
"""Run and summarize bounded FMI/EMT validation cases.

All generated runs, metrics, and figures are written below build-fmi. The 1 us
weak-grid run is a DPsim numerical reference for the same FMU/network. When the
regenerated MATLAB CSVs exist, their matching long converter cases are also
compared as independent physical references.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = pathlib.Path(__file__).resolve().parents[3]
CASE = ROOT / "dpsim-fmi/examples/inverter"
BUILD = ROOT / "build-fmi/fmi-validation/gate6-weak-grid"
EXE = ROOT / "build-fmi/dpsim-fmi/FMI_EMT_Ph3_Mixed"
FMU = CASE / "sources/GridTiedInverterPortedV3.x86_64-linux.fmu"
BASE_CONFIG = CASE / "sources/weak_grid_sag_50us.json"
MATLAB_REFERENCE = (
    ROOT
    / "build-fmi/fmi-validation/matlab-reference-attempt/voltage_sag_reference_trap_1us.csv"
)
NONLINEAR_TEST = ROOT / "build-fmi/dpsim-fmi/tests/FMI_NonlinearTest"
NONLINEAR_FMU = ROOT / "build-fmi/dpsim-fmi/tests/nonlinear_threephase.fmu"
SAMPLED_TEST = ROOT / "build-fmi/dpsim-fmi/tests/FMI_SampledTest"
SAMPLED_FMU = ROOT / "build-fmi/dpsim-fmi/tests/sampled_threephase.fmu"


def load_case(step_us: int, rerun: bool) -> tuple[pd.DataFrame, float]:
    folder = BUILD / f"{step_us}us"
    folder.mkdir(parents=True, exist_ok=True)
    config = json.loads(BASE_CONFIG.read_text())
    config.update(
        step_s=step_us * 1e-6, log_downsampling=1, multirate_mode="single_pass_explicit"
    )
    config_path = folder / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    csv_path = folder / "logs/fmi_mixed.csv"
    run_log = folder / "run.log"
    runtime_manifest = BUILD / "runtime_manifest.json"
    saved_runtimes = (
        json.loads(runtime_manifest.read_text()) if runtime_manifest.is_file() else {}
    )
    elapsed = float("nan")
    if rerun or not csv_path.is_file():
        started = time.perf_counter()
        result = subprocess.run(
            [str(EXE), str(FMU), "config.json"],
            cwd=folder,
            capture_output=True,
            text=True,
        )
        elapsed = time.perf_counter() - started
        run_log.write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(
                f"{step_us} us run failed:\n{result.stdout}\n{result.stderr}"
            )
        saved_runtimes[str(step_us)] = elapsed
        runtime_manifest.write_text(json.dumps(saved_runtimes, indent=2) + "\n")
    elif str(step_us) in saved_runtimes:
        elapsed = float(saved_runtimes[str(step_us)])
    data = pd.read_csv(csv_path, skipinitialspace=True)
    data.columns = data.columns.str.strip()
    expected = round(config["stop_time_s"] / config["step_s"]) + 1
    if len(data) != expected:
        raise RuntimeError(
            f"{step_us} us: expected {expected} samples, found {len(data)}"
        )
    return data, elapsed


def alpha_beta(frame: pd.DataFrame, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    a, b, c = (frame[f"{prefix}_{k}"].to_numpy() for k in range(3))
    return (2 / 3) * (a - 0.5 * b - 0.5 * c), (1 / np.sqrt(3)) * (b - c)


def power(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    va, vb = alpha_beta(frame, "pcc_voltage")
    ia, ib = alpha_beta(frame, "converter_current")
    return 1.5 * (va * ia + vb * ib), 1.5 * (vb * ia - va * ib)


def run_fixture_studies(out: pathlib.Path) -> None:
    work = out / "fixture_runs"
    work.mkdir(parents=True, exist_ok=True)
    nonlinear = {}
    for label, args in (("baseline", []), ("tight_newton_and_jacobian", ["tight"])):
        result = subprocess.run(
            [str(NONLINEAR_TEST), str(NONLINEAR_FMU), *args],
            cwd=work,
            capture_output=True,
            text=True,
        )
        (work / f"nonlinear_{label}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        match = re.search(
            r"maximum relative Norton W error=([^,]+), history error=([^;]+); "
            r"FMI evaluation requests=([^,]+), evaluation wall time=([^ ]+) s; "
            r"max iterations=([^,]+), max scaled residual=([^;]+)",
            result.stdout,
        )
        if not match:
            raise RuntimeError("Could not parse nonlinear fixture diagnostics")
        fields = (
            "max_relative_W_error",
            "max_relative_history_error",
            "evaluation_requests",
            "evaluation_wall_time_s",
            "max_newton_iterations",
            "max_scaled_residual",
        )
        nonlinear[label] = {
            key: float(value) for key, value in zip(fields, match.groups())
        }

    sampled = {}
    for label, mode in (
        ("rollback_newton", "multirate"),
        ("single_pass_explicit", "single_pass"),
    ):
        result = subprocess.run(
            [str(SAMPLED_TEST), str(SAMPLED_FMU), mode],
            cwd=work,
            capture_output=True,
            text=True,
        )
        (work / f"sampled_{label}.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        match = re.search(
            r"max explicit lag=([0-9.eE+-]+) A, max post-step KCL defect="
            r"([0-9.eE+-]+) A; independent recurrence state error < ([0-9.eE+-]+)",
            result.stdout,
        )
        if not match:
            raise RuntimeError("Could not parse sampled fixture diagnostics")
        fields = (
            "max_current_lag_A",
            "max_poststep_kcl_defect_A",
            "independent_state_error_bound",
        )
        sampled[label] = {
            key: float(value) for key, value in zip(fields, match.groups())
        }
    (out / "fixture_studies.json").write_text(
        json.dumps(
            {
                "nonlinear_sensitivity": nonlinear,
                "same_network_multirate_coupling": sampled,
            },
            indent=2,
        )
        + "\n"
    )


def compare_matlab_voltage_sag(out: pathlib.Path) -> None:
    if not MATLAB_REFERENCE.is_file():
        print(
            "No regenerated MATLAB voltage-sag CSV; skipping independent reference metrics."
        )
        return
    d = pd.read_csv(
        ROOT / "build-fmi/fmi-validation/voltage_sag/logs/fmi_mixed.csv",
        skipinitialspace=True,
    )
    d.columns = d.columns.str.strip()
    reference_columns = [
        "time",
        *(f"i_grid_{p}" for p in "abc"),
        *(f"v_pcc_{p}" for p in "abc"),
        "v_d_pu",
        "v_q_pu",
        "pi_d",
        "pi_q",
        "pll_frequency",
        "pll_theta",
        "pll_magnitude",
    ]
    ref = pd.read_csv(MATLAB_REFERENCE, usecols=reference_columns)
    # The native run is sampled every 100 us; the physical reference is 1 us.
    ref_index = np.rint(ref.time.to_numpy() * 1e6).astype(np.int64)
    ref_lookup = pd.Series(np.arange(len(ref)), index=ref_index)
    native_index = np.rint(d.time.to_numpy() * 1e6).astype(np.int64)
    row_indices = ref_lookup.reindex(native_index)
    if row_indices.isna().any():
        raise RuntimeError("MATLAB and DPsim voltage-sag times are not coincident")
    ref = ref.iloc[row_indices.to_numpy(dtype=np.int64)].reset_index(drop=True)

    voltage_base = 25000 * np.sqrt(2 / 3)
    current_base = 500000 / 25000 * np.sqrt(2 / 3)
    bases = {
        "voltage": voltage_base,
        "current": current_base,
        "pu": 1.0,
        "frequency": 60.0,
        "angle": 2 * np.pi,
    }
    pairs = {}
    for k, phase in enumerate("abc"):
        pairs[f"v_pcc_{phase}"] = (
            d[f"pcc_voltage_{k}"].to_numpy(),
            ref[f"v_pcc_{phase}"].to_numpy(),
            bases["voltage"],
        )
        pairs[f"i_grid_{phase}"] = (
            -d[f"line_current_{k}"].to_numpy(),
            ref[f"i_grid_{phase}"].to_numpy(),
            bases["current"],
        )
        pairs[f"i_converter_{phase}"] = (
            d[f"converter_current_{k}"].to_numpy(),
            ref[f"i_grid_{phase}"].to_numpy()
            - ref[f"v_pcc_{phase}"].to_numpy() / 12500,
            bases["current"],
        )
        pairs[f"i_load_{phase}"] = (
            d[f"load_current_{k}"].to_numpy(),
            ref[f"v_pcc_{phase}"].to_numpy() / 12500,
            bases["current"],
        )
    for signal in (
        "v_d_pu",
        "v_q_pu",
        "pi_d",
        "pi_q",
        "pll_frequency",
        "pll_theta",
        "pll_magnitude",
    ):
        base = (
            bases["frequency"]
            if signal == "pll_frequency"
            else bases["angle"] if signal == "pll_theta" else bases["pu"]
        )
        pairs[signal] = (d[signal].to_numpy(), ref[signal].to_numpy(), base)

    windows = {
        "full_overlap": (0.0, 0.199901),
        "sag_170_190ms": (0.17, 0.19),
        "recovery_190_1999ms": (0.19, 0.199901),
    }
    rows = []
    for window, (start, end) in windows.items():
        mask = (d.time.to_numpy() >= start - 1e-12) & (d.time.to_numpy() <= end + 1e-12)
        if window != "full_overlap":
            mask &= d.time.to_numpy() < end - 1e-12
        for signal, (actual, expected, base) in pairs.items():
            error = actual[mask] - expected[mask]
            raw_max = float(np.max(np.abs(error)))
            if signal == "pll_theta":
                error = (error + np.pi) % (2 * np.pi) - np.pi
            rows.append(
                {
                    "window": window,
                    "signal": signal,
                    "max_abs_error": float(np.max(np.abs(error))),
                    "rms_error": float(np.sqrt(np.mean(error * error))),
                    "max_error_percent_base": 100 * float(np.max(np.abs(error))) / base,
                    "rms_error_percent_base": 100
                    * float(np.sqrt(np.mean(error * error)))
                    / base,
                    "raw_max_abs_error": raw_max,
                }
            )
    pd.DataFrame(rows).to_csv(out / "matlab_voltage_sag_metrics.csv", index=False)
    (out / "matlab_voltage_sag_metrics.json").write_text(
        json.dumps(
            {
                "reference": str(MATLAB_REFERENCE),
                "reference_timestep_s": 1e-6,
                "native_timestep_s": 1e-6,
                "native_log_spacing_s": 1e-4,
                "coincident_samples": int(len(d)),
                "overlap_end_s": float(d.time.iloc[-1]),
                "metrics": rows,
            },
            indent=2,
        )
        + "\n"
    )

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, signal, label in zip(
        axes,
        ("v_pcc_a", "i_converter_a", "pi_q"),
        (
            "PCC phase-A voltage [V]",
            "Converter phase-A current [A]",
            "q-axis controller output",
        ),
    ):
        actual, expected, _ = pairs[signal]
        ax.plot(
            d.time * 1e3, expected, label="Simulink physical reference", linewidth=0.75
        )
        ax.plot(d.time * 1e3, actual, "--", label="DPsim imported FMU", linewidth=0.75)
        ax.axvline(150, color="black", alpha=0.35, linewidth=0.7)
        ax.axvspan(170, 190, color="tab:red", alpha=0.09)
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[0].legend(ncol=2)
    axes[-1].set_xlabel("Time [ms]")
    fig.tight_layout()
    fig.savefig(out / "matlab_voltage_sag_overlay.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for ax, signal, base, label in zip(
        axes,
        ("v_pcc_a", "i_converter_a", "pll_frequency"),
        (voltage_base, current_base, bases["frequency"]),
        (
            "PCC voltage A [% phase-peak base]",
            "Converter current A [% rated-current base]",
            "PLL frequency [% 60 Hz base]",
        ),
    ):
        actual, expected, _ = pairs[signal]
        error = actual - expected
        ax.plot(d.time * 1e3, 100 * error / base, linewidth=0.8)
        ax.axvline(150, color="black", alpha=0.35, linewidth=0.7)
        ax.axvspan(170, 190, color="tab:red", alpha=0.09)
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("Time [ms]")
    fig.tight_layout()
    fig.savefig(out / "matlab_voltage_sag_errors.png", dpi=180)
    plt.close(fig)


def compare_matlab_load_steps(out: pathlib.Path) -> None:
    """Compare independent physical Simulink load-step cases with DPsim logs."""
    reference_dir = ROOT / "build-fmi/fmi-validation/matlab-reference-attempt"
    reference_columns = [
        "time",
        *(f"i_grid_{p}" for p in "abc"),
        *(f"v_pcc_{p}" for p in "abc"),
        "v_d_pu",
        "v_q_pu",
        "pi_d",
        "pi_q",
        "pll_frequency",
        "pll_theta",
        "pll_magnitude",
    ]
    voltage_base = 25000 * np.sqrt(2 / 3)
    current_base = 500000 / 25000 * np.sqrt(2 / 3)
    rows = []
    plot_data = {}
    for scenario, step_resistances in (
        ("balanced_load", (25000.0, 25000.0, 25000.0)),
        ("unbalanced_load", (25000.0, 12500.0, 50000.0)),
    ):
        reference_path = reference_dir / f"{scenario}_reference_trap_1us.csv"
        if not reference_path.is_file():
            print(f"No regenerated MATLAB {scenario} CSV; skipping this comparison.")
            continue
        d = pd.read_csv(
            ROOT / f"build-fmi/fmi-validation/{scenario}/logs/fmi_mixed.csv",
            skipinitialspace=True,
        )
        d.columns = d.columns.str.strip()
        ref = pd.read_csv(reference_path, usecols=reference_columns)
        ref_index = np.rint(ref.time.to_numpy() * 1e6).astype(np.int64)
        lookup = pd.Series(np.arange(len(ref)), index=ref_index)
        native_index = np.rint(d.time.to_numpy() * 1e6).astype(np.int64)
        selected = lookup.reindex(native_index)
        if selected.isna().any():
            raise RuntimeError(f"MATLAB and DPsim {scenario} samples do not coincide")
        ref = ref.iloc[selected.to_numpy(dtype=np.int64)].reset_index(drop=True)

        pairs = {}
        step_enabled = ref.time.to_numpy() >= 0.17 - 1e-12
        native_voltage_abc, native_current_abc = [], []
        reference_voltage_abc, reference_current_abc = [], []
        for k, phase in enumerate("abc"):
            ref_grid = ref[f"i_grid_{phase}"].to_numpy()
            ref_voltage = ref[f"v_pcc_{phase}"].to_numpy()
            ref_step_load = np.where(
                step_enabled, ref_voltage / step_resistances[k], 0.0
            )
            pairs[f"pcc_voltage_{phase}"] = (
                d[f"pcc_voltage_{k}"].to_numpy(),
                ref_voltage,
                voltage_base,
            )
            pairs[f"grid_current_{phase}"] = (
                -d[f"line_current_{k}"].to_numpy(),
                ref_grid,
                current_base,
            )
            pairs[f"converter_current_{phase}"] = (
                d[f"converter_current_{k}"].to_numpy(),
                ref_grid - ref_voltage / 12500.0 - ref_step_load,
                current_base,
            )
            native_voltage_abc.append(d[f"pcc_voltage_{k}"].to_numpy())
            native_current_abc.append(d[f"converter_current_{k}"].to_numpy())
            reference_voltage_abc.append(ref_voltage)
            reference_current_abc.append(
                ref_grid - ref_voltage / 12500.0 - ref_step_load
            )
            pairs[f"base_load_current_{phase}"] = (
                d[f"load_current_{k}"].to_numpy(),
                ref_voltage / 12500.0,
                current_base,
            )
            pairs[f"step_load_current_{phase}"] = (
                d[f"step_load_current_{k}"].to_numpy(),
                ref_step_load,
                current_base,
            )
        for signal in (
            "v_d_pu",
            "v_q_pu",
            "pi_d",
            "pi_q",
            "pll_frequency",
            "pll_theta",
            "pll_magnitude",
        ):
            base = (
                60.0
                if signal == "pll_frequency"
                else 2 * np.pi if signal == "pll_theta" else 1.0
            )
            pairs[signal] = (d[signal].to_numpy(), ref[signal].to_numpy(), base)

        def instantaneous_power(
            voltage: list[np.ndarray], current: list[np.ndarray]
        ) -> tuple[np.ndarray, np.ndarray]:
            va, vb, vc = voltage
            ia, ib, ic = current
            valpha = (2 / 3) * (va - 0.5 * vb - 0.5 * vc)
            vbeta = (1 / np.sqrt(3)) * (vb - vc)
            ialpha = (2 / 3) * (ia - 0.5 * ib - 0.5 * ic)
            ibeta = (1 / np.sqrt(3)) * (ib - ic)
            return (
                1.5 * (valpha * ialpha + vbeta * ibeta),
                1.5 * (vbeta * ialpha - valpha * ibeta),
            )

        native_p, native_q = instantaneous_power(native_voltage_abc, native_current_abc)
        reference_p, reference_q = instantaneous_power(
            reference_voltage_abc, reference_current_abc
        )
        pairs["instantaneous_active_power"] = (native_p, reference_p, 500000.0)
        pairs["instantaneous_reactive_power"] = (native_q, reference_q, 500000.0)
        plot_data[scenario] = (d, ref, pairs)

        # The breaker closes at 150 ms; report the established loaded interval and
        # the post-170 ms load-step response separately from the startup transient.
        windows = {
            "unenergized_0_150ms": (0.0, 0.15),
            "breaker_startup_150_160ms": (0.15, 0.16),
            "pre_load_step_160_170ms": (0.16, 0.17),
            "load_step_170_1999ms": (0.17, 0.1999),
        }
        for window, (start, end) in windows.items():
            mask = (d.time.to_numpy() >= start - 1e-12) & (
                d.time.to_numpy() < end - 1e-12
            )
            for signal, (actual, expected, base) in pairs.items():
                error = actual[mask] - expected[mask]
                if signal == "pll_theta":
                    error = (error + np.pi) % (2 * np.pi) - np.pi
                rows.append(
                    {
                        "scenario": scenario,
                        "window": window,
                        "signal": signal,
                        "max_abs_error": float(np.max(np.abs(error))),
                        "rms_error": float(np.sqrt(np.mean(error * error))),
                        "max_error_percent_base": 100
                        * float(np.max(np.abs(error)))
                        / base,
                        "rms_error_percent_base": 100
                        * float(np.sqrt(np.mean(error * error)))
                        / base,
                    }
                )
    if rows:
        pd.DataFrame(rows).to_csv(out / "matlab_load_step_metrics.csv", index=False)
        (out / "matlab_load_step_metrics.json").write_text(
            json.dumps(
                {
                    "reference_timestep_s": 1e-6,
                    "native_timestep_s": 1e-6,
                    "native_log_spacing_s": 1e-4,
                    "breaker_close_s": 0.15,
                    "load_step_s": 0.17,
                    "converter_current_reference": "MATLAB grid current minus base and switched resistive-load currents, from PCC KCL",
                    "power_definition": "instantaneous alpha-beta power from PCC voltage and converter current",
                    "metrics": rows,
                },
                indent=2,
            )
            + "\n"
        )
        fig, axes = plt.subplots(3, 2, figsize=(12, 8), sharex=True)
        for col, scenario in enumerate(("balanced_load", "unbalanced_load")):
            if scenario not in plot_data:
                continue
            d, _, pairs = plot_data[scenario]
            for ax, signal, label in zip(
                axes[:, col],
                ("pcc_voltage_a", "converter_current_a", "pll_frequency"),
                (
                    "PCC phase-A voltage [V]",
                    "Converter phase-A current [A]",
                    "PLL frequency [Hz]",
                ),
            ):
                actual, expected, _ = pairs[signal]
                ax.plot(d.time * 1e3, expected, label="Simulink", linewidth=0.75)
                ax.plot(d.time * 1e3, actual, "--", label="DPsim", linewidth=0.75)
                ax.axvline(150, color="black", alpha=0.35, linewidth=0.7)
                ax.axvline(170, color="tab:red", alpha=0.45, linewidth=0.7)
                ax.set_ylabel(label)
                ax.grid(alpha=0.25)
            axes[0, col].set_title(scenario.replace("_", " ").title())
        axes[0, 0].legend(ncol=2)
        axes[-1, 0].set_xlabel("Time [ms]")
        axes[-1, 1].set_xlabel("Time [ms]")
        fig.tight_layout()
        fig.savefig(out / "matlab_load_step_overlay.png", dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true", help="reuse completed CSV runs")
    args = parser.parse_args()
    if not EXE.is_file() or not FMU.is_file():
        raise SystemExit(
            "Build FMI_EMT_Ph3_Mixed and confirm the Linux FMU is present."
        )
    runs, runtime = {}, {}
    for step_us in (1, 25, 50):
        runs[step_us], runtime[step_us] = load_case(step_us, not args.reuse)

    reference = runs[1].copy()
    reference["sample_us"] = np.rint(reference.time * 1e6).astype(int)
    reference = reference.set_index("sample_us")
    windows = {
        "full": (0.0, 0.020),
        "sag_5_10ms": (0.005, 0.010),
        "recovery_10_20ms": (0.010, 0.020),
    }
    signals = [
        *(f"pcc_voltage_{k}" for k in range(3)),
        *(f"converter_current_{k}" for k in range(3)),
        "pll_frequency",
        "pi_d",
        "pi_q",
        "v_d_pu",
        "v_q_pu",
    ]
    metrics = []
    runtime_rows = []
    for step_us, data in runs.items():
        data["sample_us"] = np.rint(data.time * 1e6).astype(int)
        aligned = reference.loc[data.sample_us.to_numpy()]
        for window, (start, end) in windows.items():
            mask = (data.time.to_numpy() >= start - 1e-12) & (
                data.time.to_numpy() <= end + 1e-12
            )
            if window == "sag_5_10ms":
                mask &= data.time.to_numpy() < end - 1e-12
            for signal in signals:
                error = data[signal].to_numpy() - aligned[signal].to_numpy()
                e = error[mask]
                metrics.append(
                    {
                        "step_us": step_us,
                        "window": window,
                        "signal": signal,
                        "max_abs_error": float(np.max(np.abs(e))),
                        "rms_error": float(np.sqrt(np.mean(e * e))),
                    }
                )
            p, q = power(data)
            pr, qr = power(aligned)
            for signal, values, ref_values in (("P_W", p, pr), ("Q_var", q, qr)):
                e = values[mask] - ref_values[mask]
                metrics.append(
                    {
                        "step_us": step_us,
                        "window": window,
                        "signal": signal,
                        "max_abs_error": float(np.max(np.abs(e))),
                        "rms_error": float(np.sqrt(np.mean(e * e))),
                    }
                )
        after_startup = data.time.to_numpy() >= 0.001
        runtime_rows.append(
            {
                "step_us": step_us,
                "wall_time_s": runtime[step_us],
                "fmi_evaluation_requests": int(data.fmi_evaluation_requests.iloc[-1]),
                "fmi_evaluation_wall_time_s": float(
                    data.fmi_evaluation_wall_time_s.iloc[-1]
                ),
                "max_newton_iterations": int(data.newton_iterations.max()),
                "max_reported_residual": float(data.newton_residual.max()),
                "max_coupling_error_A_after_1ms": float(
                    data.loc[after_startup, "fmi_coupling_error_A"].max()
                ),
                "internal_steps_min": int(data.fmi_internal_steps.iloc[1:].min()),
                "internal_steps_max": int(data.fmi_internal_steps.iloc[1:].max()),
            }
        )

    out = BUILD / "paper"
    out.mkdir(parents=True, exist_ok=True)
    run_fixture_studies(out)
    compare_matlab_voltage_sag(out)
    compare_matlab_load_steps(out)
    pd.DataFrame(metrics).to_csv(out / "weak_grid_refinement_metrics.csv", index=False)
    pd.DataFrame(runtime_rows).to_csv(
        out / "weak_grid_runtime_metrics.csv", index=False
    )
    (out / "weak_grid_refinement_metrics.json").write_text(
        json.dumps(
            {
                "reference": "same DPsim model at 1 us",
                "metrics": metrics,
                "runtime": runtime_rows,
            },
            indent=2,
        )
        + "\n"
    )

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    names = (
        ("pcc_voltage_0", "PCC phase-A voltage [V]"),
        ("converter_current_0", "Converter phase-A current [A]"),
        ("pi_q", "q-axis current controller output"),
    )
    for ax, (signal, label) in zip(axes, names):
        for step_us, frame in runs.items():
            stride = max(1, int(50 / step_us))
            ax.plot(
                frame.time.iloc[::stride] * 1e3,
                frame[signal].iloc[::stride],
                label=f"{step_us} µs",
                linewidth=0.8,
            )
        ax.axvspan(5, 10, color="tab:red", alpha=0.09)
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[0].legend(ncol=3)
    axes[-1].set_xlabel("Time [ms]")
    fig.tight_layout()
    fig.savefig(out / "weak_grid_refinement_waveforms.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for step_us in (25, 50):
        frame = runs[step_us]
        aligned = reference.loc[frame.sample_us.to_numpy()]
        for ax, signal, base in (
            (axes[0], "pcc_voltage_0", 20412.4),
            (axes[1], "converter_current_0", 2.5),
        ):
            delta = frame[signal].to_numpy() - aligned[signal].to_numpy()
            ax.plot(
                frame.time * 1e3,
                delta / base * 100,
                label=f"{step_us} µs",
                linewidth=0.8,
            )
    for ax, label in zip(
        axes,
        (
            "PCC phase-A error [% voltage peak base]",
            "Converter phase-A error [% 2.5 A base]",
        ),
    ):
        ax.axvspan(5, 10, color="tab:red", alpha=0.09)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("Time [ms]")
    fig.tight_layout()
    fig.savefig(out / "weak_grid_refinement_errors.png", dpi=180)
    plt.close(fig)

    case_rows = []
    case_frames = {}
    gates = ROOT / "build-fmi/fmi-validation/gate3"
    gate3_cases = {
        "loaded_start": (0.010, 0.020),
        "weak_grid_sag": (0.005, 0.010),
        "unbalanced_sag": (0.005, 0.010),
    }
    for name, (start, end) in gate3_cases.items():
        frame = pd.read_csv(gates / name / "logs/fmi_mixed.csv", skipinitialspace=True)
        frame.columns = frame.columns.str.strip()
        case_frames[name] = frame
        mask = (frame.time >= start - 1e-12) & (frame.time <= end + 1e-12)
        if name != "loaded_start":
            mask &= frame.time < end - 1e-12
        v = frame.loc[mask, [f"pcc_voltage_{k}" for k in range(3)]].to_numpy()
        i = frame.loc[mask, [f"converter_current_{k}" for k in range(3)]].to_numpy()
        p = np.sum(v * i, axis=1)
        va, vb = alpha_beta(frame.loc[mask], "pcc_voltage")
        ia, ib = alpha_beta(frame.loc[mask], "converter_current")
        q = 1.5 * (vb * ia - va * ib)
        lag_mask = (frame.time >= max(start, 0.001)) & (frame.time < end + 1e-12)
        case_rows.append(
            {
                "case": name,
                "window": f"{start*1e3:g}-{end*1e3:g} ms",
                "samples": int(mask.sum()),
                "pcc_voltage_rms_V": float(np.sqrt(np.mean(v * v))),
                "converter_current_rms_A": float(np.sqrt(np.mean(i * i))),
                "mean_P_W": float(np.mean(p)),
                "mean_Q_var": float(np.mean(q)),
                "max_coupling_lag_A": float(
                    frame.loc[lag_mask, "fmi_coupling_error_A"].max()
                ),
                "max_network_passes": int(frame.newton_iterations.max()),
                "internal_steps_per_macro": int(
                    frame.fmi_internal_steps.iloc[1:].max()
                ),
            }
        )
    pd.DataFrame(case_rows).to_csv(out / "converter_case_metrics.csv", index=False)
    (out / "converter_case_metrics.json").write_text(
        json.dumps(case_rows, indent=2) + "\n"
    )
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    for row, (name, frame) in enumerate(case_frames.items()):
        axes[row, 0].plot(frame.time * 1e3, frame.pcc_voltage_0, linewidth=0.8)
        axes[row, 1].plot(frame.time * 1e3, frame.converter_current_0, linewidth=0.8)
        axes[row, 0].set_ylabel(name.replace("_", " ") + "\nPCC A [V]")
        axes[row, 1].set_ylabel("converter A [A]")
        axes[row, 0].grid(alpha=0.25)
        axes[row, 1].grid(alpha=0.25)
        if name != "loaded_start":
            for ax in axes[row, :]:
                ax.axvspan(5, 10, color="tab:red", alpha=0.09)
    axes[-1, 0].set_xlabel("Time [ms]")
    axes[-1, 1].set_xlabel("Time [ms]")
    fig.tight_layout()
    fig.savefig(out / "converter_case_waveforms.png", dpi=180)
    plt.close(fig)
    print(pd.DataFrame(runtime_rows).to_string(index=False))
    print(f"Wrote weak-grid metrics and figures under {out}")


if __name__ == "__main__":
    main()
