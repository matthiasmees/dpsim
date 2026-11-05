# SPDX-License-Identifier: MPL-2.0
"""Validate nonlinear SG runs against MATLAB and timestep/tolerance checks.

Usage: python compare_nonlinear_sg.py dpsim-fmi/examples/synchronous_generator/results
Requires numpy/pandas. Engineering bases come from MATLAB's saved configuration;
no generator equations or fit parameters are used here.
"""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def read(path):
    frame = pd.read_csv(path, skipinitialspace=True)
    frame.columns = frame.columns.str.strip()
    if "input_vref_pu" in frame and "vref_pu" not in frame:
        frame["vref_pu"] = frame["input_vref_pu"]
    assert len(frame) > 1 and np.isfinite(frame.to_numpy()).all(), path
    return frame


def signal_base(name, config):
    bases = config["bases"]
    if name in ("i_a", "i_b", "i_c") or name.startswith("i_line_"):
        return bases["phase_current_A"]
    if name.startswith("v_bus_"):
        return bases["phase_voltage_V"]
    if name.startswith("p_"):
        return bases["power_W"]
    if name == "field_current_A":
        return bases["field_current_A"]
    if name == "psi_f_Wb":
        return bases["field_flux_Wb"]
    if name.startswith("psi_"):
        return bases["stator_flux_Wb"]
    if name in ("theta_e_rad", "theta_state_rad"):
        return bases["angle_rad"]
    if name == "g_load_S":
        return bases["conductance_S"]
    return bases["per_unit"]


def compare(actual, reference, config):
    ta = actual["time" if "time" in actual else "time_s"].to_numpy()
    tr = reference["time" if "time" in reference else "time_s"].to_numpy()
    assert len(ta) == len(tr) and np.max(np.abs(ta - tr)) < 1e-9
    result = {}
    for name in reference.columns:
        if (
            name
            in (
                "time",
                "time_s",
                "newton_iterations",
                "newton_residual",
                "scaled_residual",
            )
            or name not in actual
        ):
            continue
        # Reference outputs comprise physical states, control diagnostics and
        # external grid signals. Importer-only bookkeeping is excluded.
        if name.startswith("network_current_") or name.startswith("input_"):
            continue
        error = actual[name].to_numpy() - reference[name].to_numpy()
        base = signal_base(name, config)
        maximum = float(np.max(np.abs(error)))
        rmse = float(np.sqrt(np.mean(error**2)))
        result[name] = dict(
            max_abs=maximum,
            rmse=rmse,
            engineering_base=base,
            max_normalized=maximum / base,
            rms_normalized=rmse / base,
        )
    assert len(result) == 34, "all 34 benchmark signals must be compared"
    return result


def validate(directory):
    directory = Path(directory)
    config = json.loads((directory / "matlab/benchmark.json").read_text())
    matlab = {
        step: read(directory / f"matlab/reference_trapezoidal_{step}us.csv")
        for step in (50, 25)
    }
    actual = {
        step: read(directory / f"dpsim-{step}us/logs/fmi_sg_nonlinear.csv")
        for step in (50, 25)
    }
    tight = read(directory / "dpsim-50us-tight/logs/fmi_sg_nonlinear.csv")
    ode = read(directory / "matlab-origin/reference_matlab_ode15s.csv")
    report = {
        "acceptance_max_and_rms_normalized": 1e-3,
        "engineering_bases": config["bases"],
        "logger_note": "DPsim CSV rounds signal values to six decimals; small errors can be quantization-limited.",
    }
    for step in (50, 25):
        report[f"dpsim_vs_trapezoidal_{step}us"] = compare(
            actual[step], matlab[step], config
        )
        for metric in report[f"dpsim_vs_trapezoidal_{step}us"].values():
            assert metric["max_normalized"] < 1e-3 and metric["rms_normalized"] < 1e-3
        frame = actual[step]
        assert len(frame) == round(config["stop_time_s"] / (step * 1e-6)) + 1
        assert frame.newton_residual.max() <= 1
        assert frame.newton_iterations.iloc[1:].min() >= 1
        residual = max(
            float(np.max(np.abs(frame[f"i_{p}"] - frame[f"network_current_{k}"])))
            for k, p in enumerate("abc")
        )
        assert residual < 1.1e-6
        report[f"solver_{step}us"] = dict(
            max_scaled_residual=float(frame.newton_residual.max()),
            max_newton_iterations=int(frame.newton_iterations.max()),
            max_kcl_residual_A=residual,
        )
    report["tight_vs_default_tolerance"] = compare(tight, actual[50], config)
    # Both steps are compared to ode15s on the common original 50-us grid,
    # so no linear interpolation of the sharp switching transient is needed.
    report["trapezoidal_vs_ode15s_50us"] = compare(matlab[50], ode, config)
    report["trapezoidal_vs_ode15s_25us"] = compare(
        matlab[25].iloc[::2].reset_index(drop=True), ode, config
    )
    for step in (50, 25):
        for metric in report[f"trapezoidal_vs_ode15s_{step}us"].values():
            assert metric["max_normalized"] < 1e-3 and metric["rms_normalized"] < 1e-3
    report["bus_voltage_refinement_ratio"] = (
        report["trapezoidal_vs_ode15s_50us"]["v_bus_a_V"]["max_abs"]
        / report["trapezoidal_vs_ode15s_25us"]["v_bus_a_V"]["max_abs"]
    )
    assert 3.5 < report["bus_voltage_refinement_ratio"] < 4.5
    (directory / "comparison.json").write_text(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    result = validate(sys.argv[1])
    for label, value in result.items():
        if isinstance(value, dict) and "i_a" in value:
            worst = max(value, key=lambda name: value[name]["max_normalized"])
            print(label, "worst:", worst, value[worst])
    print(
        "PASS: nonlinear MNA, MATLAB agreement, KCL, convergence and second-order refinement"
    )
