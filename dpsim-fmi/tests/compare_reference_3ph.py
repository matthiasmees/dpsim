# SPDX-License-Identifier: MPL-2.0
"""Compare the three-phase FMU/SSN run with MATLAB's saved references.

Usage: python3 compare_reference_3ph.py logs/fmi_rlc_3ph.csv MATLAB-output-directory
Uses only the Python standard library. Assertions check the time grid, all nine
channels against MATLAB trapezoidal integration, and per-phase network KCL.
"""
import bisect
import csv
import json
import math
from pathlib import Path
import sys


def read(path):
    with open(path, newline="") as stream:
        rows = [
            {key.strip(): float(value) for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]
    if not rows or not all(math.isfinite(v) for row in rows for v in row.values()):
        raise ValueError("Empty or non-finite results: " + str(path))
    return rows


def compare(actual, reference):
    times = [row["time_s"] for row in reference]
    assert all(b > a for a, b in zip(times, times[1:]))
    errors = {f"{name}_{p}": [] for name in ("u", "i", "v_C") for p in "abc"}
    for row in actual:
        t = row["time"]
        assert times[0] - 1e-12 <= t <= times[-1] + 1e-12
        j = min(max(bisect.bisect_right(times, t), 1), len(times) - 1)
        fraction = (t - times[j - 1]) / (times[j] - times[j - 1])
        for key, values in errors.items():
            ref_key = key + ("_A" if key.startswith("i_") else "_V")
            value = reference[j - 1][ref_key] + fraction * (
                reference[j][ref_key] - reference[j - 1][ref_key]
            )
            values.append(row[key] - value)
    return {
        key: {
            "max_abs": max(map(abs, values)),
            "rmse": math.sqrt(sum(v * v for v in values) / len(values)),
        }
        for key, values in errors.items()
    }


def validate(actual_path, reference_directory):
    actual = read(actual_path)
    reference_directory = Path(reference_directory)
    assert len(actual) == 10001, "expected 10001 samples"
    assert abs(actual[0]["time"]) < 1e-12
    assert abs(actual[-1]["time"] - 0.1) < 1e-12
    assert all(
        abs(b["time"] - a["time"] - 1e-5) < 1e-12 for a, b in zip(actual, actual[1:])
    ), "unexpected time grid"
    report = {
        name: compare(actual, read(reference_directory / f"reference_{name}.csv"))
        for name in ("trapezoidal", "simulink")
    }
    for metric in report["trapezoidal"].values():
        assert metric["max_abs"] < 6e-7, "difference exceeds CSV rounding"
    report["max_kcl_residual_A"] = max(
        abs(row[f"i_{p}"] + row[f"source_current_{k}"])
        for row in actual
        for k, p in enumerate("abc")
    )
    assert report["max_kcl_residual_A"] < 1.1e-6
    report["max_current_sum_A"] = max(
        abs(sum(row[f"i_{p}"] for p in "abc")) for row in actual
    )
    assert report["max_current_sum_A"] < 1.1e-6
    return report


if __name__ == "__main__":
    print(json.dumps(validate(*sys.argv[1:]), indent=2))
