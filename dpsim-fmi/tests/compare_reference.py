# SPDX-License-Identifier: MPL-2.0
"""Compare the MATLAB RLC benchmark with the imported SSN's logged outputs.

Usage: python3 compare_reference.py <logs/fmi_rlc.csv> <MATLAB-output-directory>
The MATLAB trapezoidal reference must agree to the logger's six-decimal
rounding. The ode45 comparison is reported as an integration-error measure.
"""
import bisect
import csv
import json
import math
import pathlib
import sys


def read(path):
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        rows = [
            {key.strip(): float(value) for key, value in row.items()} for row in reader
        ]
    if not rows or not all(
        math.isfinite(value) for row in rows for value in row.values()
    ):
        raise ValueError("Empty or non-finite results: " + str(path))
    return rows


def compare(actual, reference):
    times = [row["time_s"] for row in reference]
    assert all(
        right > left for left, right in zip(times, times[1:])
    ), "reference times must increase"
    errors = {"i_A": [], "v_C_V": [], "u_V": []}
    for row in actual:
        time = row["time"]
        assert (
            times[0] - 1e-12 <= time <= times[-1] + 1e-12
        ), "time outside reference interval"
        index = min(max(bisect.bisect_right(times, time), 1), len(times) - 1)
        fraction = (time - times[index - 1]) / (times[index] - times[index - 1])
        for key, actual_key in (("i_A", "i"), ("v_C_V", "v_C"), ("u_V", "u_V")):
            value = reference[index - 1][key] + fraction * (
                reference[index][key] - reference[index - 1][key]
            )
            errors[key].append(row[actual_key] - value)
    return {
        key: {
            "max_abs": max(map(abs, values)),
            "rmse": math.sqrt(sum(v * v for v in values) / len(values)),
        }
        for key, values in errors.items()
    }


actual_path, reference_directory = map(pathlib.Path, sys.argv[1:])
actual = read(actual_path)
assert len(actual) == 10001, "expected 0..0.1 s at 10 microseconds"
assert abs(actual[0]["time"]) < 1e-12 and abs(actual[-1]["time"] - 0.1) < 1e-12
assert all(
    abs((b["time"] - a["time"]) - 1e-5) < 1e-12 for a, b in zip(actual, actual[1:])
)
report = {}
for name in ("trapezoidal", "simulink"):
    report[name] = compare(
        actual, read(reference_directory / ("reference_" + name + ".csv"))
    )
for metric in report["trapezoidal"].values():
    assert (
        metric["max_abs"] < 6e-7
    ), "SSN differs from MATLAB trapezoidal reference beyond CSV rounding"
print(json.dumps(report, indent=2))
