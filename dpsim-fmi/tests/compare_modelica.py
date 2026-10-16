# SPDX-License-Identifier: MPL-2.0
"""Validate a port adaptation against unmodified MSL example instances."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def compare(root):
    root = Path(root)
    reference = pd.read_csv(root / "NonlinearThreePhase.Reference_res.csv")
    reference = reference.drop_duplicates("time", keep="last")
    metrics = {}
    for name in ["dpsim-10ms", "dpsim-5ms", "dpsim-200us", "dpsim-100us"]:
        data = pd.read_csv(
            root / name / "logs/fmi_prescribed.csv", skipinitialspace=True
        )
        data.columns = data.columns.str.strip()
        assert data.time.iloc[0] == 0 and data.time.iloc[-1] >= 6.28
        error = np.array(
            [
                data[f"i_{p}"]
                - np.interp(data.time, reference.time, reference[f"i_{p}"])
                for p in "abc"
            ]
        )
        kcl = np.array(
            [data[f"i_{p}"] + data[f"source_current_{k}"] for k, p in enumerate("abc")]
        )
        assert np.isfinite(error).all()
        assert np.max(abs(kcl)) < 1.1e-6
        assert data.newton_residual.max() <= 1.000001
        metrics[name] = {
            "peak_current_error_A": float(np.max(abs(error))),
            "rms_current_error_A": float(np.sqrt(np.mean(error**2))),
            "peak_kcl_A": float(np.max(abs(kcl))),
            "max_newton_iterations": int(data.newton_iterations.max()),
        }
    ratio = (
        metrics["dpsim-10ms"]["rms_current_error_A"]
        / metrics["dpsim-5ms"]["rms_current_error_A"]
    )
    assert 3.5 < ratio < 4.5
    assert metrics["dpsim-200us"]["peak_current_error_A"] < 1e-6
    assert metrics["dpsim-100us"]["peak_current_error_A"] < 1e-6
    result = {"runs": metrics, "rms_refinement_ratio": ratio, "passed": True}
    (root / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    print(json.dumps(compare(parser.parse_args().root), indent=2))
