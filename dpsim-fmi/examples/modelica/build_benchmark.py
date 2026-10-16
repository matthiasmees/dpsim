# SPDX-License-Identifier: MPL-2.0
"""Export the nonlinear port and simulate original MSL examples on this host."""
import argparse
import json
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--omc", type=Path, required=True)
parser.add_argument("--msl", type=Path, required=True, help="MSL 4.0.0 checkout root")
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
(output / "tmp").mkdir(exist_ok=True)
msl = args.msl.resolve()
model = Path(__file__).with_name("NonlinearThreePhase.mo").resolve()
# Load known local dependencies explicitly; do not use the online package index.
loads = [
    msl / "Complex.mo",
    msl / "ModelicaServices/package.mo",
    msl / "Modelica/package.mo",
    model,
]
script = "\n".join(f"loadFile({json.dumps(str(p))}, uses=false);" for p in loads)
script += """
getErrorString();
buildModelFMU(NonlinearThreePhase.Ported, version="3.0", fmuType="me");
getErrorString();
simulate(NonlinearThreePhase.Reference, startTime=0, stopTime=6.2832, numberOfIntervals=62832, tolerance=1e-10, method="rungekutta", outputFormat="csv", variableFilter="i_a|i_b|i_c");
getErrorString();
"""
(output / "build.mos").write_text(script)
env = os.environ.copy()
env["TMPDIR"] = str(output / "tmp")
for name in ["Ported.fmu", "NonlinearThreePhase.Reference_res.csv"]:
    (output / name).unlink(missing_ok=True)
with (output / "build.log").open("w") as log:
    subprocess.run(
        [str(args.omc.resolve()), "build.mos"],
        cwd=output,
        env=env,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
    )
# OMC may return zero even when compilation fails. Check artifacts as well.
for name in ["Ported.fmu", "NonlinearThreePhase.Reference_res.csv"]:
    if not (output / name).is_file():
        raise RuntimeError(f"Missing {name}; inspect {output / 'build.log'}")
config = dict(
    voltage_inputs=["v_a", "v_b", "v_c"],
    current_outputs=["i_a", "i_b", "i_c"],
    angular_frequency_rad_s=1,
    voltage_peak_V=1.25,
    current_sign_to_passive=1,
    step_s=0.01,
    stop_time_s=6.28,
)
(output / "benchmark.json").write_text(json.dumps(config, indent=2) + "\n")
print(f"FMU and native reference saved in {output}")
