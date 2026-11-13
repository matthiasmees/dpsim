#!/usr/bin/env python3
"""Build reproducible FMI 2.0 three-channel ME and CS FMUs from C++."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEADER = HERE.parents[2] / "third_party" / "fmi2"
CONFIGS = (("ME", "LinearRL3PhME", []), ("CS", "LinearRL3PhCS", ["-DFMI_RL_COSIM=1"]))


def build(mode: str, identifier: str, defines: list[str]) -> None:
    destination = HERE / f"LinearRL3Ph_{mode}.fmu"
    xml = HERE / f"modelDescription_{mode}.xml"
    with tempfile.TemporaryDirectory(prefix=f"linear_rl_3ph_{mode}_") as temp:
        root = Path(temp)
        binary = root / "binaries" / "linux64" / f"{identifier}.so"
        binary.parent.mkdir(parents=True)
        subprocess.run(
            [
                os.environ.get("CXX", "c++"),
                "-std=c++17",
                "-O2",
                "-fPIC",
                "-shared",
                f"-I{HEADER}",
                *defines,
                str(HERE / "LinearRL3Ph.cpp"),
                "-o",
                str(binary),
            ],
            check=True,
        )
        shutil.copy2(xml, root / "modelDescription.xml")
        with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
            archive.write(root / "modelDescription.xml", "modelDescription.xml")
            archive.write(binary, f"binaries/linux64/{identifier}.so")
    print(f"Created {destination}")


if __name__ == "__main__":
    for config in CONFIGS:
        build(*config)
