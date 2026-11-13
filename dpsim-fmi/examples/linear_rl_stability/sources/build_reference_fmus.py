#!/usr/bin/env python3
"""Build self-contained FMI 2.0 Linux FMUs from the checked-in C++ source."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEADER = HERE.parents[2] / "third_party" / "fmi2"
CONFIGS = (("ME", "LinearRL_ME", []), ("CS", "LinearRL_CS", ["-DFMI_RL_COSIM=1"]))


def build(mode: str, identifier: str, defines: list[str]) -> None:
    destination = HERE / f"LinearRL_{mode}.fmu"
    xml = HERE / f"modelDescription_{mode}.xml"
    with tempfile.TemporaryDirectory(prefix=f"linear_rl_{mode}_") as temp:
        root = Path(temp)
        binary = root / "binaries" / "linux64" / f"{identifier}.so"
        binary.parent.mkdir(parents=True)
        compiler = os.environ.get("CXX", "c++")
        subprocess.run(
            [
                compiler,
                "-std=c++17",
                "-O2",
                "-fPIC",
                "-shared",
                f"-I{HEADER}",
                *defines,
                str(HERE / "LinearRL.cpp"),
                "-o",
                str(binary),
            ],
            check=True,
        )
        shutil.copy2(xml, root / "modelDescription.xml")
        with zipfile.ZipFile(
            destination, "w", compression=zipfile.ZIP_DEFLATED
        ) as archive:
            archive.write(root / "modelDescription.xml", "modelDescription.xml")
            archive.write(binary, f"binaries/linux64/{identifier}.so")
    print(f"Created {destination}")


if __name__ == "__main__":
    for mode, identifier, defines in CONFIGS:
        build(mode, identifier, defines)
