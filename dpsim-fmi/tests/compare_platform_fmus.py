#!/usr/bin/env python3
"""Check that the retained macOS and Linux FMUs describe the same models."""

from __future__ import annotations

import argparse
import hashlib
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


PAIRS = (
    ("linear_scalar", "dpsim_rlc_me"),
    ("linear_three_phase", "dpsim_rlc_3ph_me"),
    ("synchronous_generator", "dpsim_sg_3ph_me"),
    ("two_generators", "first_generator"),
    ("two_generators", "dpsim_sg_3ph_second"),
    ("inverter", "GridTiedInverterPortedV3"),
)


def normalized_description(data: bytes) -> bytes:
    root = ET.fromstring(data)
    root.attrib.pop("generationDateAndTime", None)
    root.attrib.pop("instantiationToken", None)
    return ET.tostring(root, encoding="utf-8")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def compare_pair(sources: Path, stem: str) -> None:
    mac_path = sources / f"{stem}.fmu"
    linux_path = sources / f"{stem}.x86_64-linux.fmu"
    for path in (mac_path, linux_path):
        if not path.is_file():
            raise AssertionError(f"missing FMU: {path}")

    with zipfile.ZipFile(mac_path) as mac, zipfile.ZipFile(linux_path) as linux:
        mac_names = set(mac.namelist())
        linux_names = set(linux.namelist())
        if not any(name.startswith("binaries/aarch64-darwin/") for name in mac_names):
            raise AssertionError(f"missing ARM64 macOS binary in {mac_path}")
        if not any(name.startswith("binaries/x86_64-linux/") for name in linux_names):
            raise AssertionError(f"missing x86-64 Linux binary in {linux_path}")

        mac_portable = {name for name in mac_names if not name.startswith("binaries/")}
        linux_portable = {
            name for name in linux_names if not name.startswith("binaries/")
        }
        if mac_portable != linux_portable:
            raise AssertionError(
                f"non-binary archive entries differ for {stem}: "
                f"macOS-only={sorted(mac_portable - linux_portable)}, "
                f"Linux-only={sorted(linux_portable - mac_portable)}"
            )

        mac_xml = normalized_description(mac.read("modelDescription.xml"))
        linux_xml = normalized_description(linux.read("modelDescription.xml"))
        if mac_xml != linux_xml:
            raise AssertionError(f"normalized modelDescription.xml differs for {stem}")

        # MATLAB redraws the PNG icon during export; it has no runtime role.
        ignored_payloads = {"modelDescription.xml", "terminalsAndIcons/icon.png"}
        for name in sorted(mac_portable - ignored_payloads):
            if digest(mac.read(name)) != digest(linux.read(name)):
                raise AssertionError(f"portable payload differs for {stem}: {name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--examples",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "examples",
        help="path to dpsim-fmi/examples",
    )
    args = parser.parse_args()

    for case, stem in PAIRS:
        compare_pair(args.examples / case / "sources", stem)
        print(f"PASS {case}/{stem}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, ET.ParseError, zipfile.BadZipFile) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
