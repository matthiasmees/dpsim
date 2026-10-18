# SPDX-License-Identifier: MPL-2.0
"""Exercise archive/metadata failures through the actual C++ importer."""
import pathlib
import subprocess
import sys
import tempfile
import uuid
import xml.etree.ElementTree as ET
import zipfile

executable, fixture = sys.argv[1:]
with zipfile.ZipFile(fixture) as archive:
    original = {name: archive.read(name) for name in archive.namelist()}


def xml_variant(change):
    files = dict(original)
    root = ET.fromstring(files["modelDescription.xml"])
    change(root)
    files["modelDescription.xml"] = ET.tostring(root)
    return files


cases = [
    ("missing_xml", {"other.txt": b"no metadata"}, "invalid modelDescription.xml"),
    (
        "version",
        xml_variant(lambda root: root.set("fmiVersion", "1.0")),
        "requires an FMI 3.0",
    ),
    (
        "missing_binary",
        {k: v for k, v in original.items() if not k.startswith("binaries/")},
        "missing host binary",
    ),
    (
        "wrong_input",
        xml_variant(
            lambda root: root.find("ModelVariables/Float64").set("name", "other")
        ),
        "unexpected input",
    ),
    (
        "array",
        xml_variant(
            lambda root: ET.SubElement(
                root.find("ModelVariables/Float64"), "Dimension", start="1"
            )
        ),
        "scalar Float64",
    ),
    (
        "dtd",
        {
            "modelDescription.xml": b'<!DOCTYPE fmiModelDescription [<!ENTITY test "x">]><fmiModelDescription/>'
        },
        "DTDs are unsupported",
    ),
]
marker = pathlib.Path(tempfile.gettempdir()) / ("dpsim-fmi-escape-" + uuid.uuid4().hex)
cases.append(("traversal", {"../" + marker.name: b"escape"}, "unsafe archive path"))
cases.append(("absolute", {str(marker): b"escape"}, "unsafe archive path"))
with tempfile.TemporaryDirectory() as directory:
    for name, files, message in cases:
        path = pathlib.Path(directory) / (name + ".fmu")
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
            for entry, content in files.items():
                archive.writestr(entry, content)
        subprocess.run([executable, str(path), "reject:" + message], check=True)
assert not marker.exists(), "ZIP extraction escaped its temporary directory"
print("PASS: eight malformed archive/metadata cases")
