# SPDX-License-Identifier: MPL-2.0
"""Check the exported nonlinear SG benchmark through the native FMI 3.0 ABI.

Reads the generated port contract and initial Simulink row, initializes the
controlled-generator FMU, verifies all 23 outputs and its ten-state periodic
operating point. This checks export/initialization, not transient integration.
Uses only the standard library; the host macOS binary must be in the FMU.
"""
from pathlib import Path
import csv
import ctypes as c
import json
import platform
import tempfile
import xml.etree.ElementTree as ET
import zipfile
import math
import sys

if len(sys.argv) != 2:
    raise SystemExit("Usage: python3 validate_sg_fmu.py MATLAB-benchmark-directory")
if sys.platform != "darwin":
    raise SystemExit("This native-export validation helper currently supports macOS.")
folder = Path(sys.argv[1]).resolve()
mapping = json.loads((folder / "port_mapping.json").read_text())
with (folder / "reference_uniform.csv").open() as f:
    expected = next(csv.DictReader(f))
with tempfile.TemporaryDirectory(prefix="dpsim-sg-fmu-check-") as temp:
    root = Path(temp)
    with zipfile.ZipFile(folder / "dpsim_sg_3ph_me.fmu") as z:
        z.extractall(root)
    xml = ET.parse(root / "modelDescription.xml").getroot()
    me = xml.find("ModelExchange")
    arch = "aarch64" if platform.machine() == "arm64" else platform.machine()
    lib = c.CDLL(
        str(
            root
            / "binaries"
            / f"{arch}-darwin"
            / (me.attrib["modelIdentifier"] + ".dylib")
        )
    )
    variables = {
        v.attrib["name"]: int(v.attrib["valueReference"])
        for v in xml.find("ModelVariables")
    }
    logtype = c.CFUNCTYPE(None, c.c_void_p, c.c_int, c.c_char_p, c.c_char_p)
    log = logtype(lambda env, status, category, message: None)

    def fn(name, args, restype=c.c_int):
        f = getattr(lib, "fmi3" + name)
        f.argtypes = args
        f.restype = restype
        return f

    vp = c.c_void_p
    sz = c.c_size_t
    dbl = c.c_double
    boolean = c.c_bool
    instantiate = fn(
        "InstantiateModelExchange",
        [c.c_char_p, c.c_char_p, c.c_char_p, boolean, boolean, vp, logtype],
        vp,
    )
    instance = instantiate(
        b"sg_validation",
        xml.attrib["instantiationToken"].encode(),
        (str(root / "resources") + "/").encode(),
        False,
        False,
        None,
        log,
    )
    assert instance, "FMI instantiation failed"

    def check(status):
        assert status in (0, 1), f"FMI status {status}"

    terminate = fn("Terminate", [vp])
    free = fn("FreeInstance", [vp], None)
    try:
        check(
            fn("EnterInitializationMode", [vp, boolean, dbl, dbl, boolean, dbl])(
                instance, False, 0, 0, False, 0
            )
        )
        names = ["u_a", "u_b", "u_c", "vref_pu", "pmref_pu"]
        refs = (c.c_uint32 * 5)(*[variables[n] for n in names])
        values = (dbl * 5)(
            *mapping["initial_voltages_V"],
            mapping["initial_vref_pu"],
            mapping["initial_pmref_pu"],
        )
        check(
            fn("SetFloat64", [vp, c.POINTER(c.c_uint32), sz, c.POINTER(dbl), sz])(
                instance, refs, 5, values, 5
            )
        )
        check(fn("ExitInitializationMode", [vp])(instance))
        update = fn(
            "UpdateDiscreteStates", [vp] + [c.POINTER(boolean)] * 5 + [c.POINTER(dbl)]
        )
        for _ in range(100):
            flags = [boolean() for _ in range(5)]
            nexttime = dbl()
            check(update(instance, *[c.byref(v) for v in flags], c.byref(nexttime)))
            assert not flags[1].value and not flags[4].value
            if not flags[0].value:
                break
        else:
            raise AssertionError("FMI initial event iteration failed")
        check(fn("EnterContinuousTimeMode", [vp])(instance))
        n = sz()
        ne = sz()
        check(
            fn("GetNumberOfContinuousStates", [vp, c.POINTER(sz)])(instance, c.byref(n))
        )
        check(
            fn("GetNumberOfEventIndicators", [vp, c.POINTER(sz)])(instance, c.byref(ne))
        )
        assert n.value == 10, n.value
        outputs = [
            v.attrib["name"]
            for v in xml.find("ModelVariables")
            if v.attrib.get("causality") == "output"
        ]
        ore = (c.c_uint32 * len(outputs))(*[variables[n] for n in outputs])
        val = (dbl * len(outputs))()
        check(
            fn("GetFloat64", [vp, c.POINTER(c.c_uint32), sz, c.POINTER(dbl), sz])(
                instance, ore, len(outputs), val, len(outputs)
            )
        )
        errors = {
            name: abs(val[k] - float(expected[name])) for k, name in enumerate(outputs)
        }
        assert all(
            errors[name] < 1e-8 * (1 + abs(float(expected[name]))) for name in outputs
        ), errors
        derivative = (dbl * n.value)()
        check(
            fn("GetContinuousStateDerivatives", [vp, c.POINTER(dbl), sz])(
                instance, derivative, n.value
            )
        )
        sorted_abs = sorted(abs(v) for v in derivative)
        assert sorted_abs[-2] < 1e-7 and abs(sorted_abs[-1] - 2 * math.pi * 50) < 1e-7
        report = {
            "continuous_states": n.value,
            "event_indicators": ne.value,
            "outputs_checked": len(outputs),
            "max_initial_output_abs_error": max(errors.values()),
            "initialized_periodic_orbit": True,
            "host_binary": f"{arch}-darwin",
            "initial_output_errors": errors,
        }
        (folder / "fmu_validation.json").write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=2))
        check(terminate(instance))
    finally:
        free(instance)
