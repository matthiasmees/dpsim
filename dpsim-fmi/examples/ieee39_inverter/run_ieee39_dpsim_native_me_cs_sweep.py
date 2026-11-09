#!/usr/bin/env python3
"""Run and compare native DPsim, FMI ME, and FMI CS IEEE39 inverter cases."""
from __future__ import annotations
import argparse, csv, json, math, subprocess, time
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CASE = Path(__file__).resolve().parent
BUILD = ROOT / "build-fmi"
EXE = BUILD / "dpsim-fmi" / "FMI_IEEE39_AvgInverter"
FMU = BUILD / "ieee39-inverter-fmu" / "IEEE39AvgInverter.fmu"
OUT = BUILD / "ieee39-inverter" / "dpsim-native-me-cs-10s-nosubcomponents"
DEFAULT_STEPS = (10, 25, 50, 100, 500, 1000, 5000)
STOP = 10.0
EVENT = 5.0
ROLL = 0.1
MAX_POINTS = 12000


def run(mode: int | str, step_us: int, out: Path) -> dict:
    dt = step_us * 1e-6
    directory = out / f"dt_{step_us}us" / str(mode)
    directory.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    p = subprocess.run(
        [str(EXE), str(mode), str(FMU), str(STOP), str(dt), str(EVENT)],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    wall = time.perf_counter() - started
    (directory / "runner.log").write_text(p.stdout)
    if p.returncode:
        raise RuntimeError(f"{step_us} us {mode} failed:\n{p.stdout[-4000:]}")
    line = next(
        (x for x in p.stdout.splitlines() if x.startswith("IEEE39_Inverter_")), None
    )
    if not line:
        raise RuntimeError(f"Missing timing line for {mode} {step_us} us")
    fields = dict(x.split("=", 1) for x in line.split() if "=" in x)
    raw = directory / "logs" / f"IEEE39_Inverter_{mode}" / f"IEEE39_Inverter_{mode}.csv"
    if not raw.is_file():
        raise FileNotFoundError(raw)
    with raw.open("r", encoding="utf-8") as stream:
        sample_count = sum(1 for _ in stream) - 1
    expected = round(STOP / dt) + 1
    if sample_count not in (expected, expected - 1):
        raise ValueError(
            f"{mode} {step_us}us produced {sample_count} rows; expected {expected} or {expected-1}"
        )
    for pattern in (
        "*_TaskTimes*.csv",
        "*_SimulationStepTimes.csv",
        "*_SolverPhaseTimes.csv",
        "*_SolverMatrixStats.csv",
    ):
        for profile in (directory / "logs").rglob(pattern):
            profile.unlink(missing_ok=True)
    newton = {}
    if mode in ("me", "cs"):
        with raw.open(newline="", encoding="utf-8") as stream:
            rows = csv.DictReader(stream, skipinitialspace=True)
            iterations = []
            residuals = []
            for row in rows:
                iterations.append(float(row["newton_iterations"]))
                residuals.append(float(row["newton_residual"]))
        newton = {
            "mean_iterations": sum(iterations) / len(iterations),
            "max_iterations": max(iterations),
            "max_scaled_residual": max(residuals),
        }
    return {
        "mode": str(mode),
        "step_us": step_us,
        "step_s": dt,
        "duration_s": STOP,
        "load_step_s": EVENT,
        "samples_expected": sample_count,
        "initialization": "same IEEE39 load-flow initialization; native inverter filter state set to the FMU initial P/Q current",
        "process_wall_s": wall,
        "init_s": float(fields["initialization_s"]),
        "simulation_s": float(fields["simulation_s"]),
        "solver_wall_s": float(fields["runtime_s"]),
        "newton_metrics": newton,
        "raw_path": str(raw),
        "directory": str(directory),
    }


def vals(row):
    return [float(row[f"N30.V_{j}"]) for j in range(3)] + [
        float(row[f"inverter.i_{j}"]) for j in range(3)
    ]


def compare(reference: dict, candidate: dict, mode: str) -> dict:
    ref = Path(reference["raw_path"])
    src = Path(candidate["raw_path"])
    directory = Path(candidate["directory"])
    dt = candidate["step_s"]
    n = candidate["samples_expected"]
    stride = max(1, math.ceil(n / MAX_POINTS))
    window = max(1, round(ROLL / dt))
    ring = deque()
    sv = si = rv = ri = 0.0
    full = {"ve": 0.0, "vr": 0.0, "ie": 0.0, "ir": 0.0, "n": 0}
    sections = {
        k: {"ve": 0.0, "vr": 0.0, "ie": 0.0, "ir": 0.0, "n": 0}
        for k in ("pre_load_step", "post_load_step")
    }
    first = None
    count = 0
    errout = directory / f"{mode}_rolling_error.csv"
    traceout = directory / f"{mode}_trace.csv"
    eventout = directory / f"{mode}_load_step_zoom.csv"
    native_trace = Path(reference["directory"]) / "native_trace.csv"
    write_native = not native_trace.exists()
    with (
        ref.open(newline="") as rf,
        src.open(newline="") as cf,
        errout.open("w", newline="") as ef,
        traceout.open("w", newline="") as tf,
        eventout.open("w", newline="") as zf,
        (
            native_trace.open("w", newline="")
            if write_native
            else open("/dev/null", "w", newline="")
        ) as nf,
    ):
        rr = csv.DictReader(rf, skipinitialspace=True)
        cr = csv.DictReader(cf, skipinitialspace=True)
        ew = csv.writer(ef)
        tw = csv.writer(tf)
        zw = csv.writer(zf)
        nw = csv.writer(nf)
        ew.writerow(
            [
                "time_s",
                "pcc_voltage_rms_error_percent_100ms",
                "inverter_current_rms_error_percent_100ms",
            ]
        )
        tw.writerow(
            [
                "time_s",
                "pcc_v_a_V",
                "pcc_v_b_V",
                "pcc_v_c_V",
                "current_a_A",
                "current_b_A",
                "current_c_A",
            ]
        )
        if write_native:
            nw.writerow(
                [
                    "time_s",
                    "pcc_v_a_V",
                    "pcc_v_b_V",
                    "pcc_v_c_V",
                    "current_a_A",
                    "current_b_A",
                    "current_c_A",
                ]
            )
        zw.writerow(
            ["time_s", "native_v_a_V", "native_i_a_A", f"{mode}_v_a_V", f"{mode}_i_a_A"]
        )
        for idx, (r, c) in enumerate(zip(rr, cr)):
            t = float(c["time"])
            tr = float(r["time"])
            a = vals(r)
            b = vals(c)
            if abs(t - tr) > max(1e-12, dt * 1e-8):
                raise ValueError(f"{mode} time grid mismatch at {idx}: {t} vs {tr}")
            if not all(math.isfinite(x) for x in (*a, *b)):
                raise ValueError(f"nonfinite value in {mode} at {idx}")
            if first is None:
                first = {
                    "native_v_abc_V": a[:3],
                    "native_i_abc_A": a[3:],
                    "candidate_v_abc_V": b[:3],
                    "candidate_i_abc_A": b[3:],
                }
            ve = sum((b[j] - a[j]) ** 2 for j in range(3))
            vr = sum(a[j] ** 2 for j in range(3))
            ie = sum((b[j + 3] - a[j + 3]) ** 2 for j in range(3))
            ir = sum(a[j + 3] ** 2 for j in range(3))
            for d in (
                full,
                sections["pre_load_step" if t < EVENT else "post_load_step"],
            ):
                d["ve"] += ve
                d["vr"] += vr
                d["ie"] += ie
                d["ir"] += ir
                d["n"] += 1
            ring.append((ve, vr, ie, ir))
            sv += ve
            rv += vr
            si += ie
            ri += ir
            if len(ring) > window:
                old = ring.popleft()
                sv -= old[0]
                rv -= old[1]
                si -= old[2]
                ri -= old[3]
            if idx % stride == 0 or idx == n - 1:
                ew.writerow(
                    [
                        t,
                        100 * math.sqrt(max(sv, 0) / rv) if rv else 0,
                        100 * math.sqrt(max(si, 0) / ri) if ri else 0,
                    ]
                )
                tw.writerow([t, *b])
                if write_native:
                    nw.writerow([t, *a])
            if abs(t - EVENT) <= 0.025 and idx % max(1, round(0.0001 / dt)) == 0:
                zw.writerow([t, a[0], a[3], b[0], b[3]])
            count += 1
    if count != n:
        raise ValueError(f'{mode} {candidate["step_us"]}us has {count}, expected {n}')

    def score(d):
        return {
            "samples": d["n"],
            "voltage_rms_error_V": math.sqrt(d["ve"] / (3 * d["n"])),
            "voltage_rms_error_percent": (
                100 * math.sqrt(d["ve"] / d["vr"]) if d["vr"] else None
            ),
            "current_rms_error_A": math.sqrt(d["ie"] / (3 * d["n"])),
            "current_rms_error_percent": (
                100 * math.sqrt(d["ie"] / d["ir"]) if d["ir"] else None
            ),
        }

    return {
        "mode": mode.upper(),
        "samples": count,
        "initial": first,
        "full_run": score(full),
        "windows": {k: score(v) for k, v in sections.items()},
        "rolling_error_trace": str(errout),
        "compact_trace": str(traceout),
        "load_step_zoom_trace": str(eventout),
    }


def main():
    global STOP, EVENT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steps-us", type=int, nargs="+", default=list(DEFAULT_STEPS))
    ap.add_argument("--duration-s", type=float, default=STOP)
    ap.add_argument("--load-step-s", type=float, default=EVENT)
    ap.add_argument("--output-name", default=OUT.name)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    STOP, EVENT = args.duration_s, args.load_step_s
    if STOP <= 0 or EVENT <= 0 or EVENT >= STOP:
        ap.error("duration and load-step times must satisfy 0 < load-step < duration")
    out = BUILD / "ieee39-inverter" / args.output_name
    out.mkdir(parents=True, exist_ok=True)
    if not EXE.is_file():
        subprocess.run(
            [
                "cmake",
                "--build",
                str(BUILD),
                "--target",
                "FMI_IEEE39_AvgInverter",
                "-j4",
            ],
            cwd=ROOT,
            check=True,
        )
    if not FMU.is_file():
        subprocess.run(
            ["bash", str(CASE / "sources" / "build_fmu.sh"), str(FMU.parent)],
            cwd=ROOT,
            check=True,
        )
    result = {
        "configuration": {
            "grid": "IEEE39_Ph3 EMT network, inverter at bus 30",
            "dynamic_network_components": "Conventional PiLines, RX loads, and transformers with subcomponents disabled; conventional VBR SG4 generators",
            "duration_s": STOP,
            "load_step_s": EVENT,
            "load_step_MW": 100,
            "steps_us": args.steps_us,
            "methods": ["native", "me", "cs"],
            "initialization": "same solved IEEE39 powerflow; native and FMU initial inverter current match at t=0",
            "rolling_error_window_s": ROLL,
            "error_reference": "native DPsim at matching timestep",
            "CS_coupling": "FMI state rollback/Newton with four midpoint voltage replay substeps",
            "ME_newton_stopping": "residual-based; current source default ceiling 20",
        },
        "cases": {},
    }
    metrics_path = out / "metrics.json"
    if args.resume and metrics_path.is_file():
        prior = json.loads(metrics_path.read_text())
        result["configuration"] = prior.get("configuration", result["configuration"])
        result["cases"].update(prior.get("cases", {}))
    for us in args.steps_us:
        print(f"\n=== IEEE39 {us} us ===", flush=True)
        folder = out / f"dt_{us}us"
        folder.mkdir(parents=True, exist_ok=True)
        case = {}
        if args.resume and str(us) in result["cases"]:
            done = result["cases"][str(us)]["comparisons"]
            if all(
                Path(done[m][key]).is_file()
                for m in ("me", "cs")
                for key in (
                    "rolling_error_trace",
                    "compact_trace",
                    "load_step_zoom_trace",
                )
            ):
                print(
                    f"{us} us: reusing completed metrics and compact traces", flush=True
                )
                continue
        for mode in ("native", "me", "cs"):
            recfile = folder / mode / "run.json"
            if args.resume and recfile.is_file():
                rec = json.loads(recfile.read_text())
                if Path(rec["raw_path"]).is_file():
                    with Path(rec["raw_path"]).open("r", encoding="utf-8") as stream:
                        rec["samples_expected"] = sum(1 for _ in stream) - 1
                    case[mode] = rec
                    print(
                        f'{mode}: reusing completed raw trace ({rec["samples_expected"]:,} samples)',
                        flush=True,
                    )
                    continue
            rec = run(mode, us, out)
            recfile.write_text(json.dumps(rec, indent=2) + "\n")
            case[mode] = rec
            print(
                f"{mode:6} {rec['samples_expected']:,} samples; solve {rec['solver_wall_s']:.2f}s; process {rec['process_wall_s']:.2f}s",
                flush=True,
            )
        comparisons = {}
        for mode in ("me", "cs"):
            comparisons[mode] = compare(case["native"], case[mode], mode)
            s = comparisons[mode]["full_run"]
            print(
                f"{mode.upper()} vs native: V {s['voltage_rms_error_percent']:.5g}% RMS; I {s['current_rms_error_percent']:.5g}% RMS",
                flush=True,
            )
            Path(case[mode]["raw_path"]).unlink(missing_ok=True)
        Path(case["native"]["raw_path"]).unlink(missing_ok=True)
        result["cases"][str(us)] = {
            "runtimes": {
                m: {k: v for k, v in case[m].items() if k not in ("raw_path",)}
                for m in case
            },
            "comparisons": comparisons,
        }
        (out / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"\nMetrics: {out}/metrics.json", flush=True)


if __name__ == "__main__":
    main()
