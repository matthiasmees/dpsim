"""Run and plot the two MMC examples for EMT_MMC_MNA_Comparison.ipynb.

The notebook keeps its editable cells short. CSV interpretation and the final
comparison are shared with EMT_MMC_MNA_Comparison.py. All subprocesses use
argument lists and run in new output directories; interrupting a cell stops
its child process as well.
"""

from dataclasses import asdict, dataclass
from datetime import datetime
from html import escape
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time

try:
    import numpy as np
    from IPython.display import HTML, Image, Markdown, display
except ImportError as error:
    raise ImportError(
        "Select the 'Python (DPsim)' kernel at the top of the notebook. "
        "It contains NumPy, Matplotlib, and IPython."
    ) from error

# The CLI module selects Agg. Figures are explicitly saved and displayed as
# images here, so no particular interactive Matplotlib backend is required.
COMPARISON_SCRIPT = Path(__file__).with_name("EMT_MMC_MNA_Comparison.py")
_spec = importlib.util.spec_from_file_location("mmc_csv_comparison", COMPARISON_SCRIPT)
_comparison = importlib.util.module_from_spec(_spec)
try:
    _spec.loader.exec_module(_comparison)
except ImportError as error:
    raise ImportError("Select the 'Python (DPsim)' notebook kernel.") from error
import matplotlib.pyplot as plt

PROGRAMS = {
    "ssn": "EMT_SSN_MMC_Matlab_P2P_Results",
    "mna": "EMT_MMC_MNA_Matlab_P2P_Results",
}
LABELS = {"ssn": "SSN MMC", "mna": "MMC MNA"}


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _duration(seconds):
    minutes, seconds = divmod(seconds, 60)
    return f"{int(minutes)} min {seconds:04.1f} s"


def _stop_process(process):
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.wait()


def _native_environment():
    """Keep the notebook's Python libraries out of the native C++ processes.

    The installed kernel adds its Conda lib directory to LD_LIBRARY_PATH.
    Those libraries can conflict with the system toolchain used for DPSim.
    Preserve other library paths and leave the kernel environment unchanged.
    """
    environment = os.environ.copy()
    python_lib = (Path(sys.prefix) / "lib").resolve()
    for key in ("LD_LIBRARY_PATH", "LIBRARY_PATH"):
        if key not in environment:
            continue
        paths = [part for part in environment[key].split(os.pathsep)
                 if Path(part or ".").resolve() != python_lib]
        if paths:
            environment[key] = os.pathsep.join(paths)
        else:
            environment.pop(key)
    return environment


def _run_process(command, cwd, log_path, label, *, environment=None):
    """Show elapsed time and the last log line, including during silent periods."""
    started = time.perf_counter()
    print(f"{label} is starting. Log: {log_path}", flush=True)
    status = display(HTML(f"<b>{escape(label)}:</b> starting …"), display_id=True)
    latest = "Initializing …"
    process = None

    def update(state):
        text = (f"<b>{escape(label)} — {escape(state)}</b><br>"
                f"Elapsed wall time: {_duration(time.perf_counter() - started)}"
                f"<br><small>{escape(latest[-400:])}</small>")
        if status is not None:
            status.update(HTML(text))

    with Path(log_path).open("w") as log:
        try:
            process = subprocess.Popen(
                [str(part) for part in command], cwd=cwd, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=(os.name == "posix"),
                env=environment,
            )
            with Path(log_path).open() as reader:
                next_update = 0.0
                while True:
                    lines = reader.read().splitlines()
                    if lines:
                        latest = lines[-1]
                    if time.perf_counter() >= next_update:
                        update("running")
                        next_update = time.perf_counter() + 2
                    if process.poll() is not None:
                        lines = reader.read().splitlines()
                        if lines:
                            latest = lines[-1]
                        break
                    time.sleep(0.2)
            elapsed = time.perf_counter() - started
            update("complete" if process.returncode == 0 else "failed")
        except BaseException:
            if process is not None:
                _stop_process(process)
            update("interrupted")
            raise
    if process.returncode:
        tail = Path(log_path).read_text(errors="replace")[-2500:]
        raise RuntimeError(
            f"{label} exited with code {process.returncode}.\n"
            f"Log: {log_path}\n\n{tail}"
        )
    print(f"{label}: completed in {_duration(elapsed)}.", flush=True)
    return elapsed


@dataclass(frozen=True)
class RunSettings:
    time_step: float
    final_time: float
    theta: float
    log_every: int


@dataclass(frozen=True)
class SimulationRun:
    kind: str
    directory: Path
    csv: Path
    elapsed_seconds: float
    settings: RunSettings
    csv_sha256: str
    binary_sha256: str


class MMCExperiment:
    """One immutable configuration and the latest successful run of each MMC."""

    def __init__(self, repository, *, time_step=20e-6, final_time=10.0,
                 theta=0.5, log_interval=200e-6):
        self.repository = Path(repository).resolve()
        for name, value in (("Time step", time_step), ("Final time", final_time),
                            ("Theta", theta), ("Logging interval", log_interval)):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        if not 0.5 <= theta <= 1:
            raise ValueError("Theta must be between 0.5 and 1.")
        if final_time < time_step:
            raise ValueError("The final time must cover at least one EMT step.")
        stride = round(log_interval / time_step)
        if stride < 1 or not math.isclose(stride * time_step, log_interval,
                                         rel_tol=1e-9, abs_tol=1e-15):
            raise ValueError(
                "The logging interval must be an integer multiple of the "
                "EMT step, e.g. 200e-6 for 20e-6 or 10e-6."
            )
        self._settings = RunSettings(float(time_step), float(final_time), float(theta), stride)
        self.build_dir = self.repository / "build-mmc-mna"
        self.bin_dir = self.build_dir / "dpsim/examples/cxx"
        parent = self.repository / "outputs/mmc_mna/notebook_runs"
        parent.mkdir(parents=True, exist_ok=True)
        prefix = datetime.now().strftime("%Y%m%d_%H%M%S_")
        self.directory = Path(tempfile.mkdtemp(prefix=prefix, dir=parent))
        self.runs = {}
        self.last_comparison = None
        self._write_json(self.directory / "settings.json", asdict(self.settings))

    @property
    def settings(self):
        return self._settings

    @staticmethod
    def _write_json(path, values):
        Path(path).write_text(json.dumps(values, indent=2, ensure_ascii=False) + "\n")

    def prepare(self, rebuild=True):
        missing = [self.bin_dir / name for name in PROGRAMS.values()
                   if not (self.bin_dir / name).is_file()]
        if rebuild or missing:
            if not (self.build_dir / "CMakeCache.txt").is_file():
                raise FileNotFoundError(
                    f"The configured build directory is missing: {self.build_dir}\n"
                    "See EMT_MMC_MNA_Comparison.md for configuration instructions."
                )
            if not shutil.which("cmake"):
                raise FileNotFoundError("CMake is not available in this environment.")
            _run_process(
                ["cmake", "--build", self.build_dir, "--target", *PROGRAMS.values(), "-j", "4"],
                self.repository, self.directory / "build.log", "Building examples",
                environment=_native_environment(),
            )
        for name in PROGRAMS.values():
            executable = self.bin_dir / name
            if not executable.is_file() or not os.access(executable, os.X_OK):
                raise FileNotFoundError(f"Example executable is missing: {executable}")
        s = self.settings
        display(Markdown(
            f"**Ready:** {s.final_time:g} s simulation time, "
            f"{s.time_step * 1e6:g} µs EMT step, theta {s.theta:g}; "
            f"Logging every {s.log_every} steps "
            f"({s.log_every * s.time_step * 1e6:g} µs).\n\n"
            f"**Output directory:** `{self.directory}`\n\n"
            "Run the SSN cell next."
        ))

    def run_and_plot(self, kind):
        if kind not in PROGRAMS:
            raise ValueError("The model must be 'ssn' or 'mna'.")
        # Re-running/interrupting a cell must never leave its older run selected.
        self.runs.pop(kind, None)
        self.last_comparison = None
        directory = Path(tempfile.mkdtemp(prefix=f"{kind}_", dir=self.directory))
        program = PROGRAMS[kind]
        executable = self.bin_dir / program
        settings = self.settings
        manifest = {
            "model": kind, "settings": asdict(settings), "status": "running",
            "started_at": datetime.now().astimezone().isoformat(),
            "binary": str(executable), "binary_sha256": _sha256(executable),
        }
        self._write_json(directory / "run.json", manifest)
        try:
            elapsed = _run_process(
                [executable, f"{settings.time_step:.17g}", f"{settings.final_time:.17g}",
                 f"{settings.theta:.17g}", str(settings.log_every)],
                directory, directory / "run.log", LABELS[kind],
                environment=_native_environment(),
            )
            csv = directory / "logs" / program / f"{program}.csv"
            data = _comparison.read_csv(csv, settings.time_step)
            if abs(float(data["time"][-1]) - settings.final_time) > 1.1 * settings.time_step:
                raise ValueError("The CSV does not reach the configured simulation end time.")
            for station in _comparison.STATIONS:
                for name, *_ in _comparison.SIGNALS:
                    if f"{station}.{name}" not in data:
                        raise ValueError(f"Required signal is missing: {station}.{name}")
            result = SimulationRun(kind, directory, csv, elapsed, settings,
                                   _sha256(csv), manifest["binary_sha256"])
            manifest.update(status="complete", elapsed_seconds=elapsed, csv=str(csv),
                            csv_sha256=result.csv_sha256, samples=len(data["time"]),
                            physical_end_time=float(data["time"][-1]))
        except BaseException as error:
            manifest.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                            error=f"{type(error).__name__}: {error}")
            self._write_json(directory / "run.json", manifest)
            raise
        self._write_json(directory / "run.json", manifest)
        self.runs[kind] = result
        print(f"CSV checked: {len(data['time']):,} samples. Creating individual plots …", flush=True)
        self._plot(result, data)
        return result

    def _plot(self, result, data):
        output = result.directory / "plots"
        output.mkdir(exist_ok=True)
        time_values = data["time"]
        color = "#1864ab" if result.kind == "ssn" else "#d9480f"
        style = {"font.size": 11, "axes.grid": True, "grid.alpha": 0.22}
        with plt.rc_context(style):
            for station in ("inverter", "rectifier"):
                fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True,
                                         layout="constrained")
                for ax, (name, label, unit, divisor, _) in zip(axes.flat, _comparison.SIGNALS):
                    ax.plot(time_values, data[f"{station}.{name}"] / divisor,
                            color=color, lw=1.2, label=LABELS[result.kind])
                    if name == "p_ac" and station == "inverter":
                        ax.plot(time_values, data["reference.p_inverter"] / divisor,
                                color="0.35", ls=":", lw=1, label="Setpoint")
                    ax.set(title=label, ylabel=unit)
                    for event in (4, 6, 8):
                        if event <= time_values[-1]:
                            ax.axvline(event, color="0.6", lw=0.6, ls=":")
                axes[0, 0].legend(fontsize=9)
                for ax in axes[-1]:
                    ax.set_xlabel("Time [s]")
                fig.suptitle(f"{LABELS[result.kind]} – {station.capitalize()}")
                self._save_display(fig, output, f"{station}_overview")
            mask = time_values >= max(0, time_values[-1] - 0.04)
            fig, axes = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True,
                                     layout="constrained")
            waveforms = (
                ("ac.rectifier_terminal.v_abc", "Rectifier: AC terminal voltage", "kV"),
                ("ac.inverter_terminal.v_abc", "Inverter: AC terminal voltage", "kV"),
                ("ac.right_grid.i_abc", "Right grid component: current", "kA"),
            )
            for ax, (prefix, title, unit) in zip(axes, waveforms):
                for phase, phase_color in zip(range(3), ("#1864ab", "#2b8a3e", "#c92a2a")):
                    ax.plot(time_values[mask], data[f"{prefix}_{phase}"][mask] / 1000,
                            color=phase_color, label=f"Phase {'ABC'[phase]}", lw=1.2)
                ax.set(title=title, ylabel=unit)
            axes[0].legend(ncol=3, fontsize=9)
            axes[-1].set_xlabel("Time [s]")
            fig.suptitle(f"{LABELS[result.kind]} – AC instantaneous values, final 40 ms")
            self._save_display(fig, output, "ac_waveforms")
        print(f"Individual plots saved as PNG and PDF: {output}", flush=True)

    @staticmethod
    def _save_display(fig, output, name):
        try:
            fig.savefig(output / f"{name}.png", dpi=120)
            fig.savefig(output / f"{name}.pdf")
        finally:
            plt.close(fig)
        display(Image(filename=str(output / f"{name}.png"), width=1000))

    def plot_run(self, kind):
        """Redisplay an existing successful run without repeating the simulation."""
        result = self.runs.get(kind)
        if result is None:
            raise RuntimeError(f"Run the {kind.upper()} cell first.")
        if _sha256(result.csv) != result.csv_sha256:
            raise RuntimeError("The CSV has changed since the run. Repeat the simulation.")
        self._plot(result, _comparison.read_csv(result.csv, result.settings.time_step))

    def compare(self):
        missing = [kind.upper() for kind in PROGRAMS if kind not in self.runs]
        if missing:
            raise RuntimeError(
                "This experiment is missing a successful run: " + ", ".join(missing)
                + ". Execute the corresponding run cells first."
            )
        for result in self.runs.values():
            if result.settings != self.settings:
                raise RuntimeError("The runs use different settings.")
            if _sha256(result.csv) != result.csv_sha256:
                raise RuntimeError(f"The CSV was modified after the run: {result.csv}")
        ssn, mna = self.runs["ssn"], self.runs["mna"]
        output = Path(tempfile.mkdtemp(prefix="comparison_", dir=self.directory))
        self.last_comparison = None
        _run_process(
            [sys.executable, COMPARISON_SCRIPT, ssn.csv, mna.csv,
             "--reference-dt", f"{ssn.settings.time_step:.17g}",
             "--mna-dt", f"{mna.settings.time_step:.17g}", "--output", output],
            self.repository, output / "comparison.log", "Comparing SSN and MNA",
        )
        summary = json.loads((output / "metrics.json").read_text())
        worst = max(summary["metrics"], key=lambda item: item["max_error_percent"])
        rows = ["| Run | Measured wall time |", "|---|---:|",
                f"| SSN | {_duration(ssn.elapsed_seconds)} |",
                f"| MNA | {_duration(mna.elapsed_seconds)} |"]
        display(Markdown("\n".join(rows)))
        display(Markdown(
            f"**Largest normalized error:** `{worst['signal']}` = "
            f"**{worst['max_error_percent']:.6g} %**.\n\n"
            "Percentages use fixed engineering bases "
            "(e.g. 1 GW and 640 kV), rather than each instantaneous value. "
            "The difference is always **MNA minus SSN**."
        ))
        table = ["| Signal | Maximum error | RMS error | Maximum / base |",
                 "|---|---:|---:|---:|"]
        for item in summary["metrics"]:
            table.append(
                f"| {item['signal']} | {item['max_error']:.6g} {item['unit']} "
                f"| {item['rms_error']:.6g} {item['unit']} "
                f"| {item['max_error_percent']:.6g} % |"
            )
        display(Markdown("\n".join(table)))
        captions = {
            "inverter_overview": "Inverter: SSN and MNA overlaid",
            "rectifier_overview": "Rectifier: SSN and MNA overlaid",
            "differences": "Differences: MNA minus SSN",
            "step_4s": "Setpoint step at 4 s",
            "step_6s": "Setpoint step at 6 s",
            "step_8s": "Setpoint step at 8 s",
            "ac_waveforms": "Three-phase AC waveforms",
        }
        for name, caption in captions.items():
            image_path = output / f"{name}.png"
            if image_path.is_file():
                display(Markdown(f"### {caption}"))
                display(Image(filename=str(image_path), width=1000))
        self.last_comparison = output
        print(f"Comparison complete. PNG, PDF, report.md, and metrics: {output}", flush=True)
        return summary
