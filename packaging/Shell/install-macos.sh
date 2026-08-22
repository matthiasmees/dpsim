#!/usr/bin/env bash
#
# SPDX-FileCopyrightText: 2026 The DPsim Authors
# SPDX-License-Identifier: MPL-2.0
#
# Native DPsim development setup for macOS.
#
# This script reproduces the tested native macOS development environment:
#
#   - Homebrew build dependencies
#   - Eigen 3 from Homebrew
#   - Graphviz from Homebrew
#   - SuiteSparse and spdlog fetched by DPsim
#   - repository-local Python 3.14 virtual environment: dpsim-python
#   - pybind11 Python bindings
#   - NumPy / pandas / SciPy / matplotlib
#   - JupyterLab / ipykernel / notebook support
#   - VILLAS dataprocessing Python package
#   - pytest tooling
#   - pre-commit tooling
#   - modern Homebrew Ruby for the markdownlint pre-commit hook
#   - native CMake/Ninja build including C++ examples
#   - independent `pip install .` package build
#   - Jupyter kernel registration
#
# Usage:
#
#   chmod +x packaging/Shell/install-macos.sh
#   ./packaging/Shell/install-macos.sh
#
# After setup:
#
#   source dpsim-python/bin/activate
#   cmake --build build --parallel "$(sysctl -n hw.ncpu)"
#
# Optional environment variables:
#
#   BUILD_DIR=build
#   BUILD_TYPE=Release
#   JOBS=<number of parallel build jobs>
#   VENV_DIR=dpsim-python
#

set -euo pipefail

log() {
	printf '\n==> %s\n' "$*"
}

die() {
	printf '\nERROR: %s\n' "$*" >&2
	exit 1
}

# ---------------------------------------------------------------------------
# Platform checks
# ---------------------------------------------------------------------------

[[ "$(uname -s)" == "Darwin" ]] || die "This setup script is intended for macOS."

ARCH="$(uname -m)"

case "${ARCH}" in
	arm64|x86_64)
		;;
	*)
		die "Unsupported macOS architecture: ${ARCH}"
		;;
esac

log "Detected macOS architecture: ${ARCH}"

# AppleClang, git and the macOS SDK are provided by the Xcode Command Line
# Tools.
if ! xcode-select -p >/dev/null 2>&1; then
	cat >&2 <<'EOF'
Xcode Command Line Tools are required.

Install them with:

  xcode-select --install

After installation has completed, run this script again.
EOF
	exit 1
fi

# ---------------------------------------------------------------------------
# Homebrew
# ---------------------------------------------------------------------------

if ! command -v brew >/dev/null 2>&1; then
	log "Homebrew not found; installing Homebrew"

	command -v curl >/dev/null 2>&1 \
		|| die "curl is required to install Homebrew."

	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

	# Homebrew's installer does not necessarily update the current shell.
	if [[ -x /opt/homebrew/bin/brew ]]; then
		eval "$(/opt/homebrew/bin/brew shellenv)"
	elif [[ -x /usr/local/bin/brew ]]; then
		eval "$(/usr/local/bin/brew shellenv)"
	else
		die "Homebrew installation completed, but brew could not be found."
	fi
fi

BREW_PREFIX="$(brew --prefix)"

log "Using Homebrew at ${BREW_PREFIX}"

# ---------------------------------------------------------------------------
# Homebrew dependencies
# ---------------------------------------------------------------------------

# Keep compiled DPsim dependencies minimal:
#
#   Eigen 3     -> Homebrew
#   Graphviz    -> Homebrew
#   SuiteSparse -> fetched by DPsim
#   spdlog      -> fetched by DPsim
#
# Python 3.14 is used explicitly because this DPsim revision is tested with
# Python 3.14 and pybind11 >= 3.
#
# Homebrew Ruby is needed because Apple's system Ruby 2.6 is too old for the
# dependency chain used by the markdownlint pre-commit hook.
BREW_PACKAGES=(
	cmake
	ninja
	eigen@3
	graphviz
	python@3.14
	ruby
	pkg-config
)

log "Installing required Homebrew packages"
brew install "${BREW_PACKAGES[@]}"

EIGEN_PREFIX="$(brew --prefix eigen@3)"
GRAPHVIZ_PREFIX="$(brew --prefix graphviz)"
PYTHON_PREFIX="$(brew --prefix python@3.14)"
RUBY_PREFIX="$(brew --prefix ruby)"

PYTHON_BASE="${PYTHON_PREFIX}/bin/python3.14"

[[ -x "${PYTHON_BASE}" ]] \
	|| die "Homebrew Python 3.14 executable not found at ${PYTHON_BASE}."

# Homebrew Ruby is keg-only. Put it ahead of Apple's /usr/bin/ruby for this
# setup process so pre-commit creates its markdownlint environment with a
# sufficiently recent Ruby.
export PATH="${RUBY_PREFIX}/bin:${PATH}"

command -v ruby >/dev/null 2>&1 \
	|| die "Ruby could not be found after installing Homebrew Ruby."

RUBY_VERSION="$(ruby -e 'print RUBY_VERSION')"
RUBY_MAJOR="${RUBY_VERSION%%.*}"
RUBY_REST="${RUBY_VERSION#*.}"
RUBY_MINOR="${RUBY_REST%%.*}"

if (( RUBY_MAJOR < 3 || (RUBY_MAJOR == 3 && RUBY_MINOR < 1) )); then
	die "Ruby >= 3.1 is required for the markdownlint pre-commit hook; found ${RUBY_VERSION}."
fi

log "Dependency versions"
printf '  CMake:    %s\n' "$(cmake --version | head -n1)"
printf '  Ninja:    %s\n' "$(ninja --version)"
printf '  Eigen 3:  %s\n' "$(brew list --versions eigen@3)"
printf '  Graphviz: %s\n' "$(dot -V 2>&1)"
printf '  Python:   %s\n' "$("${PYTHON_BASE}" --version 2>&1)"
printf '  Ruby:     %s\n' "$(ruby --version)"

# ---------------------------------------------------------------------------
# Repository
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel 2>/dev/null)"; then
	:
else
	die "Could not determine the DPsim repository root. Run this script from a cloned DPsim repository."
fi

cd "${REPO_ROOT}"

log "DPsim repository: ${REPO_ROOT}"

# ---------------------------------------------------------------------------
# Verify required portability fixes
# ---------------------------------------------------------------------------

# These checks make failures explicit if this script is used with a DPsim
# revision that does not yet contain the portability changes required by the
# tested macOS build.

grep -q 'GRAPHVIZ_RENDERDATA_USES_SIZE_T' cmake/FindGraphviz.cmake \
	|| die "Graphviz gvRenderData API compatibility detection is missing from cmake/FindGraphviz.cmake."

grep -q 'cgraph.h' cmake/FindGraphviz.cmake \
	|| die "Graphviz-specific header discovery is missing from cmake/FindGraphviz.cmake."

grep -q 'FILESYSTEM_HAS_NATIVE_SUPPORT' cmake/FindFilesystem.cmake \
	|| die "Native std::filesystem detection is missing from cmake/FindFilesystem.cmake."

grep -q 'template lpNorm<Eigen::Infinity>' dpsim/src/MNASolver.cpp \
	|| die "The dependent-template lpNorm fix is missing from dpsim/src/MNASolver.cpp."

grep -q 'DPSIM_PYTHON_WITH_VILLAS' setup.py \
	|| die "The platform-neutral optional VILLAS Python packaging logic is missing from setup.py."

grep -q 'FETCH_SPDLOG=ON' setup.py \
	|| die "The reproducible fetched-spdlog Python packaging configuration is missing from setup.py."

grep -q 'shlex.split' setup.py \
	|| die "Robust CMAKE_ARGS/CMAKE_OPTS parsing is missing from setup.py."

# ---------------------------------------------------------------------------
# Build / virtual-environment settings
# ---------------------------------------------------------------------------

BUILD_DIR="${BUILD_DIR:-build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
JOBS="${JOBS:-$(sysctl -n hw.ncpu)}"
VENV_DIR="${VENV_DIR:-dpsim-python}"

if [[ "${VENV_DIR}" = /* ]]; then
	VENV_PATH="${VENV_DIR}"
else
	VENV_PATH="${REPO_ROOT}/${VENV_DIR}"
fi

# ---------------------------------------------------------------------------
# Python virtual environment
# ---------------------------------------------------------------------------

if [[ ! -x "${VENV_PATH}/bin/python" ]]; then
	log "Creating Python virtual environment: ${VENV_PATH}"
	"${PYTHON_BASE}" -m venv "${VENV_PATH}"
else
	log "Using existing Python virtual environment: ${VENV_PATH}"
fi

PYTHON="${VENV_PATH}/bin/python"
PIP="${PYTHON} -m pip"

# Verify that an existing environment is actually based on Python 3.14.
PYTHON_MM="$("${PYTHON}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

[[ "${PYTHON_MM}" == "3.14" ]] \
	|| die "${VENV_PATH} uses Python ${PYTHON_MM}; Python 3.14 is required. Remove the environment and rerun this script."

log "Upgrading Python packaging tools"
"${PYTHON}" -m pip install --upgrade \
	pip \
	setuptools \
	wheel

log "Installing pybind11 build dependencies"
"${PYTHON}" -m pip install \
	"pybind11>=3.0.0" \
	"pybind11-stubgen>=2.5"

log "Installing scientific Python dependencies"
"${PYTHON}" -m pip install \
	"numpy>=2.0.0" \
	"pandas>=2.0.0" \
	"scipy>=1.10.0" \
	matplotlib

log "Installing Jupyter notebook dependencies"
"${PYTHON}" -m pip install \
	jupyterlab \
	ipykernel \
	ipywidgets \
	nbformat \
	"nbconvert==7.17.1"

log "Installing DPsim notebook/data-processing dependencies"
"${PYTHON}" -m pip install \
	"villas-dataprocessing>=0.2.6"

log "Installing Python test dependencies"
"${PYTHON}" -m pip install \
	pytest \
	pytest-xdist \
	pytest-cov \
	pyyaml

log "Installing pre-commit"
"${PYTHON}" -m pip install pre-commit

PYBIND11_DIR="$("${PYTHON}" -m pybind11 --cmakedir)"

[[ -d "${PYBIND11_DIR}" ]] \
	|| die "pybind11 CMake directory not found: ${PYBIND11_DIR}"

log "Python development environment"
printf '  Virtual env: %s\n' "${VENV_PATH}"
printf '  Python:      %s\n' "$("${PYTHON}" --version)"
printf '  pybind11:    %s\n' "$("${PYTHON}" -c 'import pybind11; print(pybind11.__version__)')"
printf '  pybind dir:  %s\n' "${PYBIND11_DIR}"
printf '  matplotlib:  %s\n' "$("${PYTHON}" -c 'import matplotlib; print(matplotlib.__version__)')"
printf '  JupyterLab:  %s\n' "$("${PYTHON}" -c 'import jupyterlab; print(jupyterlab.__version__)')"

# ---------------------------------------------------------------------------
# CMake dependency paths
# ---------------------------------------------------------------------------

# Put Eigen 3 first deliberately. Do not allow another Eigen installation to
# win package discovery.
CMAKE_PREFIX_PATH_VALUE="${EIGEN_PREFIX};${GRAPHVIZ_PREFIX};${BREW_PREFIX}"

EIGEN3_DIR="${EIGEN_PREFIX}/share/eigen3/cmake"

if [[ ! -d "${EIGEN3_DIR}" ]]; then
	# Keep this robust against possible Homebrew layout changes.
	EIGEN3_DIR="$(find "${EIGEN_PREFIX}" -type d -path '*/eigen3/cmake' -print -quit 2>/dev/null || true)"
fi

[[ -n "${EIGEN3_DIR}" && -d "${EIGEN3_DIR}" ]] \
	|| die "Could not locate Eigen3Config.cmake below ${EIGEN_PREFIX}."

# ---------------------------------------------------------------------------
# Configure native developer build
# ---------------------------------------------------------------------------

log "Configuring native DPsim build"
printf '  Build directory: %s\n' "${BUILD_DIR}"
printf '  Build type:      %s\n' "${BUILD_TYPE}"
printf '  Architecture:    %s\n' "${ARCH}"
printf '  Parallel jobs:   %s\n' "${JOBS}"
printf '  Eigen prefix:    %s\n' "${EIGEN_PREFIX}"
printf '  Graphviz prefix: %s\n' "${GRAPHVIZ_PREFIX}"
printf '  Python:          %s\n' "${PYTHON}"

cmake \
	-S . \
	-B "${BUILD_DIR}" \
	-G Ninja \
	-DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
	-DCMAKE_OSX_ARCHITECTURES="${ARCH}" \
	-DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH_VALUE}" \
	-DEigen3_DIR="${EIGEN3_DIR}" \
	-DGraphviz_ROOT="${GRAPHVIZ_PREFIX}" \
	-DPython3_EXECUTABLE="${PYTHON}" \
	-DPYTHON_EXECUTABLE="${PYTHON}" \
	-Dpybind11_DIR="${PYBIND11_DIR}" \
	-DFETCH_EIGEN=OFF \
	-DFETCH_SUITESPARSE=ON \
	-DFETCH_SPDLOG=ON \
	-DFETCH_PYBIND=OFF \
	-DWITH_PYBIND=ON \
	-DWITH_GRAPHVIZ=ON \
	-DWITH_CIM=OFF \
	-DWITH_VILLAS=OFF \
	-DDPSIM_BUILD_EXAMPLES=ON \
	-DDPSIM_BUILD_DOC=OFF

# Verify that Graphviz discovery did not regress to the broad Homebrew include
# directory, which previously allowed stale headers below /opt/homebrew/include
# to shadow dependencies from the Python virtual environment.
if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
	GRAPHVIZ_INCLUDE_CACHE="$(grep '^GRAPHVIZ_INCLUDE_DIR:PATH=' "${BUILD_DIR}/CMakeCache.txt" | cut -d= -f2- || true)"

	if [[ -n "${GRAPHVIZ_INCLUDE_CACHE}" ]]; then
		case "${GRAPHVIZ_INCLUDE_CACHE}" in
			*/include/graphviz)
				;;
			*)
				die "Graphviz include directory is unexpectedly broad: ${GRAPHVIZ_INCLUDE_CACHE}. Expected a graphviz-specific include directory."
				;;
		esac
	fi
fi

# ---------------------------------------------------------------------------
# Build native developer tree
# ---------------------------------------------------------------------------

log "Building DPsim, Python bindings and C++ examples"
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"

# Make sure the Python module was generated by the normal developer build.
if ! find "${BUILD_DIR}" -maxdepth 2 -type f -name 'dpsimpy*.so' -print -quit | grep -q .; then
	die "The native CMake build completed but dpsimpy was not found below ${BUILD_DIR}."
fi

# ---------------------------------------------------------------------------
# Independent pip package build/install
# ---------------------------------------------------------------------------

# `pip install .` performs its own isolated CMake build. It therefore needs the
# same package-discovery information as the normal developer build.
#
# Clean only setuptools-generated macOS build directories. Do NOT remove the
# normal CMake build tree.
log "Preparing independent pip package build"

rm -rf \
	"${BUILD_DIR}"/temp.macosx-* \
	"${BUILD_DIR}"/lib.macosx-* \
	"${BUILD_DIR}"/bdist.macosx-*

PIP_CMAKE_ARGS="-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_VALUE} -DEigen3_DIR=${EIGEN3_DIR} -DGraphviz_ROOT=${GRAPHVIZ_PREFIX}"

log "Installing DPsim into ${VENV_PATH} with pip"

CMAKE_GENERATOR=Ninja \
ARCHFLAGS="-arch ${ARCH}" \
CMAKE_ARGS="${PIP_CMAKE_ARGS}" \
DPSIM_PYTHON_WITH_VILLAS=0 \
"${PYTHON}" -m pip install . -v

# ---------------------------------------------------------------------------
# Python import verification
# ---------------------------------------------------------------------------

log "Verifying installed Python environment"

"${PYTHON}" - <<'PY'
import sys

import dpsim
import dpsimpy
import ipykernel
import jupyterlab
import matplotlib
import matplotlib.pyplot as plt
import nbconvert
import nbformat
import numpy
import pandas
import scipy

from villas.dataprocessing.readtools import read_timeseries_csv

print("Python:", sys.executable)
print("dpsimpy:", dpsimpy.__file__)
print("dpsim:", dpsim.__file__)
print("NumPy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("SciPy:", scipy.__version__)
print("matplotlib:", matplotlib.__version__)
print("JupyterLab:", jupyterlab.__version__)
print("ipykernel:", ipykernel.__version__)
print("nbformat:", nbformat.__version__)
print("nbconvert:", nbconvert.__version__)
print("VILLAS dataprocessing import: OK")
print("DPsim Python installation: OK")
PY

"${PYTHON}" -m pip check

# ---------------------------------------------------------------------------
# Jupyter kernel
# ---------------------------------------------------------------------------

log "Registering Jupyter kernel: DPsim Python"

KERNEL_DIR="${HOME}/Library/Jupyter/kernels/dpsim-python"

# Remove only this project's previously registered kernel so an old absolute
# interpreter path cannot survive a recreated virtual environment.
rm -rf "${KERNEL_DIR}"

"${PYTHON}" -m ipykernel install \
	--user \
	--name dpsim-python \
	--display-name "DPsim Python"

KERNEL_JSON="${KERNEL_DIR}/kernel.json"

[[ -f "${KERNEL_JSON}" ]] \
	|| die "Jupyter kernelspec was not created at ${KERNEL_JSON}."

if ! grep -q "${VENV_PATH}/bin/python" "${KERNEL_JSON}"; then
	die "The DPsim Jupyter kernelspec does not point to ${VENV_PATH}/bin/python."
fi

# ---------------------------------------------------------------------------
# Pre-commit
# ---------------------------------------------------------------------------

log "Installing pre-commit hooks"

# PATH still contains Homebrew Ruby here. This is required when pre-commit
# creates the Ruby environment for the markdownlint hook.
printf '  Ruby used for hook installation: %s\n' "$(command -v ruby)"
printf '  Ruby version:                    %s\n' "$(ruby --version)"

"${VENV_PATH}/bin/pre-commit" validate-config
"${VENV_PATH}/bin/pre-commit" install --install-hooks

# We intentionally do not run `pre-commit run --all-files` here. Several hooks
# are formatters and may modify source files. Developers can run it explicitly
# after setup and inspect the resulting diff.

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

log "DPsim macOS setup completed successfully"

printf '\nRepository:\n  %s\n' "${REPO_ROOT}"
printf '\nBuild directory:\n  %s/%s\n' "${REPO_ROOT}" "${BUILD_DIR}"
printf '\nVirtual environment:\n  %s\n' "${VENV_PATH}"
printf '\nJupyter kernel:\n  DPsim Python\n'
printf '\nExample executables are located below:\n  %s/%s/dpsim/examples/cxx\n' "${REPO_ROOT}" "${BUILD_DIR}"

cat <<EOF

Activate the environment with:

  source "${VENV_PATH}/bin/activate"

Rebuild DPsim with:

  cmake --build "${BUILD_DIR}" --parallel "\$(sysctl -n hw.ncpu)"

Reinstall the Python package after binding/package changes with:

  CMAKE_GENERATOR=Ninja \\
  ARCHFLAGS="-arch ${ARCH}" \\
  CMAKE_ARGS="-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_VALUE} -DEigen3_DIR=${EIGEN3_DIR} -DGraphviz_ROOT=${GRAPHVIZ_PREFIX}" \\
  DPSIM_PYTHON_WITH_VILLAS=0 \\
  python -m pip install . -v

Run the repository pre-commit checks with:

  pre-commit run --all-files

In VS Code/Jupyter select the kernel:

  DPsim Python

EOF
