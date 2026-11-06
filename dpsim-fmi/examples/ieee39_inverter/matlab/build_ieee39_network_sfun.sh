#!/usr/bin/env bash
set -euo pipefail

MATLAB_ROOT="${MATLAB_ROOT:-/usr/local/MATLAB/R2025a}"
REPO_ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
MATLAB_DIR="$(cd "$(dirname "$0")" && pwd)"
NETWORK_LIB="${REPO_ROOT}/build-fmi/dpsim-fmi"
MEX_EXT="mexa64"

g++ -std=c++17 -O2 -fPIC -shared -DMATLAB_MEX_FILE -DMX_COMPAT_64 \
  -I"${MATLAB_ROOT}/extern/include" \
  -I"${MATLAB_ROOT}/simulink/include" \
  "${MATLAB_DIR}/ieee39_network_sfun.cpp" \
  -L"${NETWORK_LIB}" -lieee39_simulink_network \
  -Wl,-rpath,"${NETWORK_LIB}" \
  -o "${MATLAB_DIR}/ieee39_network_sfun.${MEX_EXT}"
