#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$case_dir/../../../../" && pwd)"
build_dir="${1:-$repo/build-fmi-cosim/dual-interface-rl-fmu}"
staging="$build_dir/ThreePhaseRL.fmutmp"
mkdir -p "$staging/binaries/x86_64-linux"
cp "$case_dir/sources/modelDescription.xml" "$staging/modelDescription.xml"
c++ -std=c++17 -O2 -fPIC -shared \
  -I"$repo/dpsim-fmi/third_party/fmi3" \
  "$case_dir/sources/ThreePhaseRL.cpp" \
  -o "$staging/binaries/x86_64-linux/ThreePhaseRL.so"
rm -f "$build_dir/ThreePhaseRL.fmu"
(cd "$staging" && zip -qr "$build_dir/ThreePhaseRL.fmu" .)
printf '%s\n' "$build_dir/ThreePhaseRL.fmu"
