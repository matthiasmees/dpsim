#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$case_dir/../../../../" && pwd)"
build_dir="${1:-$repo/build-fmi-cosim/ieee9-inverter-fmu}"
mkdir -p "$build_dir"
build_dir="$(cd "$build_dir" && pwd)"
staging="$build_dir/IEEE9AvgInverter.fmutmp"
mkdir -p "$staging/binaries/x86_64-linux"
cp "$case_dir/modelDescription.xml" "$staging/modelDescription.xml"
c++ -std=c++17 -O2 -fPIC -shared -I"$repo/dpsim-fmi/third_party/fmi3" \
  "$case_dir/IEEE9AvgInverter.cpp" \
  -o "$staging/binaries/x86_64-linux/IEEE9AvgInverter.so"
rm -f "$build_dir/IEEE9AvgInverter.fmu"
(cd "$staging" && zip -qr "$build_dir/IEEE9AvgInverter.fmu" .)
printf '%s\n' "$build_dir/IEEE9AvgInverter.fmu"
