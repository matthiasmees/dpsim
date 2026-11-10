#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$case_dir/../../../../" && pwd)"
if [[ "$(uname -s)" == "Darwin" && -z "${SDKROOT:-}" ]]; then
  for sdk in /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX*.sdk \
             /Volumes/SSDMatthias/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX*.sdk; do
    if [[ -d "$sdk" ]]; then SDKROOT="$sdk"; export SDKROOT; break; fi
  done
fi
build_dir="${1:-$repo/build-fmi-cosim/ieee39-inverter-fmu}"
if [[ "$build_dir" != /* ]]; then
  build_dir="$repo/$build_dir"
fi
mkdir -p "$build_dir"
build_dir="$(cd "$build_dir" && pwd)"
staging="$build_dir/IEEE39AvgInverter.fmutmp"
rm -rf "$staging"
if [[ "$(uname -s)" == "Darwin" ]]; then
  platform_dir="aarch64-darwin"
  binary_name="IEEE39AvgInverter.dylib"
  linker_flags=(-dynamiclib)
else
  platform_dir="x86_64-linux"
  binary_name="IEEE39AvgInverter.so"
  linker_flags=(-fPIC -shared)
fi
mkdir -p "$staging/binaries/$platform_dir"
cp "$case_dir/modelDescription.xml" "$staging/modelDescription.xml"
c++ -std=c++17 -O2 "${linker_flags[@]}" -I"$repo/dpsim-fmi/third_party/fmi3" \
  "$case_dir/IEEE39AvgInverter.cpp" \
  -o "$staging/binaries/$platform_dir/$binary_name"
rm -f "$build_dir/IEEE39AvgInverter.fmu"
(cd "$staging" && zip -qr "$build_dir/IEEE39AvgInverter.fmu" .)
printf '%s\n' "$build_dir/IEEE39AvgInverter.fmu"
cp "$build_dir/IEEE39AvgInverter.fmu" "$case_dir/IEEE39AvgInverter.fmu"
printf '%s\n' "$case_dir/IEEE39AvgInverter.fmu"
