#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
EXAMPLE="$ROOT/dpsim-fmi/examples/mpc_grid"
OUT="$ROOT/build-fmi/mpc-grid-fmu"
STAGE="$OUT/stage"
mkdir -p "$STAGE/binaries/x86_64-linux"
cp "$EXAMPLE/sources/modelDescription.xml" "$STAGE/modelDescription.xml"
cc -std=c11 -O2 -fPIC -shared   -I "$ROOT/dpsim-fmi/third_party/fmi3"   "$EXAMPLE/sources/MpcGridController.c" -lm   -o "$STAGE/binaries/x86_64-linux/mpc_grid_mpc.so"
rm -f "$OUT/mpc_grid_mpc.fmu"
(cd "$STAGE" && zip -qr "$OUT/mpc_grid_mpc.fmu" .)
echo "Built $OUT/mpc_grid_mpc.fmu"
