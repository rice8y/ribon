#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE BAND_KCAL [DANGLES] [VIENNA_SOURCE]" >&2
  exit 2
fi

build=$1
sequence=$2
band=$3
dangles=${4:-2}
source=${5:-artifacts/ViennaRNA}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-suboptimal"
cc -O2 \
  -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_suboptimal_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -o "$binary"
"$binary" "$sequence" "$band" "$dangles"
