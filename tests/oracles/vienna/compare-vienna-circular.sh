#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE [VIENNA_SOURCE]" >&2
  exit 2
fi

build=$1
sequence=$2
source=${3:-artifacts/ViennaRNA}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-circular"
cc -O2 -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_circular_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -lstdc++ -o "$binary"
"$binary" "$sequence"
