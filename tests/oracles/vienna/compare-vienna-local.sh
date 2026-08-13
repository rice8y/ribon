#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 || $# -gt 6 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE WINDOW SPAN UNPAIRED [VIENNA_SOURCE]" >&2
  exit 2
fi

build=$1
sequence=$2
window=$3
span=$4
unpaired=$5
source=${6:-artifacts/ViennaRNA}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-local"
cc -O2 -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_local_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -lstdc++ -o "$binary"
"$binary" "$sequence" "$window" "$span" "$unpaired"
