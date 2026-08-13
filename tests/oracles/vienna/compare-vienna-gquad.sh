#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 || $# -gt 6 ]]; then
  echo "usage: $0 VIENNA_BUILD STACK LINKER1 LINKER2 LINKER3 [TEMPERATURE_C]" >&2
  exit 2
fi

build=$1
stack=$2
linker1=$3
linker2=$4
linker3=$5
temperature=${6:-37}
source=${VIENNA_SOURCE:-artifacts/ViennaRNA}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-gquad"
cc -O2 -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_gquad_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -lstdc++ -o "$binary"
"$binary" "$stack" "$linker1" "$linker2" "$linker3" "$temperature"
