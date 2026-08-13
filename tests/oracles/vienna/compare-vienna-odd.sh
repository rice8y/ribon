#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE DANGLES [VIENNA_SOURCE] [STRUCTURE]" >&2
  exit 2
fi

build_dir="$1"
sequence="$2"
dangles="$3"
source_dir="${4:-${VIENNA_SOURCE:-artifacts/ViennaRNA}}"
structure="${5:-}"
archive="$build_dir/src/ViennaRNA/.libs/libRNA.a"
binary="$build_dir/ribon-vienna-odd"
object="$build_dir/ribon-vienna-odd.o"
root=$(cd "$(dirname "$0")/../../.." && pwd)

if [[ ! -f "$archive" ]]; then
  echo "missing $archive; build ViennaRNA's src/ViennaRNA target first" >&2
  exit 2
fi

cc \
  -I"$build_dir" \
  -I"$build_dir/src" \
  -I"$source_dir/src" \
  -c "$root/tests/oracles/vienna/vienna_odd_harness.c" \
  -o "$object"

c++ "$object" "$archive" -lm -o "$binary"

if [[ -n "$structure" ]]; then
  "$binary" "$sequence" "$dangles" "$structure"
else
  "$binary" "$sequence" "$dangles"
fi
