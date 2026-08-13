#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 VIENNA_BUILD STRUCTURE [VIENNA_SOURCE]" >&2
  exit 2
fi

build_dir="$1"
structure="$2"
source_dir="${3:-${VIENNA_SOURCE:-artifacts/ViennaRNA}}"
archive="$build_dir/src/ViennaRNA/.libs/libRNA.a"
binary="$build_dir/ribon-vienna-simple"
root=$(cd "$(dirname "$0")/../../.." && pwd)

if [[ ! -f "$archive" ]]; then
  echo "missing $archive; build ViennaRNA's src/ViennaRNA target first" >&2
  exit 2
fi

cc \
  -I"$build_dir" \
  -I"$build_dir/src" \
  -I"$source_dir/src" \
  "$root/tests/oracles/vienna/vienna_simple_harness.c" \
  "$archive" \
  -lm \
  -lstdc++ \
  -o "$binary"

"$binary" "$structure"
