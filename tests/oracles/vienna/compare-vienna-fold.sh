#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE [VIENNA_SOURCE] [DANGLES]" >&2
  exit 2
fi

build_dir="$1"
sequence="$2"
source_dir="${3:-${VIENNA_SOURCE:-artifacts/ViennaRNA}}"
dangles="${4:-2}"
archive="$build_dir/src/ViennaRNA/.libs/libRNA.a"
binary="$build_dir/ribon-vienna-fold"
object="$build_dir/ribon-vienna-fold.o"
root=$(cd "$(dirname "$0")/../../.." && pwd)

if [[ ! -f "$archive" ]]; then
  echo "missing $archive; build ViennaRNA's src/ViennaRNA target first" >&2
  exit 2
fi

cc \
  -I"$build_dir" \
  -I"$build_dir/src" \
  -I"$source_dir/src" \
  -c "$root/tests/oracles/vienna/vienna_fold_harness.c" \
  -o "$object"

c++ \
  "$object" \
  "$archive" \
  -lm \
  -o "$binary"

"$binary" "$sequence" "$dangles"
