#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 VIENNA_BUILD [STRUCTURE]" >&2
  exit 2
fi

build_dir="$1"
structure="${2:-((..((....)).(((....))).))}"
source_dir="${VIENNA_SOURCE:-artifacts/ViennaRNA}"
object="$build_dir/src/ViennaRNA/plotting/naview/.libs/naview.o"
binary="$build_dir/ribon-vienna-naview"
root=$(cd "$(dirname "$0")/../../.." && pwd)

if [[ ! -f "$object" ]]; then
  echo "missing $object; configure and build ViennaRNA's naview target first" >&2
  exit 2
fi

cc \
  -I"$build_dir" \
  -I"$source_dir/src" \
  "$root/tests/oracles/vienna/vienna_naview_harness.c" \
  "$object" \
  -lm \
  -o "$binary"

"$binary" "$structure"
