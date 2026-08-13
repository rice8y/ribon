#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "usage: $0 VIENNA_BUILD turtle|puzzler STRUCTURE [VIENNA_SOURCE]" >&2
  exit 2
fi

build_dir="$1"
method="$2"
structure="$3"
source_dir="${4:-${VIENNA_SOURCE:-artifacts/ViennaRNA}}"
archive="$build_dir/src/ViennaRNA/.libs/libRNA.a"
binary="$build_dir/ribon-vienna-modern-layout"
root=$(cd "$(dirname "$0")/../../.." && pwd)

cc -I"$build_dir" -I"$build_dir/src" -I"$source_dir/src" \
  "$root/tests/oracles/vienna/vienna_modern_layout_harness.c" "$archive" -lm -lstdc++ -o "$binary"
"$binary" "$method" "$structure"
