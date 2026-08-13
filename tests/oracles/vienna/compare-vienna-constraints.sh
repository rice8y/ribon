#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 8 ]]; then
  echo "usage: $0 VIENNA_BUILD VIENNA_SOURCE SEQ DANGLES HC NOLP NOGU MAXSPAN [MODE [MODE_ARGS...]]" >&2
  exit 2
fi

build_dir="$1"
source_dir="$2"
shift 2
archive="$build_dir/src/ViennaRNA/.libs/libRNA.a"
binary="$build_dir/ribon-vienna-constraints"
object="$build_dir/ribon-vienna-constraints.o"
root=$(cd "$(dirname "$0")/../../.." && pwd)

cc \
  -I"$build_dir" \
  -I"$build_dir/src" \
  -I"$source_dir/src" \
  -c "$root/tests/oracles/vienna/vienna_constraints_harness.c" \
  -o "$object"
c++ "$object" "$archive" -lm -o "$binary"

if [[ $# -eq 6 ]]; then
  "$binary" "$@" none
else
  "$binary" "$@"
fi
