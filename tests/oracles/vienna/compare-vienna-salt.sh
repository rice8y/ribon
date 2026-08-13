#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 6 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE SALT_MOLAR [VIENNA_SOURCE] [DANGLES] [STRUCTURE]" >&2
  exit 2
fi

build=$1
sequence=$2
salt=$3
source=${4:-artifacts/ViennaRNA}
dangles=${5:-2}
structure=${6:-}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-salt"
cc -O2 \
  -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_salt_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -o "$binary"
if [[ -n "$structure" ]]; then
  "$binary" "$sequence" "$salt" "$dangles" "$structure"
else
  "$binary" "$sequence" "$salt" "$dangles"
fi
