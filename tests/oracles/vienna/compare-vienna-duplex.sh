#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
  echo "usage: $0 VIENNA_BUILD SEQUENCE_A SEQUENCE_B [SALT_MOLAR] [VIENNA_SOURCE]" >&2
  exit 2
fi

build=$1
sequence_a=$2
sequence_b=$3
salt=${4:-1.021}
source=${5:-artifacts/ViennaRNA}
root=$(cd "$(dirname "$0")/../../.." && pwd)
binary="$build/ribon-vienna-duplex"
cc -O2 \
  -I"$build/src" -I"$source/src" -I"$source" \
  "$root/tests/oracles/vienna/vienna_duplex_harness.c" \
  "$build/src/ViennaRNA/.libs/libRNA.a" -lm -o "$binary"
"$binary" "$sequence_a" "$sequence_b" "$salt"
