#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
built="$root/target/wasm32-unknown-unknown/release/ribon_plugin.wasm"
distributed="$root/package/ribon_plugin.wasm"

cargo build --release --locked --target wasm32-unknown-unknown -p ribon-plugin --manifest-path "$root/Cargo.toml"

if ! cmp -s "$built" "$distributed"; then
  echo "package/ribon_plugin.wasm is not synchronized; run 'just plugin'" >&2
  exit 1
fi

echo "package/ribon_plugin.wasm matches the release build"
