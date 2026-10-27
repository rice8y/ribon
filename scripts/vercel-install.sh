#!/usr/bin/env bash
set -euo pipefail

typage_version="${TYPAGE_VERSION:-0.2.1}"
typage_rust_version="${TYPAGE_RUST_VERSION:-1.92.0}"
starter_version="${TYPST_DOCS_STARTER_VERSION:-0.1.4}"
starter="${TYPST_DOCS_STARTER:-github:rice8y/typage-starter-typst-docs#v${starter_version}}"
docs_root="${TYPAGE_DOCS_ROOT:-ribon-docs}"

case "$docs_root" in
  "" | "." | "/" | /* | *".."*)
    echo "Unsafe TYPAGE_DOCS_ROOT: ${docs_root}" >&2
    exit 1
    ;;
esac

export PATH="$HOME/.cargo/bin:/rust/bin:$PATH"

echo "[setup] installing Rust ${typage_rust_version} for Typage"
rustup toolchain install "$typage_rust_version" --profile minimal
# Override the runtime's rust-toolchain.toml for documentation tooling.
export RUSTUP_TOOLCHAIN="$typage_rust_version"

echo "[setup] installing Typage ${typage_version} with $(rustc --version)"
cargo install typage --version "$typage_version" --locked

echo "[setup] initializing Typage docs starter from ${starter}"
rm -rf "$docs_root"
typage init "$docs_root" \
  --starter "$starter" \
  --install-plugins
