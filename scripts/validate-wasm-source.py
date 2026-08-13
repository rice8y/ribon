#!/usr/bin/env python3
"""Verify that the distributed WASM records the current build-input fingerprint."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA = "ribon.wasm-source-fingerprint/1"


def input_paths(root: Path) -> list[Path]:
    paths = [
        root / "Cargo.lock",
        root / "Cargo.toml",
        root / "rust-toolchain.toml",
        root / "crates/ribon-core/Cargo.toml",
        root / "crates/ribon-core/build.rs",
        root / "crates/ribon-plugin/Cargo.toml",
    ]
    paths.extend((root / "crates/ribon-core/src").rglob("*.rs"))
    paths.extend(
        path for path in (root / "crates/ribon-core/data").rglob("*") if path.is_file()
    )
    paths.extend((root / "crates/ribon-plugin/src").rglob("*.rs"))
    return sorted(set(paths), key=lambda path: path.relative_to(root).as_posix())


def fingerprint(root: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    paths = input_paths(root)
    for path in paths:
        relative = path.relative_to(root).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return {
        "algorithm": "sha256(path + NUL + content + NUL)",
        "file_count": len(paths),
        "schema": SCHEMA,
        "sha256": digest.hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--update",
        action="store_true",
        help="record the current fingerprint after rebuilding the distributed WASM",
    )
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    report_path = root / "tests/reports/wasm_source_validation.json"
    current = fingerprint(root)
    if arguments.update:
        report_path.write_text(
            json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        if not report_path.is_file():
            raise AssertionError("WASM source fingerprint is missing; run 'just plugin'")
        recorded = json.loads(report_path.read_text(encoding="utf-8"))
        if recorded != current:
            raise AssertionError(
                "package/ribon_plugin.wasm is not synchronized with its build inputs; "
                "run 'just plugin'"
            )
    print(json.dumps(current, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
