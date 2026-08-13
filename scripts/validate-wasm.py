#!/usr/bin/env python3
"""Validate distributable WASM modules and their complete host-import surface."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path


EXPECTED = {
    ("typst_env", "wasm_minimal_protocol_write_args_to_buffer"),
    ("typst_env", "wasm_minimal_protocol_send_result_to_host"),
}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    rows = []
    for relative in ("package/ribon_plugin.wasm",):
        path = root / relative
        subprocess.run(["wasm-tools", "validate", str(path)], cwd=root, check=True)
        wat = subprocess.run(
            ["wasm-tools", "print", str(path)],
            cwd=root,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
        imports = set(re.findall(r'\(import "([^"]+)" "([^"]+)"', wat))
        if imports != EXPECTED:
            raise AssertionError(f"unexpected host imports in {relative}: {sorted(imports)}")
        exports = set(re.findall(r'\(export "([^"]+)"', wat))
        if "run" not in exports:
            raise AssertionError(f"stable public exports missing from {relative}: {sorted(exports)}")
        if "version" in exports:
            raise AssertionError(f"obsolete version export remains in {relative}")
        rows.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "imports": [list(entry) for entry in sorted(imports)],
                "public_exports": sorted({"run"} & exports),
                "valid": True,
            }
        )
    report = {"wasm_tools_validation": True, "modules": rows}
    output = root / "tests/reports/wasm_validation.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
