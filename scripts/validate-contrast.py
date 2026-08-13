#!/usr/bin/env python3
"""Validate strict WCAG node-text contrast failure contracts."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASES = (
    (
        "tests/typst/contrast_aaa_failure.typ",
        "node 1 cannot attain WCAG AAA text contrast",
    ),
    (
        "tests/typst/contrast_transparency_failure.typ",
        "a transparent node fill requires contrast-background",
    ),
    (
        "tests/typst/contrast_paint_failure.typ",
        "contrast background must be a solid color",
    ),
)


def main() -> int:
    validated = []
    with tempfile.TemporaryDirectory(prefix="ribon-contrast-") as directory:
        for index, (source, expected) in enumerate(CASES):
            output = Path(directory) / f"case-{index}.pdf"
            process = subprocess.run(
                ["typst", "compile", "--root", str(ROOT), source, str(output)],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            if process.returncode == 0:
                raise AssertionError(f"contrast failure fixture unexpectedly compiled: {source}")
            if expected not in process.stderr:
                raise AssertionError(
                    f"contrast failure mismatch for {source}: expected {expected!r}\n{process.stderr}"
                )
            validated.append({"source": source, "message": expected})
    print(json.dumps({"schema": "ribon.contrast-validation/1", "validated": validated}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
