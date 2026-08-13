#!/usr/bin/env python3
"""Pixel-exact regression test for the 72-page real-data render matrix.

The compact golden file stores SHA-256 digests of Poppler-rendered pages rather
than raster images.  A deliberate visual change therefore requires an explicit
``--update-golden`` run, while ordinary release checks fail on a one-pixel drift.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from pixel_golden import validate_pixel_golden


def run(command: list[str], root: Path) -> str:
    return subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def version(command: list[str], root: Path) -> str:
    output = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout
    return output.splitlines()[0].strip()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pdf",
        nargs="?",
        type=Path,
        default=Path("target/qa/ribon-real-data-validation.pdf"),
    )
    parser.add_argument(
        "--golden",
        type=Path,
        default=Path("tests/golden/render_golden_sha256.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/reports/render_golden_validation.json"),
    )
    parser.add_argument("--update-golden", action="store_true")
    parser.add_argument("--dpi", type=int, default=110)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    pdf = (root / arguments.pdf).resolve()
    if not pdf.is_file():
        raise FileNotFoundError(pdf)

    with tempfile.TemporaryDirectory(prefix="ribon-golden-") as directory:
        prefix = Path(directory) / "page"
        run(
            ["pdftoppm", "-png", "-r", str(arguments.dpi), str(pdf), str(prefix)],
            root,
        )
        pages = sorted(Path(directory).glob("page-*.png"))
        if not pages:
            raise AssertionError("Poppler produced no rendered pages")
        hashes = [digest(page) for page in pages]

    manifest = {
        "schema": 1,
        "source": "tests/typst/real_data_render.typ",
        "dpi": arguments.dpi,
        "page_count": len(hashes),
        "sha256": hashes,
    }
    golden_path = root / arguments.golden
    if arguments.update_golden:
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        if not golden_path.is_file():
            raise FileNotFoundError(
                f"visual golden is missing: run {Path(__file__).name} --update-golden"
            )
        golden = json.loads(golden_path.read_text())
        validate_pixel_golden(manifest, golden, label="real-data render")

    report = {
        "pdf": arguments.pdf.as_posix(),
        "golden": arguments.golden.as_posix(),
        "pages": len(hashes),
        "dpi": arguments.dpi,
        "pixel_exact": True,
        "changed_pages": [],
        "rasterizer": version(["pdftoppm", "-v"], root),
        "typst": version(["typst", "--version"], root),
        "golden_updated": arguments.update_golden,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
