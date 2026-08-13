#!/usr/bin/env python3
"""Numeric, semantic, vector, and pixel-exact checks for quantitative figures."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

from imagemagick import convert_command


def run(command: list[str], root: Path) -> str:
    return subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def first_line(command: list[str], root: Path) -> str:
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
    parser.add_argument("pdf", nargs="?", type=Path, default=Path("target/qa/typst/plot-quality.pdf"))
    parser.add_argument("--update-golden", action="store_true")
    parser.add_argument(
        "--golden",
        type=Path,
        default=Path("tests/golden/plot_quality_golden_sha256.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/reports/plot_quality_validation.json"),
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pdf = root / arguments.pdf

    run(["qpdf", "--check", str(pdf)], root)
    information = run(["pdfinfo", str(pdf)], root)
    if not re.search(r"^Pages:\s+2$", information, re.MULTILINE):
        raise AssertionError("plot quality contract must be exactly two pages")
    images = run(["pdfimages", "-list", str(pdf)], root)
    if any(re.match(r"^\s*\d+\s+\d+\s+", line) for line in images.splitlines()):
        raise AssertionError("plot quality contract contains a raster image XObject")

    metadata = json.loads(
        run(
            [
                "typst",
                "query",
                "--root",
                ".",
                "tests/typst/plot_quality.typ",
                "<ribon-plot-qa>",
                "--field",
                "value",
                "--one",
            ],
            root,
        )
    )
    expected_metadata = {
        "schema": "ribon.plot-qa/1",
        "exact-expected": [0.5, 0.75, 0.5, 0.0],
        "exact-reference": [1, 2, 1, 0],
        "exact-length": 4,
        "comparison-pair-count": 2,
    }
    if metadata != expected_metadata:
        raise AssertionError(f"plot numeric metadata drift: {metadata}")

    extracted = run(["pdftotext", str(pdf), "-"], root)
    required_labels = [
        "Two-ensemble probability matrix",
        "Pair probability",
        "Comparison probability",
        "Enclosing base pairs",
        "Sequence position",
        "Reference pair",
        "Normalized position",
    ]
    missing = [label for label in required_labels if label not in extracted]
    if missing:
        raise AssertionError(f"plot labels missing from vector text: {missing}")

    pages_report: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ribon-plot-quality-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "150", str(pdf), str(prefix)], root)
        pages = sorted(Path(directory).glob("page-*.png"))
        if len(pages) != 2:
            raise AssertionError(f"Poppler rendered {len(pages)} plot pages")
        hashes = [digest(page) for page in pages]
        for page_number, page in enumerate(pages, 1):
            metric = run(
                convert_command(
                    str(page),
                    "-fuzz",
                    "2%",
                    "-trim",
                    "-format",
                    "%w %h %[fx:page.x] %[fx:page.y] %[fx:page.width] "
                    "%[fx:page.height] %[fx:mean] %[fx:standard_deviation]",
                    "info:",
                ),
                root,
            )
            width, height, x, y, canvas_width, canvas_height, mean, deviation = map(
                float, metric.split()
            )
            margins = [x, y, canvas_width - x - width, canvas_height - y - height]
            if min(margins) < 35:
                raise AssertionError(f"plot page {page_number} approaches the crop edge: {margins}")
            if 1.0 - mean < 0.004 or deviation < 0.035:
                raise AssertionError(f"plot page {page_number} is blank or lacks contrast")
            pages_report.append(
                {
                    "page": page_number,
                    "margins_pixels": margins,
                    "ink_fraction": 1.0 - mean,
                    "standard_deviation": deviation,
                }
            )

    manifest = {
        "schema": "ribon.plot-pixel-golden/1",
        "source": "tests/typst/plot_quality.typ",
        "renderer": first_line(["pdftoppm", "-v"], root),
        "resolution_ppi": 150,
        "page_count": 2,
        "sha256": hashes,
    }
    golden = root / arguments.golden
    if arguments.update_golden:
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        expected = json.loads(golden.read_text())
        stable_keys = ("schema", "source", "renderer", "resolution_ppi", "page_count")
        for key in stable_keys:
            if manifest[key] != expected.get(key):
                raise AssertionError(f"plot golden environment differs for {key}")
        if manifest["sha256"] != expected.get("sha256"):
            changed = [
                index
                for index, (actual, wanted) in enumerate(
                    zip(manifest["sha256"], expected.get("sha256", []), strict=False),
                    1,
                )
                if actual != wanted
            ]
            raise AssertionError(f"plot pixel golden drift: pages={changed}")

    report = {
        "schema": "ribon.plot-validation/1",
        "pdf": arguments.pdf.as_posix(),
        "numeric_exact": True,
        "vector_only": True,
        "pixel_exact": True,
        "golden_updated": arguments.update_golden,
        "metadata": metadata,
        "pages": pages_report,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
