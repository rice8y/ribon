#!/usr/bin/env python3
"""Validate generic axis, viewport, and legend placement contracts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

from imagemagick import convert_command


def run(command: list[str], root: Path, *, stderr_to_stdout: bool = False) -> str:
    result = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT if stderr_to_stdout else subprocess.PIPE,
    )
    return result.stdout


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pdf", nargs="?", type=Path, default=Path("target/qa/typst/plot-layout-quality.pdf")
    )
    parser.add_argument("--update-golden", action="store_true")
    parser.add_argument(
        "--golden",
        type=Path,
        default=Path("tests/golden/plot_layout_golden_sha256.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/reports/plot_layout_validation.json"),
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pdf = root / arguments.pdf

    run(["qpdf", "--check", str(pdf)], root)
    information = run(["pdfinfo", str(pdf)], root)
    if not re.search(r"^Pages:\s+3$", information, re.MULTILINE):
        raise AssertionError("plot layout contract must be exactly three pages")
    images = run(["pdfimages", "-list", str(pdf)], root)
    if any(re.match(r"^\s*\d+\s+\d+\s+", line) for line in images.splitlines()):
        raise AssertionError("plot layout contract contains a raster image XObject")

    metadata = json.loads(
        run(
            [
                "typst", "query", "--root", ".", "tests/typst/plot_layout_quality.typ",
                "<ribon-plot-layout-qa>", "--field", "value", "--one",
            ],
            root,
        )
    )
    expected_metadata = {
        "schema": "ribon.plot-layout-qa/1",
        "outer-positions": ["top", "bottom", "left", "right"],
        "inner-position-count": 9,
        "explicit-position": [0.28, 0.32],
        "axis-modes": ["linear", "log"],
        "secondary-axes": ["x2", "y2"],
        "exact-profile": [0.5, 0.75, 0.5, 0.0],
        "exact-reference-profiles": [[1, 2, 1, 0], [1, 1, 1, 0]],
        "reference-axes": [["x", "y"], ["x", "y"]],
        "secondary-series-axes": ["x", "y2"],
    }
    if metadata != expected_metadata:
        raise AssertionError(f"plot layout metadata drift: {metadata}")

    extracted = run(["pdftotext", str(pdf), "-"], root)
    labels = [
        "Outer legend placement",
        "One shared legend for two independent plots",
        "inner-north-west",
        "inner-center",
        "inner-south-east",
        "Explicit 28%, 32% coordinate",
        "Custom domain, tick labels, minor grid, and aspect",
        "Reverse position",
        "Log scale",
        "Top index",
        "Secondary",
        "Vertical reversed continuous legend",
    ]
    missing = [label for label in labels if label not in extracted]
    if missing:
        raise AssertionError(f"plot layout labels missing from vector text: {missing}")

    pages_report: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ribon-plot-layout-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "150", str(pdf), str(prefix)], root)
        pages = sorted(Path(directory).glob("page-*.png"))
        if len(pages) != 3:
            raise AssertionError(f"Poppler rendered {len(pages)} layout pages")
        hashes = [digest(page) for page in pages]
        for page_number, page in enumerate(pages, 1):
            metric = run(
                convert_command(
                    str(page), "-fuzz", "2%", "-trim", "-format",
                    "%w %h %[fx:page.x] %[fx:page.y] %[fx:page.width] "
                    "%[fx:page.height] %[fx:mean] %[fx:standard_deviation]", "info:",
                ),
                root,
            )
            width, height, x, y, canvas_width, canvas_height, mean, deviation = map(
                float, metric.split()
            )
            margins = [x, y, canvas_width - x - width, canvas_height - y - height]
            if min(margins) < 35:
                raise AssertionError(f"layout page {page_number} approaches crop edge: {margins}")
            if 1.0 - mean < 0.003 or deviation < 0.03:
                raise AssertionError(f"layout page {page_number} is blank or too faint")
            pages_report.append(
                {
                    "page": page_number,
                    "margins_pixels": margins,
                    "ink_fraction": 1.0 - mean,
                    "standard_deviation": deviation,
                }
            )

    renderer = run(["pdftoppm", "-v"], root, stderr_to_stdout=True).splitlines()[0]
    manifest = {
        "schema": "ribon.plot-layout-pixel-golden/1",
        "source": "tests/typst/plot_layout_quality.typ",
        "renderer": renderer,
        "resolution_ppi": 150,
        "page_count": 3,
        "sha256": hashes,
    }
    golden = root / arguments.golden
    if arguments.update_golden:
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        expected = json.loads(golden.read_text())
        if manifest != expected:
            changed = [
                index
                for index, (actual, wanted) in enumerate(
                    zip(manifest["sha256"], expected.get("sha256", []), strict=False), 1
                )
                if actual != wanted
            ]
            raise AssertionError(f"plot layout pixel golden drift: pages={changed}")

    report = {
        "schema": "ribon.plot-layout-validation/1",
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
