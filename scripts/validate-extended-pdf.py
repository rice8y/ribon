#!/usr/bin/env python3
"""Render-level checks for the compact extended-feature Typst matrix."""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], root: Path) -> str:
    return subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    pdf = root / "target/qa/typst/extended.pdf"
    information = run(["pdfinfo", str(pdf)], root)
    if not re.search(r"^Pages:\s+1$", information, re.MULTILINE):
        raise AssertionError("extended matrix must be exactly one page")
    images = run(["pdfimages", "-list", str(pdf)], root)
    if any(re.match(r"^\s*\d+\s+\d+\s+", line) for line in images.splitlines()):
        raise AssertionError("extended matrix contains a raster image XObject")
    extracted = run(["pdftotext", str(pdf), "-"], root)
    labels = [
        "Cofold",
        "Circular RNA",
        "Modified base",
        "G-quadruplex",
        "Pseudoknot",
        "Local accessibility",
        "Exact energy landscape",
        "Exact inverse design",
        "Ligand microstate ensemble",
        "aptamer",
        "p=0.",
        "m6A",
    ]
    missing = [label for label in labels if label not in extracted]
    if missing:
        raise AssertionError(f"missing extended visual labels: {missing}")
    with tempfile.TemporaryDirectory(prefix="ribon-extended-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "140", "-singlefile", str(pdf), str(prefix)], root)
        metric = run(
            [
                "magick",
                str(prefix.with_suffix(".png")),
                "-fuzz",
                "2%",
                "-trim",
                "-format",
                "%w %h %[fx:page.x] %[fx:page.y] %[fx:page.width] "
                "%[fx:page.height] %[fx:mean] %[fx:standard_deviation]",
                "info:",
            ],
            root,
        )
    width, height, x, y, canvas_width, canvas_height, mean, deviation = map(
        float, metric.split()
    )
    margins = [x, y, canvas_width - x - width, canvas_height - y - height]
    ink_fraction = 1.0 - mean
    if min(margins) < 8 or ink_fraction < 0.005 or deviation < 0.04:
        raise AssertionError(
            f"extended matrix visual metrics failed: margins={margins}, ink={ink_fraction}, sd={deviation}"
        )
    report = {
        "pdf": "target/qa/typst/extended.pdf",
        "pages": 1,
        "feature_labels": labels,
        "embedded_raster_images": 0,
        "content_margins_pixels": margins,
        "ink_fraction": ink_fraction,
        "grayscale_standard_deviation": deviation,
    }
    output = root / "tests/reports/extended_render_validation.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
