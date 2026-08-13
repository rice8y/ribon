#!/usr/bin/env python3
"""Semantic and pixel-exact checks for the public rendering contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

from imagemagick import convert_command
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


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", nargs="?", type=Path, default=Path("target/qa/typst/publication.pdf"))
    parser.add_argument("--update-golden", action="store_true")
    parser.add_argument(
        "--golden",
        type=Path,
        default=Path("tests/golden/publication_render_golden_sha256.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/reports/publication_render_validation.json"),
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pdf = root / arguments.pdf

    run(["qpdf", "--check", str(pdf)], root)
    information = run(["pdfinfo", str(pdf)], root)
    if not re.search(r"^Pages:\s+4$", information, re.MULTILINE):
        raise AssertionError("publication contract must be exactly four pages")
    images = run(["pdfimages", "-list", str(pdf)], root)
    if any(re.match(r"^\s*\d+\s+\d+\s+", line) for line in images.splitlines()):
        raise AssertionError("publication contract contains a raster image XObject")

    extracted = run(["pdftotext", str(pdf), "-"], root)
    labels = [
        "Aspect-preserving circular layout",
        "Rotation, reflection, and direction",
        "apical loop",
        "reactive site",
        "tertiary contact",
        "Reference/alternative pair classification",
        "Multi-strand identifiers and termini",
        "guide",
        "target",
        "SHAPE reactivity",
        "Two-ensemble probability dot plot",
        "Expected and discrete mountain profiles",
        "Fixed-structure evaluation",
        "Boltzmann sample",
        "Suboptimal structure",
        "Accessibility overlay",
        "Local accessibility overlay",
        "MFE, centroid, and MEA comparison",
        "600 nt linear overview",
        "Hand-edited scene coordinates",
    ]
    missing = [label for label in labels if label not in extracted]
    if missing:
        raise AssertionError(f"publication labels missing from vector text: {missing}")

    metadata = json.loads(
        run(
            [
                "typst",
                "query",
                "--root",
                ".",
                "tests/typst/publication.typ",
                "<ribon-publication-qa>",
                "--field",
                "value",
                "--one",
            ],
            root,
        )
    )
    expected_metadata = {
        "schema": "ribon.publication-qa/1",
        "sequence-length": 18,
        "stems": 2,
        "hairpins": 2,
        "shared-pairs": 3,
        "reference-only-pairs": 3,
        "alternative-only-pairs": 2,
        "annotation-controls": [
            "nucleotide-label",
            "strand-label",
            "interaction-label",
        ],
        "resource-errors": ["resource_limit", "resource_limit", "resource_limit"],
    }
    if metadata != expected_metadata:
        raise AssertionError(f"publication semantic metadata drift: {metadata}")

    page_metrics: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ribon-publication-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "120", str(pdf), str(prefix)], root)
        pages = sorted(Path(directory).glob("page-*.png"))
        if len(pages) != 4:
            raise AssertionError(f"Poppler rendered {len(pages)} publication pages")
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
            if min(margins) < 12:
                raise AssertionError(f"page {page_number} approaches the crop edge: {margins}")
            if 1.0 - mean < 0.002 or deviation < 0.025:
                raise AssertionError(f"page {page_number} appears blank or too faint")
            page_metrics.append(
                {
                    "page": page_number,
                    "margins": margins,
                    "ink_fraction": 1.0 - mean,
                    "grayscale_standard_deviation": deviation,
                }
            )

        # The first panel is deliberately placed on a wide, short canvas. A
        # saturation mask isolates its colored nucleotide nodes; a near-square
        # box proves that the renderer did not stretch circular coordinates.
        circle = run(
            convert_command(
                str(pages[0]),
                "-crop",
                "400x300+80+50",
                "-colorspace",
                "HSL",
                "-channel",
                "G",
                "-separate",
                "+channel",
                "-threshold",
                "12%",
                "-trim",
                "-format",
                "%w %h",
                "info:",
            ),
            root,
        )
        circle_width, circle_height = map(float, circle.split())
        circle_ratio = circle_width / circle_height
        if not 0.97 <= circle_ratio <= 1.03:
            raise AssertionError(
                f"aspect-preserving circular panel is distorted: ratio={circle_ratio}"
            )

    manifest = {
        "schema": 1,
        "source": "tests/typst/publication.typ",
        "dpi": 120,
        "page_count": 4,
        "sha256": hashes,
    }
    golden = root / arguments.golden
    if arguments.update_golden:
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        expected = json.loads(golden.read_text())
        validate_pixel_golden(manifest, expected, label="publication")

    report = {
        "pdf": arguments.pdf.as_posix(),
        "pages": 4,
        "embedded_raster_images": 0,
        "vector_text_labels": labels,
        "semantic_metadata": metadata,
        "aspect_preserving_circle_ratio": circle_ratio,
        "pixel_exact": True,
        "golden_updated": arguments.update_golden,
        "pages_detail": page_metrics,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "pages_detail"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
