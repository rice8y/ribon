#!/usr/bin/env python3
"""Render and inspect the 72-page native-vector Rfam validation document."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pdf", type=Path, default=Path("target/qa/ribon-real-data-validation.pdf"), nargs="?"
    )
    parser.add_argument(
        "--corpus", type=Path, default=Path("tests/data/rfam_real_24.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("tests/reports/render_validation.json")
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pdf = root / arguments.pdf
    corpus = json.loads((root / arguments.corpus).read_text())
    expected_pages = 3 * len(corpus["cases"])

    run(["qpdf", "--check", str(pdf)], root)

    information = run(["pdfinfo", str(pdf)], root).stdout
    pages_match = re.search(r"^Pages:\s+(\d+)$", information, re.MULTILINE)
    if not pages_match or int(pages_match.group(1)) != expected_pages:
        raise AssertionError(f"expected {expected_pages} PDF pages")

    images = run(["pdfimages", "-list", str(pdf)], root).stdout
    image_rows = [line for line in images.splitlines() if re.match(r"^\s*\d+\s+\d+\s+", line)]
    if image_rows:
        raise AssertionError(f"expected native vectors, found {len(image_rows)} raster images")

    extracted = run(["pdftotext", str(pdf), "-"], root).stdout
    missing_accessions = [
        case["accession"] for case in corpus["cases"] if case["accession"] not in extracted
    ]
    if missing_accessions:
        raise AssertionError(f"missing page labels: {missing_accessions}")
    required_feature_labels = [
        "NAView",
        "Affine loop geometry",
        "Collision-reduced loop geometry",
        "AA/AB cofold macrostate",
        "Local pair/accessibility",
        "Circular ensemble MFE",
        "Modified-base profile",
        "G-quadruplex states",
        "Multi-component H-type pseudoknot",
        "Gap-aware covariation consensus",
        "Decoder source probabilities",
        "RNAstructure 6.6 DNA model",
        "Custom normalized model",
        "Fatgraph topology annotation",
        "MFE/MEA structure comparison",
        "Exact minimum-saddle landscape",
        "Exact inverse folding",
        "Exact ligand microstates",
    ]
    missing_features = [
        label for label in required_feature_labels if extracted.count(label) < len(corpus["cases"])
    ]
    if missing_features:
        raise AssertionError(f"feature panels missing from one or more real-data cases: {missing_features}")

    metrics = []
    with tempfile.TemporaryDirectory(prefix="ribon-pdf-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "110", str(pdf), str(prefix)], root)
        rendered = sorted(Path(directory).glob("page-*.png"))
        if len(rendered) != expected_pages:
            raise AssertionError(f"rendered {len(rendered)} of {expected_pages} pages")
        for page, image in enumerate(rendered, 1):
            output = run(
                [
                    "magick",
                    str(image),
                    "-fuzz",
                    "2%",
                    "-trim",
                    "-format",
                    "%w %h %[fx:page.x] %[fx:page.y] %[fx:page.width] "
                    "%[fx:page.height] %[fx:mean] %[fx:standard_deviation]",
                    "info:",
                ],
                root,
            ).stdout
            width, height, x, y, canvas_width, canvas_height, mean, deviation = map(
                float, output.split()
            )
            margins = [x, y, canvas_width - x - width, canvas_height - y - height]
            if min(margins) < 8:
                raise AssertionError(f"page {page} content approaches the crop edge: {margins}")
            if deviation < 0.02 or mean > 0.999:
                raise AssertionError(f"page {page} appears blank: mean={mean}, sd={deviation}")
            metrics.append(
                {
                    "page": page,
                    "content_margins_pixels": margins,
                    "ink_fraction": 1.0 - mean,
                    "grayscale_standard_deviation": deviation,
                }
            )

    report = {
        "pdf": arguments.pdf.as_posix(),
        "pages": expected_pages,
        "embedded_raster_images": 0,
        "qpdf_syntax_check": True,
        "all_accession_labels_extractable": True,
        "feature_panels_per_case": len(required_feature_labels),
        "all_feature_labels_extractable_per_case": True,
        "minimum_content_margin_pixels": min(
            min(page["content_margins_pixels"]) for page in metrics
        ),
        "minimum_ink_fraction": min(page["ink_fraction"] for page in metrics),
        "maximum_ink_fraction": max(page["ink_fraction"] for page in metrics),
        "minimum_grayscale_standard_deviation": min(
            page["grayscale_standard_deviation"] for page in metrics
        ),
        "pages_detail": metrics,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "pages_detail"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
