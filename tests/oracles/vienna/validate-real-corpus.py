#!/usr/bin/env python3
"""Numerical and coordinate validation on the pinned 24-family Rfam corpus."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import tempfile
from pathlib import Path


def run_json(command: list[str], root: Path) -> dict:
    completed = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return json.loads(completed.stdout)


def pair_set(structure: str) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    stacks: dict[str, list[int]] = {symbol: [] for symbol in "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
    closing = {")": "(", "]": "[", "}": "{", ">": "<"}
    closing.update({chr(ord("a") + i): chr(ord("A") + i) for i in range(26)})
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in closing and stacks[closing[symbol]]:
            left = stacks[closing[symbol]].pop()
            pairs.add((left, position))
    return pairs


def planar_structure(structure: str) -> str:
    result = ["."] * len(structure)
    stack: list[int] = []
    for index, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")" and stack:
            left = stack.pop()
            result[left] = "("
            result[index] = ")"
    return "".join(result)


def distance_signature(points: list[dict]) -> list[float]:
    backbone = statistics.fmean(
        math.hypot(b["x"] - a["x"], b["y"] - a["y"])
        for a, b in zip(points, points[1:])
    )
    return [
        math.hypot(points[j]["x"] - points[i]["x"], points[j]["y"] - points[i]["y"])
        / backbone
        for i in range(len(points))
        for j in range(i + 1, len(points))
    ]


def layout_rms(reference: dict, result: dict) -> float:
    restored = [
        {"x": point["x"] * result["aspect_ratio"], "y": point["y"]}
        for point in result["points"]
    ]
    left = distance_signature(reference["points"])
    right = distance_signature(restored)
    return math.sqrt(statistics.fmean((a - b) ** 2 for a, b in zip(left, right)))


def normalized_points(scene: dict, ribon: bool, size: int = 512) -> list[tuple[float, float]]:
    points = scene["points"]
    if ribon:
        points = [
            {"x": point["x"] * scene["aspect_ratio"], "y": point["y"]}
            for point in points
        ]
    minimum_x = min(point["x"] for point in points)
    maximum_x = max(point["x"] for point in points)
    minimum_y = min(point["y"] for point in points)
    maximum_y = max(point["y"] for point in points)
    span_x = max(maximum_x - minimum_x, 1.0e-12)
    span_y = max(maximum_y - minimum_y, 1.0e-12)
    scale = (size - 48) / max(span_x, span_y)
    offset_x = (size - span_x * scale) / 2
    offset_y = (size - span_y * scale) / 2
    return [
        (
            offset_x + (point["x"] - minimum_x) * scale,
            offset_y + (point["y"] - minimum_y) * scale,
        )
        for point in points
    ]


def scene_svg(points: list[tuple[float, float]], pairs: set[tuple[int, int]]) -> str:
    elements = ['<rect width="512" height="512" fill="white"/>']
    for first, second in zip(points, points[1:]):
        elements.append(
            f'<line x1="{first[0]:.8f}" y1="{first[1]:.8f}" '
            f'x2="{second[0]:.8f}" y2="{second[1]:.8f}" '
            'stroke="#666" stroke-width="1.25" stroke-linecap="round"/>'
        )
    for left, right in sorted(pairs):
        first = points[left - 1]
        second = points[right - 1]
        elements.append(
            f'<line x1="{first[0]:.8f}" y1="{first[1]:.8f}" '
            f'x2="{second[0]:.8f}" y2="{second[1]:.8f}" '
            'stroke="#315eaa" stroke-width="1.5" stroke-linecap="round"/>'
        )
    for x, y in points:
        elements.append(
            f'<circle cx="{x:.8f}" cy="{y:.8f}" r="1.45" '
            'fill="white" stroke="#333" stroke-width="0.75"/>'
        )
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" width="512" height="512" '
        'viewBox="0 0 512 512">' + "".join(elements) + "</svg>\n"
    )


def image_rmse(
    reference: dict,
    result: dict,
    pairs: set[tuple[int, int]],
    root: Path,
) -> float:
    with tempfile.TemporaryDirectory(prefix="ribon-layout-") as directory:
        temporary = Path(directory)
        reference_svg = temporary / "reference.svg"
        result_svg = temporary / "result.svg"
        reference_png = temporary / "reference.png"
        result_png = temporary / "result.png"
        reference_svg.write_text(scene_svg(normalized_points(reference, False), pairs))
        result_svg.write_text(scene_svg(normalized_points(result, True), pairs))
        for source, destination in (
            (reference_svg, reference_png),
            (result_svg, result_png),
        ):
            subprocess.run(
                ["rsvg-convert", "-o", str(destination), str(source)],
                cwd=root,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
        comparison = subprocess.run(
            [
                "magick",
                "compare",
                "-metric",
                "RMSE",
                str(reference_png),
                str(result_png),
                "null:",
            ],
            cwd=root,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if comparison.returncode not in (0, 1):
            raise RuntimeError(comparison.stderr)
        metric = comparison.stderr.strip()
        return float(metric.rsplit("(", 1)[1].rstrip(")")) if "(" in metric else 0.0


def prediction_accuracy(predicted: str, reference: str) -> dict:
    predicted_pairs = pair_set(predicted)
    reference_pairs = pair_set(reference)
    true_positive = len(predicted_pairs & reference_pairs)
    precision = true_positive / len(predicted_pairs) if predicted_pairs else float(not reference_pairs)
    sensitivity = true_positive / len(reference_pairs) if reference_pairs else float(not predicted_pairs)
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    return {"precision": precision, "sensitivity": sensitivity, "f1": f1}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vienna_build", type=Path)
    parser.add_argument("--vienna-source", type=Path, default=Path("artifacts/ViennaRNA"))
    parser.add_argument("--corpus", type=Path, default=Path("tests/data/rfam_real_24.json"))
    parser.add_argument("--output", type=Path, default=Path("tests/reports/rfam_validation.json"))
    # The public MFE API and vrna_ep_t probability entries are float-valued;
    # the default partition recurrence itself uses double precision. Ribon
    # keeps all corresponding recurrences in f64.
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--energy-tolerance",
        type=float,
        default=1.0e-5,
        help="tolerance for ViennaRNA public MFE/free-energy values",
    )
    parser.add_argument(
        "--aggregate-tolerance",
        type=float,
        default=1.0e-6,
        help="tolerance for aggregate ensemble observables",
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    corpus = json.loads((root / arguments.corpus).read_text())

    subprocess.run(
        ["cargo", "build", "--release", "--offline", "-p", "ribon-core", "--examples"],
        cwd=root,
        check=True,
    )
    first_sequence = corpus["cases"][0]["sequence"]
    subprocess.run(
        [
            "bash",
            "tests/oracles/vienna/compare-vienna-fold.sh",
            str(arguments.vienna_build),
            first_sequence,
            str(arguments.vienna_source),
            "2",
        ],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "bash",
            "tests/oracles/vienna/compare-vienna-modern-layout.sh",
            str(arguments.vienna_build),
            "turtle",
            "(((...)))",
            str(arguments.vienna_source),
        ],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    ours_analyze = root / "target/release/examples/analyze"
    ours_layout = root / "target/release/examples/layout"
    vienna_fold = arguments.vienna_build / "ribon-vienna-fold"
    vienna_layout = arguments.vienna_build / "ribon-vienna-modern-layout"
    results = []
    probability_checks = 0
    for case in corpus["cases"]:
        sequence = case["sequence"]
        ours = run_json([str(ours_analyze), sequence], root)
        reference = run_json([str(vienna_fold), sequence, "2"], root)
        mfe_difference = abs(ours["mfe_energy_kcal_mol"] - reference["mfe_energy"])
        ensemble_difference = abs(
            ours["ensemble_free_energy_kcal_mol"] - reference["ensemble_free_energy"]
        )
        if (
            mfe_difference > arguments.energy_tolerance
            or ensemble_difference > arguments.energy_tolerance
        ):
            raise AssertionError(
                f"{case['accession']}: energy differences {mfe_difference}, {ensemble_difference}"
            )
        reference_pairs = {
            (entry["i"], entry["j"]): entry["p"]
            for entry in reference["pair_probabilities"]
            if entry["p"] > 1.0e-8
        }
        ours_pairs = {
            (entry["i"], entry["j"]): entry["probability"]
            for entry in ours["pair_probabilities"]
            if entry["probability"] > 1.0e-8
        }
        probability_difference = 0.0
        for pair in reference_pairs.keys() | ours_pairs.keys():
            probability_checks += 1
            difference = abs(reference_pairs.get(pair, 0.0) - ours_pairs.get(pair, 0.0))
            probability_difference = max(probability_difference, difference)
            if difference > arguments.tolerance:
                raise AssertionError(f"{case['accession']}: pair {pair} differs by {difference}")

        mean_distance_difference = abs(
            ours["ensemble"]["mean_base_pair_distance"] - reference["mean_base_pair_distance"]
        )
        entropy_difference = max(
            abs(a - b)
            for a, b in zip(
                ours["ensemble"]["positional_entropy_bits"],
                reference["positional_entropy_bits"],
            )
        )
        if (
            mean_distance_difference > arguments.aggregate_tolerance
            or entropy_difference > arguments.tolerance
        ):
            raise AssertionError(
                f"{case['accession']}: ensemble observable differences "
                f"{mean_distance_difference}, {entropy_difference}"
            )

        layout_metrics = {}
        scaffold = planar_structure(case["reference_structure"])
        scaffold_pairs = pair_set(scaffold)
        for method in ("turtle", "puzzler"):
            vienna_scene = run_json([str(vienna_layout), method, scaffold], root)
            ribon_scene = run_json(
                [str(ours_layout), sequence, case["reference_structure"], method], root
            )
            reference_finite = vienna_scene.get("finite", True)
            layout_metrics[method] = {
                "vienna_reference_finite": reference_finite,
                "distance_signature_rms": (
                    layout_rms(vienna_scene, ribon_scene) if reference_finite else None
                ),
                "standardized_image_rmse": (
                    image_rmse(vienna_scene, ribon_scene, scaffold_pairs, root)
                    if reference_finite
                    else None
                ),
                "ribon_crossings": ribon_scene["crossings"],
            }
            if (
                method == "turtle"
                and reference_finite
                and layout_metrics[method]["distance_signature_rms"] > 5.0e-6
            ):
                raise AssertionError(
                    f"{case['accession']}: Turtle geometry RMS "
                    f"{layout_metrics[method]['distance_signature_rms']} exceeds 5e-6"
                )
            if (
                method == "turtle"
                and reference_finite
                and layout_metrics[method]["standardized_image_rmse"] > 1.0e-4
            ):
                raise AssertionError(
                    f"{case['accession']}: Turtle standardized image RMSE "
                    f"{layout_metrics[method]['standardized_image_rmse']} exceeds 1e-4"
                )

        results.append(
            {
                "accession": case["accession"],
                "family_id": case["family_id"],
                "length": case["length"],
                "mfe_energy_abs_difference": mfe_difference,
                "ensemble_energy_abs_difference": ensemble_difference,
                "max_pair_probability_abs_difference": probability_difference,
                "mean_base_pair_distance_abs_difference": mean_distance_difference,
                "max_positional_entropy_abs_difference": entropy_difference,
                "mfe_structure_tie": ours["mfe_structure"] != reference["mfe_structure"],
                "prediction_against_rfam": prediction_accuracy(
                    ours["mfe_structure"], case["reference_structure"]
                ),
                "layouts": layout_metrics,
            }
        )

    report = {
        "corpus": arguments.corpus.as_posix(),
        "cases": len(results),
        "probability_checks": probability_checks,
        "numeric_tolerance": arguments.tolerance,
        "energy_numeric_tolerance": arguments.energy_tolerance,
        "aggregate_numeric_tolerance": arguments.aggregate_tolerance,
        "reference_numeric_precision": (
            "ViennaRNA MFE/public probability list=float, partition recurrence=double; Ribon=f64"
        ),
        "max_mfe_energy_abs_difference": max(r["mfe_energy_abs_difference"] for r in results),
        "max_ensemble_energy_abs_difference": max(
            r["ensemble_energy_abs_difference"] for r in results
        ),
        "max_pair_probability_abs_difference": max(
            r["max_pair_probability_abs_difference"] for r in results
        ),
        "max_mean_base_pair_distance_abs_difference": max(
            r["mean_base_pair_distance_abs_difference"] for r in results
        ),
        "max_positional_entropy_abs_difference": max(
            r["max_positional_entropy_abs_difference"] for r in results
        ),
        "mfe_structure_ties": sum(r["mfe_structure_tie"] for r in results),
        "median_prediction_f1_against_rfam": statistics.median(
            r["prediction_against_rfam"]["f1"] for r in results
        ),
        "max_turtle_rms": max(
            r["layouts"]["turtle"]["distance_signature_rms"]
            for r in results
            if r["layouts"]["turtle"]["distance_signature_rms"] is not None
        ),
        "max_puzzler_rms": max(
            r["layouts"]["puzzler"]["distance_signature_rms"]
            for r in results
            if r["layouts"]["puzzler"]["distance_signature_rms"] is not None
        ),
        "max_turtle_image_rmse": max(
            r["layouts"]["turtle"]["standardized_image_rmse"]
            for r in results
            if r["layouts"]["turtle"]["standardized_image_rmse"] is not None
        ),
        "max_puzzler_image_rmse": max(
            r["layouts"]["puzzler"]["standardized_image_rmse"]
            for r in results
            if r["layouts"]["puzzler"]["standardized_image_rmse"] is not None
        ),
        "vienna_nonfinite_layouts": sum(
            not metrics["vienna_reference_finite"]
            for result in results
            for metrics in result["layouts"].values()
        ),
        "results": results,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "results"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
