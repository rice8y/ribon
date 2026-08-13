#!/usr/bin/env python3
"""Validate the independent Rust pseudoknot decoder on the 24-record ShPK set."""

from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path


OPEN = {"(": ")", "[": "]", "{": "}", "<": ">"}


def pairs(structure: str) -> set[tuple[int, int]]:
    stacks = {symbol: [] for symbol in OPEN}
    closing = {value: key for key, value in OPEN.items()}
    result: set[tuple[int, int]] = set()
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in closing:
            opening = closing[symbol]
            if not stacks[opening]:
                raise AssertionError(f"unbalanced structure at {position}")
            result.add((stacks[opening].pop(), position))
    if any(stacks.values()):
        raise AssertionError("unbalanced structure")
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    output = subprocess.run(
        [
            "typst",
            "query",
            "--root",
            ".",
            "tests/typst/pseudoknot_query.typ",
            "<ribon-pseudoknot-validation>",
            "--field",
            "value",
            "--one",
        ],
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    records = json.loads(output)
    rows = []
    for record in records:
        prediction = record["prediction"]
        arbitrary = record["arbitrary"]
        expected = pairs(record["expected-structure"])
        method_structures = {
            "probknot": prediction["structure"],
            "hybrid": prediction["hybrid_structure"],
            "matching_centroid": prediction["matching_centroid_structure"],
            "matching_mea": prediction["matching_mea_structure"],
            "restricted_mfe": prediction["restricted_mfe_structure"],
            "restricted_centroid": prediction["restricted_centroid_structure"],
            "restricted_mea": prediction["restricted_mea_structure"],
        }
        if prediction["restricted_state_count"] < 1:
            raise AssertionError(f"{record['id']}: restricted ensemble has no state")
        exact_state_count = prediction["restricted_state_count_exact"]
        if not isinstance(exact_state_count, str) or not exact_state_count.isdecimal():
            raise AssertionError(f"{record['id']}: exact state count is not decimal")
        if int(exact_state_count) < prediction["restricted_state_count"]:
            raise AssertionError(f"{record['id']}: saturated state count exceeds exact count")
        for value in (
            prediction["restricted_mfe_energy_kcal_mol"],
            prediction["restricted_ensemble_free_energy_kcal_mol"],
            prediction["restricted_log_partition_function"],
        ):
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise AssertionError(f"{record['id']}: non-finite restricted ensemble value")
        pair_mass = [0.0] * len(record["sequence"])
        for pair in prediction["restricted_pair_probabilities"]:
            pair_mass[pair["i"] - 1] += pair["probability"]
            pair_mass[pair["j"] - 1] += pair["probability"]
        for index, (paired, unpaired) in enumerate(
            zip(pair_mass, prediction["restricted_unpaired_probabilities"], strict=True), 1
        ):
            if abs(paired + unpaired - 1.0) > 1.0e-10:
                raise AssertionError(f"{record['id']}: probability mass drift at {index}")
        if not arbitrary["state_space_complete"]:
            raise AssertionError(f"{record['id']}: arbitrary matching state space is incomplete")
        if arbitrary["state_count"] < 1 or arbitrary["state_count_exact"] != str(
            arbitrary["state_count"]
        ):
            raise AssertionError(f"{record['id']}: arbitrary matching state count is invalid")
        if not math.isfinite(arbitrary["log_partition_function"]):
            raise AssertionError(f"{record['id']}: arbitrary matching log PF is non-finite")
        arbitrary_pair_mass = [0.0] * len(record["arbitrary-sequence"])
        for pair in arbitrary["pair_probabilities"]:
            arbitrary_pair_mass[pair["i"] - 1] += pair["probability"]
            arbitrary_pair_mass[pair["j"] - 1] += pair["probability"]
        for index, (paired, unpaired) in enumerate(
            zip(
                arbitrary_pair_mass,
                arbitrary["unpaired_probabilities"],
                strict=True,
            ),
            1,
        ):
            if abs(paired + unpaired - 1.0) > 1.0e-10:
                raise AssertionError(
                    f"{record['id']}: arbitrary probability mass drift at {index}"
                )
        method_metrics = {}
        for name, structure in method_structures.items():
            predicted = pairs(structure)
            true_positive = len(expected & predicted)
            method_metrics[name] = {
                "predicted_pair_count": len(predicted),
                "sensitivity": true_positive / len(expected) if expected else 1.0,
                "precision": true_positive / len(predicted) if predicted else float(not expected),
            }
        rows.append(
            {
                "id": record["id"],
                "length": len(record["sequence"]),
                "reference_pair_count": len(expected),
                "methods": method_metrics,
                "probknot_crossing_count": prediction["crossing_count"],
                "hybrid_crossing_count": prediction["hybrid_crossing_count"],
                "restricted_state_count": prediction["restricted_state_count"],
                "restricted_state_count_exact": exact_state_count,
                "arbitrary_prefix_length": len(record["arbitrary-sequence"]),
                "arbitrary_state_count": arbitrary["state_count"],
            }
        )
    method_summary = {}
    for name in rows[0]["methods"]:
        method_summary[name] = {
            "mean_sensitivity": sum(row["methods"][name]["sensitivity"] for row in rows)
            / len(rows),
            "mean_precision": sum(row["methods"][name]["precision"] for row in rows)
            / len(rows),
        }
    report = {
        "source": "Andronescu-Pop-Condon S-Test ShPK",
        "engine": "independent Rust probability-directed pseudoknot decoder",
        "cases": len(rows),
        "methods": method_summary,
        "probknot_predictions_with_crossing_pairs": sum(
            row["probknot_crossing_count"] > 0 for row in rows
        ),
        "hybrid_predictions_with_crossing_pairs": sum(
            row["hybrid_crossing_count"] > 0 for row in rows
        ),
        "restricted_probability_mass_normalized": True,
        "restricted_values_finite": True,
        "arbitrary_prefix_state_spaces_complete": True,
        "arbitrary_prefix_probability_mass_normalized": True,
        "arbitrary_prefix_minimum_state_count": min(
            row["arbitrary_state_count"] for row in rows
        ),
        "arbitrary_prefix_maximum_state_count": max(
            row["arbitrary_state_count"] for row in rows
        ),
        "results": rows,
    }
    if report["cases"] < 24:
        raise AssertionError("pseudoknot corpus must contain at least 24 cases")
    if method_summary["hybrid"]["mean_sensitivity"] < 0.75:
        raise AssertionError("hybrid pseudoknot mean sensitivity fell below 0.75")
    if method_summary["hybrid"]["mean_precision"] < 0.74:
        raise AssertionError("hybrid pseudoknot mean precision fell below 0.74")
    if report["hybrid_predictions_with_crossing_pairs"] < 20:
        raise AssertionError("hybrid pseudoknot crossing recovery fell below 20/24")
    report_path = root / "tests/reports/pseudoknot_validation.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "results"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
