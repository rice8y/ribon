#!/usr/bin/env python3
"""Strict 24-Rfam WASM/Typst validation for exact exponential features."""

from __future__ import annotations

import json
import math
import subprocess
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


def pair_set(structure: str) -> set[tuple[int, int]]:
    stack: list[int] = []
    result: set[tuple[int, int]] = set()
    for position, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(position)
        elif symbol == ")":
            if not stack:
                raise AssertionError(f"unbalanced path structure: {structure}")
            result.add((stack.pop(), position))
    if stack:
        raise AssertionError(f"unbalanced path structure: {structure}")
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    output = run(
        [
            "typst",
            "query",
            "--root",
            ".",
            "tests/typst/exact_features_real_query.typ",
            "<ribon-exact-features-real>",
            "--field",
            "value",
            "--one",
        ],
        root,
    )
    rows = json.loads(output)
    if len(rows) != 24 or len({row["accession"] for row in rows}) != 24:
        raise AssertionError("exact-feature query must contain 24 distinct Rfam accessions")

    minimum_states = math.inf
    maximum_states = 0
    maximum_mass_error = 0.0
    maximum_saddle_error = 0.0
    for row in rows:
        sequence = row["sequence"]
        landscape = row["landscape"]
        if not landscape["state-space-complete"]:
            raise AssertionError(f"{row['accession']}: incomplete landscape")
        if landscape["state-count"] <= 0:
            raise AssertionError(f"{row['accession']}: empty landscape")
        minimum_states = min(minimum_states, landscape["state-count"])
        maximum_states = max(maximum_states, landscape["state-count"])
        path = landscape["path"]
        for left, right in zip(path, path[1:], strict=False):
            symmetric_difference = pair_set(left["structure"]) ^ pair_set(right["structure"])
            if len(symmetric_difference) != 1:
                raise AssertionError(
                    f"{row['accession']}: path step is not one pair move: {left} -> {right}"
                )
        saddle_error = abs(
            max(step["energy_kcal_mol"] for step in path) - landscape["saddle-energy"]
        )
        maximum_saddle_error = max(maximum_saddle_error, saddle_error)
        if saddle_error > 1.0e-10:
            raise AssertionError(f"{row['accession']}: saddle mismatch {saddle_error}")

        design = row["inverse-design"]
        if not design["search-complete"]:
            raise AssertionError(f"{row['accession']}: incomplete inverse design")
        if design["candidate-sequence-count"] != 1 or design["evaluated-sequence-count"] != 1:
            raise AssertionError(f"{row['accession']}: fixed-template enumeration was not exact")
        if design["sequence"] != sequence:
            raise AssertionError(f"{row['accession']}: inverse design changed a fixed template")
        if not 0.0 < design["target-probability"] <= 1.0:
            raise AssertionError(f"{row['accession']}: invalid target probability")

        ligand = row["ligand"]
        if not ligand["state-space-complete"] or not 0.0 <= ligand["occupancy"] <= 1.0:
            raise AssertionError(f"{row['accession']}: invalid ligand ensemble")
        paired = [0.0] * len(sequence)
        for pair in ligand["pair-probabilities"]:
            paired[pair["i"] - 1] += pair["probability"]
            paired[pair["j"] - 1] += pair["probability"]
        for index, unpaired in enumerate(ligand["unpaired-probabilities"]):
            maximum_mass_error = max(maximum_mass_error, abs(paired[index] + unpaired - 1.0))
    if maximum_mass_error > 1.0e-9:
        raise AssertionError(f"ligand probability mass error {maximum_mass_error}")

    report = {
        "schema": 1,
        "case_count": len(rows),
        "distinct_rfam_accessions": len({row["accession"] for row in rows}),
        "prefix_length": 10,
        "landscape_state_space_complete": True,
        "minimum_landscape_state_count": minimum_states,
        "maximum_landscape_state_count": maximum_states,
        "maximum_landscape_saddle_error_kcal_mol": maximum_saddle_error,
        "inverse_design_search_complete": True,
        "inverse_design_fixed_template_count_exact": True,
        "ligand_state_space_complete": True,
        "maximum_ligand_probability_mass_error": maximum_mass_error,
        "rows": rows,
    }
    path = root / "tests/reports/exact_features_real_validation.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in report if key != "rows"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
