#!/usr/bin/env python3
"""Black-box differential validation against the public CParty 1.0 CLI.

The executable is treated strictly as an external oracle: this script only
invokes its documented command-line interface and parses its textual output.
No CParty source or object is linked, copied, or inspected by Ribon.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
import subprocess
from pathlib import Path


OPEN = {"(": ")", "[": "]", "{": "}", "<": ">"}
RESULT = re.compile(
    r"^(?P<structure>\S+)\s+\((?P<energy>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\)$"
)
MFE_FREQUENCY = re.compile(
    r"frequency of MFE structure in ensemble:\s*"
    r"(?P<frequency>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
    re.IGNORECASE,
)
GAS_CONSTANT_KCAL = 0.001_987_17
REFERENCE_TEMPERATURE_KELVIN = 310.15


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], cwd: Path, timeout: int = 300) -> str:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
    ).stdout


def pairs(structure: str) -> set[tuple[int, int]]:
    stacks = {opening: [] for opening in OPEN}
    closing = {closing: opening for opening, closing in OPEN.items()}
    result: set[tuple[int, int]] = set()
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in closing:
            opening = closing[symbol]
            if not stacks[opening]:
                raise AssertionError(f"unbalanced black-box structure at {position}: {structure}")
            result.add((stacks[opening].pop(), position))
    if any(stacks.values()):
        raise AssertionError(f"unbalanced black-box structure: {structure}")
    return result


def pair_scores(reference: set[tuple[int, int]], observed: set[tuple[int, int]]) -> dict[str, float]:
    common = len(reference & observed)
    precision = common / len(observed) if observed else float(not reference)
    recall = common / len(reference) if reference else float(not observed)
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"precision": precision, "recall": recall, "f1": f1}


def parse_cparty(output: str) -> tuple[str, float, str, float, float]:
    matches = [RESULT.fullmatch(line.strip()) for line in output.splitlines()]
    results = [match for match in matches if match is not None]
    if len(results) < 2:
        raise AssertionError(f"CParty output did not contain MFE and ensemble rows:\n{output}")
    frequency_match = MFE_FREQUENCY.search(output)
    if frequency_match is None:
        raise AssertionError(f"CParty output did not contain the MFE frequency:\n{output}")
    mfe, ensemble = results[0], results[1]
    frequency = float(frequency_match.group("frequency"))
    if not 0.0 < frequency <= 1.0:
        raise AssertionError(f"invalid CParty MFE frequency {frequency}")
    return (
        mfe.group("structure"),
        float(mfe.group("energy")),
        ensemble.group("structure"),
        float(ensemble.group("energy")),
        frequency,
    )


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise AssertionError("correlation requires equal nontrivial samples")
    mean_left = statistics.fmean(left)
    mean_right = statistics.fmean(right)
    numerator = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right, strict=True))
    denominator = math.sqrt(
        sum((a - mean_left) ** 2 for a in left)
        * sum((b - mean_right) ** 2 for b in right)
    )
    return numerator / denominator if denominator else 1.0


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    mfe_deltas = [float(row["mfe_energy_delta_kcal_mol"]) for row in rows]
    ensemble_deltas = [float(row["ensemble_energy_delta_kcal_mol"]) for row in rows]
    f1 = [float(row["mfe_pair_scores"]["f1"]) for row in rows]  # type: ignore[index]
    ribon_mfe = [float(row["ribon_mfe_energy_kcal_mol"]) for row in rows]
    cparty_mfe = [float(row["cparty_mfe_energy_kcal_mol"]) for row in rows]
    ribon_ensemble = [float(row["ribon_ensemble_energy_kcal_mol"]) for row in rows]
    cparty_ensemble = [float(row["cparty_ensemble_energy_kcal_mol"]) for row in rows]
    frequency_ensemble_deltas = [
        float(row["frequency_implied_ensemble_energy_delta_kcal_mol"]) for row in rows
    ]
    cparty_frequency_ensemble = [
        float(row["cparty_frequency_implied_ensemble_energy_kcal_mol"]) for row in rows
    ]
    cparty_self_residuals = [
        float(row["cparty_reported_ensemble_identity_residual_kcal_mol"]) for row in rows
    ]
    mfe_frequency_deltas = [float(row["mfe_frequency_delta"]) for row in rows]
    return {
        "count": len(rows),
        "exact_mfe_pair_set_count": sum(bool(row["mfe_pair_set_exact"]) for row in rows),
        "mean_mfe_pair_f1": statistics.fmean(f1),
        "minimum_mfe_pair_f1": min(f1),
        "mfe_energy_mae_kcal_mol": statistics.fmean(abs(value) for value in mfe_deltas),
        "mfe_energy_max_abs_kcal_mol": max(abs(value) for value in mfe_deltas),
        "ensemble_energy_mae_kcal_mol": statistics.fmean(abs(value) for value in ensemble_deltas),
        "ensemble_energy_max_abs_kcal_mol": max(abs(value) for value in ensemble_deltas),
        "mfe_energy_pearson": pearson(ribon_mfe, cparty_mfe),
        "ensemble_energy_pearson": pearson(ribon_ensemble, cparty_ensemble),
        "frequency_implied_ensemble_energy_mae_kcal_mol": statistics.fmean(
            abs(value) for value in frequency_ensemble_deltas
        ),
        "frequency_implied_ensemble_energy_max_abs_kcal_mol": max(
            abs(value) for value in frequency_ensemble_deltas
        ),
        "frequency_implied_ensemble_energy_pearson": pearson(
            ribon_ensemble, cparty_frequency_ensemble
        ),
        "cparty_reported_ensemble_identity_residual_mae_kcal_mol": statistics.fmean(
            abs(value) for value in cparty_self_residuals
        ),
        "cparty_reported_ensemble_identity_residual_max_abs_kcal_mol": max(
            abs(value) for value in cparty_self_residuals
        ),
        "mfe_frequency_mae": statistics.fmean(abs(value) for value in mfe_frequency_deltas),
        "mfe_frequency_max_abs": max(abs(value) for value in mfe_frequency_deltas),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cparty-root", type=Path, default=Path("/private/tmp/ribon-cparty"))
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    cparty_root = arguments.cparty_root.resolve()
    executable = cparty_root / "build/CParty"
    rna_parameters = cparty_root / "params/rna_Turner04.par"
    if not executable.is_file() or not rna_parameters.is_file():
        raise FileNotFoundError(
            "CParty 1.0 black-box oracle not found; pass --cparty-root containing "
            "build/CParty and params/rna_Turner04.par"
        )

    query = run(
        [
            "typst",
            "query",
            "--root",
            ".",
            "tests/typst/cparty_blackbox_query.typ",
            "<ribon-cparty-blackbox>",
            "--field",
            "value",
            "--one",
        ],
        root,
    )
    records = json.loads(query)
    if len(records) != 64:
        raise AssertionError(f"expected 64 RNA/DNA/dangle black-box cases, found {len(records)}")

    rows: list[dict[str, object]] = []
    for record in records:
        command = [str(executable), "--noPS", "-r", record["seed"], "-d", str(record["dangles"])]
        if record["family"] == "rna":
            command.extend(["-P", "params/rna_Turner04.par"])
        else:
            command.append("--noConv")
        command.append(record["sequence"])
        output = run(command, cparty_root)
        try:
            (
                cparty_structure,
                cparty_mfe,
                cparty_probability_structure,
                cparty_ensemble,
                cparty_mfe_frequency,
            ) = parse_cparty(output)
        except AssertionError as error:
            raise AssertionError(
                f"{record['id']} ({record['family']}, dangles={record['dangles']}): {error}"
            ) from error
        ribon_pairs = pairs(record["mfe-structure"])
        cparty_pairs = pairs(cparty_structure)
        seed_pairs = pairs(record["seed"])
        if not seed_pairs <= ribon_pairs or not seed_pairs <= cparty_pairs:
            raise AssertionError(f"{record['id']}: a fixed seed pair was lost")
        for value in (
            record["mfe-energy"],
            record["ensemble-energy"],
            cparty_mfe,
            cparty_ensemble,
        ):
            if not math.isfinite(float(value)):
                raise AssertionError(f"{record['id']}: non-finite black-box comparison value")
        if float(record["ensemble-energy"]) > float(record["mfe-energy"]) + 1.0e-9:
            raise AssertionError(f"{record['id']}: Ribon ensemble free energy exceeds MFE")
        if cparty_ensemble > cparty_mfe + 1.0e-6:
            raise AssertionError(f"{record['id']}: CParty ensemble free energy exceeds MFE")
        scores = pair_scores(cparty_pairs, ribon_pairs)
        rt = GAS_CONSTANT_KCAL * REFERENCE_TEMPERATURE_KELVIN
        cparty_frequency_ensemble = cparty_mfe + rt * math.log(cparty_mfe_frequency)
        ribon_mfe_frequency = math.exp(
            (float(record["ensemble-energy"]) - float(record["mfe-energy"])) / rt
        )
        rows.append(
            {
                "id": record["id"],
                "family": record["family"],
                "length": len(record["sequence"]),
                "dangles": record["dangles"],
                "seed_pair_count": len(seed_pairs),
                "ribon_mfe_structure": record["mfe-structure"],
                "cparty_mfe_structure": cparty_structure,
                "cparty_probability_structure": cparty_probability_structure,
                "mfe_pair_set_exact": ribon_pairs == cparty_pairs,
                "mfe_pair_scores": scores,
                "ribon_mfe_energy_kcal_mol": record["mfe-energy"],
                "cparty_mfe_energy_kcal_mol": cparty_mfe,
                "mfe_energy_delta_kcal_mol": float(record["mfe-energy"]) - cparty_mfe,
                "ribon_ensemble_energy_kcal_mol": record["ensemble-energy"],
                "cparty_ensemble_energy_kcal_mol": cparty_ensemble,
                "ensemble_energy_delta_kcal_mol": float(record["ensemble-energy"]) - cparty_ensemble,
                "cparty_mfe_frequency": cparty_mfe_frequency,
                "ribon_mfe_frequency": ribon_mfe_frequency,
                "mfe_frequency_delta": ribon_mfe_frequency - cparty_mfe_frequency,
                "cparty_frequency_implied_ensemble_energy_kcal_mol": cparty_frequency_ensemble,
                "frequency_implied_ensemble_energy_delta_kcal_mol": (
                    float(record["ensemble-energy"]) - cparty_frequency_ensemble
                ),
                "cparty_reported_ensemble_identity_residual_kcal_mol": (
                    cparty_ensemble - cparty_frequency_ensemble
                ),
            }
        )

    groups: dict[str, dict[str, object]] = {}
    for family in ("rna", "dna"):
        for dangles in (0, 2):
            key = f"{family}-dangles-{dangles}"
            groups[key] = summarize(
                [row for row in rows if row["family"] == family and row["dangles"] == dangles]
            )

    # These gates intentionally test correlation and structural agreement, not
    # byte identity: CParty translates Vienna parameter files internally while
    # Ribon consumes independently generated RNAstructure 6.6 tables. Bounds
    # are tight enough to flag a single-pair or centi-kcal regression in the
    # standard RNA dangles=2 path while preserving known oracle/model deltas.
    gates = {
        "rna-dangles-0": {
            "count": 24,
            "exact": 24,
            "mean_f1": 1.0,
            "mfe_mae": 0.25,
            "mfe_max": 4.30,
            "ensemble_mae": 0.40,
            "ensemble_max": 4.30,
            "mfe_r": 0.996,
            "ensemble_r": 0.996,
            "frequency_ensemble_mae": 0.21,
            "frequency_ensemble_max": 4.20,
            "frequency_ensemble_r": 0.997,
            "mfe_frequency_mae": 0.016,
            "mfe_frequency_max": 0.060,
        },
        "rna-dangles-2": {
            "count": 24,
            "exact": 24,
            "mean_f1": 1.0,
            "mfe_mae": 0.06,
            "mfe_max": 0.50,
            "ensemble_mae": 0.25,
            "ensemble_max": 1.00,
            "mfe_r": 0.999,
            "ensemble_r": 0.999,
            "frequency_ensemble_mae": 0.055,
            "frequency_ensemble_max": 0.40,
            "frequency_ensemble_r": 0.9999,
            "mfe_frequency_mae": 0.016,
            "mfe_frequency_max": 0.060,
        },
        "dna-dangles-0": {
            "count": 8,
            "exact": 5,
            "mean_f1": 0.90,
            "mfe_mae": 0.25,
            "mfe_max": 1.70,
            "ensemble_mae": 0.20,
            "ensemble_max": 0.60,
            "mfe_r": 0.997,
            "ensemble_r": 0.999,
            "frequency_ensemble_mae": 0.060,
            "frequency_ensemble_max": 0.11,
            "frequency_ensemble_r": 0.9999,
            "mfe_frequency_mae": 0.016,
            "mfe_frequency_max": 0.060,
        },
        "dna-dangles-2": {
            "count": 8,
            "exact": 5,
            "mean_f1": 0.88,
            "mfe_mae": 0.07,
            "mfe_max": 0.50,
            "ensemble_mae": 0.30,
            "ensemble_max": 0.80,
            "mfe_r": 0.999,
            "ensemble_r": 0.999,
            "frequency_ensemble_mae": 0.23,
            "frequency_ensemble_max": 0.82,
            "frequency_ensemble_r": 0.9995,
            "mfe_frequency_mae": 0.016,
            "mfe_frequency_max": 0.060,
        },
    }
    for key, gate in gates.items():
        group = groups[key]
        checks = {
            "count": int(group["count"]) == gate["count"],
            "exact pair-set count": int(group["exact_mfe_pair_set_count"]) >= gate["exact"],
            "mean pair F1": float(group["mean_mfe_pair_f1"]) + 1.0e-12 >= gate["mean_f1"],
            "MFE MAE": float(group["mfe_energy_mae_kcal_mol"]) <= gate["mfe_mae"],
            "MFE max abs": float(group["mfe_energy_max_abs_kcal_mol"]) <= gate["mfe_max"],
            "ensemble MAE": float(group["ensemble_energy_mae_kcal_mol"])
            <= gate["ensemble_mae"],
            "ensemble max abs": float(group["ensemble_energy_max_abs_kcal_mol"])
            <= gate["ensemble_max"],
            "MFE correlation": float(group["mfe_energy_pearson"]) >= gate["mfe_r"],
            "ensemble correlation": float(group["ensemble_energy_pearson"])
            >= gate["ensemble_r"],
            "frequency-implied ensemble MAE": float(
                group["frequency_implied_ensemble_energy_mae_kcal_mol"]
            )
            <= gate["frequency_ensemble_mae"],
            "frequency-implied ensemble max abs": float(
                group["frequency_implied_ensemble_energy_max_abs_kcal_mol"]
            )
            <= gate["frequency_ensemble_max"],
            "frequency-implied ensemble correlation": float(
                group["frequency_implied_ensemble_energy_pearson"]
            )
            >= gate["frequency_ensemble_r"],
            "MFE frequency MAE": float(group["mfe_frequency_mae"])
            <= gate["mfe_frequency_mae"],
            "MFE frequency max abs": float(group["mfe_frequency_max_abs"])
            <= gate["mfe_frequency_max"],
        }
        failed = [name for name, passed in checks.items() if not passed]
        if failed:
            raise AssertionError(f"{key}: black-box regression gates failed: {failed}; {group}")

    report = {
        "schema": 1,
        "oracle": "CParty 1.0 documented CLI black box",
        "oracle_version": run([str(executable), "--version"], cparty_root).strip(),
        "oracle_executable_sha256": digest(executable),
        "oracle_rna_parameter_sha256": digest(rna_parameters),
        "linked_or_copied_into_release": False,
        "thermodynamic_identity": "G_ensemble = E_MFE + RT ln(p_MFE)",
        "reference_temperature_kelvin": REFERENCE_TEMPERATURE_KELVIN,
        "case_count": len(rows),
        "distinct_biological_sequences": len({row["id"] for row in rows}),
        "groups": groups,
        "gates": gates,
        "rows": rows,
    }
    output = root / "tests/reports/cparty_blackbox_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"case_count": len(rows), "groups": groups}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
