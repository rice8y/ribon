#!/usr/bin/env python3
"""Compare Ribon's independent Rust engine with RNAstructure 6.6.

This is a differential validation harness, not a runtime dependency.  It uses
the official command-line release as an external oracle and records both
structure-level and numerical ensemble metrics on 24 diverse Rfam sequences.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests" / "data" / "rfam_real_24.json"
DEFAULT_REFERENCE = Path("/private/tmp/ribon-rnastructure-6.6")
DATAPATH_RELATIVE = Path("share/rnastructure/data_tables")
EXPECTED_ARCHIVE_SHA256 = (
    "8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8"
)


def run(arguments: list[str], *, env: dict[str, str], stdin: str | None = None) -> str:
    completed = subprocess.run(
        arguments,
        input=stdin,
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(arguments)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed.stdout


def dot_pairs(structure: str) -> set[tuple[int, int]]:
    stacks: dict[str, list[int]] = {"(": [], "[": [], "{": [], "<": []}
    close = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs: set[tuple[int, int]] = set()
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in close:
            opening = close[symbol]
            if not stacks[opening]:
                raise ValueError(f"unbalanced structure: {structure}")
            pairs.add((stacks[opening].pop(), position))
        elif "A" <= symbol <= "Z":
            stacks.setdefault(symbol, []).append(position)
        elif "a" <= symbol <= "z":
            opening = symbol.upper()
            if not stacks.get(opening):
                raise ValueError(f"unbalanced structure: {structure}")
            pairs.add((stacks[opening].pop(), position))
    if any(stacks.values()):
        raise ValueError(f"unbalanced structure: {structure}")
    return pairs


def pair_scores(predicted: str, reference: str) -> tuple[float, float]:
    observed = dot_pairs(predicted)
    expected = dot_pairs(reference)
    hits = len(observed & expected)
    sensitivity = 1.0 if not expected else hits / len(expected)
    precision = 1.0 if not observed else hits / len(observed)
    return sensitivity, precision


def dbn_from_output(output: str) -> tuple[str, float | None]:
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    structure = next(
        line.split()[0]
        for line in reversed(lines)
        if set(line.split()[0]) <= set(".()[]{}<>")
    )
    match = re.search(r"ENERGY\s*=\s*([+-]?[0-9.]+)", output)
    return structure, float(match.group(1)) if match else None


def reference_case(
    case: dict[str, object], reference: Path, env: dict[str, str], work: Path
) -> dict[str, object]:
    accession = str(case["accession"])
    sequence = str(case["sequence"])[:80]
    fasta = work / f"{accession}.fa"
    pfs = work / f"{accession}.pfs"
    dotplot = work / f"{accession}.dp"
    mea_ct = work / f"{accession}-mea.ct"
    centroid_ct = work / f"{accession}-centroid.ct"
    fasta.write_text(f">{accession}\n{sequence}\n", encoding="ascii")

    bin_dir = reference / "bin"
    mfe_output = run(
        [
            str(bin_dir / "Fold"),
            str(fasta),
            "-",
            "--bracket",
            "--MFE",
            "--isolated",
        ],
        env=env,
    )
    mfe_structure, mfe_energy = dbn_from_output(mfe_output)
    if mfe_energy is None:
        raise RuntimeError(f"Fold did not report an energy for {accession}")

    run(
        [str(bin_dir / "partition"), str(fasta), str(pfs), "--isolated", "--quiet"],
        env=env,
    )
    ensemble_output = run(
        [str(bin_dir / "EnsembleEnergy"), str(pfs), "--silent"], env=env
    )
    match = re.search(r":\s*([+-]?[0-9.]+)\s+kcal/mol", ensemble_output)
    if not match:
        raise RuntimeError(f"EnsembleEnergy did not report an energy for {accession}")
    ensemble_energy = float(match.group(1))

    run(
        [str(bin_dir / "ProbabilityPlot"), str(pfs), str(dotplot), "--text"],
        env=env,
    )
    probabilities: dict[tuple[int, int], float] = {}
    for line in dotplot.read_text(encoding="utf-8").splitlines()[2:]:
        fields = line.split()
        if len(fields) == 3:
            probabilities[(int(fields[0]), int(fields[1]))] = 10.0 ** -float(fields[2])

    run(
        [
            str(bin_dir / "MaxExpect"),
            str(pfs),
            str(mea_ct),
            "--structures",
            "1",
            "--window",
            "0",
            "--gamma",
            "1.0",
        ],
        env=env,
    )
    mea_output = run(
        [str(bin_dir / "ct2dot"), str(mea_ct), "1", "-"], env=env
    )
    mea_structure, _ = dbn_from_output(mea_output)

    run(
        [
            str(bin_dir / "ProbablePair"),
            str(pfs),
            str(centroid_ct),
            "--threshold",
            "0.5",
        ],
        env=env,
    )
    centroid_output = run(
        [str(bin_dir / "ct2dot"), str(centroid_ct), "1", "-"], env=env
    )
    centroid_structure, _ = dbn_from_output(centroid_output)
    return {
        "accession": accession,
        "length": len(sequence),
        "sequence": sequence,
        "mfe_structure": mfe_structure,
        "mfe_energy_kcal_mol": mfe_energy,
        "ensemble_free_energy_kcal_mol": ensemble_energy,
        "centroid_structure": centroid_structure,
        "mea_structure": mea_structure,
        "pair_probabilities": [
            {"i": i, "j": j, "probability": probability}
            for (i, j), probability in sorted(probabilities.items())
        ],
    }


def ribon_case(case: dict[str, object]) -> dict[str, object]:
    sequence = str(case["sequence"])
    analyze = ROOT / "target" / "debug" / "examples" / "analyze"
    fold = ROOT / "target" / "debug" / "examples" / "fold"
    analysis = json.loads(run([str(analyze), sequence, "2", "1.021"], env=os.environ.copy()))
    mfe_by_dangles = {
        str(dangles): json.loads(
            run([str(fold), sequence, str(dangles)], env=os.environ.copy())
        )
        for dangles in range(4)
    }
    mfe = mfe_by_dangles["3"]
    analysis["mfe_structure"] = mfe["structure"]
    analysis["mfe_energy_kcal_mol"] = mfe["energy_kcal_mol"]
    analysis["mfe_by_dangles"] = mfe_by_dangles
    return analysis


def compare(ours: dict[str, object], reference: dict[str, object]) -> dict[str, object]:
    mfe_sensitivity, mfe_precision = pair_scores(
        str(ours["mfe_structure"]), str(reference["mfe_structure"])
    )
    centroid_sensitivity, centroid_precision = pair_scores(
        str(ours["centroid_structure"]), str(reference["centroid_structure"])
    )
    mea_sensitivity, mea_precision = pair_scores(
        str(ours["mea_structure"]), str(reference["mea_structure"])
    )
    observed = {
        (int(entry["i"]), int(entry["j"])): float(entry["probability"])
        for entry in ours["pair_probabilities"]
    }
    expected = {
        (int(entry["i"]), int(entry["j"])): float(entry["probability"])
        for entry in reference["pair_probabilities"]
    }
    coordinates = set(observed) | set(expected)
    pair_errors = sorted(
        (
            (
                abs(observed.get(pair, 0.0) - expected.get(pair, 0.0)),
                pair,
                observed.get(pair, 0.0),
                expected.get(pair, 0.0),
            )
            for pair in coordinates
        ),
        reverse=True,
    )
    errors = [entry[0] for entry in pair_errors]
    compact_reference = {
        key: value
        for key, value in reference.items()
        if key != "pair_probabilities"
    }
    dangle_metrics = {}
    for dangles, result in ours["mfe_by_dangles"].items():
        sensitivity, precision = pair_scores(
            str(result["structure"]), str(reference["mfe_structure"])
        )
        dangle_metrics[dangles] = {
            "energy_absolute_error_kcal_mol": abs(
                float(result["energy_kcal_mol"])
                - float(reference["mfe_energy_kcal_mol"])
            ),
            "pair_sensitivity": sensitivity,
            "pair_precision": precision,
        }
    return {
        "accession": reference["accession"],
        "length": reference["length"],
        "mfe_energy_absolute_error_kcal_mol": abs(
            float(ours["mfe_energy_kcal_mol"])
            - float(reference["mfe_energy_kcal_mol"])
        ),
        "ensemble_energy_absolute_error_kcal_mol": abs(
            float(ours["ensemble_free_energy_kcal_mol"])
            - float(reference["ensemble_free_energy_kcal_mol"])
        ),
        "mfe_pair_sensitivity": mfe_sensitivity,
        "mfe_pair_precision": mfe_precision,
        "centroid_pair_sensitivity": centroid_sensitivity,
        "centroid_pair_precision": centroid_precision,
        "mea_pair_sensitivity": mea_sensitivity,
        "mea_pair_precision": mea_precision,
        "pair_probability_mean_absolute_error": sum(errors) / max(1, len(errors)),
        "pair_probability_maximum_absolute_error": max(errors, default=0.0),
        "largest_pair_probability_errors": [
            {
                "i": pair[0],
                "j": pair[1],
                "absolute_error": error,
                "ribon_probability": observed_probability,
                "rnastructure_probability": reference_probability,
            }
            for error, pair, observed_probability, reference_probability in pair_errors[:20]
        ],
        "rnastructure_significant_pair_probabilities": [
            {"i": pair[0], "j": pair[1], "probability": probability}
            for pair, probability in sorted(expected.items())
            if probability >= 0.01
        ],
        "mfe_by_dangles": dangle_metrics,
        "ribon": {
            "mfe_structure": ours["mfe_structure"],
            "mfe_energy_kcal_mol": ours["mfe_energy_kcal_mol"],
            "ensemble_free_energy_kcal_mol": ours["ensemble_free_energy_kcal_mol"],
            "centroid_structure": ours["centroid_structure"],
            "mea_structure": ours["mea_structure"],
        },
        "rnastructure": compact_reference,
    }


def mean(cases: list[dict[str, object]], key: str) -> float:
    return sum(float(case[key]) for case in cases) / len(cases)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "tests" / "reports" / "rnastructure_6_6_validation.json",
    )
    parser.add_argument("--enforce", action="store_true")
    arguments = parser.parse_args()

    required = ["Fold", "partition", "EnsembleEnergy", "ProbabilityPlot", "MaxExpect"]
    missing = [name for name in required if not (arguments.reference / "bin" / name).is_file()]
    if missing:
        raise SystemExit(
            f"RNAstructure 6.6 reference is missing {', '.join(missing)} under "
            f"{arguments.reference}; pass --reference"
        )
    subprocess.run(
        [
            "cargo",
            "build",
            "-p",
            "ribon-core",
            "--example",
            "analyze",
            "--example",
            "fold",
            "--offline",
        ],
        cwd=ROOT,
        check=True,
    )
    env = os.environ.copy()
    env["DATAPATH"] = str(arguments.reference / DATAPATH_RELATIVE)
    corpus = json.loads(CORPUS.read_text(encoding="utf-8"))["cases"]
    comparisons: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ribon-rnastructure-") as temporary:
        work = Path(temporary)
        for source in corpus:
            query = dict(source)
            query["sequence"] = str(source["sequence"])[:80]
            reference = reference_case(query, arguments.reference, env, work)
            comparisons.append(compare(ribon_case(query), reference))

    metrics = {
        "case_count": len(comparisons),
        "mean_mfe_energy_absolute_error_kcal_mol": mean(
            comparisons, "mfe_energy_absolute_error_kcal_mol"
        ),
        "maximum_mfe_energy_absolute_error_kcal_mol": max(
            float(case["mfe_energy_absolute_error_kcal_mol"]) for case in comparisons
        ),
        "mean_ensemble_energy_absolute_error_kcal_mol": mean(
            comparisons, "ensemble_energy_absolute_error_kcal_mol"
        ),
        "mean_mfe_pair_sensitivity": mean(comparisons, "mfe_pair_sensitivity"),
        "mean_mfe_pair_precision": mean(comparisons, "mfe_pair_precision"),
        "mean_centroid_pair_sensitivity": mean(comparisons, "centroid_pair_sensitivity"),
        "mean_centroid_pair_precision": mean(comparisons, "centroid_pair_precision"),
        "mean_mea_pair_sensitivity": mean(comparisons, "mea_pair_sensitivity"),
        "mean_mea_pair_precision": mean(comparisons, "mea_pair_precision"),
        "mean_pair_probability_absolute_error": mean(
            comparisons, "pair_probability_mean_absolute_error"
        ),
        "maximum_pair_probability_absolute_error": max(
            float(case["pair_probability_maximum_absolute_error"])
            for case in comparisons
        ),
    }
    metrics["mfe_dangle_comparison"] = {
        dangles: {
            "mean_energy_absolute_error_kcal_mol": sum(
                float(case["mfe_by_dangles"][dangles]["energy_absolute_error_kcal_mol"])
                for case in comparisons
            )
            / len(comparisons),
            "mean_pair_sensitivity": sum(
                float(case["mfe_by_dangles"][dangles]["pair_sensitivity"])
                for case in comparisons
            )
            / len(comparisons),
            "mean_pair_precision": sum(
                float(case["mfe_by_dangles"][dangles]["pair_precision"])
                for case in comparisons
            )
            / len(comparisons),
        }
        for dangles in ("0", "1", "2", "3")
    }
    report = {
        "schema_version": 1,
        "reference": "RNAstructure 6.6 official osx-arm64 conda package",
        "reference_archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "corpus": "Rfam CC0 24-family corpus; first 80 nt per family",
        "ribon_mfe_dangles": 3,
        "ribon_partition_dangles": 2,
        "rnastructure_isolated_pairs": True,
        "metrics": metrics,
        "cases": comparisons,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))

    if arguments.enforce:
        gates = [
            (metrics["case_count"] >= 24, "fewer than 24 real-data cases"),
            (metrics["mean_mfe_pair_sensitivity"] >= 0.79, "MFE sensitivity below 0.79"),
            (metrics["mean_mfe_pair_precision"] >= 0.81, "MFE precision below 0.81"),
            (
                metrics["mean_centroid_pair_sensitivity"] >= 0.88,
                "centroid sensitivity below 0.88",
            ),
            (
                metrics["mean_centroid_pair_precision"] >= 0.90,
                "centroid precision below 0.90",
            ),
            (metrics["mean_mea_pair_sensitivity"] >= 0.87, "MEA sensitivity below 0.87"),
            (metrics["mean_mea_pair_precision"] >= 0.88, "MEA precision below 0.88"),
            (
                metrics["mean_mfe_energy_absolute_error_kcal_mol"] <= 1.1,
                "mean MFE energy error above 1.1 kcal/mol",
            ),
            (
                metrics["maximum_mfe_energy_absolute_error_kcal_mol"] <= 3.0,
                "maximum MFE energy error above 3.0 kcal/mol",
            ),
            (
                metrics["mean_ensemble_energy_absolute_error_kcal_mol"] <= 1.6,
                "mean ensemble energy error above 1.6 kcal/mol",
            ),
            (
                metrics["mean_pair_probability_absolute_error"] <= 0.01,
                "pair-probability MAE above 0.01",
            ),
        ]
        failures = [message for passed, message in gates if not passed]
        if failures:
            raise SystemExit("reference gates failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
