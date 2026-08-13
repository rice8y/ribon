#!/usr/bin/env python3
"""Differential validation of exact k-best enumeration against ViennaRNA."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


CASES = [
    ("GGGAAACCC", 5.0, 2),
    ("GGGAAACCC", 5.0, 0),
    ("GCAUACGUC", 4.0, 2),
    ("AUGGCUACGAU", 3.0, 2),
    ("GCGCAAAAGCGC", 3.0, 0),
    ("GGAUCCAAAGGAUCC", 2.0, 2),
]


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vienna_build", type=Path)
    parser.add_argument("--vienna-source", type=Path, default=Path("artifacts/ViennaRNA"))
    parser.add_argument(
        "--output", type=Path, default=Path("tests/reports/suboptimal_validation.json")
    )
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    subprocess.run(
        ["cargo", "build", "--release", "--offline", "-p", "ribon-core", "--example", "suboptimal-json"],
        cwd=root,
        check=True,
    )
    first = CASES[0]
    subprocess.run(
        [
            "bash",
            "tests/oracles/vienna/compare-vienna-suboptimal.sh",
            str(arguments.vienna_build),
            first[0],
            str(first[1]),
            str(first[2]),
            str(arguments.vienna_source),
        ],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    vienna_binary = arguments.vienna_build / "ribon-vienna-suboptimal"
    ours_binary = root / "target/release/examples/suboptimal-json"

    results = []
    maximum_energy_difference = 0.0
    total_structures = 0
    for sequence, band, dangles in CASES:
        reference = run_json(
            [str(vienna_binary), sequence, str(band), str(dangles)], root
        )
        ours = run_json(
            [str(ours_binary), sequence, str(band), str(dangles), "200"], root
        )
        if ours["truncated"]:
            raise AssertionError(f"{sequence}: Ribon result truncated")
        reference_map = {
            entry["structure"]: entry["energy"] for entry in reference["structures"]
        }
        ours_map = {
            entry["structure"]: entry["energy_kcal_mol"] for entry in ours["structures"]
        }
        if reference_map.keys() != ours_map.keys():
            missing = sorted(reference_map.keys() - ours_map.keys())
            extra = sorted(ours_map.keys() - reference_map.keys())
            raise AssertionError(f"{sequence}: missing={missing}, extra={extra}")
        case_maximum = max(
            abs(reference_map[structure] - ours_map[structure])
            for structure in reference_map
        )
        if case_maximum > 1.0e-6:
            raise AssertionError(f"{sequence}: energy difference {case_maximum}")
        maximum_energy_difference = max(maximum_energy_difference, case_maximum)
        total_structures += len(reference_map)
        results.append(
            {
                "sequence": sequence,
                "energy_band_kcal_mol": band,
                "dangles": dangles,
                "structures": len(reference_map),
                "max_energy_abs_difference": case_maximum,
            }
        )

    report = {
        "reference": "ViennaRNA 2.7.2 vrna_subopt (Wuchty), uniq_ML=true",
        "cases": len(results),
        "structures_compared": total_structures,
        "max_energy_abs_difference": maximum_energy_difference,
        "results": results,
    }
    output = root / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "results"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
