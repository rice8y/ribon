#!/usr/bin/env python3
"""Differential fuzzing for ViennaRNA/Ribon odd-dangle MFE and evaluation."""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path


def run_json(command: list[str], cwd: Path) -> dict:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return json.loads(completed.stdout)


def close(actual: float, expected: float, tolerance: float) -> bool:
    return abs(actual - expected) <= tolerance


def loop_map(entries: list[dict], energy_key: str) -> dict[tuple[int, int], float]:
    return {(entry["i"], entry["j"]): entry[energy_key] for entry in entries}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vienna_build", type=Path)
    parser.add_argument("--vienna-source", type=Path, default=Path("artifacts/ViennaRNA"))
    parser.add_argument("--cases", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--min-length", type=int, default=18)
    parser.add_argument("--max-length", type=int, default=80)
    # ViennaRNA's public API returns `float`; 5e-6 covers its serialization
    # round-off while remaining below one Turner centi-kcal unit by 2000x.
    parser.add_argument("--tolerance", type=float, default=5.0e-6)
    arguments = parser.parse_args()

    root = Path(__file__).resolve().parents[3]
    subprocess.run(
        [
            "bash",
            "tests/oracles/vienna/compare-vienna-odd.sh",
            str(arguments.vienna_build),
            "GGGAAACCC",
            "1",
            str(arguments.vienna_source),
        ],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        ["cargo", "build", "--offline", "-p", "ribon-core", "--example", "fold"],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    vienna = arguments.vienna_build / "ribon-vienna-odd"
    ribon = root / "target/debug/examples/fold"
    randomizer = random.Random(arguments.seed)
    structure_ties = 0
    comparisons = 0

    for case_index in range(arguments.cases):
        length = randomizer.randint(arguments.min_length, arguments.max_length)
        sequence = "".join(randomizer.choice("GGGCCCAAUU") for _ in range(length))
        for dangles in (1, 3):
            reference = run_json([str(vienna), sequence, str(dangles)], root)
            folded = run_json([str(ribon), sequence, str(dangles)], root)
            comparisons += 1

            if not close(
                folded["energy_kcal_mol"],
                reference["mfe_energy"],
                arguments.tolerance,
            ):
                raise AssertionError(
                    f"MFE energy mismatch case={case_index} d={dangles}: "
                    f"sequence={sequence} Vienna={reference['mfe_energy']} "
                    f"Ribon={folded['energy_kcal_mol']}"
                )
            if folded["structure"] != reference["mfe_structure"]:
                structure_ties += 1

            # Evaluate Vienna's chosen structure with Ribon and compare every
            # loop, which remains meaningful even when an MFE tie chooses a
            # different dot-bracket string.
            evaluated = run_json(
                [str(ribon), sequence, str(dangles), reference["mfe_structure"]], root
            )
            if not close(
                evaluated["total_kcal_mol"],
                reference["evaluated_energy"],
                arguments.tolerance,
            ):
                raise AssertionError(
                    f"evaluation mismatch case={case_index} d={dangles}: "
                    f"sequence={sequence} structure={reference['mfe_structure']} "
                    f"Vienna={reference['evaluated_energy']} "
                    f"Ribon={evaluated['total_kcal_mol']}"
                )

            reference_loops = loop_map(reference["loop_energies"], "energy")
            ribon_loops = loop_map(evaluated["loop_energies"], "energy_kcal_mol")
            if reference_loops.keys() != ribon_loops.keys():
                raise AssertionError(
                    f"loop key mismatch case={case_index} d={dangles}: "
                    f"Vienna={sorted(reference_loops)} Ribon={sorted(ribon_loops)}"
                )
            for key, expected in reference_loops.items():
                if not close(ribon_loops[key], expected, arguments.tolerance):
                    raise AssertionError(
                        f"loop mismatch case={case_index} d={dangles} loop={key}: "
                        f"sequence={sequence} structure={reference['mfe_structure']} "
                        f"Vienna={expected} Ribon={ribon_loops[key]}"
                    )

            # Cross-evaluate Ribon's chosen structure in ViennaRNA too.
            reference_for_ribon = run_json(
                [str(vienna), sequence, str(dangles), folded["structure"]], root
            )
            if not close(
                folded["evaluated_energy_kcal_mol"],
                reference_for_ribon["evaluated_energy"],
                arguments.tolerance,
            ):
                raise AssertionError(
                    f"cross-evaluation mismatch case={case_index} d={dangles}: "
                    f"sequence={sequence} structure={folded['structure']} "
                    f"Vienna={reference_for_ribon['evaluated_energy']} "
                    f"Ribon={folded['evaluated_energy_kcal_mol']}"
                )

        if (case_index + 1) % 100 == 0:
            print(f"checked {case_index + 1}/{arguments.cases} sequences", file=sys.stderr)

    print(
        json.dumps(
            {
                "seed": arguments.seed,
                "sequences": arguments.cases,
                "dangle_comparisons": comparisons,
                "structure_ties": structure_ties,
                "status": "ok",
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
