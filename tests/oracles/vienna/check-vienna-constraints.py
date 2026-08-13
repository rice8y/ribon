#!/usr/bin/env python3
"""Deterministic ViennaRNA/Ribon differential checks for constraints and probing."""

from __future__ import annotations

import argparse
import json
import random
import subprocess
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


def close(left: float, right: float, tolerance: float) -> bool:
    return abs(left - right) <= tolerance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vienna_build", type=Path)
    parser.add_argument("--vienna-source", type=Path, default=Path("artifacts/ViennaRNA"))
    parser.add_argument("--tolerance", type=float, default=5.0e-6)
    parser.add_argument("--random-cases", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260810)
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[3]

    subprocess.run(
        [
            "bash",
            "tests/oracles/vienna/compare-vienna-constraints.sh",
            str(arguments.vienna_build),
            str(arguments.vienna_source),
            "GGGAAACCC",
            "2",
            "-",
            "0",
            "0",
            "-1",
        ],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        ["cargo", "build", "--offline", "-p", "ribon-core", "--example", "constrained"],
        cwd=root,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    reactivities = "0,0.1,0.2,0.4,0.8,1.2,0.7,0.3,null"
    values = [0.0, 0.1, 0.2, 0.4, 0.8, 1.2, 0.7, 0.3, None]
    cases = [
        (
            "hard-pair-unpaired",
            "GGGAAACCC",
            2,
            "(..xxx..)",
            0,
            0,
            -1,
            ["none"],
            {"force-pairs": [{"i": 1, "j": 9}], "force-unpaired": [4, 5, 6]},
        ),
        ("no-lonely-pairs", "GGGAAACCC", 2, "-", 1, 0, -1, ["none"], {"no-lonely-pairs": True}),
        ("no-gu", "GGGAAAUCC", 2, "-", 0, 1, -1, ["none"], {"no-gu": True}),
        ("max-span", "GGGAAACCC", 2, "-", 0, 0, 8, ["none"], {"max-span": 8}),
        (
            "unpaired-energy",
            "GGGAAACCC",
            0,
            "-",
            0,
            0,
            -1,
            ["up", "5", "-2"],
            {"soft": {"unpaired": [{"position": 5, "energy-kcal-mol": -2.0}]}},
        ),
        (
            "pair-energy",
            "GAAAC",
            0,
            "-",
            0,
            0,
            -1,
            ["pair", "1", "5", "-10"],
            {"soft": {"pairs": [{"i": 1, "j": 5, "energy-kcal-mol": -10.0}]}},
        ),
        (
            "shape-deigan",
            "GGGAAACCC",
            2,
            "-",
            0,
            0,
            -1,
            ["deigan", "1.8", "-0.6", reactivities],
            {"probing": {"method": "deigan", "reactivities": values}},
        ),
        (
            "shape-zarringhalam",
            "GGGAAACCC",
            2,
            "-",
            0,
            0,
            -1,
            ["zarringhalam", "0.89", "O", "0.5", reactivities],
            {"probing": {"method": "zarringhalam", "reactivities": values}},
        ),
    ]

    binary = arguments.vienna_build / "ribon-vienna-constraints"
    ours = root / "target/debug/examples/constrained"
    probability_checks = 0
    structure_ties = 0

    rng = random.Random(arguments.seed)
    random_cases = []
    for index in range(arguments.random_cases):
        length = rng.randint(12, 42)
        sequence = "".join(rng.choice("ACGU") for _ in range(length))
        dangles = rng.choice((0, 2))
        no_lp = int(rng.random() < 0.3)
        no_gu = int(rng.random() < 0.2)
        span = -1 if rng.random() < 0.65 else rng.randint(5, length)
        forced_unpaired = [
            position
            for position in range(1, length + 1)
            if rng.random() < 0.08
        ]
        hard = "".join(
            "x" if position in forced_unpaired else "."
            for position in range(1, length + 1)
        )
        config = {
            "force-unpaired": forced_unpaired,
            "no-lonely-pairs": bool(no_lp),
            "no-gu": bool(no_gu),
        }
        if span >= 0:
            config["max-span"] = span
        random_cases.append(
            (
                f"random-{index}",
                sequence,
                dangles,
                hard,
                no_lp,
                no_gu,
                span,
                ["none"],
                config,
            )
        )

    for name, sequence, dangles, hard, no_lp, no_gu, span, mode, config in [
        *cases,
        *random_cases,
    ]:
        reference = run_json(
            [str(binary), sequence, str(dangles), hard, str(no_lp), str(no_gu), str(span), *mode],
            root,
        )
        result = run_json([str(ours), sequence, str(dangles), json.dumps(config)], root)
        if not close(reference["mfe"], result["mfe_energy_kcal_mol"], arguments.tolerance):
            raise AssertionError(f"{name}: MFE {reference['mfe']} != {result['mfe_energy_kcal_mol']}")
        if not close(reference["ensemble"], result["ensemble_free_energy_kcal_mol"], arguments.tolerance):
            raise AssertionError(
                f"{name}: ensemble {reference['ensemble']} != {result['ensemble_free_energy_kcal_mol']}"
            )
        if reference["structure"] != result["mfe_structure"]:
            structure_ties += 1

        # Vienna's probability-list helper drops entries below 1e-12. Compare
        # both implementations above a slightly more conservative threshold so
        # an otherwise irrelevant sparse-output convention cannot fail a run.
        cutoff = 1.0e-10
        reference_pairs = {
            (pair["i"], pair["j"]): pair["p"]
            for pair in reference["pairs"]
            if pair["p"] > cutoff
        }
        result_pairs = {
            (pair["i"], pair["j"]): pair["probability"]
            for pair in result["pair_probabilities"]
            if pair["probability"] > cutoff
        }
        if reference_pairs.keys() != result_pairs.keys():
            raise AssertionError(f"{name}: pair-probability support differs")
        for pair, probability in reference_pairs.items():
            probability_checks += 1
            if not close(probability, result_pairs[pair], arguments.tolerance):
                raise AssertionError(
                    f"{name}: pair {pair} probability {probability} != {result_pairs[pair]}"
                )

    print(
        json.dumps(
            {
                "cases": len(cases),
                "random_cases": len(random_cases),
                "probability_checks": probability_checks,
                "seed": arguments.seed,
                "structure_ties": structure_ties,
                "status": "ok",
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
