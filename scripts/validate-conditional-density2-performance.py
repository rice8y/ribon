#!/usr/bin/env python3
"""Validate the production density-2 engine's cubic scaling contract."""

from __future__ import annotations

import json
import math
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "tests/reports/conditional_density2_performance.json"
MEASURED_LENGTHS = (240, 320, 400)
REPETITIONS = 3


def main() -> None:
    # The shorter first case warms allocator/code pages and is deliberately
    # excluded from the regression.  Repeated sizes make the result robust to
    # ordinary CI scheduling noise without hiding a systematic exponent rise.
    arguments = ["120"] + [str(n) for _ in range(REPETITIONS) for n in MEASURED_LENGTHS]
    command = [
        "cargo",
        "run",
        "--release",
        "--offline",
        "-q",
        "-p",
        "ribon-core",
        "--example",
        "conditional-density2-bench",
        "--",
        *arguments,
    ]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    rows = json.loads(completed.stdout)
    samples: dict[int, list[dict[str, float | int]]] = defaultdict(list)
    for row in rows[1:]:
        length = int(row["length"])
        if length not in MEASURED_LENGTHS:
            raise AssertionError(f"unexpected benchmark length {length}")
        elapsed = float(row["elapsed_seconds"])
        log_z = float(row["log_partition_function"])
        entries = int(row["pair_probability_entries"])
        peak_heap = int(row["peak_heap_bytes"])
        if not math.isfinite(elapsed) or elapsed <= 0.0:
            raise AssertionError(f"invalid elapsed time at n={length}: {elapsed}")
        if not math.isfinite(log_z):
            raise AssertionError(f"non-finite log partition function at n={length}")
        if entries <= 0:
            raise AssertionError(f"empty pair-probability output at n={length}")
        if peak_heap <= 0:
            raise AssertionError(f"empty heap measurement at n={length}")
        samples[length].append(row)

    medians = {}
    for length in MEASURED_LENGTHS:
        if len(samples[length]) != REPETITIONS:
            raise AssertionError(
                f"n={length} has {len(samples[length])} samples, expected {REPETITIONS}"
            )
        elapsed = statistics.median(
            float(row["elapsed_seconds"]) for row in samples[length]
        )
        medians[length] = {
            "elapsed_seconds": elapsed,
            "seconds_per_n_cubed": elapsed / length**3,
            "peak_heap_bytes": statistics.median(
                int(row["peak_heap_bytes"]) for row in samples[length]
            ),
            "pair_probability_entries": int(samples[length][0]["pair_probability_entries"]),
            "log_partition_function": float(samples[length][0]["log_partition_function"]),
        }

    exponents = []
    for left, right in zip(MEASURED_LENGTHS, MEASURED_LENGTHS[1:]):
        exponent = math.log(
            medians[right]["elapsed_seconds"] / medians[left]["elapsed_seconds"]
        ) / math.log(right / left)
        exponents.append({"from": left, "to": right, "exponent": exponent})
        if not 2.0 <= exponent <= 3.8:
            raise AssertionError(
                f"observed scaling exponent {exponent:.3f} for {left}->{right} "
                "is inconsistent with the cubic implementation contract"
            )

    normalized = [medians[n]["seconds_per_n_cubed"] for n in MEASURED_LENGTHS]
    normalized_ratio = max(normalized) / min(normalized)
    if normalized_ratio > 1.5:
        raise AssertionError(
            f"normalized cubic runtime spread {normalized_ratio:.3f} exceeds 1.5"
        )

    memory_exponents = []
    for left, right in zip(MEASURED_LENGTHS, MEASURED_LENGTHS[1:]):
        exponent = math.log(
            medians[right]["peak_heap_bytes"] / medians[left]["peak_heap_bytes"]
        ) / math.log(right / left)
        memory_exponents.append({"from": left, "to": right, "exponent": exponent})
        if not 1.5 <= exponent <= 2.5:
            raise AssertionError(
                f"observed heap scaling exponent {exponent:.3f} for {left}->{right} "
                "is inconsistent with the quadratic storage contract"
            )

    report = {
        "algorithm": "conditional-density2 polynomial interval hypergraph",
        "build": "cargo release --offline",
        "repetitions": REPETITIONS,
        "medians": {str(length): medians[length] for length in MEASURED_LENGTHS},
        "pairwise_scaling_exponents": exponents,
        "pairwise_heap_scaling_exponents": memory_exponents,
        "normalized_cubic_spread": normalized_ratio,
        "criteria": {
            "pairwise_time_exponent_range": [2.0, 3.8],
            "pairwise_heap_exponent_range": [1.5, 2.5],
            "maximum_normalized_cubic_spread": 1.5,
            "finite_partition_outputs": True,
        },
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
