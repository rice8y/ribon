#!/usr/bin/env python3
"""Bounded public-WASM performance and output-shape checks for local ensembles."""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    rows = []
    for length in (120, 240, 480):
        started = time.perf_counter()
        process = subprocess.run(
            [
                "typst",
                "query",
                "--root",
                ".",
                "--input",
                f"length={length}",
                "tests/typst/performance_query.typ",
                "<ribon-performance>",
                "--field",
                "value",
                "--one",
            ],
            cwd=root,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30.0,
        )
        elapsed = time.perf_counter() - started
        result = json.loads(process.stdout)
        if result["length"] != length:
            raise AssertionError("performance query returned the wrong sequence length")
        if "exact banded sliding-window Turner inside/outside ensembles" not in result["method"]:
            raise AssertionError("performance query did not use the Rust local-ensemble backend")
        expected_windows = max(1, length - min(50, length) + 1)
        if result["window-count"] != expected_windows:
            raise AssertionError(
                f"unexpected local window count at {length}: {result['window-count']}"
            )
        if result["pair-count"] > length * min(35, length - 1):
            raise AssertionError("local result is not band-limited")
        if elapsed > 30.0:
            raise AssertionError(f"{length}-nt local analysis exceeded 30 s: {elapsed:.3f} s")
        rows.append(result | {"elapsed_seconds": elapsed})

    ratio = rows[-1]["elapsed_seconds"] / max(rows[0]["elapsed_seconds"], 1.0e-9)
    if ratio > 10.0:
        raise AssertionError(f"120→480 nt wall-time ratio exceeded 10: {ratio:.3f}")
    report = {
        "backend": "public Typst/WASM exact banded sliding-window Turner inside/outside ensembles",
        "source": "longest record in tests/data/rfam_real_24.json",
        "window_size": 50,
        "maximum_pair_span": 35,
        "maximum_unpaired": 1,
        "per_case_ceiling_seconds": 30.0,
        "maximum_scaling_ratio_120_to_480": 10.0,
        "observed_scaling_ratio_120_to_480": ratio,
        "results": rows,
    }
    output = root / "tests/reports/performance_validation.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
