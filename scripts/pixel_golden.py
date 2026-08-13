"""Shared comparison logic for deterministic pixel-golden manifests."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


def _changed_pages(actual: Sequence[str], expected: Sequence[str]) -> list[int]:
    changed = [
        index
        for index, (actual_hash, expected_hash) in enumerate(
            zip(actual, expected, strict=False), 1
        )
        if actual_hash != expected_hash
    ]
    if len(actual) != len(expected):
        changed.extend(
            range(
                min(len(actual), len(expected)) + 1,
                max(len(actual), len(expected)) + 1,
            )
        )
    return changed


def validate_pixel_golden(
    actual: Mapping[str, Any], expected: Mapping[str, Any], *, label: str
) -> None:
    """Compare stable manifest metadata and page hashes with useful failures."""
    actual_metadata = {key: value for key, value in actual.items() if key != "sha256"}
    expected_metadata = {key: value for key, value in expected.items() if key != "sha256"}
    if actual_metadata != expected_metadata:
        keys = sorted(
            key
            for key in actual_metadata.keys() | expected_metadata.keys()
            if actual_metadata.get(key) != expected_metadata.get(key)
        )
        raise AssertionError(f"{label} pixel-golden metadata drift: keys={keys}")

    actual_hashes = actual.get("sha256")
    expected_hashes = expected.get("sha256")
    if (
        not isinstance(actual_hashes, list)
        or not isinstance(expected_hashes, list)
        or not all(isinstance(value, str) for value in actual_hashes + expected_hashes)
    ):
        raise AssertionError(f"{label} pixel golden must contain a sha256 list")
    changed = _changed_pages(actual_hashes, expected_hashes)
    if changed:
        raise AssertionError(f"{label} pixel-golden image drift: pages={changed}")
