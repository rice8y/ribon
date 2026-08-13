#!/usr/bin/env python3
"""Build the pinned 24-family real-RNA validation corpus from Rfam SEED data."""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from datetime import date
from pathlib import Path


FAMILIES = [
    "RF00001",  # 5S rRNA
    "RF00002",  # 5.8S rRNA
    "RF00003",  # U2 snRNA
    "RF00004",  # U1 snRNA
    "RF00005",  # tRNA
    "RF00008",  # Hammerhead ribozyme III
    "RF00010",  # bacterial RNase P
    "RF00012",  # U3 snoRNA
    "RF00017",  # bacterial SRP RNA
    "RF00023",  # tmRNA
    "RF00026",  # U6 snRNA
    "RF00029",  # Intron GP I
    "RF00050",  # FMN riboswitch
    "RF00059",  # TPP riboswitch
    "RF00080",  # yybP-ykoY riboswitch
    "RF00094",  # HDV ribozyme
    "RF00162",  # SAM riboswitch
    "RF00167",  # Purine riboswitch
    "RF00168",  # Lysine riboswitch
    "RF00174",  # Cobalamin riboswitch
    "RF00234",  # glmS ribozyme
    "RF00504",  # Glycine riboswitch
    "RF01051",  # cyclic-di-GMP-I riboswitch
    "RF01734",  # Fluoride riboswitch
]

OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
CLOSE_TO_OPEN = {value: key for key, value in OPEN_TO_CLOSE.items()}


def fetch(accession: str) -> str:
    request = urllib.request.Request(
        f"https://rfam.org/family/{accession}/alignment",
        headers={"User-Agent": "ribon-validation/1.0 (+https://github.com/yoneyama/ribon)"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def parse_stockholm(text: str, accession: str) -> dict:
    sequences: dict[str, str] = {}
    structure_parts: list[str] = []
    family_id = accession
    description = accession
    structure_source = ""
    for line in text.splitlines():
        if line.startswith("#=GF ID "):
            family_id = line[8:].strip()
        elif line.startswith("#=GF DE "):
            description = line[8:].strip()
        elif line.startswith("#=GF SS "):
            structure_source = line[8:].strip()
        elif line.startswith("#=GC SS_cons "):
            structure_parts.append(line.split(maxsplit=2)[2])
        elif line and not line.startswith("#") and line != "//":
            fields = line.split()
            if len(fields) == 2:
                sequences[fields[0]] = sequences.get(fields[0], "") + fields[1]
    consensus = "".join(structure_parts)
    if not consensus or not sequences:
        raise ValueError(f"{accession}: missing SEED sequence or SS_cons")

    candidates = []
    for sequence_id, aligned in sequences.items():
        ungapped = re.sub(r"[-.~]", "", aligned).upper().replace("T", "U")
        if 20 <= len(ungapped) <= 500 and re.fullmatch(r"[ACGURYSWKMBDHVN]+", ungapped):
            gap_fraction = 1.0 - len(ungapped) / len(aligned)
            candidates.append((gap_fraction, -len(ungapped), sequence_id, aligned, ungapped))
    if not candidates:
        raise ValueError(f"{accession}: no suitable 20-500 nt IUPAC sequence")
    _, _, sequence_id, aligned, sequence = min(candidates)
    if len(aligned) != len(consensus):
        raise ValueError(f"{accession}: alignment/SS_cons length mismatch")
    structure = project_structure(consensus, aligned)
    if len(structure) != len(sequence):
        raise AssertionError(f"{accession}: projected structure length mismatch")
    return {
        "accession": accession,
        "family_id": family_id,
        "description": description,
        "structure_source": structure_source,
        "sequence_id": sequence_id,
        "sequence": sequence,
        "reference_structure": structure,
        "length": len(sequence),
    }


def project_structure(consensus: str, aligned_sequence: str) -> str:
    partners: dict[int, int] = {}
    stacks: dict[str, list[int]] = {symbol: [] for symbol in OPEN_TO_CLOSE}
    letter_stacks: dict[str, list[int]] = {}
    for index, symbol in enumerate(consensus):
        if symbol in OPEN_TO_CLOSE:
            stacks[symbol].append(index)
        elif symbol in CLOSE_TO_OPEN:
            opening = CLOSE_TO_OPEN[symbol]
            if stacks[opening]:
                left = stacks[opening].pop()
                partners[left] = index
                partners[index] = left
        elif symbol.isupper():
            letter_stacks.setdefault(symbol, []).append(index)
        elif symbol.islower():
            opening = symbol.upper()
            if letter_stacks.get(opening):
                left = letter_stacks[opening].pop()
                partners[left] = index
                partners[index] = left

    kept_columns = [index for index, base in enumerate(aligned_sequence) if base not in "-.~"]
    projected_index = {column: index for index, column in enumerate(kept_columns)}
    result = ["."] * len(kept_columns)
    levels = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
    accepted: list[tuple[int, int, int]] = []
    for left in kept_columns:
        right = partners.get(left)
        if right is None or right not in projected_index or left >= right:
            continue
        i = projected_index[left]
        j = projected_index[right]
        level = 0
        while any(
            other_level == level and ((a < i < b < j) or (i < a < j < b))
            for a, b, other_level in accepted
        ):
            level += 1
        accepted.append((i, j, level))
        if level < len(levels):
            opening, closing = levels[level]
        else:
            opening = chr(ord("A") + level - len(levels))
            closing = opening.lower()
        result[i] = opening
        result[j] = closing
    return "".join(result)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("tests/data/rfam_real_24.json"))
    arguments = parser.parse_args()
    cases = [parse_stockholm(fetch(accession), accession) for accession in FAMILIES]
    payload = {
        "source": "Rfam SEED alignments",
        "source_url": "https://rfam.org/",
        "license": "CC0-1.0",
        "retrieved": date.today().isoformat(),
        "selection": "one minimum-gap 20-500 nt IUPAC sequence from each curated SEED family",
        "cases": cases,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(arguments.output), "cases": len(cases)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
