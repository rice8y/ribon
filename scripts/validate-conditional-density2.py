#!/usr/bin/env python3
"""Numeric and pixel-level validation of the conditional density-2 engine."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path

from imagemagick import convert_command


OPEN = {"(": ")", "[": "]", "{": "}", "<": ">"}
ENERGY_TERMS = (
    "seed_layer_kcal_mol",
    "added_layer_kcal_mol",
    "spanning_stack_correction_kcal_mol",
    "spanning_internal_correction_kcal_mol",
    "spanning_multiloop_correction_kcal_mol",
    "pseudoloop_initiation_kcal_mol",
    "multiloop_pseudoknot_kcal_mol",
    "nested_pseudoknot_kcal_mol",
    "band_kcal_mol",
    "pseudoloop_unpaired_kcal_mol",
    "closed_subregion_kcal_mol",
    "constraint_kcal_mol",
    "decomposition_alignment_kcal_mol",
)


def run(command: list[str], root: Path) -> str:
    return subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def pairs(structure: str) -> set[tuple[int, int]]:
    stacks = {symbol: [] for symbol in OPEN}
    closing = {value: key for key, value in OPEN.items()}
    result: set[tuple[int, int]] = set()
    for position, symbol in enumerate(structure, 1):
        if symbol in stacks:
            stacks[symbol].append(position)
        elif symbol in closing:
            opening = closing[symbol]
            if not stacks[opening]:
                raise AssertionError(f"unbalanced structure at {position}")
            result.add((stacks[opening].pop(), position))
    if any(stacks.values()):
        raise AssertionError("unbalanced structure")
    return result


def pair_map(entries: list[dict[str, object]]) -> dict[tuple[int, int], float]:
    return {
        (int(entry["i"]), int(entry["j"])): float(entry["probability"])
        for entry in entries
    }


def validate_full_production_ensemble(
    record: dict[str, object],
    ensemble: dict[str, object],
    mfe_evaluation: dict[str, object],
    supplied_evaluation: dict[str, object],
    label: str,
) -> None:
    case_id = str(record["id"])
    if not ensemble["state_space_complete"]:
        raise AssertionError(f"{case_id}: {label} full production ensemble incomplete")
    if ensemble["time_complexity"] != "O(n^3)" or ensemble["space_complexity"] != "O(n^2)":
        raise AssertionError(f"{case_id}: {label} polynomial complexity contract changed")
    if abs(
        float(ensemble["mfe_energy_kcal_mol"])
        - float(mfe_evaluation["energy_kcal_mol"])
    ) > 1.0e-9:
        raise AssertionError(f"{case_id}: {label} MFE energy evaluator round-trip mismatch")
    if not mfe_evaluation["derivation_unique"]:
        raise AssertionError(f"{case_id}: {label} MFE derivation is ambiguous")
    if float(ensemble["ensemble_free_energy_kcal_mol"]) > float(
        ensemble["mfe_energy_kcal_mol"]
    ) + 1.0e-9:
        raise AssertionError(f"{case_id}: {label} ensemble free energy exceeds MFE")
    if float(ensemble["mfe_energy_kcal_mol"]) > float(
        supplied_evaluation["energy_kcal_mol"]
    ) + 1.0e-9:
        raise AssertionError(f"{case_id}: {label} MFE is above supplied-state energy")
    seed_pairs = pairs(str(ensemble["seed_structure"]))
    for decoder in ("mfe_structure", "centroid_structure", "mea_structure"):
        decoded = str(ensemble[decoder])
        if len(decoded) != len(str(record["sequence"])):
            raise AssertionError(f"{case_id}: {label} {decoder} length changed")
        if not seed_pairs <= pairs(decoded):
            raise AssertionError(f"{case_id}: {label} {decoder} lost a fixed seed pair")
    ensemble_pairs = pair_map(ensemble["pair_probabilities"])
    for pair in seed_pairs:
        if abs(ensemble_pairs.get(pair, 0.0) - 1.0) > 1.0e-12:
            raise AssertionError(f"{case_id}: {label} fixed seed probability changed at {pair}")
    full_mass = [0.0] * len(str(record["sequence"]))
    for (i, j), probability in ensemble_pairs.items():
        if not 0.0 <= probability <= 1.0:
            raise AssertionError(f"{case_id}: {label} invalid pair probability at {(i, j)}")
        full_mass[i - 1] += probability
        full_mass[j - 1] += probability
    for position, (paired, unpaired) in enumerate(
        zip(full_mass, ensemble["unpaired_probabilities"], strict=True), 1
    ):
        if abs(paired + float(unpaired) - 1.0) > 1.0e-9:
            raise AssertionError(
                f"{case_id}: {label} full-ensemble probability mass drift at {position}"
            )


def validate_planar_reduction(
    record: dict[str, object],
    conditional: dict[str, object],
    oracle: dict[str, object] | None,
    planar: dict[str, object],
    label: str,
) -> None:
    case_id = str(record["id"])
    if not conditional["state_space_complete"]:
        raise AssertionError(f"{case_id}: {label} conditional prefix state space incomplete")
    if abs(
        float(conditional["log_partition_function"])
        - float(planar["log_partition_function"])
    ) > 1.0e-10:
        raise AssertionError(f"{case_id}: {label} empty-seed log PF mismatch")
    if oracle is not None and abs(
        float(conditional["log_partition_function"])
        - float(oracle["log_partition_function"])
    ) > 1.0e-10:
        raise AssertionError(f"{case_id}: {label} polynomial/oracle log PF mismatch")
    conditional_pairs = pair_map(conditional["pair_probabilities"])
    planar_pairs = pair_map(planar["pair_probabilities"])
    oracle_pairs = pair_map(oracle["pair_probabilities"]) if oracle is not None else {}
    for pair in conditional_pairs.keys() | planar_pairs.keys():
        if abs(conditional_pairs.get(pair, 0.0) - planar_pairs.get(pair, 0.0)) > 1.0e-10:
            raise AssertionError(f"{case_id}: {label} planar marginal mismatch at {pair}")
    if oracle is not None:
        for pair in conditional_pairs.keys() | oracle_pairs.keys():
            if abs(conditional_pairs.get(pair, 0.0) - oracle_pairs.get(pair, 0.0)) > 1.0e-10:
                raise AssertionError(f"{case_id}: {label} oracle marginal mismatch at {pair}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def version(command: list[str], root: Path) -> str:
    output = subprocess.run(
        command,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout
    return output.splitlines()[0].strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--update-golden", action="store_true")
    arguments = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    query = run(
        [
            "typst",
            "query",
            "--root",
            ".",
            "tests/typst/conditional_density2_query.typ",
            "<ribon-conditional-density2-validation>",
            "--field",
            "value",
            "--one",
        ],
        root,
    )
    records = json.loads(query)
    if len(records) < 24:
        raise AssertionError("conditional density-2 corpus must contain at least 24 cases")

    numeric_rows = []
    for record in records:
        evaluation = record["evaluation"]
        oracle_evaluation = record["oracle-evaluation"]
        energy = oracle_evaluation["energy"]
        if pairs(evaluation["structure"]) != pairs(record["expected-structure"]):
            raise AssertionError(f"{record['id']}: production evaluated pair table changed")
        if pairs(oracle_evaluation["structure"]) != pairs(record["expected-structure"]):
            raise AssertionError(f"{record['id']}: oracle evaluated pair table changed")
        if not math.isfinite(float(evaluation["energy_kcal_mol"])):
            raise AssertionError(f"{record['id']}: non-finite production energy")
        if not evaluation["derivation_unique"]:
            raise AssertionError(f"{record['id']}: production derivation is ambiguous")
        values = [float(energy[name]) for name in ENERGY_TERMS]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError(f"{record['id']}: non-finite energy component")
        if abs(sum(values) - float(energy["total_kcal_mol"])) > 1.0e-9:
            raise AssertionError(f"{record['id']}: energy components do not sum to total")
        if int(energy["maximum_band_density"]) > 2:
            raise AssertionError(f"{record['id']}: density-2 reference was rejected")
        if int(energy["pseudoknot_band_count"]) < 2:
            raise AssertionError(f"{record['id']}: pseudoknot bands were not recovered")

        ensemble = record["ensemble"]
        mfe_evaluation = record["mfe-evaluation"]
        if not ensemble["state_space_complete"]:
            raise AssertionError(f"{record['id']}: full production ensemble incomplete")
        if ensemble["time_complexity"] != "O(n^3)" or ensemble["space_complexity"] != "O(n^2)":
            raise AssertionError(f"{record['id']}: polynomial complexity contract changed")
        if abs(
            float(ensemble["mfe_energy_kcal_mol"])
            - float(mfe_evaluation["energy_kcal_mol"])
        ) > 1.0e-9:
            raise AssertionError(f"{record['id']}: MFE energy evaluator round-trip mismatch")
        if not mfe_evaluation["derivation_unique"]:
            raise AssertionError(f"{record['id']}: MFE derivation is ambiguous")
        if float(ensemble["ensemble_free_energy_kcal_mol"]) > float(
            ensemble["mfe_energy_kcal_mol"]
        ) + 1.0e-9:
            raise AssertionError(f"{record['id']}: ensemble free energy exceeds MFE")
        if float(ensemble["mfe_energy_kcal_mol"]) > float(
            evaluation["energy_kcal_mol"]
        ) + 1.0e-9:
            raise AssertionError(f"{record['id']}: MFE is above the supplied reference state")
        seed_pairs = pairs(ensemble["seed_structure"])
        for decoder in ("mfe_structure", "centroid_structure", "mea_structure"):
            decoded = ensemble[decoder]
            if len(decoded) != len(record["sequence"]):
                raise AssertionError(f"{record['id']}: {decoder} length changed")
            if not seed_pairs <= pairs(decoded):
                raise AssertionError(f"{record['id']}: {decoder} lost a fixed seed pair")
        ensemble_pairs = pair_map(ensemble["pair_probabilities"])
        for pair in seed_pairs:
            if abs(ensemble_pairs.get(pair, 0.0) - 1.0) > 1.0e-12:
                raise AssertionError(f"{record['id']}: fixed seed probability changed at {pair}")
        full_mass = [0.0] * len(record["sequence"])
        for (i, j), probability in ensemble_pairs.items():
            if not 0.0 <= probability <= 1.0:
                raise AssertionError(f"{record['id']}: invalid pair probability at {(i, j)}")
            full_mass[i - 1] += probability
            full_mass[j - 1] += probability
        for position, (paired, unpaired) in enumerate(
            zip(full_mass, ensemble["unpaired_probabilities"], strict=True), 1
        ):
            if abs(paired + float(unpaired) - 1.0) > 1.0e-9:
                raise AssertionError(
                    f"{record['id']}: full-ensemble probability mass drift at {position}"
                )

        evaluation_d2 = record["evaluation-d2"]
        if pairs(evaluation_d2["structure"]) != pairs(record["expected-structure"]):
            raise AssertionError(f"{record['id']}: dangles=2 evaluated pair table changed")
        if not math.isfinite(float(evaluation_d2["energy_kcal_mol"])):
            raise AssertionError(f"{record['id']}: non-finite dangles=2 production energy")
        if not evaluation_d2["derivation_unique"]:
            raise AssertionError(f"{record['id']}: dangles=2 production derivation is ambiguous")
        validate_full_production_ensemble(
            record,
            record["ensemble-d2"],
            record["mfe-evaluation-d2"],
            evaluation_d2,
            "dangles=2",
        )

        conditional = record["conditional"]
        conditional_oracle = record["conditional-oracle"]
        planar = record["planar"]
        if not conditional["state_space_complete"]:
            raise AssertionError(f"{record['id']}: conditional prefix state space incomplete")
        if abs(
            float(conditional["log_partition_function"])
            - float(planar["log_partition_function"])
        ) > 1.0e-10:
            raise AssertionError(f"{record['id']}: empty-seed log PF mismatch")
        if abs(
            float(conditional["log_partition_function"])
            - float(conditional_oracle["log_partition_function"])
        ) > 1.0e-10:
            raise AssertionError(f"{record['id']}: polynomial/oracle log PF mismatch")
        conditional_pairs = pair_map(conditional["pair_probabilities"])
        planar_pairs = pair_map(planar["pair_probabilities"])
        oracle_pairs = pair_map(conditional_oracle["pair_probabilities"])
        for pair in conditional_pairs.keys() | planar_pairs.keys():
            if abs(conditional_pairs.get(pair, 0.0) - planar_pairs.get(pair, 0.0)) > 1.0e-10:
                raise AssertionError(f"{record['id']}: pair marginal mismatch at {pair}")
        for pair in conditional_pairs.keys() | oracle_pairs.keys():
            if abs(conditional_pairs.get(pair, 0.0) - oracle_pairs.get(pair, 0.0)) > 1.0e-10:
                raise AssertionError(
                    f"{record['id']}: polynomial/oracle pair marginal mismatch at {pair}"
                )
        mass = [0.0] * len(record["prefix"])
        for (i, j), probability in conditional_pairs.items():
            mass[i - 1] += probability
            mass[j - 1] += probability
        for position, (paired, unpaired) in enumerate(
            zip(mass, conditional["unpaired_probabilities"], strict=True), 1
        ):
            if abs(paired + float(unpaired) - 1.0) > 1.0e-10:
                raise AssertionError(f"{record['id']}: probability mass drift at {position}")
        validate_planar_reduction(
            record,
            record["conditional-d2"],
            None,
            record["planar-d2"],
            "dangles=2",
        )
        numeric_rows.append(
            {
                "id": record["id"],
                "length": len(record["sequence"]),
                "density": energy["maximum_band_density"],
                "pseudoknot_bands": energy["pseudoknot_band_count"],
                "conditional_energy_kcal_mol": evaluation["energy_kcal_mol"],
                "oracle_diagnostic_energy_kcal_mol": energy["total_kcal_mol"],
                "mfe_energy_kcal_mol": ensemble["mfe_energy_kcal_mol"],
                "ensemble_free_energy_kcal_mol": ensemble[
                    "ensemble_free_energy_kcal_mol"
                ],
                "full_pair_probability_entries": len(ensemble["pair_probabilities"]),
                "prefix_state_count": conditional_oracle["state_count"],
            }
        )

    output_pdf = root / "target/qa/ribon-conditional-density2-validation.pdf"
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "typst",
            "compile",
            "--root",
            ".",
            "tests/typst/conditional_density2_render.typ",
            str(output_pdf),
        ],
        root,
    )
    run(["qpdf", "--check", str(output_pdf)], root)
    information = run(["pdfinfo", str(output_pdf)], root)
    pages_match = re.search(r"^Pages:\s+(\d+)$", information, re.MULTILINE)
    if not pages_match or int(pages_match.group(1)) != len(records):
        raise AssertionError("conditional PDF page count mismatch")
    images = run(["pdfimages", "-list", str(output_pdf)], root)
    if any(re.match(r"^\s*\d+\s+\d+\s+", line) for line in images.splitlines()):
        raise AssertionError("conditional PDF contains an embedded raster image")
    extracted = run(["pdftotext", str(output_pdf), "-"], root)
    missing = [record["id"] for record in records if record["id"] not in extracted]
    if missing:
        raise AssertionError(f"conditional PDF labels missing: {missing}")

    page_metrics = []
    with tempfile.TemporaryDirectory(prefix="ribon-density2-") as directory:
        prefix = Path(directory) / "page"
        run(["pdftoppm", "-png", "-r", "110", str(output_pdf), str(prefix)], root)
        rendered = sorted(Path(directory).glob("page-*.png"))
        if len(rendered) != len(records):
            raise AssertionError("conditional PDF raster page count mismatch")
        hashes = [digest(page) for page in rendered]
        for page_number, image in enumerate(rendered, 1):
            metrics = run(
                convert_command(
                    str(image),
                    "-fuzz",
                    "2%",
                    "-trim",
                    "-format",
                    "%w %h %[fx:page.x] %[fx:page.y] %[fx:page.width] "
                    "%[fx:page.height] %[fx:mean] %[fx:standard_deviation]",
                    "info:",
                ),
                root,
            )
            width, height, x, y, canvas_width, canvas_height, mean, deviation = map(
                float, metrics.split()
            )
            margins = [x, y, canvas_width - x - width, canvas_height - y - height]
            if min(margins) < 8:
                raise AssertionError(f"page {page_number} approaches crop edge: {margins}")
            if deviation < 0.015 or mean > 0.999:
                raise AssertionError(f"page {page_number} appears blank")
            page_metrics.append(
                {
                    "page": page_number,
                    "content_margins_pixels": margins,
                    "ink_fraction": 1.0 - mean,
                    "grayscale_standard_deviation": deviation,
                }
            )

    manifest = {
        "schema": 1,
        "source": "tests/typst/conditional_density2_render.typ",
        "dpi": 110,
        "page_count": len(hashes),
        "sha256": hashes,
        "rasterizer": version(["pdftoppm", "-v"], root),
        "typst": version(["typst", "--version"], root),
    }
    golden_path = root / "tests/golden/conditional_density2_golden_sha256.json"
    if arguments.update_golden:
        golden_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        golden = json.loads(golden_path.read_text())
        for key in ("schema", "dpi", "page_count", "rasterizer", "typst", "sha256"):
            if manifest[key] != golden.get(key):
                raise AssertionError(f"conditional visual golden differs for {key}")

    report = {
        "source": "Andronescu-Pop-Condon S-Test ShPK",
        "cases": len(records),
        "all_production_and_oracle_reference_pair_tables_exact": True,
        "all_production_derivations_unique": True,
        "all_full_ensembles_complete": True,
        "all_dangles_0_and_2_full_ensembles_complete": True,
        "all_mfe_energy_round_trips_exact": True,
        "all_decoders_preserve_fixed_seed": True,
        "all_energy_components_finite_and_additive": True,
        "maximum_observed_density": max(row["density"] for row in numeric_rows),
        "empty_seed_prefix_pf_and_marginals_match_planar_dp_and_exhaustive_oracle": True,
        "probability_mass_normalized": True,
        "pdf": output_pdf.relative_to(root).as_posix(),
        "pages": len(page_metrics),
        "native_vector_only": True,
        "pixel_exact": True,
        "minimum_content_margin_pixels": min(
            min(page["content_margins_pixels"]) for page in page_metrics
        ),
        "minimum_ink_fraction": min(page["ink_fraction"] for page in page_metrics),
        "results": numeric_rows,
        "pages_detail": page_metrics,
    }
    report_path = root / "tests/reports/conditional_density2_validation.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key not in {"results", "pages_detail"}}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
