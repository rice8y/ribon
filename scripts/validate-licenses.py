#!/usr/bin/env python3
"""Validate the independently implemented release license boundary."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess


MODEL_BUNDLE_SHA256 = "0c00a31400f1dedbe9a3e161b2f9b1b74cde54941144ee988f48173d33bbcd7b"
DNA_MODEL_BUNDLE_SHA256 = "019ad1d5c3dac421df37e0a5aeded6d3da50da03deecc23ba0ae5a6d5d06b977"
REFERENCE_ARCHIVE_SHA256 = "8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8"


def run(arguments: list[str], root: Path) -> str:
    return subprocess.run(
        arguments,
        cwd=root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def parameter_fingerprint(directory: Path, prefix: str) -> tuple[int, str]:
    files = sorted(path for path in directory.glob(f"{prefix}.*") if path.is_file())
    digest = hashlib.sha256()
    for path in files:
        digest.update(path.name.encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return len(files), digest.hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = (root / "package/typst.toml").read_text(encoding="utf-8")
    cargo = (root / "Cargo.toml").read_text(encoding="utf-8")
    if 'license = "GPL-2.0-only"' not in manifest or 'license = "GPL-2.0-only"' not in cargo:
        raise AssertionError("Cargo and Typst manifests must both use GPL-2.0-only")
    for license_file in (root / "LICENSE", root / "package/LICENSE"):
        text = license_file.read_text(encoding="utf-8")
        if "GNU GENERAL PUBLIC LICENSE" not in text or "Version 2, June 1991" not in text:
            raise AssertionError(f"GPL-2.0 text is missing from {license_file}")

    data = root / "crates/ribon-core/data/rnastructure-6.6"
    count, fingerprint = parameter_fingerprint(data, "rna")
    if count != 34 or fingerprint != MODEL_BUNDLE_SHA256:
        raise AssertionError(f"parameter bundle drift: count={count}, sha256={fingerprint}")
    dna_count, dna_fingerprint = parameter_fingerprint(data, "dna")
    if dna_count != 33 or dna_fingerprint != DNA_MODEL_BUNDLE_SHA256:
        raise AssertionError(
            f"DNA parameter bundle drift: count={dna_count}, sha256={dna_fingerprint}"
        )
    data_license = (data / "GPL-2.0.txt").read_text(encoding="utf-8")
    if "Version 2, June 1991" not in data_license:
        raise AssertionError("RNAstructure parameter license text is missing")
    workspace_members = re.search(r"members\s*=\s*\[([^]]+)\]", cargo, re.DOTALL)
    if not workspace_members or "cparty" in workspace_members.group(1).lower():
        raise AssertionError("linked CParty backend remains a Cargo workspace member")
    package_source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((root / "package").rglob("*.typ"))
    )
    plugin_source = (root / "crates/ribon-plugin/src/lib.rs").read_text(encoding="utf-8")
    core_source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (root / "crates/ribon-core/src").glob("*.rs")
    )
    release_source = "\n".join((package_source, plugin_source, core_source)).lower()
    for forbidden in ("ribon_cparty_plugin", "mod vienna_backend"):
        if forbidden in release_source:
            raise AssertionError(f"obsolete linked backend reference remains: {forbidden}")

    wasm = root / "package/ribon_plugin.wasm"
    run(["wasm-tools", "validate", str(wasm)], root)
    wat = run(["wasm-tools", "print", str(wasm)], root)
    imports = re.findall(r'^\s*\(import\s+"([^"]+)"\s+"([^"]+)"', wat, re.MULTILINE)
    exports = re.findall(r'^\s*\(export\s+"([^"]+)"', wat, re.MULTILINE)
    expected_imports = {
        ("typst_env", "wasm_minimal_protocol_write_args_to_buffer"),
        ("typst_env", "wasm_minimal_protocol_send_result_to_host"),
    }
    if set(imports) != expected_imports:
        raise AssertionError(f"unexpected WASM imports: {imports}")
    if "run" not in exports:
        raise AssertionError(f"stable protocol exports missing: {exports}")
    obsolete_exports = {
        "version",
        "predict",
        "predict_cparty",
        "exact_analyze",
        "parameter_profile",
    }
    if obsolete_exports & set(exports):
        raise AssertionError(f"obsolete multi-entry ABI remains: {obsolete_exports & set(exports)}")

    stale_linked_paths = [
        "crates/ribon-cparty-plugin",
        "crates/ribon-core/vendor/vienna",
        "crates/ribon-core/src/vienna_backend.rs",
        "package/ribon_cparty_plugin.wasm",
        "package/VIENNARNA-LICENSE",
    ]
    for relative in stale_linked_paths:
        path = root / relative
        if path.exists():
            raise AssertionError(f"stale linked-backend path remains: {relative}")
        if relative.startswith("package/") and Path(relative).name in manifest:
            raise AssertionError(f"stale linked-backend manifest reference remains: {relative}")
    package_files = {
        path.relative_to(root / "package").as_posix()
        for path in (root / "package").rglob("*")
        if path.is_file()
    }
    expected_package_files = {
        "LICENSE",
        "NOTICE.md",
        "README.md",
        "THIRD_PARTY.md",
        "docs/documentation.pdf",
        "docs/documentation.typ",
        "docs/references.bib",
        "justfile",
        "lib.typ",
        "ribon_plugin.wasm",
        "src/analysis.typ",
        "src/annotations.typ",
        "src/chart.typ",
        "src/constraints.typ",
        "src/plots.typ",
        "src/protocol.typ",
        "src/render.typ",
        "typst.toml",
    }
    for name in (
        "annotations",
        "comparison-dot-plot",
        "continuous-annotation",
        "dot-plot",
        "edited-scene",
        "decoder-comparison",
        "mountain-plot",
        "multi-strand",
        "predicted-structure",
        "pseudoknot",
        "secondary-structure",
        "structure-comparison",
    ):
        expected_package_files.add(f"examples/{name}.typ")
        expected_package_files.add(f"examples/{name}.png")
    if package_files != expected_package_files:
        raise AssertionError(
            "runtime package file set drift: "
            f"unexpected={sorted(package_files - expected_package_files)}, "
            f"missing={sorted(expected_package_files - package_files)}"
        )
    for relative in ("artifacts/ViennaRNA", "artifacts/VARNA"):
        path = root / relative
        if path.exists():
            ignored = subprocess.run(
                ["git", "check-ignore", "-q", relative], cwd=root, check=False
            ).returncode == 0
            if not ignored:
                raise AssertionError(f"external oracle is not release-excluded: {relative}")

    notice = (root / "package/NOTICE.md").read_text(encoding="utf-8")
    third_party = (root / "THIRD_PARTY.md").read_text(encoding="utf-8")
    package_third_party = (root / "package/THIRD_PARTY.md").read_text(encoding="utf-8")
    for value in (
        MODEL_BUNDLE_SHA256,
        DNA_MODEL_BUNDLE_SHA256,
        REFERENCE_ARCHIVE_SHA256,
        "GPL-2.0-only",
    ):
        if value not in notice + third_party + package_third_party + manifest:
            raise AssertionError(f"release provenance is missing {value}")
    for required in ("RNAstructure 6.6", "GPL-2.0-only", "wasm-minimal-protocol"):
        if required not in package_third_party:
            raise AssertionError(f"package THIRD_PARTY.md is missing {required}")

    report = {
        "project_license": "GPL-2.0-only",
        "parameter_source": "RNAstructure 6.6 data_tables/rna.*",
        "parameter_file_count": count,
        "parameter_bundle_sha256": fingerprint,
        "dna_parameter_source": "RNAstructure 6.6 data_tables/dna.*",
        "dna_parameter_file_count": dna_count,
        "dna_parameter_bundle_sha256": dna_fingerprint,
        "reference_archive_sha256": REFERENCE_ARCHIVE_SHA256,
        "linked_viennarna": False,
        "linked_cparty": False,
        "linked_varna": False,
        "linked_original_naview": False,
        "wasm_imports": [f"{module}.{name}" for module, name in imports],
        "wasm_public_exports": [name for name in exports if name == "run"],
        "stale_linked_backend_paths_absent": True,
        "external_oracles_release_excluded": True,
        "runtime_package_file_set_exact": True,
        "typst_universe_license_blocker": False,
    }
    output = root / "tests/reports/license_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
