# Third-party notices

This document is a technical notice and does not replace the license text of any component.

## RNAstructure 6.6 thermodynamic tables

Ribon includes 34 `data_tables/rna.*` files and 33 `data_tables/dna.*` files from the official RNAstructure 6.6 distribution in its source tree. `ribon-core/build.rs` parses and validates the free-energy and enthalpy data independently, then converts them into separate Rust constants for RNA and DNA. No RNAstructure C++ program code is compiled, linked, copied, or distributed with Ribon.

- Upstream: <https://rna2.urmc.rochester.edu/RNAstructureDownload.html>
- Version: 6.6
- Official macOS arm64 Conda archive SHA-256: `8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8`
- Normalized 34-file RNA bundle SHA-256: `0c00a31400f1dedbe9a3e161b2f9b1b74cde54941144ee988f48173d33bbcd7b`
- Normalized 33-file DNA bundle SHA-256: `019ad1d5c3dac421df37e0a5aeded6d3da50da03deecc23ba0ae5a6d5d06b977`
- License: `GPL-2.0-only`
- License text: `crates/ribon-core/data/rnastructure-6.6/GPL-2.0.txt`

Because these tables are included, the Ribon Rust and Typst sources and the distributed WASM module are provided under `GPL-2.0-only`.

## Rust dependencies

`Cargo.lock` is the authoritative record of exact dependency versions.

| Component | License |
|---|---|
| serde / serde_json / proc-macro2 / quote / syn / itoa | MIT OR Apache-2.0 |
| unicode-ident | (MIT OR Apache-2.0) AND Unicode-3.0 |
| memchr | Unlicense OR MIT |
| wasm-minimal-protocol | Unlicense |

## Modified-nucleotide thermodynamic data

`modified_parameters.rs` contains independently transcribed numerical facts from published nearest-neighbor tables for m6A, pseudouridine, inosine-C/U, 7-deazaadenosine, and purine/nebularine. Dihydrouridine uses the separately identified published model correction rather than presenting it as an experimental nearest-neighbor table. Every modified-analysis result reports DOI provenance and whether enthalpy-based temperature scaling is available. The sources are indexed in `docs/REFERENCES.md`; no article text, figures, or upstream program code are included.

## Validation data

- `tests/data/rfam_real_24.json` is derived from Rfam SEED alignments distributed under CC0-1.0. It is used only for validation and is not part of the Typst package.
- `tests/data/pseudoknot_real_24.json` contains 24 published records from the Andronescu-Pop-Condon S-Test ShPK set. It is used only for validation and is not part of the Typst package. RNA STRAND, PseudoBase, and literature provenance for each record are documented in `tests/data/README.md`.

## External development oracles

The following programs and publications are used only for algorithm research or differential validation. They are not linked into or embedded in the distributed WASM module.

- The official RNAstructure 6.6 CLI provides a 24-case oracle for MFE structures and energies, ensemble energies, pair probabilities, centroid structures, and MEA structures.
- ViennaRNA and RNAplot 2.7.2 provide historical dangle and layout fixtures.
- VARNA informs palette and annotation-expression comparisons.
- The NAView paper and existing output provide geometric regression references for the independent Rust layout.
- CParty, ProbKnot, and the Dirks-Pierce publications provide pseudoknot-model comparisons. The documented CParty 1.0 CLI is also used for a 64-case black-box numerical comparison.

The Ps, Psm, Psp, Pb, Pup, and Pps states, ordinary and spanning loop coefficients in `conditional_density2.rs`, and the 15 states and border recurrences in `conditional_density2/polynomial.rs` are independently implemented from the CParty supplementary information, the published Jabbari-Condon-Zhao recurrences, and the cited DP09 tables. CParty, HFold, and HotKnots source code, headers, objects, and generated tables are not copied, compiled, linked, or distributed. A CParty executable may be invoked only as an optional black-box development oracle and is not a release-test or runtime dependency.

External validation adapters and source checkouts are local development tools and are not distributed in this repository, release archives, or the Typst package. The ignored `artifacts/` directory is reserved for local external checkouts.

## Runtime boundary

The release WASM module imports only the two `typst_env` functions required by the minimal protocol. Its inspected exports are `memory`, `run`, and toolchain data/heap globals. The distributed package contains no third-party native archive, C/C++ runtime, filesystem bridge, process bridge, or network bridge.
