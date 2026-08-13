# Development scripts

This directory contains repository-owned validation and maintenance commands that operate only on Ribon sources and generated artifacts. External comparison programs and their adapters live under `tests/oracles/` so they cannot be confused with release inputs.

## Release validation

- `validate-wasm.py` checks the distributable WASM module, imports, exports, size, and digest.
- `check-wasm-sync.sh` rebuilds the plugin and verifies that the distributed WASM is byte-for-byte current.
- `validate-licenses.py` checks the parameter provenance, license boundary, and exact Typst package file set.
- `validate-contrast.py`, `validate-extended-pdf.py`, `validate-publication-render.py`, `validate-plot-quality.py`, and `validate-plot-layout.py` check Typst rendering contracts.
- `validate-rendered-pdf.py` and `validate-render-golden.py` check the 24-family real-data rendering matrix.
- `validate-pseudoknot.py`, `validate-conditional-density2.py`, `validate-exact-features-real.py`, and the performance scripts check advanced numerical behavior.

The root `justfile` is the supported entry point for these commands. Generated PDFs and rasterized pages are written below `target/qa/`; machine-readable summaries are written to `tests/reports/`.

## Maintenance

- `format-markdown.py` enforces one physical source line per prose paragraph and requires alternative text for images.
- `fetch-rfam-corpus.py` refreshes the pinned real-data corpus and is not run by ordinary CI.

Run the complete local gate with:

```sh
just release-check
```
