# Repository validation

This directory contains the commands used by repository and release validation. They operate only on Ribon sources and generated artifacts.

## Required checks

Run the same bounded checks as pull-request CI with:

```sh
just ci-check
```

This gate covers Rust formatting, the bounded unit/integration/reference/protocol test set, lints, byte-for-byte synchronization of the distributed WASM module, and compilation of the Typst smoke and complete public-API fixtures. It deliberately avoids installing PDF rasterizers, the long-sequence log-domain stress test, and performance-sensitive or exhaustive release-scale validation.

## Release validation

- `validate-wasm.py` checks the distributable WASM module, imports, exports, size, and digest.
- `check-wasm-sync.sh` rebuilds the plugin and verifies that the distributed WASM is byte-for-byte current.
- `validate-licenses.py` checks the parameter provenance, license boundary, and exact Typst package file set.
- `validate-contrast.py`, `validate-extended-pdf.py`, `validate-publication-render.py`, `validate-plot-quality.py`, and `validate-plot-layout.py` check Typst rendering contracts.
- `validate-rendered-pdf.py` and `validate-render-golden.py` check the 24-family real-data rendering matrix.
- `validate-pseudoknot.py`, `validate-conditional-density2.py`, `validate-exact-features-real.py`, and the performance scripts check advanced numerical behavior.

The root `justfile` is the supported entry point for these commands. Generated PDFs and rasterized pages are written below `target/qa/`; machine-readable summaries are written to `tests/reports/`.

Run the complete release gate locally with:

```sh
just release-check
```
