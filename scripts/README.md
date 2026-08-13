# Repository validation

This directory contains the commands used by repository and release validation. They operate only on Ribon sources and generated artifacts.

## Required checks

Run the same bounded checks as pull-request CI with:

```sh
just ci-check
```

This gate covers Rust formatting, the bounded unit/integration/reference/protocol test set, lints, host-independent source synchronization of the distributed WASM module, and compilation of the Typst smoke and complete public-API fixtures. It validates the committed distribution without modifying it and deliberately avoids installing PDF rasterizers, the long-sequence log-domain stress test, and performance-sensitive or exhaustive release-scale validation.

## Release validation

- `validate-wasm.py` checks the distributable WASM module, imports, exports, size, and digest.
- `validate-wasm-source.py` fingerprints every build input and verifies that the distributed WASM was generated from the current source, manifests, lockfile, toolchain, and parameter data.
- `validate-licenses.py` checks the parameter provenance, license boundary, and exact Typst package file set.
- `validate-contrast.py`, `validate-extended-pdf.py`, `validate-publication-render.py`, `validate-plot-quality.py`, and `validate-plot-layout.py` check Typst rendering contracts.
- `validate-rendered-pdf.py` and `validate-render-golden.py` check the 24-family real-data rendering matrix.
- `validate-pseudoknot.py`, `validate-conditional-density2.py`, `validate-exact-features-real.py`, and the performance scripts check advanced numerical behavior.

The root `justfile` is the supported entry point for these commands. Generated PDFs and rasterized pages are written below `target/qa/`; machine-readable summaries are written to `tests/reports/`.

Run the complete release gate locally with:

```sh
just release-check
```
