# Tests

The test tree separates executable fixtures, immutable inputs, image baselines, external oracles, and generated reports.

- `typst/` contains public-API, error-contract, rendering, and query fixtures.
- `data/` contains the pinned Rfam and pseudoknot corpora with provenance.
- `golden/` contains compact SHA-256 manifests for pixel-exact rendering checks.
- `oracles/` contains optional adapters for independent external implementations; none are release dependencies.
- `reports/` contains machine-readable validation summaries referenced by `docs/VALIDATION.md`.

Rust unit and integration tests remain beside their implementation in `crates/ribon-core/` and `crates/ribon-plugin/`. Temporary PDFs and rasterized pages are generated below `target/qa/` and are never source-controlled.
