# External oracles

These adapters compare Ribon with independent command-line tools or public library APIs. They are optional development tests, are not compiled into Ribon, and are excluded from the Typst package assembled from `package/`.

- `rnastructure.py` performs the pinned RNAstructure 6.6 numerical differential.
- `cparty.py` performs a black-box comparison through the documented CParty CLI.
- `vienna/` contains small public-API harnesses and drivers for layout and thermodynamic differential testing.

External source trees and build products belong under the ignored `artifacts/` directory or an explicitly supplied path. Exact provenance, hashes, scope, and tolerances are documented in `docs/VALIDATION.md` and `THIRD_PARTY.md`.
