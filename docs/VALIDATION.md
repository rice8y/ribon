# Validation

Ribon validates numerical outputs, algorithmic invariants, the public WASM and Typst APIs, and rendered images as separate contracts. Machine-readable results are stored in `tests/reports/*.json`.

## RNAstructure 6.6 differential validation

The differential validation recorded in `tests/reports/rnastructure_6_6_validation.json` uses the official macOS ARM64 conda package as an external oracle; its archive SHA-256 is `8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8`. The corpus contains the first 80 nucleotides from 24 distinct Rfam families. Both programs run at 37 °C with isolated pairs permitted.

The comparison covers MFE structure and energy, ensemble free energy, every base-pair probability, the 0.5-threshold centroid, and MEA with gamma set to 1. Ribon uses dangles=3, including coaxial stacking, for MFE and dangles=2 for the integrated ensemble. The programs use the same parameter family but not an identical grammar, so the report records both structural statistics and numerical errors instead of requiring bitwise equality.

| Metric | 24-case result | Release gate |
|---|---:|---:|
| Mean absolute MFE-energy error | `1.00000 kcal/mol` | At most `1.1` |
| Maximum absolute MFE-energy error | `2.7000 kcal/mol` | At most `3.0` |
| Mean absolute ensemble-energy error | `1.48940 kcal/mol` | At most `1.6` |
| Mean MFE pair sensitivity / precision | `0.79816 / 0.82318` | At least `0.79 / 0.81` |
| Mean centroid sensitivity / precision | `0.89657 / 0.91438` | At least `0.88 / 0.90` |
| Mean MEA sensitivity / precision | `0.88267 / 0.89590` | At least `0.87 / 0.88` |
| All-pair probability MAE | `0.0041549` | At most `0.01` |
| Maximum single-pair probability error | `0.58989` | Report only |

To keep large pairwise discrepancies visible, `rnastructure_6_6_validation.json` stores the 20 largest pair-probability errors for every case. A Rust integration test recomputes the reference structures, energies, and all probabilities at or above `0.01` from that report.

## Exact internal validation

Short sequences are compared with complete structure enumeration to validate the MFE dynamic program, partition inside/outside values, sampling frequencies, and k-best energy order. A 44 nt structure with a forced 36 nt internal loop verifies that the default dynamic program does not impose the historical 30 nt limit and that the loop is excluded only when the caller explicitly sets `Some(30)`. Cofold tests cover reduction to a bare duplex, reduction to the product of independent monomer partition functions, and MFE re-evaluation across a strand break. Fixed-structure tests require the total energy to equal the sum of its loop breakdown and verify constraints and coaxial selection. Probability tests enforce per-position mass balance at approximately `1e-10`, and long GC-rich inputs verify finite log-domain results.

The 24-family integration corpus exercises the following paths for every Rfam case:

- Global analysis, fold and evaluation with dangle models 0–3, partition functions, and ensemble defect
- Sampling, suboptimal folding, and accessibility
- Duplex, cofold, local, circular, modified-base, and G-quadruplex analysis
- Multicomponent pseudoknot ensembles, comparative covariation, and parameter profiles
- Simple, NAView, circular, Turtle, Puzzler, and linear layouts

## Pseudoknot benchmark

The 24 published ShPK records in `tests/data/pseudoknot_real_24.json` validate extended-dot-bracket round trips, crossing-bracket coloring, pair/unpaired normalization, and finite multicomponent H-type partition functions. On small systems, weighted-interval dynamic programming agrees with independent complete-state enumeration within `1e-12`. Exact decimal state counts remain correct even at:

```math
2^{64}.
```

The opt-in arbitrary-topology matching ensemble is also completed for the first 10 nt of all 24 real sequences, covering 4–116 states per input and checking the complete-state-space flag, finite log partition functions, and probability mass balance at every position. Mean sensitivity and precision against experimental structures are pinned in `tests/reports/pseudoknot_validation.json`. Because the fast ensemble is restricted to span-disjoint H-type component states, Ribon does not claim energy identity with an unrestricted pseudoknot model.

The current baseline has mean sensitivity/precision `0.47405 / 0.57568` for pure ProbKnot under its published conditions and `0.76569 / 0.74715` for the thermodynamic-core hybrid; the hybrid recovers a crossing pair in 24/24 cases. Arbitrary-topology matching centroid and MEA achieve `0.37916 / 0.53414` and `0.40157 / 0.56257`, respectively. Restricted MFE, centroid, and MEA achieve `0.70568 / 0.71803`, `0.68512 / 0.72468`, and `0.70123 / 0.72063`. Restricted pair/unpaired mass is normalized and MFE and partition values are finite in every case. This small corpus contains only pseudoknotted ShPK structures and is not the same population as the full multi-family benchmark in the ProbKnot publication.

## Conditional density-2 benchmark

`scripts/validate-conditional-density2.py` decomposes the same 24 real ShPK structures into two planar layers. Across lengths 28–91 nt, the production evaluator and an independent complete enumerator reconstruct the supplied pair table in 24/24 cases, and the production derivation is unique. Every fixed-seed production ensemble completes at full length. Re-evaluated MFE energies agree within `1e-9 kcal/mol`; MFE, centroid, and MEA preserve every seed pair; and per-position pair/unpaired mass agrees within `1e-9`. The energy ordering is:

```math
G_{\mathrm{ensemble}} \leq G_{\mathrm{MFE}} \leq G_{\mathrm{supplied\ state}}.
```

Every independently enumerated component energy is finite and additive, and maximum band density is 2.

For an independent short-input reduction, the first 10 nt of each sequence are evaluated with an empty seed by both the production polynomial engine and the ordinary planar inside/outside dynamic program under dangle models 0 and 2. Log partition functions and every pair marginal agree within `1e-10`; model 0 adds a third comparison with the exponential density-2 oracle. Rust unit tests compare structure counts with complete enumeration for 14 varied fixed seeds and every small planar seed generated at 12 nt. For a separate nonempty seed, every allowed G-prime layer under models 0 and 2 is independently enumerated and weighted by the production evaluator; log partition function, all pair marginals, MFE, centroid objective, and MEA objective agree within `1e-10`. Outside pair marginals also agree with central finite differences of the log partition function with respect to pair energy within `2e-8`. The polynomial function rejects models 1 and 3, while Rust and WASM tests confirm that the public API dispatches those requests to the complete exhaustive engine.

`scripts/validate-conditional-density2-performance.py` measures release builds with a tracking allocator three times each at sequence lengths `240`, `320`, and `400`. The CI gate requires an adjacent-size time exponent of `2.0–3.8`, a heap exponent of `1.5–2.5`, a maximum-to-minimum cubic-normalized time ratio no greater than `1.5`, a finite partition function, and nonempty pair output. Median measurements, fitted exponents, peak heap use, and normalized spread are written to `tests/reports/conditional_density2_performance.json` on every run.

The same 24 cases are rendered one per page as production MFE, centroid, and MEA structures in `target/qa/ribon-conditional-density2-validation.pdf`. Validation requires a valid qpdf document, 24 pages, every case label, no embedded raster image, and no blank or clipped PNG rendering. The minimum content margin is 38 px and the minimum ink fraction is `0.040481`. Per-page SHA-256 hashes in `tests/golden/conditional_density2_golden_sha256.json` reject any one-pixel drift at the release gate.

## CParty CLI black-box differential validation

The black-box differential validation recorded in `tests/reports/cparty_blackbox_validation.json` invokes only the documented CParty 1.0 command-line interface as an external process. It evaluates the 24 published ShPK RNAs with fixed seeds under dangle models 0 and 2 and DNA homologs of the first eight sequences with an empty seed under both models, for 64 conditions in total. The report pins the external binary, RNA parameter file, and version by SHA-256. None of these external artifacts is part of the distributed WASM module.

For RNA, MFE pair sets agree in 24/24 cases under both dangle conditions. With the standard dangles=2 setting, MFE-energy MAE is `0.04329 kcal/mol`, with a maximum of `0.395`. The MAE against the CLI's displayed ensemble energy is `0.20403`. Ensemble energy reconstructed from the same CLI's MFE frequency uses:

```math
G_{\mathrm{ensemble}} = E_{\mathrm{MFE}} + RT\ln p_{\mathrm{MFE}}.
```

The reconstructed value reduces the MAE to `0.04931 kcal/mol`, with a maximum of `0.35697` and Pearson correlation `0.999970`. The CLI's own displayed values have an identity residual with MAE `0.20810` and maximum `0.87997 kcal/mol` for that group, so Ribon does not imitate the display offset. Separate gates record pair F1, exact-pair-set count, and the MAE, maximum error, and correlation for MFE energy, displayed ensemble energy, reconstructed ensemble energy, and MFE frequency in every RNA and DNA group. All rows are stored in `tests/reports/cparty_blackbox_validation.json`.

## Public API and WASM

`tests/typst/api_all.typ` runs 19 standard and eight advanced operations through the public wrappers, asserts each result kind and principal shape, then renders analysis, circular, modified-base, G-quadruplex, pseudoknot, conditional density-2, exact landscape saddle, inverse-design target, ligand MEA and occupancy annotation, dot plot, and mountain plot outputs. Crossing structures round-trip directly through `evaluate-pseudoknot` and `evaluate-conditional-density2`. An independent Rust test requires landscape search to select a longer path with a lower saddle over a shorter path with a higher saddle, every returned step to be one pair move, and the path maximum to equal the reported saddle. Inverse design for a 5 nt target evaluates all six canonical-pair choices and all 64 assignments of its three unpaired positions, or `384` sequences; the GC=1 constraint evaluates exactly the 16 explicitly allowed sequences, and tests cover probability ranges and total ranking order. A one-site ligand partition function agrees within `1e-10` with the independent identity:

```math
Z_{\mathrm{joint}} = Z_{\mathrm{RNA}} + Z_{\mathrm{target}}\exp\left(-\frac{\Delta G_{\mathrm{eff}}}{RT}\right).
```

A separate test prevents simultaneous occupancy of two completely overlapping sites.

`scripts/validate-exact-features-real.py` repeats these three exact features through the distributed WASM and Typst wrapper on the first 10 nt of all 24 Rfam cases. All 24 landscapes use complete state spaces, every adjacent path state differs by one pair move, and the maximum saddle discrepancy is zero. Fixed-template inverse design evaluates 1/1 sequence in every case, and the maximum ligand pair/unpaired mass error is `2.220446049250313e-16`. Every row is stored in `tests/reports/exact_features_real_validation.json`.

In addition to `wasm-tools validate`, the suite inspects the binary interface:

```text
imports: typst minimal protocol write_args_to_buffer / send_result_to_host
exports: memory, run, toolchain data/heap globals
```

The byte length and SHA-256 of the release WASM are pinned on every build in `tests/reports/wasm_validation.json`.

## Image-level validation

`tests/typst/plot_quality.typ` validates WCAG 2.2 sRGB relative luminance, including a 21:1 black/white ratio, optimal text colors over dark and light fills, the AA 4.5:1 and AAA 7:1 thresholds, and compositing of transparent colors over a background. `scripts/validate-contrast.py` checks stable explicit errors for unattainable AAA contrast on mid-luminance colors, transparent fills without a background, and gradient fills. `tests/typst/publication.typ` and the README continuous-annotation image pin white text over dark nodes with pixel goldens.

`tests/typst/plot_quality.typ` fixes structure comparison, two-ensemble dot plots, reference-pair dot plots, continuous scales, and expected and discrete mountain profiles across two pages. `scripts/validate-plot-quality.py` requires exact agreement with the hand-calculated four-nucleotide expected profile `[0.5, 0.75, 0.5, 0.0]`, reference profile `[1, 2, 1, 0]`, and comparison pair counts. The release gate also requires zero PDF `/Image` XObjects, vector text, acceptable margins and ink coverage, and pixel-golden equality on every page.

`tests/typst/plot_layout_quality.typ` fixes all four outer sides, nine inner anchors, arbitrary viewport coordinates, offsets, row and column flow, column counts, legends shared by multiple plots, custom/reversed/logarithmic axes, major and minor grids, aspect and padding controls, open frames, `x2` and `y2` secondary axes, and vertical, reversed, and custom-formatted continuous legends across three pages. `scripts/validate-plot-layout.py` requires exact expected and discrete profiles and exact series-to-axis assignments, then verifies native vectors, semantic metadata, margins, ink coverage, and pixel SHA-256 on every page.

`tests/typst/real_data_render.typ` renders three pages for each of the 24 Rfam cases, for 72 pages total:

- Page 1: full-length NAView, Turtle, Puzzler, and linear layouts
- Page 2: cofold, local accessibility, circular, modified-base, G-quadruplex, pseudoknot, comparative, dot plot, DNA model, custom normalized model, fatgraph-topology annotation, and MFE/MEA comparison
- Page 3: complete-state-space minimum-saddle landscape, complete fixed-template inverse folding, and complete ligand microstate ensemble

`scripts/validate-rendered-pdf.py` checks qpdf syntax, page count, all 24 accessions, 18 labels per case, PNG rendering of every page, blank pages, clipping, content margins, and zero PDF `/Image` XObjects. The current baseline has a minimum margin of 20 px, an ink fraction range of `0.005544–0.024680`, and a minimum grayscale standard deviation of `0.054526`.

`scripts/validate-render-golden.py` compares all 72 pages with SHA-256 hashes after 110 dpi Poppler rasterization. A one-pixel change fails the gate; only an explicit `--update-golden` operation updates the baseline. The manifest also pins Typst 0.15.0 and Poppler 25.11.0.

Representative pages 1, 2, 47, and 48 receive an additional visual inspection for clipping or overlap in nucleotide labels, numbering, pair edges, annotations, and captions. The generated real-data QA artifact is `target/qa/ribon-real-data-validation.pdf`.

`tests/typst/publication.typ` is a four-page public rendering contract. It jointly checks circular geometry on a nonsquare canvas, rotation and reflection, measured labels and leaders, multistrand termini, structure differences, continuous-value legends, two-ensemble dot plots, expected and discrete mountain plots, direct rendering from analysis results, level of detail at 600 nt, and manually edited scenes. Semantic metadata pins stem, hairpin, difference, and resource-limit values, while PDF text, native vectors, margins, ink coverage, and pixel goldens verify the final pages. The circular panel has a saturation-mask bounding-box ratio of `0.99038`, within the accepted interval `0.97–1.03`.

## Release gate

Run the same bounded checks as pull-request CI with:

```sh
just ci-check
```

This checks Rust formatting, the bounded unit/integration/reference/protocol test set, Clippy, byte-for-byte WASM synchronization, and Typst smoke/public-API compilation. It excludes the long-sequence log-domain stress test, platform-sensitive performance measurements, exhaustive feature matrices, and release-scale image validation.

Run the complete reproducible release gate with:

```sh
just release-check
```

This task includes the bounded CI gate and adds WASM validation, license boundaries, complete Typst renderer tests, pseudoknots, conditional density-2 exactness and performance, exact-feature validation, general performance, and real-data image validation. It runs locally before publication and in the tag-triggered release workflow rather than on every pull request.

External differential adapters and their reference installations are maintained as local development tools rather than distributed with Ribon. Their pinned results remain available in `tests/reports/`. Differences against an external oracle may reflect either an implementation defect or a declared model boundary. Passing the release gate therefore does not imply identity with RNAstructure, ViennaRNA, or another package. Ribon's reproducibility unit is the model identifier, parameter-bundle hash, API schema, and package version.
