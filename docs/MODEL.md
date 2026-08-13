# Analysis Models

This document defines the thermodynamic, decoding, pseudoknot, and layout models exposed by Ribon. Citations for the underlying methods are maintained in the [reference bibliography](../package/docs/references.bib) and attached to the corresponding claims in the [complete manual](../package/docs/documentation.pdf).

## Standard RNA parameter family

The public RNA model identifier is `ribon-rnastructure-6.6-rna`. Its stack, hairpin, bulge, internal-loop, 1-by-1, 1-by-2, and 2-by-2 special internal-loop, terminal-mismatch, dangle, special-hairpin, multiloop, and duplex-initiation tables are generated from the official RNAstructure 6.6 `rna.*` free-energy and enthalpy parameters. This is the RNAstructure 6.6 standard RNA family: a Turner 2004 lineage with subsequent RNAstructure revisions. Values away from 37 °C are interpolated as follows:

```math
\Delta G(T) = \Delta H - T\frac{\Delta H - \Delta G_{37}}{310.15\,\mathrm{K}}
```

Public energies use kcal/mol; generated internal tables use centi-kcal/mol.

MFE folding and fixed-structure evaluation accept dangle models 0–3. Model 3 includes flush coaxial stacking between directly adjacent helices in a multiloop; mismatch-mediated coaxial stacking across one intervening nucleotide is not part of this model. Partition functions for models 0 and 2 use log-domain inside/outside dynamic programming. Models 1 and 3 enumerate every planar structure and apply the requested exclusive single-dangle or coaxial-stacking optimization to each fixed-structure energy before forming the Boltzmann sum. Ribon does not substitute model 2 for an odd dangle model, so these paths preserve the requested energy definition at the cost of exponential running time. Monovalent salt corrections use the published polyelectrolyte expression for stacks, loops, multiloops, and duplex initiation.

## DNA and custom parameter families

`ribon-rnastructure-6.6-dna` is generated from the 33 official RNAstructure 6.6 `dna.*` files and identifies the standard DNA family compiled by the Mathews group. It has distinct stack, mismatch, dangle, 1-by-1, 1-by-2, and 2-by-2 internal-loop, loop-initiation, multiloop, duplex-initiation, special-loop, free-energy, and enthalpy tables. Input thymine uses the uracil-shaped thermodynamic index only for table lookup and remains `T` in public results and drawings. The RNA salt correction is not applied to DNA; `dna-model` therefore has no `salt` parameter, and low-level DNA requests reject any `salt_molar` other than the neutral internal value `1.021`. The same `EnergyModel` is used by MFE, partition functions, evaluation, sampling, suboptimal folding, accessibility, duplex and cofold, local and circular folding, pseudoknot, conditional density-2, and comparative operations.

`ribon-custom-thermodynamic-v1` starts from an RNA or DNA family and replaces complete normalized centi-kcal/mol tables field by field. Before analysis, Ribon validates the schema version, every table dimension, finite values, ranges, special-loop alphabet and length constraints, duplicate entries, and the 64-character SHA-256 provenance field. Only omitted fields inherit from the base family. Modified-base and G-quadruplex models reject a DNA base family because their calibrations are RNA-specific; a custom RNA overlay composes its canonical tables with those corrections.

## MFE and ensemble calculations

MFE folding uses a pseudoknot-free interval dynamic program. Internal loops are not truncated at 30 nt: every geometrically possible size is considered, and loop initiations beyond the table range use the published Turner logarithmic extrapolation. The Rust API limits this search only when `EnergyModel.max_internal_loop = Some(n)` is set explicitly.

Partition functions use log-domain inside/outside recurrences:

```math
\begin{aligned}
Z &= \sum_S \exp\left(-\frac{E(S)}{RT}\right), \\
G_{\mathrm{ensemble}} &= -RT\log Z, \\
p(i,j) &= \frac{Z_{i,j\,\mathrm{paired}}}{Z}.
\end{aligned}
```

Marginals obey the following identity at every position:

```math
p(i\ \mathrm{unpaired}) + \sum_j p(i,j) = 1.
```

`log_partition_function` remains finite when the ordinary `partition_function` exceeds the floating-point range.

The centroid decoder minimizes expected base-pair distance. The MEA decoder uses the following pair and unpaired-position scores, then finds the maximum-scoring noncrossing structure:

```math
s_{\mathrm{pair}}(i,j) = 2\gamma p(i,j), \qquad s_{\mathrm{unpaired}}(i) = p(i\ \mathrm{unpaired}).
```

Ensemble summaries report the expected number of pairs, mean base-pair distance, and positional entropy.

## Constraints and probing

The shared constraint model supports forced-unpaired positions, forced-paired positions, forced and forbidden pairs, maximum pair span, GU-pair exclusion, lonely-pair exclusion, position-specific paired/unpaired/stack pseudo-energies, and pair-specific energies. Contradictory constraints, multiple forced partners, and crossing forced pairs are rejected before dynamic programming begins.

SHAPE and DMS observations are converted to generic pseudo-energies with the Deigan or Zarringhalam method. MFE, partition, decoding, and fixed-structure evaluation consume the same compiled values.

## Specialized analyses

- `sample` performs seeded Boltzmann backtracking from the same inside dynamic program.
- `suboptimal` returns deterministic k-best structures within a specified kcal/mol band above the MFE.
- `accessibility` computes the exact probability that a contiguous region is jointly unpaired and its opening energy.
- `local` aggregates pair and accessibility tracks from overlapping window ensembles.
- `duplex` uses a connected antiparallel grammar in which every pair is intermolecular.
- `cofold` considers all noncrossing structures across a strand break, including intramolecular hairpins, internal loops, and multiloops within associated states, and solves concentration-dependent mass balance for arbitrary A, B, AA, AB, and BB species. Dangle models 0 and 2 use polynomial inside/outside dynamic programming; models 1 and 3 use complete fixed-structure enumeration.
- `circular` uses a circular root grammar with open, exterior-hairpin, exterior-internal-loop, and exterior-multiloop cases. It accepts dangle models 0–3 and the shared hard, soft, and probing constraints.
- `modified` adds measured sparse nearest-neighbor corrections for m6A, pseudouridine, inosine, 7-deazaadenosine, and purine/nebularine, plus the published model correction for dihydrouridine, to the shared MFE and partition grammar. Unmeasured contexts use the model's canonical reference state, and explicit paired, unpaired, and stack pseudo-energies may be applied concurrently.
- `gquad` integrates G-quadruplexes as exterior or multiloop terminals in the MFE and partition grammar and permits multiple or nested states.
- `comparative` averages ungapped loop energies across alignment rows by column and adds covariation within an MFE/McCaskill grammar.
- `landscape` enumerates every planar structure allowed by the constraints and runs minimax Dijkstra search on the graph of single-base-pair insertions and deletions. The saddle path and local minima are exact for the complete state space; no beam, indirect-path width, or state cap is used.
- `inverse-design` enumerates every sequence that satisfies an IUPAC template, the target's canonical pairs, and an explicit GC interval. It computes exact target probability from the ratio of a target-constrained partition function to the full partition function. The result limit truncates only the ranked output, not the search.
- `ligand` enumerates the product of all planar RNA structures and all independent sets of compatible motif sites. Overlapping sites cannot be occupied together, while nonoverlapping sites may be. RNA pair probabilities and site occupancies come from the same joint partition function. Each occupied site contributes an effective free energy defined by:

```math
\Delta G_{\mathrm{eff}} = \Delta G^\circ - RT\ln\left(\frac{c}{1\,\mathrm{M}}\right).
```

## Pseudoknots

`pseudoknot` reports two explicitly distinct result families:

1. A ProbKnot structure selected from standard pseudoknot-free pair probabilities by the mutual-maximum rule, using the published defaults of threshold 0, one iteration, a minimum helix length of 3, and one-nucleotide bulges treated as helix continuations.
2. An ensemble formed from density-2 H-type components, each consisting of two crossing helices, with span-disjoint components allowed in the same state.

`structure` and `pairs` contain the pure ProbKnot result. Fields prefixed with `hybrid_` contain an additional prediction seeded by the most favorable thermodynamic H-type core. `matching_centroid_*` and `matching_mea_*` exploit the purine–pyrimidine bipartition of the canonical-pair graph and use the Hungarian algorithm to obtain a cubic-time global optimum while allowing arbitrary crossing topology. These decoder families are not presented as the same scientific object.

The H-type state energy contains Turner stack and terminal contributions, the DP09 exterior-pseudoloop initiation term, per-band penalties, and enclosed-unpaired penalties; the open state is always present. By default, there is no limit on components, helices, or states. Weighted-interval inside/outside dynamic programming returns exact MFE, log partition function, pair and unpaired marginals, centroid, and MEA over all sets of span-disjoint H-type components. The potentially exponential state count is returned exactly as a decimal string in `restricted_state_count_exact`, while the numeric field saturates at the target runtime's `usize::MAX`. A caller-supplied `max-components` or `max-ensemble-states` limit sets `restricted_state_space_complete` to `false`. Two cores with overlapping spans represent one higher-order pseudoknot topology and are not treated as independent components.

With `exact-arbitrary-ensemble: true`, Ribon instead enumerates every canonical base-pair matching that satisfies the minimum-loop condition, without restricting crossing graphs or topology. MFE, log partition function, pair and unpaired marginals, centroid, and MEA are exact over that finite state space. Because no beam hides the NP-hard general-pseudoknot problem, this path uses exponential time and space and runs only when requested explicitly. Its energy is Ribon's published generalized DP09 diagnostic model. This arbitrary-topology ensemble is distinct from both the fast span-disjoint H-type ensemble and a fixed-seed density-2 conditional partition function.

`conditional-density2` handles the fixed-seed case. Given a fixed pseudoknot-free seed called G, it considers a second pseudoknot-free layer called G-prime that shares no nucleotide with G. The published 15-state decomposition `W, V, P, Pg, Pgw, Vp, Vpr, Vpl, Vm, Wm, Wm1, Wmp, Wi, Wip, BE` and its primed and unprimed borders are represented as an acyclic Rust interval hypergraph. Every derivation is unique. For dangle models 0 and 2, the complete allowed density-2 state space is evaluated in cubic time and quadratic memory without a beam, enumeration cap, or state-count cap.

One edge generator is evaluated with a log-sum-exp semiring for inside values, min-plus for MFE, and max-plus for centroid and MEA. Pair marginals are derivatives obtained by a reverse outside pass over the same hypergraph. Seed pairs therefore have probability 1, and variable-pair plus unpaired probability is normalized at every position. The outside calculation extends the published decomposition with base-pair probabilities without introducing a different approximate grammar. `evaluate-conditional-density2` restricts allowed pairs to a supplied G-prime layer and uses the same maximum-pair and minimum-energy semirings and backtrace to reproduce the input structure and energy. An independent exponential enumerator is used only for short-input validation and is not part of the public API.

Exclusive single dangles in model 1 and coaxial stacking in model 3 create nonlocal competition between adjacent stems. Reusing model-2 local factors in the 15-state partition function would therefore be an approximation. The public API instead dispatches models 1 and 3 to complete fixed-structure enumeration. That path has no state or energy beam and returns exact MFE, partition, marginal, centroid, and MEA values, but reports `exponential` time and space complexity. Calling the polynomial model-0/2 function directly with model 1 or 3 returns an explicit error.

The conditional energy model stores the following published coefficients as separate fields:

```math
\begin{aligned}
a &= 3.39, & b &= 0.03, & c &= 0.02, \\
P_s &= -1.38, & P_{sm} &= 10.07, & P_{sp} &= 15.00, \\
P_b &= 2.46, & P_{up} &= 0.06, & P_{ps} &= 0.96\ \mathrm{kcal\,mol^{-1}}, \\
s_{\mathrm{stack}} &= 0.89, & s_{\mathrm{internal}} &= 0.74, \\
a' &= 3.41, & b' &= 0.56, & c' &= 0.12.
\end{aligned}
```

Ordinary loops use the shared Turner evaluator. Fixed-seed energy cancels from conditional probabilities but remains part of the partition function and ensemble free energy.

The upper bound in supplementary Equation (iii) is written as:

```math
b(i,l) < l.
```

Under the accompanying border definition, that bound would make the sum empty. Ribon instead uses the following bound, consistent with the interval semantics of the decomposition and the recurrence published by the same authors for HFold:

```math
b(i,j).
```

The term is positive infinity when the required border does not exist, and the base case is:

```math
BE(i,i) = 0.
```

This base case corresponds to the multiplicative identity 1 in the Boltzmann semiring. Independent exhaustive enumeration over small planar seeds fixes these interpretations as testable contracts.

## Layout models

- `naview` is a Rust implementation of classic modified-radial loop and helix geometry.
- `simple` uses an RNAplot-style regular-bond radial polygon.
- `circular` places nucleotides equidistantly on a circle.
- `linear` draws a single strand as a backbone with pair arcs. For multiple strands, it alternates strand directions across parallel rows, draws intermolecular pairs as rungs, and draws intramolecular pairs as arcs outside each row.
- `turtle` performs affine loop and helix traversal using circle radii solved from loop-chord constraints and emits backbone arcs suitable for Typst vector primitives.
- `puzzler` begins from a Turtle layout and applies deterministic monotone collision reduction over nucleotide–nucleotide, nucleotide–edge, and segment-crossing terms. It terminates by numerical convergence of the objective rather than an iteration cap and guarantees zero crossings for planar input.

NAView and simple layouts are regressed against pairwise-distance signatures from external coordinate fixtures. Higher bracket levels do not enter the planar scaffold; they are drawn as additional pair edges.

## Declared boundaries

- The canonical global partition function is pseudoknot-free.
- The fixed-seed density-2 conditional partition function for dangle models 0 and 2 is a polynomial implementation of the published decomposition. Models 1 and 3 dispatch to complete enumeration so that nonlocal dangle and coaxial choices remain exact. The release contract is defined by the published equations, independent exhaustive enumeration, planar reductions, finite differences, and real-data backtraces rather than bitwise compatibility with another executable.
- Partition functions and consumers requested with odd dangle models apply fixed-structure single-dangle or coaxial energy to every state. This includes sampling, suboptimal, accessibility, local, circular, and cofold operations and is exponential for general long inputs.
- Modified-base parameters are sparse corrections relative to a canonical reference state. An unmeasured context contributes a zero correction and is counted in `canonical_reference_stacks`. Dihydrouridine uses a published model correction rather than a measured nearest-neighbor table.
- The comparative model evaluates loop energy on each ungapped alignment row and adds covariance, but does not claim bitwise compatibility with historical covariance scaling in other programs.
- Puzzler is an independent collision-reduction implementation and does not claim bitwise coordinate compatibility with another implementation.
