#import "../../package/lib.typ": *

#set page(
  width: 297mm,
  height: 210mm,
  margin: (x: 10mm, y: 9mm),
  footer: context align(center, text(size: 6.5pt, fill: luma(45%), [
    Ribon - Rfam real-data vector rendering - page #counter(page).display()
  ])),
)
#set text(font: "New Computer Modern", size: 8pt)

#let corpus = json("../data/rfam_real_24.json")
#let cases = corpus.cases
#let custom-profile = custom-model(
  thermodynamic-parameter-overrides(
    "render-regression-terminal-au",
    "abababababababababababababababababababababababababababababababab",
    tables: ("misc": (410, 360, 500, 370)),
  ),
)
#let reverse-complement(sequence) = sequence.clusters().rev().map(base => (
  A: "U",
  C: "G",
  G: "C",
  U: "A",
  R: "N",
  Y: "N",
  S: "N",
  W: "N",
  K: "N",
  M: "N",
  B: "N",
  D: "N",
  H: "N",
  V: "N",
  N: "N",
).at(base, default: "N")).join("")

#for (case-index, case) in cases.enumerate() {
  let n = case.length
  let radius = if n <= 160 { 2.8pt } else if n <= 300 { 1.8pt } else { 1.05pt }
  let label-size = if n <= 160 { 4.2pt } else if n <= 300 { 2.6pt } else { 1.35pt }
  let interval = if n <= 180 { 20 } else if n <= 320 { 50 } else { 100 }

  align(center, text(size: 13pt, weight: "semibold", [
    #case.accession - #case.description
  ]))
  v(2pt)
  align(center, text(size: 6.5pt, fill: luma(42%), [
    #case.family_id | #case.sequence_id | #n nt | #case.structure_source
  ]))
  v(5pt)

  grid(
    columns: (1fr, 1fr, 1fr),
    gutter: 3mm,
    figure(
      draw(
        case.sequence,
        structure: case.reference_structure,
        method: "naview",
        width: 88mm,
        height: 67mm,
        node-radius: radius,
        font-size: label-size,
        numbering: interval,
      ),
      caption: [NAView],
    ),
    figure(
      draw(
        case.sequence,
        structure: case.reference_structure,
        method: "turtle",
        width: 88mm,
        height: 67mm,
        node-radius: radius,
        font-size: label-size,
        numbering: interval,
      ),
      caption: [Affine loop geometry],
    ),
    figure(
      draw(
        case.sequence,
        structure: case.reference_structure,
        method: "puzzler",
        width: 88mm,
        height: 67mm,
        node-radius: radius,
        font-size: label-size,
        numbering: interval,
      ),
      caption: [Collision-reduced loop geometry],
    ),
  )
  v(2pt)
  figure(
    draw(
      case.sequence,
      structure: case.reference_structure,
      method: "linear",
      width: 274mm,
      height: 55mm,
      node-radius: radius,
      font-size: label-size,
      numbering: interval,
    ),
    caption: [Linear arc representation],
  )

  pagebreak()

  let fragment = case.sequence.slice(0, calc.min(24, n))
  let interaction-a = fragment.slice(0, calc.min(12, fragment.len()))
  let interaction-b = reverse-complement(interaction-a)
  let analysis = analyze(fragment)
  let dna-analysis = analyze(fragment.replace("U", "T"), model: dna-model())
  let custom-analysis = analyze(fragment, model: custom-profile)
  let exact-sequence = fragment.slice(0, calc.min(10, fragment.len()))
  let exact-open = range(exact-sequence.len()).map(_ => ".").join()
  let exact-model = analysis-model(dangles: 0)
  let exact-target = data(fold(exact-sequence, model: exact-model)).structure
  let exact-path = landscape(exact-sequence, exact-open, exact-target, model: exact-model)
  let exact-design = inverse-design(
    exact-target,
    template: exact-sequence,
    return-count: 1,
    model: exact-model,
  )
  let exact-ligand = ligand(
    exact-sequence,
    (ligand-motif(
      case.accession + "-site",
      1,
      exact-sequence,
      exact-target,
      -1.0,
      concentration: 0.001,
    ),),
    model: exact-model,
  )
  let local = local(
    fragment,
    window-size: calc.min(16, fragment.len()),
    max-pair-span: calc.min(12, fragment.len() - 1),
    max-unpaired: 1,
  )
  let circular = circular(fragment)
  let modified = modified(fragment, (
    modified-base(1, "mod", fragment.clusters().first(), stack-energy: -0.1),
  ))
  let gquad = gquad(fragment)
  let pseudoknot = pseudoknot(fragment, options: pseudoknot-options(
    max-components: 4,
    max-ensemble-states: 512,
  ))
  let pseudoknot-structure = data(pseudoknot).hybrid_structure
  let topology = fatgraph-topology(fragment, pseudoknot-structure)
  let comparative = comparative((fragment, fragment, fragment))
  let complex = cofold(interaction-a, interaction-b)

  align(center, text(size: 13pt, weight: "semibold", [
    #case.accession - extended analysis/render matrix
  ]))
  v(2pt)
  align(center, text(size: 6.5pt, fill: luma(42%), [
    real #fragment.len()-nt prefix | cofold - local - circular - modified - G-quad - pseudoknot - comparative
  ]))
  v(5pt)

  grid(
    columns: (1fr, 1fr, 1fr, 1fr),
    gutter: 2.5mm,
    row-gutter: 3mm,
    figure(
      render(complex, width: 65mm, height: 40mm, numbering: none),
      caption: [AA/AB cofold macrostate],
    ),
    figure(
      draw(
        fragment,
        structure: data(analysis).mfe_structure,
        width: 65mm,
        height: 40mm,
        numbering: none,
        probabilities: data(local).pair_probabilities,
        annotations: local-accessibility-annotations(local),
        legend: false,
      ),
      caption: [Local pair/accessibility],
    ),
    figure(
      render(circular, width: 65mm, height: 40mm, numbering: none),
      caption: [Circular ensemble MFE],
    ),
    figure(
      render(modified, width: 65mm, height: 40mm, numbering: none),
      caption: [Modified-base profile],
    ),
    figure(
      render(gquad, width: 65mm, height: 40mm, numbering: none),
      caption: [G-quadruplex states],
    ),
    figure(
      render(pseudoknot, method: "circular", width: 65mm, height: 40mm, numbering: none),
      caption: [Multi-component H-type pseudoknot ensemble],
    ),
    figure(
      render(comparative, width: 65mm, height: 40mm, numbering: none),
      caption: [Gap-aware covariation consensus],
    ),
    figure(
      dot-plot(fragment, probabilities: data(pseudoknot).source_pair_probabilities, width: 40mm, threshold: 0.02),
      caption: [Decoder source probabilities],
    ),
    figure(
      render(dna-analysis, width: 65mm, height: 40mm, numbering: none),
      caption: [Mathews DNA model],
    ),
    figure(
      render(
        custom-analysis,
        width: 65mm,
        height: 40mm,
        numbering: none,
        annotations: entropy-annotations(custom-analysis),
        legend: false,
      ),
      caption: [Custom normalized model],
    ),
    figure(
      draw(
        fragment,
        structure: pseudoknot-structure,
        method: "circular",
        width: 65mm,
        height: 40mm,
        numbering: none,
        annotations: topology-annotations(topology),
      ),
      caption: [Fatgraph topology annotation],
    ),
    figure(
      compare-structures(
        fragment,
        data(analysis).mfe_structure,
        data(analysis).mea_structure,
        width: 65mm,
        height: 40mm,
        numbering: none,
      ),
      caption: [MFE/MEA structure comparison],
    ),
  )

  pagebreak()

  align(center, text(size: 13pt, weight: "semibold", [
    #case.accession - exact state-space analyses
  ]))
  v(2pt)
  align(center, text(size: 6.5pt, fill: luma(42%), [
    real #exact-sequence.len()-nt prefix | complete landscape - complete inverse search - complete ligand microstates
  ]))
  v(5pt)

  grid(
    columns: (1fr, 1fr, 1fr),
    gutter: 3mm,
    figure(
      render(exact-path, which: "saddle", method: "circular", width: 88mm, height: 70mm, numbering: none),
      caption: [Exact minimum-saddle landscape],
    ),
    figure(
      render(exact-design, which: "target", method: "circular", width: 88mm, height: 70mm, numbering: none),
      caption: [Exact inverse folding],
    ),
    figure(
      render(exact-ligand, which: "mea", method: "circular", width: 88mm, height: 70mm, numbering: none),
      caption: [Exact ligand microstates],
    ),
  )

  if case-index + 1 < cases.len() { pagebreak() }
}
