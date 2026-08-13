#import "../../package/lib.typ": *

#set page(width: 21cm, height: 29.7cm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 7.5pt)
#set par(leading: 0.65em)

#let panel(title, body) = figure(
  body,
  caption: text(size: 6pt, weight: "semibold", title),
)

#let sequence = "GGGAAACCCGGGAAACCC"
#let reference = "(((...)))(((...)))"
#let alternative = "((....)).(((...)))"
#let model = analysis-model(dangles: 2)
#let prediction = analyze(sequence, model: model)
#let perturbed = analyze(
  sequence,
  model: model,
  constraints: folding-constraints(force-unpaired: (1, 2)),
)
#let samples = sample(sequence, count: 4, seed: 91, model: model)
#let alternatives = suboptimal(sequence, energy-band: 3.0, limit: 4, model: model)
#let evaluated = evaluate(sequence, reference, model: model)
#let access = accessibility(
  sequence,
  range(1, sequence.len() + 1).map(position => accessibility-window(position, position)),
  model: model,
)
#let local-result = local(sequence, window-size: 12, max-pair-span: 10, max-unpaired: 1, model: model)
#let elements = structure-elements(sequence, reference)
#let difference = structure-difference(sequence, reference, alternative)
#let hairpin = elements.loops.find(loop => loop.kind == "hairpin")
#let reactivity-scale = color-scale(
  minimum: 0.0,
  maximum: 1.0,
  colors: (rgb("#2166ac"), rgb("#f7f7f7"), rgb("#b2182b")),
  label: [SHAPE reactivity],
)
#let reactivities = (
  0.02, 0.08, 0.12, 0.86, 0.91, 0.72, 0.15, 0.09, 0.03,
  0.04, 0.07, 0.14, 0.74, 0.88, 0.81, 0.16, 0.06, 0.02,
)
#let label-contract = label-annotation(
  5,
  [contract label],
  dx: 12pt,
  dy: -9pt,
  leader-stroke: (paint: red, thickness: 0.6pt),
  leader-bend: 0.12,
  leader-start-gap: 2pt,
  leader-end-gap: 3pt,
  width: 32pt,
  text-align: center,
  box-fill: white,
  box-stroke: luma(60%) + 0.3pt,
  box-inset: 2pt,
  box-radius: 1.5pt,
)
#let strand-contract = strand-label(2, [target], at: "end", dx: 8pt, dy: 6pt, leader: false)
#let interaction-contract = interaction-annotation(
  5,
  14,
  label: [contact],
  label-position: 0.65,
  label-dx: 7pt,
  label-dy: -4pt,
  label-width: 28pt,
)

#assert.eq(elements.stems.len(), 2)
#assert.eq(elements.loops.filter(loop => loop.kind == "hairpin").len(), 2)
#assert.eq(difference.common.len(), 3)
#assert.eq(difference.reference-only.len(), 3)
#assert.eq(difference.alternative-only.len(), 2)
#assert.eq(data(samples).samples.len(), 4)
#assert(data(alternatives).structures.len() > 0)
#assert.eq(evaluated.result.kind, "energy")
#assert.eq(access.result.kind, "accessibility")
#assert.eq(local-result.result.kind, "local")
#assert.eq(label-contract.dx, 12pt)
#assert.eq(label-contract.dy, -9pt)
#assert.eq(label-contract.leader-bend, 0.12)
#assert.eq(label-contract.leader-start-gap, 2pt)
#assert.eq(label-contract.leader-end-gap, 3pt)
#assert.eq(label-contract.width, 32pt)
#assert.eq(label-contract.text-align, center)
#assert.eq(label-contract.box-inset, 2pt)
#assert.eq(strand-contract.at, "end")
#assert.eq(strand-contract.dx, 8pt)
#assert.eq(strand-contract.leader, false)
#assert.eq(interaction-contract.label-position, 0.65)
#assert.eq(interaction-contract.label-dx, 7pt)
#assert.eq(interaction-contract.label-width, 28pt)

#let oversized = range(0, 501).map(_ => "A").join()
#let rejected = try-request("analyze", (sequence: oversized))
#let excessive-sampling = try-request(
  "sample",
  (sequence: sequence),
  options: (count: 2001, seed: 0, unique: false),
)
#let excessive-landscape = try-request(
  "landscape",
  (
    sequence: "GGGAAACCCGGGAAA",
    start_structure: "...............",
    target_structure: "...............",
  ),
)
#assert.eq(rejected.error.code, "resource_limit")
#assert.eq(excessive-sampling.error.code, "resource_limit")
#assert.eq(excessive-landscape.error.code, "resource_limit")

= Publication rendering contract

#grid(
  columns: (1fr, 1fr),
  gutter: 5mm,
  panel(
    [Aspect-preserving circular layout],
    draw(
      sequence,
      structure: reference,
      method: "circular",
      width: 8.8cm,
      height: 5.2cm,
      fit: "contain",
      theme: varna-theme,
      numbering: numbering-style(every: 5),
    ),
  ),
  panel(
    [Rotation, reflection, and direction],
    draw(
      sequence,
      structure: reference,
      width: 8.8cm,
      height: 5.2cm,
      rotation: 24deg,
      mirror-x: true,
      show-direction: true,
      numbering: numbering-style(positions: (5, 10, 15)),
    ),
  ),
)

#panel(
  [Measured labels, motif highlight, and tertiary interaction],
  draw(
    sequence,
    structure: reference,
    width: 17.5cm,
    height: 6.2cm,
    theme: varna-theme,
    annotations: element-annotations(hairpin, fill: rgb("#ffe082").transparentize(25%)) + (
      label-annotation(4, [apical loop], dx: -35pt, dy: -21pt, leader-bend: 0.08),
      label-annotation(13, [reactive site], dx: 34pt, dy: -18pt),
      interaction-annotation(5, 14, label: [tertiary contact], label-position: 0.62, label-dx: 12pt, label-dy: 14pt),
    ),
  ),
)

#pagebreak()
= Comparison, strands, and quantitative tracks

#grid(
  columns: (1fr, 1fr),
  gutter: 5mm,
  panel(
    [Reference/alternative pair classification],
    compare-structures(
      sequence,
      reference,
      alternative,
      width: 8.7cm,
      height: 5.4cm,
      numbering: none,
    ),
  ),
  panel(
    [Multi-strand identifiers and termini],
    draw(
      "GGGG&CCCC",
      structure: "((((&))))",
      method: "linear",
      width: 8.7cm,
      height: 5.4cm,
      show-direction: true,
      numbering: numbering-style(every: 2, per-strand: true),
      annotations: (
        strand-label(1, [guide], dx: -20pt, dy: -16pt, leader: false),
        strand-label(2, [target], dx: 20pt, dy: 16pt, leader: false),
      ),
    ),
  ),
)

#grid(
  columns: (1fr, 1fr),
  gutter: 5mm,
  panel(
    [Continuous annotation with exact legend],
    draw(
      sequence,
      structure: reference,
      width: 8.7cm,
      height: 5.2cm,
      annotations: value-annotations(
        reactivities,
        scale: reactivity-scale,
        legend: (width: 3.6cm),
      ),
      numbering: none,
    ),
  ),
  panel(
    [Two-ensemble probability dot plot],
    dot-plot(
      sequence,
      probabilities: prediction,
      comparison: perturbed,
      width: 7.4cm,
      threshold: 0.005,
    ),
  ),
)

#panel(
  [Expected and discrete mountain profiles],
  mountain-plot(
    sequence,
    probabilities: prediction,
    reference-structures: (
      (label: [Reference], structure: reference),
      (label: [Alternative], structure: alternative),
    ),
    width: 17cm,
    height: 4.2cm,
  ),
)

#pagebreak()
= Analysis-to-render integration

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 3mm,
  panel([Fixed-structure evaluation], render(evaluated, width: 5.7cm, height: 4.3cm, numbering: none)),
  panel([Boltzmann sample], render(samples, item: 1, width: 5.7cm, height: 4.3cm, numbering: none)),
  panel([Suboptimal structure], render(alternatives, item: 1, width: 5.7cm, height: 4.3cm, numbering: none)),
  panel([Accessibility overlay], render(access, width: 5.7cm, height: 4.3cm, numbering: none)),
  panel([Local accessibility overlay], render(local-result, width: 5.7cm, height: 4.3cm, numbering: none)),
  panel([Integrated MFE], render(prediction, width: 5.7cm, height: 4.3cm, numbering: none)),
)

#panel(
  [MFE, centroid, and MEA comparison],
  grid(
    columns: 3,
    gutter: 3mm,
    ..("mfe", "centroid", "mea").map(which => render(
      prediction,
      which: which,
      width: 5.4cm,
      height: 4.2cm,
      numbering: none,
    )),
  ),
)

#pagebreak()
= Long-RNA level of detail

#let long-sequence = range(0, 150).map(_ => "ACGU").join()
#let long-structure = range(0, 600).map(_ => ".").join()
#assert.eq(long-sequence.len(), 600)

#panel(
  [600 nt linear overview; automatic backbone-only rendering],
  draw(
    long-sequence,
    structure: long-structure,
    method: "linear",
    width: 18cm,
    height: 4cm,
    detail: auto,
    numbering: numbering-style(every: 100),
  ),
)

#let edited-scene = data(layout(sequence, reference, method: "naview"))
#let moved = edited-scene.points.at(3)
#let revised-points = edited-scene.points.enumerate().map(((index, point)) => if index == 3 {
  (x: point.x - 0.18, y: point.y - 0.12)
} else { point })
#edited-scene.insert("points", revised-points)

#panel(
  [Hand-edited scene coordinates],
  render-scene(
    edited-scene,
    width: 12cm,
    height: 7cm,
    annotations: (label-annotation(4, [manually positioned]),),
  ),
)

#metadata((
  schema: "ribon.publication-qa/1",
  sequence-length: sequence.len(),
  stems: elements.stems.len(),
  hairpins: elements.loops.filter(loop => loop.kind == "hairpin").len(),
  shared-pairs: difference.common.len(),
  reference-only-pairs: difference.reference-only.len(),
  alternative-only-pairs: difference.alternative-only.len(),
  annotation-controls: ("nucleotide-label", "strand-label", "interaction-label"),
  resource-errors: (rejected.error.code, excessive-sampling.error.code, excessive-landscape.error.code),
)) <ribon-publication-qa>
