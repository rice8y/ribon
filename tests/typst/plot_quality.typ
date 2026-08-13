#import "../../package/lib.typ": *

#set page(width: 210mm, height: 297mm, margin: 14mm)
#set text(font: "Libertinus Serif", size: 8pt)
#set par(leading: 0.7em)

#let panel(title, body) = block(
  width: 100%,
  inset: 4mm,
  stroke: luma(86%) + 0.4pt,
  radius: 2mm,
  grid(
    columns: (1fr,),
    row-gutter: 3mm,
    text(size: 8pt, weight: "semibold", title),
    align(center, body),
  ),
)

// This hand-checkable fixture proves that plot refactoring does not alter the
// numerical definition of an expected or discrete mountain profile.
#let exact-profile = mountain-profile(
  "GCGC",
  probabilities: (
    (i: 1, j: 4, probability: 0.5),
    (i: 2, j: 3, probability: 0.25),
  ),
  reference-structures: ("(())",),
)
#assert.eq(exact-profile.expected, (0.5, 0.75, 0.5, 0.0))
#assert.eq(exact-profile.references.first().values, (1, 2, 1, 0))
#assert.eq(
  plot-theme(legend-stroke: none),
  plot-theme() + (legend-stroke: none),
)
#assert.eq(type(color-legend(color-scale(), stroke: none)), content)
#let annotation-track = value-annotations((0.0, 1.0), scale: color-scale(label: [Test scale]))
#assert.eq(annotation-track.kind, "annotation-track")
#assert.eq(annotation-track.annotations.len(), 2)
#assert.eq(annotation-track.legends.len(), 1)
#assert.eq(annotation-track.legends.first().scale.label, [Test scale])
#assert(calc.abs(wcag-contrast-ratio(black, white) - 21.0) < 0.0000001)
#let mid-blue = rgb("#4575b4")
#let mid-blue-best = wcag-contrast-ratio(white, mid-blue)
#assert(calc.abs(mid-blue-best - 4.712096296624817) < 0.0000001)
#assert(calc.abs(wcag-contrast-ratio(black, mid-blue) - 4.456615204371331) < 0.0000001)
#assert.eq(wcag-text-fill(mid-blue), white)
#assert.eq(wcag-text-fill(rgb("#313695"), level: "aaa"), white)
#assert.eq(wcag-text-fill(rgb("#74add1")), black)
#assert.eq(wcag-text-fill(mid-blue, level: "aaa", on-failure: "best"), white)
#let translucent-blue = rgb("#313695").transparentize(20%)
#let translucent-text = wcag-text-fill(translucent-blue, background-behind: white)
#assert(wcag-contrast-ratio(translucent-text, translucent-blue, background-behind: white) >= 4.5)

#let sequence = "GGGAAACCCGGGAAACCC"
#let reference = "(((...)))(((...)))"
#let alternative = "((....)).(((...)))"
#let prediction = analyze(sequence)
#let perturbed = analyze(
  sequence,
  constraints: folding-constraints(force-unpaired: (1, 2)),
)

= Quantitative figure contract

#grid(
  columns: (1fr, 1fr),
  gutter: 5mm,
  panel(
    [Categorical structure legend],
    compare-structures(
      sequence,
      reference,
      alternative,
      width: 78mm,
      height: 50mm,
      numbering: none,
    ),
  ),
  panel(
    [Two-ensemble probability matrix],
    dot-plot(
      sequence,
      probabilities: prediction,
      comparison: perturbed,
      width: 76mm,
      threshold: 0.005,
    ),
  ),
)

#v(5mm)

#panel(
  [Expected and discrete mountain profiles],
  mountain-plot(
    sequence,
    probabilities: prediction,
    reference-structures: (
      (label: [Reference], structure: reference),
      (label: [Alternative], structure: alternative),
    ),
    width: 166mm,
    height: 58mm,
  ),
)

#pagebreak()

= Scale and theme contract

#let compact-theme = plot-theme(
  text-size: 7pt,
  label-size: 7.5pt,
  grid-stroke: (paint: rgb("#dbe4ee"), thickness: 0.35pt),
  frame-stroke: (paint: rgb("#34495e"), thickness: 0.65pt),
  tick-stroke: (paint: rgb("#34495e"), thickness: 0.65pt),
  legend-fill: rgb("#f8fafc"),
  legend-stroke: (paint: rgb("#cbd5e1"), thickness: 0.4pt),
)

#grid(
  columns: (1fr, 1fr),
  gutter: 5mm,
  panel(
    [Reference-pair lower triangle],
    dot-plot(
      sequence,
      probabilities: prediction,
      reference-structure: reference,
      width: 75mm,
      threshold: 0.005,
      theme: compact-theme,
    ),
  ),
  panel(
    [Continuous color legend],
    draw(
      sequence,
      structure: reference,
      width: 75mm,
      height: 52mm,
      numbering: none,
      annotations: value-annotations(
        range(0, sequence.len()).map(index => index / (sequence.len() - 1)),
        scale: color-scale(label: [Normalized position]),
        legend: (width: 55mm, ticks: 5),
      ),
      legend-theme: compact-theme,
    ),
  ),
)

#v(5mm)

#panel(
  [Compact scientific axes without a legend],
  mountain-plot(
    "GCGC",
    probabilities: (
      (i: 1, j: 4, probability: 0.5),
      (i: 2, j: 3, probability: 0.25),
    ),
    reference-structures: ("(())",),
    width: 166mm,
    height: 55mm,
    legend: false,
    theme: compact-theme,
  ),
)

#metadata((
  schema: "ribon.plot-qa/1",
  exact-expected: exact-profile.expected,
  exact-reference: exact-profile.references.first().values,
  exact-length: exact-profile.length,
  comparison-pair-count: structure-difference(
    sequence,
    reference,
    alternative,
  ).alternative-only.len(),
)) <ribon-plot-qa>
