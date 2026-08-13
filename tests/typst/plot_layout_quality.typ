#import "../../package/lib.typ": *

#set page(width: 210mm, height: 297mm, margin: 10mm)
#set text(font: "Libertinus Serif", size: 7pt)

#let sequence = "GCGC"
#let probabilities = (
  (i: 1, j: 4, probability: 0.5),
  (i: 2, j: 3, probability: 0.25),
)
#let references = (
  (label: [Nested], structure: "(())"),
  (label: [Outer], structure: "(..)"),
)
#let exact-profile = mountain-profile(
  sequence,
  probabilities: probabilities,
  reference-structures: references,
)
#let compact-theme = plot-theme(
  text-size: 5.5pt,
  label-size: 6pt,
  legend-fill: white.transparentize(8%),
  legend-stroke: luma(78%) + 0.3pt,
)

#let specimen(title, body) = block(
  width: 100%,
  inset: 2.5mm,
  stroke: luma(86%) + 0.35pt,
  radius: 1.5mm,
  grid(columns: (1fr,), row-gutter: 1.5mm,
    text(size: 6.5pt, weight: "semibold", title),
    align(center, body),
  ),
)

#let mountain(legend, width: 61mm, height: 35mm, ..args) = mountain-plot(
  sequence,
  probabilities: probabilities,
  reference-structures: references,
  width: width,
  height: height,
  x-label: none,
  y-label: none,
  legend: legend,
  theme: compact-theme,
  ..args,
)

= Outer legend placement

#grid(
  columns: (1fr, 1fr),
  gutter: 4mm,
  specimen([Top; wrapped to two columns], mountain(legend-style(
    position: "top", columns: 2, row-gap: 2pt,
  ))),
  specimen([Bottom; one row], mountain(legend-style(
    position: "bottom", direction: "row",
  ))),
  specimen([Left; vertical flow], mountain(legend-style(
    position: "left", direction: "column", width: 24mm,
  ), width: 55mm)),
  specimen([Right; offset and explicit anchor], mountain(legend-style(
    position: "right", direction: "column", anchor: top + left,
    offset: (2pt, 4pt), width: 24mm,
  ), width: 55mm)),
)

#v(4mm)

#specimen(
  [One shared legend for two independent plots],
  place-legend(
    grid(
      columns: (1fr, 1fr),
      gutter: 5mm,
      mountain(false, width: 72mm, height: 40mm),
      mountain(false, width: 72mm, height: 40mm, x-axis: axis-style(domain: (4, 1))),
    ),
    plot-legend((
      legend-item([Expected], stroke: rgb("#315eaa") + 1pt),
      legend-item([Nested], stroke: (paint: rgb("#d73027"), thickness: 0.8pt, dash: "dashed")),
      legend-item([Outer], stroke: (paint: rgb("#1a9850"), thickness: 0.8pt, dash: "dotted")),
    ), style: legend-style(columns: 2, max-columns: 2)),
    style: legend-style(position: "bottom"),
  ),
)

#pagebreak()

= Inner legend anchors

#let inner-positions = (
  "inner-north-west", "inner-north", "inner-north-east",
  "inner-west", "inner-center", "inner-east",
  "inner-south-west", "inner-south", "inner-south-east",
)

#grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 3mm,
  ..inner-positions.map(position => specimen(
    position,
    mountain(
      legend-style(
        position: position,
        direction: "column",
        columns: 1,
        inset: 2pt,
        item-gap: 3pt,
        row-gap: 1pt,
        fill: white.transparentize(12%),
        stroke: luma(75%) + 0.25pt,
      ),
      width: 51mm,
      height: 38mm,
      layout: plot-layout(padding: (left: 18pt, right: 5pt, top: 5pt, bottom: 15pt)),
    ),
  )),
)

#v(4mm)

#grid(
  columns: (1fr, 1fr),
  gutter: 4mm,
  specimen(
    [Explicit 28%, 32% coordinate; center anchor],
    mountain(legend-style(
      position: (0.28, 0.32),
      anchor: center + horizon,
      direction: "column",
      columns: 1,
      inset: 2pt,
    ), width: 76mm, height: 45mm),
  ),
  specimen(
    [Inner north-east with negative offset],
    mountain(legend-style(
      position: "inner-north-east",
      direction: "column",
      offset: (-4pt, 4pt),
      inset: 2pt,
    ), width: 76mm, height: 45mm),
  ),
)

#pagebreak()

= Axis and plot-area configuration

#grid(
  columns: (1fr, 1fr),
  gutter: 4mm,
  specimen(
    [Custom domain, tick labels, minor grid, and aspect],
    mountain(
      false,
      width: 78mm,
      height: 48mm,
      x-axis: axis-style(
        domain: (1, 4),
        ticks: ((1, [5′]), (2, [two]), (3, [three]), (4, [3′])),
        minor-tick-step: 0.5,
        grid: "both",
      ),
      y-axis: axis-style(domain: (0, 2.5), tick-step: 0.5, grid: "major"),
      layout: plot-layout(aspect: 1.8, frame: true),
    ),
  ),
  specimen(
    [Reversed sequence axis and open L-frame],
    mountain(
      false,
      width: 78mm,
      height: 48mm,
      x-axis: axis-style(domain: (4, 1), tick-step: 1, label: [Reverse position]),
      y-axis: axis-style(domain: (0, 2), tick-step: 1),
      layout: plot-layout(frame: false, padding: (left: 28pt, right: 8pt, top: 8pt, bottom: 28pt)),
    ),
  ),
  specimen(
    [Logarithmic y axis and scientific formatter],
    mountain(
      false,
      width: 78mm,
      height: 48mm,
      y-axis: axis-style(
        domain: (0.01, 10),
        mode: "log",
        base: 10,
        format: "scientific",
        label: [Log scale],
      ),
      layout: plot-layout(aspect: 1.6),
    ),
  ),
  specimen(
    [Top and right secondary axes; one series on y2],
    mountain-plot(
      sequence,
      probabilities: probabilities,
      reference-structures: (
        (label: [Primary], structure: "(..)"),
        (label: [Secondary], structure: "(())", axes: ("x", "y2")),
      ),
      width: 78mm,
      height: 48mm,
      legend: legend-style(position: "inner-south", max-columns: 2, inset: 2pt),
      x-label: [Bottom],
      y-label: [Primary],
      x2-axis: axis-style(
        domain: (1, 4),
        ticks: ((1, [A]), (2, [B]), (3, [C]), (4, [D])),
        label: [Top index],
        grid: none,
      ),
      y2-axis: axis-style(domain: (0, 4), tick-step: 1, label: [Secondary], grid: none),
      theme: compact-theme,
    ),
  ),
)

#v(4mm)

#grid(
  columns: (1fr, 1fr),
  gutter: 4mm,
  specimen(
    [Vertical reversed continuous legend],
    color-legend(
      color-scale(label: [Probability]),
      width: 45mm,
      ticks: 5,
      orientation: "vertical",
      reverse: true,
      format: value => [#calc.round(value * 100)%],
      theme: compact-theme,
    ),
  ),
  specimen(
    [Dot plot with explicit padding and inner legend],
    dot-plot(
      sequence,
      probabilities: probabilities,
      reference-structure: "(())",
      width: 67mm,
      height: 58mm,
      x-axis: axis-style(ticks: ((1, [one]), (4, [four])), grid: "major"),
      y-axis: axis-style(ticks: (1, 4), grid: "major"),
      layout: plot-layout(padding: (left: 28pt, right: 7pt, top: 7pt, bottom: 25pt)),
      legend: legend-style(position: "inner-south-east", direction: "column", width: 25mm, inset: 2pt),
      theme: compact-theme,
    ),
  ),
)

#metadata((
  schema: "ribon.plot-layout-qa/1",
  outer-positions: ("top", "bottom", "left", "right"),
  inner-position-count: inner-positions.len(),
  explicit-position: (0.28, 0.32),
  axis-modes: ("linear", "log"),
  secondary-axes: ("x2", "y2"),
  exact-profile: exact-profile.expected,
  exact-reference-profiles: exact-profile.references.map(series => series.values),
  reference-axes: exact-profile.references.map(series => series.axes),
  secondary-series-axes: mountain-profile(
    sequence,
    probabilities: probabilities,
    reference-structures: ((structure: "(())", axes: ("x", "y2")),),
  ).references.first().axes,
)) <ribon-plot-layout-qa>
