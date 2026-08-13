#import "../../package/lib.typ": *

#set page(width: 18cm, height: auto, margin: 8mm)
#set text(size: 7pt)

#let sequence = "GGGAAACCC"
#let structure = "(((...)))"

#grid(
  columns: (1fr, 1fr),
  gutter: 6pt,
  ..("naview", "simple", "circular", "linear", "turtle", "puzzler").map(method => [
    #align(center)[*#method*]
    #draw(sequence, structure: structure, method: method, width: 7cm, height: 4cm, numbering: none)
  ]),
)

#let pk-sequence = "GCGCAAAAGCGC"
#let pk-structure = "(([[..))..]]"
#assert.eq(data(validate(pk-sequence, pk-structure)).pairs.len(), 4)

#draw(
  pk-sequence,
  structure: pk-structure,
  method: "circular",
  width: 8cm,
  theme: varna-theme,
  annotations: (
    pair-annotation(1, 8, stroke: (paint: red, thickness: 1.2pt)),
    base-annotation(5, fill: yellow),
    label-annotation(5, [pseudoknot]),
  ),
)

#let split = data(validate("GGG&CCC", "(((&)))"))
#assert.eq(split.strand_breaks, (3,))
#draw("GGG&CCC", structure: "(((&)))", method: "circular", width: 6cm)

#let linear-duplex = data(layout("GGGG&CCCC", "((((&))))", method: "linear"))
#assert.eq(linear-duplex.algorithm, "strand-aware antiparallel linear diagram")
#assert.eq(linear-duplex.strand_breaks, (4,))
#assert(linear-duplex.points.at(0).x < linear-duplex.points.at(3).x)
#assert(linear-duplex.points.at(4).x > linear-duplex.points.at(7).x)
#assert(linear-duplex.points.at(0).y < linear-duplex.points.at(4).y)
#assert(linear-duplex.pairs.all(pair => pair.interstrand))
#assert(linear-duplex.pairs.all(pair => calc.abs(
  linear-duplex.points.at(pair.i - 1).x - linear-duplex.points.at(pair.j - 1).x,
) < 0.000001))

#draw(
  "GGGG&CCCC",
  structure: "((((&))))",
  method: "linear",
  width: 11.6cm,
  height: 4.1cm,
  theme: varna-theme,
  numbering: none,
  show-ends: false,
  show-direction: true,
)

#draw(
  "NNNN&NNNN",
  structure: "([.)&..].",
  method: "linear",
  width: 8cm,
  height: 4cm,
  numbering: none,
  show-ends: false,
  show-direction: true,
)
