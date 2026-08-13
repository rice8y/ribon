#import "../../package/lib.typ": analysis-model, conditional-density2, data, render

#set page(width: 18cm, height: 18cm, margin: 9mm)
#set text(font: "Libertinus Serif", size: 8pt)

#let corpus = json("../data/pseudoknot_real_24.json")
#let model = analysis-model(min-loop: 0, dangles: 0)
#let layer(structure, opening, closing) = structure.clusters().map(symbol => {
  if symbol == opening or symbol == closing { symbol } else { "." }
}).join()

#for (index, case) in corpus.cases.enumerate() {
  let response = conditional-density2(
    case.sequence,
    layer(case.structure, "(", ")"),
    model: model,
  )
  let result = data(response)
  let energy = result.mfe_energy_kcal_mol
  let display-energy = calc.round(energy, digits: 3)
  [
    = #case.id — conditional density-2 reference state

    #table(
      columns: (1fr, 1fr),
      inset: 3pt,
      [*length*], [*conditional ΔG*],
      [#case.sequence.len() nt],
      [#display-energy kcal/mol],
    )

    #v(3mm)
    #grid(
      columns: (1fr, 1fr, 1fr),
      gutter: 2mm,
      [*MFE* #render(response, which: "mfe", method: "circular", width: 4.9cm, height: 12.2cm, numbering: 10)],
      [*centroid* #render(response, which: "centroid", method: "circular", width: 4.9cm, height: 12.2cm, numbering: 10)],
      [*MEA* #render(response, which: "mea", method: "circular", width: 4.9cm, height: 12.2cm, numbering: 10)],
    )
  ]
  if index + 1 < corpus.cases.len() { pagebreak() }
}
