#import "../../package/lib.typ": analysis-model, analyze, capabilities, data, render

#set page(width: 16cm, height: auto, margin: 8mm)
#set text(font: "Libertinus Serif", size: 9pt)

#let model = analysis-model(temperature: 37.0, min-loop: 3, dangles: 2, mea-gamma: 1.0)
#let result = analyze("GGGAAACCC", model: model, id: "api-hairpin")
#let result-data = data(result)
#let caps = data(capabilities())

#assert.eq(result.schema_version, 1)
#assert.eq(result.engine.backend, "pure-rust")
#assert.eq(result.result.kind, "analysis")
#assert.eq(result-data.sequence, "GGGAAACCC")
#assert.eq(result-data.length, 9)
#assert.eq(result-data.pair_probabilities.len() > 0, true)
#assert.eq(caps.native_vector_renderer, true)
#assert.eq(caps.operations.contains("analyze"), true)

= Ribon analysis/1

MFE: #result-data.mfe_structure (#result-data.mfe_energy_kcal_mol kcal/mol)

#render(result, which: "mfe", method: "naview", width: 7cm)

#render(result, which: "centroid", method: "simple", width: 7cm)

#render(result, which: "mea", method: "circular", width: 7cm)
