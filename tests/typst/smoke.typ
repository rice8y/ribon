#import "../../package/lib.typ": *

#let sequence = "GGGAAACCC"
#let model = analysis-model(dangles: 2)
#let result = analyze(sequence, model: model)
#let analysis = data(result)

#assert.eq(result.schema_version, 1)
#assert.eq(result.engine.backend, "pure-rust")
#assert.eq(analysis.model.dangles, 2)
#assert.eq(analysis.mfe_structure, "(((...)))")
#assert.eq(data(validate(sequence, analysis.mfe_structure)).length, 9)
#assert.eq(data(layout(sequence, analysis.mfe_structure, method: "circular")).points.len(), 9)
#assert(analysis.partition_function >= 1.0)
#assert(analysis.pair_probabilities.len() > 0)
#assert(data(evaluate(sequence, analysis.mfe_structure, model: model)).total_kcal_mol < 0.0)
#assert.eq(data(fold(sequence, model: analysis-model(dangles: 3))).dangles, 3)
#let dna = data(analyze("GGGTTTCCC", model: dna-model(dangles: 0)))
#assert.eq(dna.mfe_structure, "(((...)))")
#assert(calc.abs(dna.mfe_energy_kcal_mol + 0.2) < 0.000000001)
#assert(dna.model.parameter_set.contains("DNA"))
#let dna-conditional = data(conditional-density2(
  "GGGTTTCCC",
  ".........",
  model: dna-model(dangles: 0),
))
#assert.eq(dna-conditional.mfe_structure, "[[[...]]]")
#assert(dna-conditional.model.contains("DNA"))
#let custom = custom-model(
  thermodynamic-parameter-overrides(
    "terminal-au-test",
    "abababababababababababababababababababababababababababababababab",
    tables: ("misc": (410, 360, 500, 370)),
  ),
  dangles: 0,
)
#let custom-energy = data(evaluate("AAAAU", "(...)", model: custom)).total_kcal_mol
#let ordinary-energy = data(evaluate("AAAAU", "(...)", model: analysis-model(dangles: 0))).total_kcal_mol
#assert(calc.abs(custom-energy - ordinary-energy - 9.0) < 0.000000001)
#assert(data(conditional-density2("AAAAU", ".....", model: custom)).model.contains("custom"))

#let hard = folding-constraints(
  force-unpaired: (4, 5, 6),
  force-pairs: (constraint-pair(1, 9),),
)
#let constrained = data(analyze(sequence, model: model, constraints: hard))
#assert(constrained.constraints.hard_constraints)
#assert(data(validate(sequence, constrained.mfe_structure)).pairs.any(pair => pair.i == 1 and pair.j == 9))

#let probing = folding-constraints(probing: probing-data(
  (0.0, 0.1, 0.2, 0.4, 0.8, 1.2, 0.7, 0.3, none),
  method: "zarringhalam",
))
#let probed = analyze(sequence, model: model, constraints: probing)
#assert.eq(data(probed).constraints.probing_method, "zarringhalam")
#assert.eq(data(probed).ensemble.positional_entropy_bits.len(), 9)

#let sampled = data(sample(sequence, count: 12, seed: 42, model: model))
#assert.eq(sampled.returned, 12)
#assert.eq(data(sample(sequence, count: 1, seed: 42, model: model)).samples.first().structure, sampled.samples.first().structure)
#assert(data(suboptimal(sequence, energy-band: 5.0, limit: 12, model: model)).structures.len() >= 1)
#let conditional-sequence = "GCGCGCGCGCGC"
#let conditional-seed = "(......)...."
#let conditional-model = analysis-model(dangles: 0)
#let conditional = data(conditional-density2(
  conditional-sequence,
  conditional-seed,
  model: conditional-model,
  constraints: folding-constraints(force-unpaired: (2,)),
))
#let conditional-samples = data(conditional-density2-sample(
  conditional-sequence,
  conditional-seed,
  count: 8,
  seed: 17,
  model: conditional-model,
))
#assert.eq(conditional-samples.returned, 8)
#assert.eq(conditional-samples.samples.first().topology.genus, data(fatgraph-topology(
  conditional-sequence,
  conditional-samples.samples.first().structure,
)).genus)
#assert.eq(topology-annotations(fatgraph-topology("GCGCGCG", "(.[.).]")).len(), 2)
#assert(data(conditional-density2-suboptimal(
  conditional-sequence,
  conditional-seed,
  energy-band: 20.0,
  limit: 8,
  model: conditional-model,
)).structures.len() >= 1)
#let opening = accessibility(sequence, (accessibility-window(1, 1), accessibility-window(4, 6)), model: model)
#assert.eq(data(opening).windows.len(), 2)
#assert(data(duplex("GGG", "CCC", model: model)).standard_state_bound_probability > 0.0)
#assert(ensemble-defect(sequence, "(((...)))", result).normalized-ensemble-defect >= 0.0)

#let hidden-track = value-annotations((0.0, 1.0), legend: false)
#assert.eq(hidden-track.kind, "annotation-track")
#assert.eq(hidden-track.legends.len(), 0)

#render(result, width: 6cm, annotations: entropy-annotations(result))
#render(probed, width: 6cm, annotations: reactivity-annotations(data(probed)))
#draw(
  "GCGC",
  structure: "....",
  width: 6cm,
  annotations: (
    value-annotations((0.0, 0.4), positions: (1, 2), scale: color-scale(label: [First])),
    value-annotations((0.6, 1.0), positions: (3, 4), scale: color-scale(label: [Second])),
  ),
  legend: legend-style(position: "right", direction: "column", stroke: none),
)
#draw(
  "AAAA",
  structure: "....",
  width: 5cm,
  numbering: none,
  annotations: (
    base-annotation(
      1,
      fill: rgb("#4575b4"),
      text-contrast: "aaa",
      contrast-on-failure: "best",
    ),
    base-annotation(
      2,
      fill: rgb("#313695").transparentize(20%),
      contrast-background: white,
    ),
    // Explicit text-fill has higher priority than even a conflicting strict
    // per-node policy, so this intentionally low-contrast pair must compile.
    base-annotation(3, fill: black, text-fill: black, text-contrast: "aaa"),
  ),
)
#draw(
  "AAAA",
  structure: "....",
  width: 5cm,
  numbering: none,
  node-text-contrast: "aaa",
  annotations: (base-annotation(1, fill: rgb("#313695")),),
)
#draw(
  "AAAA",
  structure: "....",
  width: 5cm,
  numbering: none,
  theme: default-theme + (text-fill: purple, node-text-contrast: "fixed"),
  annotations: (base-annotation(1, fill: gradient.linear(red, blue)),),
)
#dot-plot(sequence, probabilities: result, width: 5cm)
#mountain-plot(sequence, probabilities: result, width: 6cm)
