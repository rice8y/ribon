#import "../../package/lib.typ": *

#set page(width: 18cm, height: auto, margin: 8mm)
#set text(font: "Libertinus Serif", size: 8pt)

#let sequence = "GGGAAACCC"
#let structure = "(((...)))"
#let model = analysis-model(dangles: 2)
#let circular-model = analysis-model(dangles: 3)

#let analysis = analyze(sequence, model: model)
#let mfe = fold(sequence, model: model)
#let energy = evaluate(sequence, structure, model: model)
#let parsed = validate(sequence, structure)
#let geometry = layout(sequence, structure, method: "naview")
#let samples = sample(sequence, count: 4, seed: 42, model: model)
#let opening = accessibility(
  sequence,
  (accessibility-window(1, 1), accessibility-window(4, 6)),
  model: model,
)
#let alternatives = suboptimal(sequence, energy-band: 4.0, limit: 8, model: model)
#let duplex-result = duplex("GGGG", "CCCC", model: model)
#let cofold-result = cofold(
  "GGGG",
  "CCCC",
  concentration-a: 1e-6,
  concentration-b: 2e-6,
  model: model,
)
#let local-result = local(
  "GGGAAACCCGGGAAACCC",
  window-size: 12,
  max-pair-span: 8,
  max-unpaired: 3,
  model: model,
)
#let circular-result = circular(sequence, model: circular-model)
#let modified-result = modified(
  "ACGAAACGU",
  (modified-base(1, "m6A", "A", kind: "m6a"),),
  model: model,
)
#let gquad-result = gquad("GGAGGAGGAGG", model: model)
#let pseudoknot-result = pseudoknot("GGGAAACCCGGGAAACCC", model: model)
#let exact-arbitrary-pseudoknot = pseudoknot(
  "GGGAAACCC",
  model: model,
  options: pseudoknot-options(exact-arbitrary-ensemble: true),
)
#let evaluated-pseudoknot = evaluate-pseudoknot(
  "GGGAAACCCGGGAAACCC",
  "([......)......]..",
  model: model,
)
#let conditional-result = conditional-density2(
  "GCGAAACGCU",
  "(........)",
  model: analysis-model(min-loop: 0, dangles: 0),
)
#let evaluated-conditional = evaluate-conditional-density2(
  "GAGAACAACU",
  "(....)....",
  "..(.....).",
  model: analysis-model(min-loop: 0, dangles: 0),
)
#let comparative-result = comparative(
  ("GGGAAACCC", "GGGAAACCC", "GAGAAACUC"),
  model: model,
)
#let landscape-result = landscape(sequence, ".........", structure, model: model)
#let design-result = inverse-design(
  "(...)",
  template: "NNNNN",
  return-count: 3,
  model: analysis-model(dangles: 0),
)
#let ligand-result = ligand(
  sequence,
  (ligand-motif("aptamer", 1, sequence, structure, -3.0),),
  model: analysis-model(dangles: 0),
)
#let parameter-data = data(parameters())
#let capability-data = data(capabilities())
#let invalid = try-request("analyze", ("sequence": "AX!"))
#let advertised-operations = capability-data.operations + capability-data.advanced_operations

#assert.eq(analysis.result.kind, "analysis")
#assert.eq(mfe.result.kind, "mfe")
#assert.eq(energy.result.kind, "energy")
#assert.eq(data(parsed).pairs.len(), 3)
#assert.eq(data(geometry).points.len(), 9)
#assert.eq(data(samples).samples.len(), 4)
#assert.eq(data(opening).windows.len(), 2)
#assert.eq(data(alternatives).structures.len() > 0, true)
#assert.eq(duplex-result.result.kind, "duplex")
#assert.eq(cofold-result.result.kind, "cofold")
#assert.eq(local-result.result.kind, "local")
#assert.eq(circular-result.result.kind, "circular")
#assert.eq(modified-result.result.kind, "modified")
#assert.eq(data(modified-result).parameter_usage.len(), 1)
#assert.eq(data(modified-result).parameter_usage.first().matched_stacks > 0, true)
#assert.eq(gquad-result.result.kind, "gquad")
#assert.eq(pseudoknot-result.result.kind, "pseudoknot")
#assert.eq(evaluated-pseudoknot.result.kind, "pseudoknot-energy")
#assert.eq(conditional-result.result.kind, "conditional-density2")
#assert.eq(evaluated-conditional.result.kind, "conditional-density2-energy")
#assert.eq(comparative-result.result.kind, "comparative")
#assert.eq(landscape-result.result.kind, "landscape")
#assert(data(landscape-result).state_space_complete)
#assert.eq(data(landscape-result).path.first().structure, ".........")
#assert.eq(data(landscape-result).path.last().structure, structure)
#assert.eq(design-result.result.kind, "inverse-design")
#assert(data(design-result).search_complete)
#assert.eq(data(design-result).candidate_sequence_count, 384)
#assert.eq(ligand-result.result.kind, "ligand")
#assert(data(ligand-result).state_space_complete)
#assert(data(ligand-result).motifs.first().occupancy_probability > 0.0)
#assert.eq(parameter-data.active.model_id, "ribon-turner-2004")
#assert.eq(capability-data.operations.len(), 19)
#assert.eq(capability-data.advanced_operations.len(), 8)
#assert.eq(advertised-operations.len(), 27)
#assert(advertised-operations.all(operation => operation not in (
  "conditional-density2-oracle",
  "evaluate-conditional-density2-oracle",
)))
#assert.eq(data(pseudoknot-result).restricted_state_count > 0, true)
#assert.eq(type(data(pseudoknot-result).restricted_state_count_exact), str)
#assert.eq(data(pseudoknot-result).hybrid_structure.len(), data(pseudoknot-result).sequence.len())
#assert(data(exact-arbitrary-pseudoknot).exact_arbitrary_ensemble.state_space_complete)
#assert(data(exact-arbitrary-pseudoknot).exact_arbitrary_ensemble.state_count > 1)
#assert.eq(data(evaluated-pseudoknot).crossing_count, 1)
#assert(data(conditional-result).state_space_complete)
#assert.eq(data(conditional-result).pair_probabilities.find(
  pair => pair.i == 1 and pair.j == 10,
).probability, 1.0)
#assert.eq(invalid.ok, false)
#assert.eq(invalid.error.code, "analysis_failed")

= Ribon analysis/1 operation matrix

#table(
  columns: (1fr, 1fr, 1fr),
  inset: 4pt,
  [*MFE*], [*Ensemble ΔG*], [*Pair entries*],
  [#data(analysis).mfe_energy_kcal_mol kcal/mol],
  [#data(analysis).ensemble_free_energy_kcal_mol kcal/mol],
  [#data(analysis).pair_probabilities.len()],
)

#grid(
  columns: (1fr, 1fr),
  gutter: 6mm,
  render(analysis, which: "mfe", method: "naview", width: 7.5cm),
  render(analysis, which: "mea", method: "simple", width: 7.5cm),
  render(circular-result, width: 7.5cm),
  render(modified-result, method: "turtle", width: 7.5cm),
  render(gquad-result, method: "linear", width: 7.5cm, height: 3.8cm),
  render(pseudoknot-result, method: "circular", width: 7.5cm),
  render(exact-arbitrary-pseudoknot, which: "arbitrary-mfe", method: "circular", width: 7.5cm),
  render(evaluated-pseudoknot, method: "circular", width: 7.5cm),
  render(conditional-result, which: "centroid", method: "circular", width: 7.5cm),
  render(evaluated-conditional, method: "circular", width: 7.5cm),
  render(landscape-result, which: "saddle", method: "circular", width: 7.5cm),
  render(design-result, which: "target", method: "circular", width: 7.5cm),
  render(ligand-result, which: "mea", method: "circular", width: 7.5cm),
)

#grid(
  columns: (1fr, 1fr),
  gutter: 6mm,
  dot-plot(sequence, probabilities: analysis, width: 7.5cm),
  mountain-plot(sequence, probabilities: analysis, width: 7.5cm),
)
