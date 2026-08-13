#import "../../package/lib.typ": *

#set page(width: 18cm, height: auto, margin: 8mm)
#set text(size: 7pt)

#let model = analysis-model()
#let manifest = data(parameters()).active
#assert.eq(manifest.model_id, "ribon-turner-2004")
#assert.eq(manifest.parameter_file_count, 34)

#let complex = cofold(
  "GGGAAACCC",
  "GGGAAACCC",
  concentration-a: 0.000001,
  concentration-b: 0.000002,
  model: model,
)
#assert(calc.abs(data(complex).equilibrium.mass_balance_a_error_molar) < 0.000000000001)

#let local-result = local("GGGAAACCC", window-size: 7, max-pair-span: 6, max-unpaired: 2, model: model)
#assert.eq(data(local-result).windows.len(), 3)
#assert(data(local-result).pair_probabilities.len() > 0)

#let circ = circular("GGGAAACCCU", model: analysis-model(dangles: 0))
#assert.eq(data(circ).unpaired_probabilities.len(), 10)
#assert(data(circ).partition_function > 0.0)

#let modified-result = modified(
  "GGGAAACCC",
  (modified-base(4, "m6A", "A", unpaired-energy: 0.2),),
  model: model,
)
#assert.eq(data(modified-result).display_symbols.at(3), "m6A")

#let gq = gquad("GGAGGAGGAGG", model: model)
#assert(data(gq).candidates.len() >= 1)
#assert(data(gq).gquad_probability >= 0.0 and data(gq).gquad_probability <= 1.0)

#let pk = pseudoknot("GGGAAACCCGGGAAACCC", model: model)
#assert.eq(data(validate(data(pk).sequence, data(pk).structure)).pairs.len(), data(pk).pair_count)
#assert(data(pk).restricted_partition_function > 0.0)
#assert(data(pk).restricted_state_count > 0)
#assert(type(data(pk).restricted_state_count_exact) == str)
#assert(data(pk).hybrid_structure.len() == data(pk).sequence.len())
#let exact-pk = pseudoknot(
  "GGGAAACCC",
  model: model,
  options: pseudoknot-options(exact-arbitrary-ensemble: true),
)
#assert(data(exact-pk).exact_arbitrary_ensemble.state_space_complete)
#let pk-energy = evaluate-pseudoknot(
  "GGGAAACCCGGGAAACCC",
  "([......)......]..",
  model: model,
)
#assert(data(pk-energy).crossing_count == 1)

#let comparison = comparative(("GGGAAACCC", "GGGAAACCC", "CCCAAAGGG"), model: model)
#assert.eq(data(comparison).sequence_count, 3)
#assert.eq(data(comparison).alignment_length, 9)
#let path = landscape("GGGAAACCC", ".........", "(((...)))", model: model)
#assert(data(path).state_space_complete)
#assert(data(path).saddle_energy_kcal_mol >= data(path).path.first().energy_kcal_mol)
#let design = inverse-design(
  "(...)",
  template: "NNNNN",
  return-count: 2,
  model: analysis-model(dangles: 0),
)
#assert(data(design).search_complete)
#let ligand-result = ligand(
  "GGGAAACCC",
  (ligand-motif("aptamer", 1, "GGGAAACCC", "(((...)))", -3.0),),
  model: analysis-model(dangles: 0),
)
#assert(data(ligand-result).state_space_complete)

#grid(
  columns: (1fr, 1fr),
  gutter: 6pt,
  figure(render(complex, width: 7cm, height: 4cm, numbering: none), caption: [Cofold]),
  figure(render(circ, width: 7cm, height: 4cm, numbering: none), caption: [Circular RNA]),
  figure(render(modified-result, width: 7cm, height: 4cm, numbering: none), caption: [Modified base]),
  figure(render(gq, width: 7cm, height: 4cm, numbering: none), caption: [G-quadruplex]),
  figure(render(pk, method: "circular", width: 7cm, height: 4cm, numbering: none), caption: [Pseudoknot]),
  figure(
    draw(
      data(local-result).sequence,
      structure: "(((...)))",
      width: 7cm,
      height: 4cm,
      numbering: none,
      annotations: local-accessibility-annotations(local-result),
    ),
    caption: [Local accessibility],
  ),
  figure(render(path, which: "saddle", width: 7cm, height: 4cm, numbering: none), caption: [Exact energy landscape]),
  figure(render(design, which: "target", width: 7cm, height: 4cm, numbering: none), caption: [Exact inverse design]),
  figure(render(ligand-result, which: "mea", width: 7cm, height: 4cm, numbering: none), caption: [Ligand microstate ensemble]),
)
