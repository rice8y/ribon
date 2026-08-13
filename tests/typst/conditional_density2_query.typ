#import "../../package/lib.typ": analysis-model, analyze, conditional-density2, data, evaluate-conditional-density2, execution-policy, request

// Exhaustive engines are intentionally absent from the publication API. They
// remain reachable here only as independent, resource-limited test oracles.
#let conditional-density2-oracle(sequence, structure, model: analysis-model()) = request(
  "conditional-density2-oracle",
  (sequence: sequence, structure: structure),
  model: model,
  execution: execution-policy(allow-expensive: true),
)
#let evaluate-conditional-density2-oracle(
  sequence,
  seed-structure,
  added-structure,
  model: analysis-model(),
) = request(
  "evaluate-conditional-density2-oracle",
  (
    sequence: sequence,
    seed_structure: seed-structure,
    added_structure: added-structure,
  ),
  model: model,
  execution: execution-policy(allow-expensive: true),
)

#let corpus = json("../data/pseudoknot_real_24.json")
#let evaluation-model = analysis-model(min-loop: 0, dangles: 0)
#let evaluation-model-d2 = analysis-model(min-loop: 0, dangles: 2)
#let planar-model = analysis-model(min-loop: 3, dangles: 0)
#let planar-model-d2 = analysis-model(min-loop: 3, dangles: 2)

#let layer(structure, opening, closing) = structure.clusters().map(symbol => {
  if symbol == opening or symbol == closing { symbol } else { "." }
}).join()

#let results = corpus.cases.map(case => {
  let seed = layer(case.structure, "(", ")")
  let added = layer(case.structure, "[", "]")
  let evaluation = data(evaluate-conditional-density2(
    case.sequence,
    seed,
    added,
    model: evaluation-model,
  ))
  let oracle-evaluation = data(evaluate-conditional-density2-oracle(
    case.sequence,
    seed,
    added,
    model: evaluation-model,
  ))
  let ensemble = data(conditional-density2(
    case.sequence,
    seed,
    model: evaluation-model,
  ))
  let mfe-evaluation = data(evaluate-conditional-density2(
    case.sequence,
    seed,
    ensemble.mfe_added_structure,
    model: evaluation-model,
  ))
  let evaluation-d2 = data(evaluate-conditional-density2(
    case.sequence,
    seed,
    added,
    model: evaluation-model-d2,
  ))
  let ensemble-d2 = data(conditional-density2(
    case.sequence,
    seed,
    model: evaluation-model-d2,
  ))
  let mfe-evaluation-d2 = data(evaluate-conditional-density2(
    case.sequence,
    seed,
    ensemble-d2.mfe_added_structure,
    model: evaluation-model-d2,
  ))
  let prefix = case.sequence.slice(0, 10)
  let conditional = data(conditional-density2(
    prefix,
    "..........",
    model: planar-model,
  ))
  let conditional-oracle = data(conditional-density2-oracle(
    prefix,
    "..........",
    model: planar-model,
  ))
  let planar = data(analyze(prefix, model: planar-model))
  let conditional-d2 = data(conditional-density2(
    prefix,
    "..........",
    model: planar-model-d2,
  ))
  let planar-d2 = data(analyze(prefix, model: planar-model-d2))
  (
    id: case.id,
    sequence: case.sequence,
    expected-structure: case.structure,
    evaluation: evaluation,
    oracle-evaluation: oracle-evaluation,
    ensemble: ensemble,
    mfe-evaluation: mfe-evaluation,
    evaluation-d2: evaluation-d2,
    ensemble-d2: ensemble-d2,
    mfe-evaluation-d2: mfe-evaluation-d2,
    prefix: prefix,
    conditional: conditional,
    conditional-oracle: conditional-oracle,
    planar: planar,
    conditional-d2: conditional-d2,
    planar-d2: planar-d2,
  )
})

#metadata(results) <ribon-conditional-density2-validation>
