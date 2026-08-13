#import "../../package/lib.typ": analysis-model, data, pseudoknot, pseudoknot-options

#let corpus = json("../data/pseudoknot_real_24.json")
#let model = analysis-model(min-loop: 3, dangles: 2)
#let results = corpus.cases.map(case => {
  let response = pseudoknot(case.sequence, model: model)
  let arbitrary-sequence = case.sequence.slice(0, calc.min(10, case.sequence.len()))
  let arbitrary = pseudoknot(
    arbitrary-sequence,
    model: model,
    options: pseudoknot-options(exact-arbitrary-ensemble: true),
  )
  (
    id: case.id,
    sequence: case.sequence,
    expected-structure: case.structure,
    prediction: data(response),
    arbitrary-sequence: arbitrary-sequence,
    arbitrary: data(arbitrary).exact_arbitrary_ensemble,
  )
})

#metadata(results) <ribon-pseudoknot-validation>
