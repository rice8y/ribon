#import "../../package/lib.typ": *

#let requested = int(sys.inputs.at("length", default: "120"))
#let corpus = json("../data/rfam_real_24.json")
#let source = corpus.cases.sorted(key: case => -case.length).first()
#assert(source.length >= requested)
#let sequence = source.sequence.slice(0, requested)
#let response = local(
  sequence,
  window-size: calc.min(50, requested),
  max-pair-span: calc.min(35, requested - 1),
  max-unpaired: 1,
)
#let result = data(response)

#metadata((
  accession: source.accession,
  length: requested,
  pair-count: result.pair_probabilities.len(),
  accessibility-count: result.accessibility.len(),
  window-count: result.windows.len(),
  method: result.method,
)) <ribon-performance>
