#import "../../package/lib.typ": *

#let corpus = json("../data/rfam_real_24.json")
#let values = corpus.cases.map(case => {
  let sequence = case.sequence.slice(0, 10)
  let open = range(sequence.len()).map(_ => ".").join()
  let model = analysis-model(dangles: 0)
  let target = data(fold(sequence, model: model)).structure
  let path = data(landscape(sequence, open, target, model: model))
  let design = data(inverse-design(
    open,
    template: sequence,
    return-count: 1,
    model: model,
  ))
  let bound = data(ligand(
    sequence,
    (ligand-motif(
      case.accession + "-domain",
      1,
      sequence,
      open,
      -1.0,
      concentration: 0.001,
    ),),
    model: model,
  ))
  (
    accession: case.accession,
    sequence: sequence,
    landscape: (
      state-space-complete: path.state_space_complete,
      state-count: path.state_count,
      edge-count: path.edge_count,
      saddle-energy: path.saddle_energy_kcal_mol,
      path: path.path,
    ),
    inverse-design: (
      search-complete: design.search_complete,
      candidate-sequence-count: design.candidate_sequence_count,
      evaluated-sequence-count: design.evaluated_sequence_count,
      sequence: design.candidates.first().sequence,
      target-probability: design.candidates.first().target_probability,
    ),
    ligand: (
      state-space-complete: bound.state_space_complete,
      occupancy: bound.motifs.first().occupancy_probability,
      pair-probabilities: bound.pair_probabilities,
      unpaired-probabilities: bound.unpaired_probabilities,
    ),
  )
})

#metadata(values) <ribon-exact-features-real>
