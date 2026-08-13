#import "../../package/lib.typ": analysis-model, conditional-density2, data, dna-model

#let corpus = json("../data/pseudoknot_real_24.json")
#let seed-layer(structure) = structure.clusters().map(symbol => {
  if symbol == "(" or symbol == ")" { symbol } else { "." }
}).join()

#let evaluate-case(case, family, dangles) = {
  let sequence = if family == "dna" { case.sequence.replace("U", "T") } else { case.sequence }
  // CParty 1.0's documented DNA mode rejects forced A-T pairs while still
  // accepting DNA sequences. Use an empty fixed layer for the DNA oracle so
  // the comparison exercises its Matthews tables without triggering that CLI
  // validation defect.
  let seed = if family == "dna" {
    range(sequence.len()).map(_ => ".").join()
  } else {
    seed-layer(case.structure)
  }
  let result = data(conditional-density2(
    sequence,
    seed,
    model: if family == "dna" {
      dna-model(min-loop: 3, dangles: dangles)
    } else {
      analysis-model(min-loop: 3, dangles: dangles)
    },
  ))
  (
    id: case.id,
    family: family,
    dangles: dangles,
    sequence: sequence,
    seed: seed,
    mfe-structure: result.mfe_structure,
    mfe-energy: result.mfe_energy_kcal_mol,
    ensemble-energy: result.ensemble_free_energy_kcal_mol,
    mea-structure: result.mea_structure,
    centroid-structure: result.centroid_structure,
  )
}

// All 24 published pseudoknot cases under both CParty-supported dangle
// settings, plus eight sequence-diverse DNA homologues under both settings.
#let rna-values = corpus.cases.map(case => (
  evaluate-case(case, "rna", 0),
  evaluate-case(case, "rna", 2),
)).flatten()
#let dna-values = corpus.cases.slice(0, 8).map(case => (
  evaluate-case(case, "dna", 0),
  evaluate-case(case, "dna", 2),
)).flatten()
#let values = rna-values + dna-values

#metadata(values) <ribon-cparty-blackbox>
