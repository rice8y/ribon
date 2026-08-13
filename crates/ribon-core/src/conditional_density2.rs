//! Conditional density-2 pseudoknot ensemble.
//!
//! A pseudoknot-free seed structure `G` is held fixed while a second,
//! pseudoknot-free base-pair set `G'` is enumerated.  The two layers may not
//! share a nucleotide and only unions whose band density is at most two are
//! retained.  This module is deliberately an exhaustive reference engine: it
//! defines the finite state space independently of the polynomial CParty
//! decomposition and is used as its correctness oracle.

use crate::constraints::{ConstraintConfig, ConstraintModel, ConstraintSummary};
use crate::energy::{EnergyBreakdown, EnergyModel, NucleicAcid, ThermodynamicParameterOverrides};
use crate::exact_enumeration::for_each_noncrossing_structure;
use crate::partition::PairProbability;
use crate::structure::{parse_structure, RnaError};
use crate::topology::fatgraph_from_pairs;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

mod polynomial;

pub use polynomial::{
    conditional_density2_polynomial, conditional_density2_polynomial_with_constraints,
    evaluate_conditional_density2_polynomial, sample_conditional_density2_polynomial,
    sample_conditional_density2_polynomial_with_constraints,
    suboptimal_conditional_density2_polynomial,
    suboptimal_conditional_density2_polynomial_with_constraints,
    ConditionalDensity2PolynomialEvaluationResult, ConditionalDensity2PolynomialResult,
    ConditionalDensity2Sample, ConditionalDensity2SamplingResult,
    ConditionalDensity2SuboptimalResult, ConditionalDensity2SuboptimalStructure,
};

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;

/// Independently stated DP09/CParty density-2 parameters at 37 C.
///
/// The Turner contributions are temperature-scaled by [`EnergyModel`].  The
/// empirical pseudoknot constants are kept explicit because the published
/// DP09 parameterization does not provide enthalpies for temperature scaling.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct ConditionalDensity2Options {
    pub gamma: f64,
    /// Thermodynamic table family for the planar loop terms. The API layer
    /// derives this from the selected model id.
    pub nucleic_acid: NucleicAcid,
    pub parameter_overrides: Option<ThermodynamicParameterOverrides>,
    /// Restrict the added layer to pairs that cross at least one seed pair.
    /// This is the `pk-only` state space used by conditional pseudoknot tools.
    pub pk_only: bool,
    /// DP09 affine multiloop initiation parameter `a`.
    pub multiloop_init_kcal_mol: f64,
    /// DP09 affine multiloop branch parameter `b`.
    pub multiloop_branch_kcal_mol: f64,
    /// DP09 affine multiloop unpaired-base parameter `c`.
    pub multiloop_unpaired_kcal_mol: f64,
    pub pseudoloop_initiation_kcal_mol: f64,
    pub multiloop_pseudoknot_kcal_mol: f64,
    pub nested_pseudoknot_kcal_mol: f64,
    pub band_kcal_mol: f64,
    pub pseudoloop_unpaired_kcal_mol: f64,
    pub closed_subregion_kcal_mol: f64,
    pub spanning_stack_factor: f64,
    pub spanning_internal_factor: f64,
    pub spanning_multiloop_init_kcal_mol: f64,
    pub spanning_multiloop_branch_kcal_mol: f64,
    pub spanning_multiloop_unpaired_kcal_mol: f64,
}

impl Default for ConditionalDensity2Options {
    fn default() -> Self {
        Self {
            gamma: 1.0,
            nucleic_acid: NucleicAcid::Rna,
            parameter_overrides: None,
            pk_only: false,
            multiloop_init_kcal_mol: 3.39,
            multiloop_branch_kcal_mol: 0.03,
            multiloop_unpaired_kcal_mol: 0.02,
            pseudoloop_initiation_kcal_mol: -1.38,
            multiloop_pseudoknot_kcal_mol: 10.07,
            nested_pseudoknot_kcal_mol: 15.00,
            band_kcal_mol: 2.46,
            pseudoloop_unpaired_kcal_mol: 0.06,
            closed_subregion_kcal_mol: 0.96,
            spanning_stack_factor: 0.89,
            spanning_internal_factor: 0.74,
            spanning_multiloop_init_kcal_mol: 3.41,
            spanning_multiloop_branch_kcal_mol: 0.56,
            spanning_multiloop_unpaired_kcal_mol: 0.12,
        }
    }
}

fn conditional_energy_model(
    temperature_celsius: f64,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
) -> Result<EnergyModel, RnaError> {
    if let Some(overrides) = &options.parameter_overrides {
        EnergyModel::with_parameter_overrides(
            temperature_celsius,
            dangles,
            salt_molar,
            options.nucleic_acid,
            overrides.clone(),
        )
    } else {
        EnergyModel::with_parameter_family(
            temperature_celsius,
            dangles,
            salt_molar,
            options.nucleic_acid,
        )
    }
}

pub(super) fn conditional_parameter_model_name(model: &EnergyModel) -> &'static str {
    if model.parameter_profile_name().is_some() {
        "Ribon conditional density-2 DP09 model with custom normalized loop parameters"
    } else {
        match model.nucleic_acid() {
            NucleicAcid::Rna => {
                "Ribon conditional density-2 DP09 model with RNAstructure 6.6 RNA loops"
            }
            NucleicAcid::Dna => {
                "Ribon conditional density-2 DP09 model with RNAstructure 6.6 DNA loops"
            }
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2Energy {
    /// Ordinary Turner energy of the fixed planar layer.
    pub seed_layer_kcal_mol: f64,
    /// Ordinary Turner energy of the variable planar layer.
    pub added_layer_kcal_mol: f64,
    pub spanning_stack_correction_kcal_mol: f64,
    pub spanning_internal_correction_kcal_mol: f64,
    pub spanning_multiloop_correction_kcal_mol: f64,
    pub pseudoloop_initiation_kcal_mol: f64,
    pub multiloop_pseudoknot_kcal_mol: f64,
    pub nested_pseudoknot_kcal_mol: f64,
    pub band_kcal_mol: f64,
    pub pseudoloop_unpaired_kcal_mol: f64,
    pub closed_subregion_kcal_mol: f64,
    /// Hard/soft-constraint pseudo-energy applied once to the union of layers.
    pub constraint_kcal_mol: f64,
    /// Difference between the independently decomposed diagnostic terms and
    /// the authoritative unambiguous interval-hypergraph energy.
    pub decomposition_alignment_kcal_mol: f64,
    pub total_kcal_mol: f64,
    pub crossing_component_count: usize,
    pub exterior_pseudoloop_count: usize,
    pub multiloop_pseudoknot_count: usize,
    pub nested_pseudoknot_count: usize,
    pub band_count: usize,
    pub pseudoknot_band_count: usize,
    pub pseudoloop_unpaired_count: usize,
    pub closed_subregion_count: usize,
    pub maximum_band_density: usize,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2Result {
    pub sequence: String,
    pub seed_structure: String,
    pub mfe_structure: String,
    pub mfe_added_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub mfe_energy: ConditionalDensity2Energy,
    pub ensemble_free_energy_kcal_mol: f64,
    pub partition_function: f64,
    pub log_partition_function: f64,
    /// Marginals for both layers.  Every seed pair has probability one.
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub centroid_structure: String,
    pub centroid_distance: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub constraints: ConstraintSummary,
    pub state_count: usize,
    pub state_count_exact: String,
    pub state_space_complete: bool,
    pub maximum_density: usize,
    pub time_complexity: &'static str,
    pub space_complexity: &'static str,
    pub algorithm: &'static str,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2EvaluationResult {
    pub sequence: String,
    pub seed_structure: String,
    pub added_structure: String,
    pub structure: String,
    pub energy_kcal_mol: f64,
    pub derivation_unique: bool,
    pub energy: ConditionalDensity2Energy,
    pub constraints: ConstraintSummary,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
#[serde(untagged)]
pub enum ConditionalDensity2Analysis {
    Polynomial(Box<ConditionalDensity2PolynomialResult>),
    Exhaustive(Box<ConditionalDensity2Result>),
}

#[derive(Clone, Debug, Serialize)]
#[serde(untagged)]
pub enum ConditionalDensity2Evaluation {
    Polynomial(ConditionalDensity2PolynomialEvaluationResult),
    Exhaustive(Box<ConditionalDensity2EvaluationResult>),
}

#[derive(Clone)]
struct State {
    added: Vec<(usize, usize)>,
    energy: ConditionalDensity2Energy,
    log_weight: f64,
}

/// Evaluate one explicitly supplied `G ∪ G'` state under the same conditional
/// model used by [`conditional_density2_ensemble`].  Both inputs must be
/// pseudoknot-free and nucleotide-disjoint; their union must have density at
/// most two.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_conditional_density2_structure(
    sequence: &str,
    seed_structure: &str,
    added_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
) -> Result<ConditionalDensity2EvaluationResult, RnaError> {
    evaluate_conditional_density2_structure_with_constraints(
        sequence,
        seed_structure,
        added_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn evaluate_conditional_density2_structure_with_constraints(
    sequence: &str,
    seed_structure: &str,
    added_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2EvaluationResult, RnaError> {
    validate_options(temperature_celsius, options)?;
    let seed = parse_structure(sequence, seed_structure)?;
    let added = parse_structure(sequence, added_structure)?;
    if !seed.strand_breaks.is_empty() || !added.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis accepts one strand".into(),
        ));
    }
    let seed_pairs = seed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    let added_pairs = added
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    validate_planar_layer("seed", &seed, &seed_pairs, min_loop)?;
    validate_planar_layer("added", &added, &added_pairs, min_loop)?;
    let seed_occupied = seed_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .collect::<HashSet<_>>();
    if let Some(position) = added_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .find(|position| seed_occupied.contains(position))
    {
        return Err(RnaError::MultiplePartners {
            position: position + 1,
        });
    }
    let mut union = seed_pairs.clone();
    union.extend_from_slice(&added_pairs);
    union.sort_unstable();
    let topology = band_topology(seed.length, &union);
    if topology.maximum_density > 2 {
        return Err(RnaError::InvalidOption(format!(
            "conditional union has band density {}, expected at most 2",
            topology.maximum_density
        )));
    }
    let model = conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    let conditional_model = conditional_parameter_model_name(&model);
    let constraints = ConstraintModel::compile(seed.length, constraint_config)?;
    let union_partner = {
        let mut partner = vec![None; seed.length];
        for &(i, j) in &union {
            partner[i] = Some(j);
            partner[j] = Some(i);
        }
        partner
    };
    constraints.validate_structure(seed.sequence.as_bytes(), &union_partner, min_loop, &model)?;
    let normalized_seed = planar_structure(seed.length, &seed_pairs, '(', ')');
    let normalized_added = planar_structure(seed.length, &added_pairs, '(', ')');
    let seed_breakdown = model.evaluate(&seed.sequence, &normalized_seed)?;
    let added_breakdown = model.evaluate(&seed.sequence, &normalized_added)?;
    let energy = conditional_energy(
        seed.length,
        &seed_pairs,
        &added_pairs,
        &seed_breakdown,
        &added_breakdown,
        &topology,
        options,
        constraints.structure_energy(&union_partner),
        conditional_model,
    );
    let energy_kcal_mol = energy.total_kcal_mol;
    Ok(ConditionalDensity2EvaluationResult {
        sequence: seed.sequence,
        seed_structure: normalized_seed,
        added_structure: normalized_added,
        structure: layered_structure(seed.length, &seed_pairs, &added_pairs),
        energy_kcal_mol,
        derivation_unique: true,
        energy,
        constraints: constraints.summary(),
        model: conditional_model,
    })
}

/// Select the polynomial production engine whenever its local factorization
/// applies, otherwise retain the complete exact constrained state space.
#[allow(clippy::too_many_arguments)]
pub fn conditional_density2(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2Analysis, RnaError> {
    let length = sequence.chars().filter(|&symbol| symbol != '&').count();
    ConstraintModel::compile(length, constraint_config)?;
    if dangles % 2 == 0 && polynomial_compatible_constraints(constraint_config) {
        Ok(ConditionalDensity2Analysis::Polynomial(Box::new(
            conditional_density2_polynomial_with_constraints(
                sequence,
                seed_structure,
                temperature_celsius,
                min_loop,
                dangles,
                salt_molar,
                options,
                constraint_config,
            )?,
        )))
    } else {
        Ok(ConditionalDensity2Analysis::Exhaustive(Box::new(
            conditional_density2_ensemble_with_constraints(
                sequence,
                seed_structure,
                temperature_celsius,
                min_loop,
                dangles,
                salt_molar,
                options,
                constraint_config,
            )?,
        )))
    }
}

#[allow(clippy::too_many_arguments)]
pub fn evaluate_conditional_density2(
    sequence: &str,
    seed_structure: &str,
    added_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2Evaluation, RnaError> {
    let length = sequence.chars().filter(|&symbol| symbol != '&').count();
    let constraints = ConstraintModel::compile(length, constraint_config)?;
    if dangles % 2 == 0 && !constraints.summary().enabled {
        Ok(ConditionalDensity2Evaluation::Polynomial(
            evaluate_conditional_density2_polynomial(
                sequence,
                seed_structure,
                added_structure,
                temperature_celsius,
                min_loop,
                dangles,
                salt_molar,
                options,
            )?,
        ))
    } else {
        Ok(ConditionalDensity2Evaluation::Exhaustive(Box::new(
            evaluate_conditional_density2_structure_with_constraints(
                sequence,
                seed_structure,
                added_structure,
                temperature_celsius,
                min_loop,
                dangles,
                salt_molar,
                options,
                constraint_config,
            )?,
        )))
    }
}

/// Compute the exact conditional ensemble by enumerating every admissible
/// second planar layer.  There is no implicit sequence-length, state-count, or
/// energy beam limit.
#[allow(clippy::too_many_arguments)]
pub fn conditional_density2_ensemble(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
) -> Result<ConditionalDensity2Result, RnaError> {
    conditional_density2_ensemble_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_density2_ensemble_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2Result, RnaError> {
    let mut enumerated = enumerate_conditional_states(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        constraint_config,
    )?;
    enumerated.states.sort_by(|left, right| {
        left.energy
            .total_kcal_mol
            .total_cmp(&right.energy.total_kcal_mol)
            .then_with(|| left.added.cmp(&right.added))
    });
    finish_ensemble(
        enumerated.sequence,
        enumerated.seed_planar,
        enumerated.seed_pairs,
        enumerated.states,
        enumerated.maximum_density,
        enumerated.rt,
        options.gamma,
        enumerated.constraints,
    )
}

struct EnumeratedConditionalStates {
    sequence: String,
    seed_planar: String,
    seed_pairs: Vec<(usize, usize)>,
    states: Vec<State>,
    maximum_density: usize,
    rt: f64,
    constraints: ConstraintSummary,
}

#[allow(clippy::too_many_arguments)]
fn enumerate_conditional_states(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<EnumeratedConditionalStates, RnaError> {
    validate_options(temperature_celsius, options)?;
    let seed = parse_structure(sequence, seed_structure)?;
    if !seed.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis accepts one strand".into(),
        ));
    }
    let seed_pairs = seed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    if !is_noncrossing(&seed_pairs) {
        return Err(RnaError::InvalidOption(
            "the conditional seed must be pseudoknot-free".into(),
        ));
    }
    if let Some(pair) = seed.pairs.iter().find(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed contains noncanonical pair {}-{}",
            pair.i, pair.j
        )));
    }
    if let Some(&(i, j)) = seed_pairs.iter().find(|&&(i, j)| j - i <= min_loop) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed pair {}-{} violates min-loop={min_loop}",
            i + 1,
            j + 1
        )));
    }

    let model = conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    let conditional_model = conditional_parameter_model_name(&model);
    let user_constraints = ConstraintModel::compile(seed.length, constraint_config)?;
    let constraint_summary = user_constraints.summary();
    // A nucleotide already paired in G cannot be reused by G'. Hard
    // requirements already fulfilled by G are removed from the added-layer
    // grammar; all soft terms are scored once on the final union below.
    let seed_occupied = seed_pairs
        .iter()
        .flat_map(|&(i, j)| [i + 1, j + 1])
        .collect::<HashSet<_>>();
    for &(i, j) in &seed_pairs {
        if !user_constraints.allows_pair(seed.sequence.as_bytes(), i, j, &model) {
            return Err(RnaError::InvalidOption(format!(
                "seed pair {}-{} violates the requested constraints",
                i + 1,
                j + 1
            )));
        }
    }
    let mut added_config = constraint_config.clone();
    added_config.soft = Default::default();
    added_config.probing = None;
    added_config
        .force_paired
        .retain(|position| !seed_occupied.contains(position));
    added_config
        .force_pairs
        .retain(|pair| !seed_pairs.contains(&(pair.i.min(pair.j) - 1, pair.i.max(pair.j) - 1)));
    if let Some(pair) = added_config
        .force_pairs
        .iter()
        .find(|pair| seed_occupied.contains(&pair.i) || seed_occupied.contains(&pair.j))
    {
        return Err(RnaError::InvalidOption(format!(
            "forced pair {}-{} reuses a nucleotide occupied by the seed",
            pair.i, pair.j
        )));
    }
    added_config
        .force_unpaired
        .extend(seed_occupied.iter().copied());
    added_config.force_unpaired.sort_unstable();
    added_config.force_unpaired.dedup();
    let constraints = ConstraintModel::compile(seed.length, &added_config)?;
    let seed_planar = planar_structure(seed.length, &seed_pairs, '(', ')');
    let seed_breakdown = model.evaluate(&seed.sequence, &seed_planar)?;
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut states = Vec::new();
    let mut maximum_density = 0usize;

    for_each_noncrossing_structure(
        seed.sequence.as_bytes(),
        min_loop,
        &model,
        &constraints,
        |added| {
            if options.pk_only
                && added
                    .iter()
                    .any(|&pair| !seed_pairs.iter().copied().any(|seed| crosses(pair, seed)))
            {
                return Ok(());
            }
            let mut union = seed_pairs.clone();
            union.extend_from_slice(added);
            union.sort_unstable();
            let topology = band_topology(seed.length, &union);
            maximum_density = maximum_density.max(topology.maximum_density);
            if topology.maximum_density > 2 {
                return Ok(());
            }
            let added_structure = planar_structure(seed.length, added, '(', ')');
            let added_breakdown = model.evaluate(&seed.sequence, &added_structure)?;
            let constraint_kcal_mol = {
                let mut union_partner = vec![None; seed.length];
                for &(i, j) in &union {
                    union_partner[i] = Some(j);
                    union_partner[j] = Some(i);
                }
                if user_constraints
                    .validate_structure(seed.sequence.as_bytes(), &union_partner, min_loop, &model)
                    .is_err()
                {
                    return Ok(());
                }
                user_constraints.structure_energy(&union_partner)
            };
            let mut energy = conditional_energy(
                seed.length,
                &seed_pairs,
                added,
                &seed_breakdown,
                &added_breakdown,
                &topology,
                options,
                constraint_kcal_mol,
                conditional_model,
            );
            if dangles % 2 == 0 && constraint_summary.enabled {
                let Ok(evaluated) = evaluate_conditional_density2_polynomial(
                    &seed.sequence,
                    &seed_planar,
                    &added_structure,
                    temperature_celsius,
                    min_loop,
                    dangles,
                    salt_molar,
                    options,
                ) else {
                    return Ok(());
                };
                let hypergraph_energy = evaluated.energy_kcal_mol + constraint_kcal_mol;
                energy.decomposition_alignment_kcal_mol = hypergraph_energy - energy.total_kcal_mol;
                energy.total_kcal_mol = hypergraph_energy;
            }
            states.push(State {
                added: added.to_vec(),
                log_weight: -energy.total_kcal_mol / rt,
                energy,
            });
            Ok(())
        },
    )?;
    if states.is_empty() {
        return Err(RnaError::Numerical(
            "conditional density-2 ensemble has no admissible state".into(),
        ));
    }
    Ok(EnumeratedConditionalStates {
        sequence: seed.sequence,
        seed_planar,
        seed_pairs,
        states,
        maximum_density,
        rt,
        constraints: constraint_summary,
    })
}

/// Public dispatcher for exact conditional Boltzmann sampling. Even dangle
/// models use stochastic hypergraph traceback; odd models use the complete
/// fixed-structure state space and therefore retain exact shared-dangle and
/// coaxial energies.
#[allow(clippy::too_many_arguments)]
pub fn sample_conditional_density2(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    count: usize,
    random_seed: u64,
    unique: bool,
) -> Result<ConditionalDensity2SamplingResult, RnaError> {
    sample_conditional_density2_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        count,
        random_seed,
        unique,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn sample_conditional_density2_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    count: usize,
    random_seed: u64,
    unique: bool,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2SamplingResult, RnaError> {
    ConstraintModel::compile(
        sequence.chars().filter(|&c| c != '&').count(),
        constraint_config,
    )?;
    if dangles % 2 == 0 && polynomial_compatible_constraints(constraint_config) {
        return sample_conditional_density2_polynomial_with_constraints(
            sequence,
            seed_structure,
            temperature_celsius,
            min_loop,
            dangles,
            salt_molar,
            options,
            count,
            random_seed,
            unique,
            constraint_config,
        );
    }
    if count == 0 {
        return Err(RnaError::InvalidOption(
            "conditional sample count must be positive".into(),
        ));
    }
    let enumerated = enumerate_conditional_states(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        constraint_config,
    )?;
    let log_z = enumerated
        .states
        .iter()
        .fold(f64::NEG_INFINITY, |total, state| {
            log_add(total, state.log_weight)
        });
    let maximum = enumerated
        .states
        .iter()
        .map(|state| state.log_weight)
        .fold(f64::NEG_INFINITY, f64::max);
    let weights = enumerated
        .states
        .iter()
        .map(|state| (state.log_weight - maximum).exp())
        .collect::<Vec<_>>();
    let total = weights.iter().sum::<f64>();
    let mut rng = ConditionalSplitMix64::new(random_seed);
    let mut seen = HashSet::new();
    let mut samples = Vec::with_capacity(count);
    for _ in 0..count {
        let mut threshold = rng.next_f64() * total;
        let mut selected = weights.len() - 1;
        for (index, &weight) in weights.iter().enumerate() {
            threshold -= weight;
            if threshold <= 0.0 {
                selected = index;
                break;
            }
        }
        let state = &enumerated.states[selected];
        let added_structure = planar_structure(enumerated.sequence.len(), &state.added, '(', ')');
        if unique && !seen.insert(added_structure.clone()) {
            continue;
        }
        let mut union = enumerated.seed_pairs.clone();
        union.extend_from_slice(&state.added);
        union.sort_unstable();
        let log_probability = state.log_weight - log_z;
        samples.push(ConditionalDensity2Sample {
            structure: layered_structure(
                enumerated.sequence.len(),
                &enumerated.seed_pairs,
                &state.added,
            ),
            added_structure,
            effective_energy_kcal_mol: state.energy.total_kcal_mol,
            probability: log_probability.exp(),
            log_probability,
            topology: fatgraph_from_pairs(&union),
        });
    }
    Ok(ConditionalDensity2SamplingResult {
        sequence: enumerated.sequence,
        seed_structure: enumerated.seed_planar,
        temperature_celsius,
        dangles,
        salt_molar,
        seed: random_seed,
        requested: count,
        returned: samples.len(),
        unique,
        log_partition_function: log_z,
        samples,
        constraints: enumerated.constraints,
        algorithm: "exact exhaustive odd-dangle conditional Boltzmann sampler",
    })
}

/// Public exact conditional k-best dispatcher.
#[allow(clippy::too_many_arguments)]
pub fn suboptimal_conditional_density2(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    energy_band_kcal_mol: f64,
    limit: usize,
) -> Result<ConditionalDensity2SuboptimalResult, RnaError> {
    suboptimal_conditional_density2_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        energy_band_kcal_mol,
        limit,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn suboptimal_conditional_density2_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    energy_band_kcal_mol: f64,
    limit: usize,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2SuboptimalResult, RnaError> {
    ConstraintModel::compile(
        sequence.chars().filter(|&c| c != '&').count(),
        constraint_config,
    )?;
    if dangles % 2 == 0 && polynomial_compatible_constraints(constraint_config) {
        return suboptimal_conditional_density2_polynomial_with_constraints(
            sequence,
            seed_structure,
            temperature_celsius,
            min_loop,
            dangles,
            salt_molar,
            options,
            energy_band_kcal_mol,
            limit,
            constraint_config,
        );
    }
    if !energy_band_kcal_mol.is_finite() || energy_band_kcal_mol < 0.0 {
        return Err(RnaError::InvalidOption(
            "conditional suboptimal energy band must be finite and non-negative".into(),
        ));
    }
    if limit == 0 {
        return Err(RnaError::InvalidOption(
            "conditional suboptimal limit must be positive".into(),
        ));
    }
    let mut enumerated = enumerate_conditional_states(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        constraint_config,
    )?;
    enumerated.states.sort_by(|left, right| {
        left.energy
            .total_kcal_mol
            .total_cmp(&right.energy.total_kcal_mol)
            .then_with(|| left.added.cmp(&right.added))
    });
    let mfe = enumerated.states[0].energy.total_kcal_mol;
    let in_band = enumerated
        .states
        .iter()
        .take_while(|state| state.energy.total_kcal_mol <= mfe + energy_band_kcal_mol + 1.0e-12)
        .count();
    let mut structures = Vec::with_capacity(in_band.min(limit));
    for (rank, state) in enumerated
        .states
        .iter()
        .take(in_band.min(limit))
        .enumerate()
    {
        let mut union = enumerated.seed_pairs.clone();
        union.extend_from_slice(&state.added);
        union.sort_unstable();
        let delta = state.energy.total_kcal_mol - mfe;
        structures.push(ConditionalDensity2SuboptimalStructure {
            rank: rank + 1,
            structure: layered_structure(
                enumerated.sequence.len(),
                &enumerated.seed_pairs,
                &state.added,
            ),
            added_structure: planar_structure(enumerated.sequence.len(), &state.added, '(', ')'),
            energy_kcal_mol: state.energy.total_kcal_mol,
            delta_energy_kcal_mol: delta,
            relative_boltzmann_weight: (-delta / enumerated.rt).exp(),
            topology: fatgraph_from_pairs(&union),
        });
    }
    Ok(ConditionalDensity2SuboptimalResult {
        sequence: enumerated.sequence,
        seed_structure: enumerated.seed_planar,
        temperature_celsius,
        dangles,
        salt_molar,
        energy_band_kcal_mol,
        requested_limit: limit,
        truncated: in_band > limit,
        structures,
        constraints: enumerated.constraints,
        algorithm: "exact exhaustive odd-dangle conditional energy ordering",
    })
}

struct ConditionalSplitMix64(u64);

impl ConditionalSplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D049BB133111EB);
        value ^= value >> 31;
        (value >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let high = a.max(b);
    high + (a.min(b) - high).exp().ln_1p()
}

fn validate_options(
    temperature_celsius: f64,
    options: &ConditionalDensity2Options,
) -> Result<(), RnaError> {
    if !temperature_celsius.is_finite() || temperature_celsius <= -273.15 {
        return Err(RnaError::InvalidOption(
            "temperature must be finite and above absolute zero".into(),
        ));
    }
    if !options.gamma.is_finite() || options.gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "conditional MEA gamma must be finite and positive".into(),
        ));
    }
    for (name, value) in [
        ("multiloop-init", options.multiloop_init_kcal_mol),
        ("multiloop-branch", options.multiloop_branch_kcal_mol),
        ("multiloop-unpaired", options.multiloop_unpaired_kcal_mol),
        (
            "pseudoloop-initiation",
            options.pseudoloop_initiation_kcal_mol,
        ),
        (
            "multiloop-pseudoknot",
            options.multiloop_pseudoknot_kcal_mol,
        ),
        ("nested-pseudoknot", options.nested_pseudoknot_kcal_mol),
        ("band", options.band_kcal_mol),
        ("pseudoloop-unpaired", options.pseudoloop_unpaired_kcal_mol),
        ("closed-subregion", options.closed_subregion_kcal_mol),
        ("spanning-stack-factor", options.spanning_stack_factor),
        ("spanning-internal-factor", options.spanning_internal_factor),
        (
            "spanning-multiloop-init",
            options.spanning_multiloop_init_kcal_mol,
        ),
        (
            "spanning-multiloop-branch",
            options.spanning_multiloop_branch_kcal_mol,
        ),
        (
            "spanning-multiloop-unpaired",
            options.spanning_multiloop_unpaired_kcal_mol,
        ),
    ] {
        if !value.is_finite() {
            return Err(RnaError::InvalidOption(format!(
                "conditional density-2 {name} must be finite"
            )));
        }
    }
    if options.spanning_stack_factor < 0.0 || options.spanning_internal_factor < 0.0 {
        return Err(RnaError::InvalidOption(
            "spanning-loop scale factors must be nonnegative".into(),
        ));
    }
    Ok(())
}

fn polynomial_compatible_constraints(config: &ConstraintConfig) -> bool {
    config.force_paired.is_empty() && config.force_pairs.is_empty() && !config.no_lonely_pairs
}

fn validate_planar_layer(
    name: &str,
    parsed: &crate::structure::ParsedStructure,
    pairs: &[(usize, usize)],
    min_loop: usize,
) -> Result<(), RnaError> {
    if !is_noncrossing(pairs) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional {name} layer must be pseudoknot-free"
        )));
    }
    if let Some(pair) = parsed.pairs.iter().find(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional {name} layer contains noncanonical pair {}-{}",
            pair.i, pair.j
        )));
    }
    if let Some(&(i, j)) = pairs.iter().find(|&&(i, j)| j - i <= min_loop) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional {name} pair {}-{} violates min-loop={min_loop}",
            i + 1,
            j + 1
        )));
    }
    Ok(())
}

#[derive(Clone, Debug)]
struct BandTopology {
    band_count: usize,
    crossing_components: Vec<Vec<usize>>,
    pseudoknot_bands: HashSet<usize>,
    pseudoloop_unpaired: HashSet<usize>,
    pseudoloops: Vec<PseudoloopTopology>,
    maximum_density: usize,
}

#[derive(Clone, Debug)]
struct PseudoloopTopology {
    region: (usize, usize),
    direct_closed_subregions: Vec<(usize, usize)>,
}

fn band_topology(length: usize, pairs: &[(usize, usize)]) -> BandTopology {
    // The published definition applies the band relation only to
    // pseudoknotted base pairs.  Ordinary nested helices therefore contribute
    // zero to pseudoloop density.
    let pseudoknotted = (0..pairs.len())
        .filter(|&index| {
            pairs
                .iter()
                .enumerate()
                .any(|(other, &pair)| other != index && crosses(pairs[index], pair))
        })
        .collect::<Vec<_>>();
    let mut parent = (0..pseudoknotted.len()).collect::<Vec<_>>();
    for outer_slot in 0..pseudoknotted.len() {
        for inner_slot in 0..pseudoknotted.len() {
            if outer_slot == inner_slot
                || !directly_banded(
                    pairs[pseudoknotted[outer_slot]],
                    pairs[pseudoknotted[inner_slot]],
                    pairs,
                )
            {
                continue;
            }
            union(&mut parent, outer_slot, inner_slot);
        }
    }
    let mut groups = HashMap::<usize, Vec<usize>>::new();
    for (slot, &pair_index) in pseudoknotted.iter().enumerate() {
        let root = find(&mut parent, slot);
        groups.entry(root).or_default().push(pair_index);
    }
    let mut bands = groups.into_values().collect::<Vec<_>>();
    bands.sort_by_key(|band| band.iter().map(|&index| pairs[index]).min());

    let mut adjacency = vec![Vec::new(); bands.len()];
    for left in 0..bands.len() {
        for right in left + 1..bands.len() {
            if bands[left]
                .iter()
                .any(|&a| bands[right].iter().any(|&b| crosses(pairs[a], pairs[b])))
            {
                adjacency[left].push(right);
                adjacency[right].push(left);
            }
        }
    }
    let mut seen = vec![false; bands.len()];
    let mut crossing_components = Vec::new();
    for start in 0..bands.len() {
        if adjacency[start].is_empty() || seen[start] {
            continue;
        }
        seen[start] = true;
        let mut stack = vec![start];
        let mut component = Vec::new();
        while let Some(band) = stack.pop() {
            component.push(band);
            for &next in &adjacency[band] {
                if !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        component.sort_unstable();
        crossing_components.push(component);
    }

    let pair_index = pairs
        .iter()
        .enumerate()
        .map(|(index, &pair)| (pair, index))
        .collect::<HashMap<_, _>>();
    let mut partner = vec![None; length];
    for &(i, j) in pairs {
        partner[i] = Some(j);
        partner[j] = Some(i);
    }
    let pseudoknotted_set = pseudoknotted.iter().copied().collect::<HashSet<_>>();
    let mut closed_regions = Vec::new();
    for i in 0..length {
        for j in i + 1..length {
            if !weakly_closed(i, j, pairs) {
                continue;
            }
            let decomposable = (i..j)
                .any(|split| weakly_closed(i, split, pairs) && weakly_closed(split + 1, j, pairs));
            if !decomposable {
                closed_regions.push((i, j));
            }
        }
    }
    let pseudoloop_regions = closed_regions
        .iter()
        .copied()
        .filter(|&(i, j)| {
            let Some(left_partner) = partner[i] else {
                return false;
            };
            let Some(right_partner) = partner[j] else {
                return false;
            };
            let left_pair = (i.min(left_partner), i.max(left_partner));
            let right_pair = (j.min(right_partner), j.max(right_partner));
            pair_index
                .get(&left_pair)
                .is_some_and(|index| pseudoknotted_set.contains(index))
                && pair_index
                    .get(&right_pair)
                    .is_some_and(|index| pseudoknotted_set.contains(index))
        })
        .collect::<Vec<_>>();

    let mut maximum_density = 0usize;
    let mut pseudoloops = Vec::new();
    let mut associated_band_set = HashSet::new();
    let paired = pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .collect::<HashSet<_>>();
    let mut pseudoloop_unpaired = HashSet::new();
    for region in pseudoloop_regions {
        let proper_closed = closed_regions
            .iter()
            .copied()
            .filter(|&candidate| strictly_contains(region, candidate))
            .collect::<Vec<_>>();
        let direct_closed_subregions = proper_closed
            .iter()
            .copied()
            .filter(|&candidate| {
                !proper_closed
                    .iter()
                    .copied()
                    .any(|other| other != candidate && strictly_contains(other, candidate))
            })
            .collect::<Vec<_>>();
        let associated_bands = bands
            .iter()
            .enumerate()
            .filter_map(|(band_index, band)| {
                let (outer, _) = band_closing_pairs(band, pairs);
                if !(region.0 <= outer.0 && outer.1 <= region.1) {
                    return None;
                }
                let nested = proper_closed
                    .iter()
                    .any(|&(i, j)| i <= outer.0 && outer.1 <= j);
                (!nested).then_some(band_index)
            })
            .collect::<Vec<_>>();
        for &band in &associated_bands {
            associated_band_set.insert(band);
        }
        for position in region.0..=region.1 {
            let density = associated_bands
                .iter()
                .filter(|&&band| {
                    let (outer, _) = band_closing_pairs(&bands[band], pairs);
                    outer.0 <= position && position <= outer.1
                })
                .count();
            maximum_density = maximum_density.max(density);
            if !paired.contains(&position)
                && !direct_closed_subregions
                    .iter()
                    .any(|&(i, j)| i <= position && position <= j)
            {
                pseudoloop_unpaired.insert(position);
            }
        }
        pseudoloops.push(PseudoloopTopology {
            region,
            direct_closed_subregions,
        });
    }
    BandTopology {
        band_count: bands.len(),
        crossing_components,
        pseudoknot_bands: associated_band_set,
        pseudoloop_unpaired,
        pseudoloops,
        maximum_density,
    }
}

fn strictly_contains(outer: (usize, usize), inner: (usize, usize)) -> bool {
    outer != inner && outer.0 <= inner.0 && inner.1 <= outer.1
}

fn band_closing_pairs(
    band: &[usize],
    pairs: &[(usize, usize)],
) -> ((usize, usize), (usize, usize)) {
    let outer = band
        .iter()
        .map(|&index| pairs[index])
        .min_by_key(|&(i, _)| i)
        .expect("a band is nonempty");
    let inner = band
        .iter()
        .map(|&index| pairs[index])
        .max_by_key(|&(i, _)| i)
        .expect("a band is nonempty");
    (outer, inner)
}

/// `inner` is directly banded by `outer` when it is strictly nested and both
/// intervening flanks are weakly closed (no base pair has exactly one endpoint
/// in the flank).  The transitive closure of this relation forms a band.
fn directly_banded(outer: (usize, usize), inner: (usize, usize), pairs: &[(usize, usize)]) -> bool {
    outer.0 < inner.0
        && inner.1 < outer.1
        && weakly_closed(outer.0 + 1, inner.0.saturating_sub(1), pairs)
        && weakly_closed(inner.1 + 1, outer.1.saturating_sub(1), pairs)
}

fn weakly_closed(start: usize, end: usize, pairs: &[(usize, usize)]) -> bool {
    if start > end {
        return true;
    }
    pairs.iter().all(|&(i, j)| {
        let left = start <= i && i <= end;
        let right = start <= j && j <= end;
        left == right
    })
}

fn find(parent: &mut [usize], index: usize) -> usize {
    if parent[index] != index {
        parent[index] = find(parent, parent[index]);
    }
    parent[index]
}

fn union(parent: &mut [usize], left: usize, right: usize) {
    let left = find(parent, left);
    let right = find(parent, right);
    if left != right {
        parent[right] = left;
    }
}

#[allow(clippy::too_many_arguments)]
fn conditional_energy(
    length: usize,
    seed_pairs: &[(usize, usize)],
    added_pairs: &[(usize, usize)],
    seed: &EnergyBreakdown,
    added: &EnergyBreakdown,
    topology: &BandTopology,
    options: &ConditionalDensity2Options,
    constraint_kcal_mol: f64,
    model_name: &'static str,
) -> ConditionalDensity2Energy {
    let stack_correction =
        spanning_loop_correction(length, seed_pairs, added_pairs, seed, options, "stack")
            + spanning_loop_correction(length, added_pairs, seed_pairs, added, options, "stack");
    let internal_correction =
        spanning_loop_correction(length, seed_pairs, added_pairs, seed, options, "internal")
            + spanning_loop_correction(length, added_pairs, seed_pairs, added, options, "internal");
    let multiloop_correction =
        spanning_loop_correction(length, seed_pairs, added_pairs, seed, options, "multiloop")
            + spanning_loop_correction(
                length,
                added_pairs,
                seed_pairs,
                added,
                options,
                "multiloop",
            );
    let nested_pseudoknot_count = topology
        .pseudoloops
        .iter()
        .filter(|pseudoloop| {
            topology
                .pseudoloops
                .iter()
                .any(|outer| strictly_contains(outer.region, pseudoloop.region))
        })
        .count();
    let multiloop_pseudoknot_count = topology
        .pseudoloops
        .iter()
        .filter(|pseudoloop| {
            !topology
                .pseudoloops
                .iter()
                .any(|outer| strictly_contains(outer.region, pseudoloop.region))
                && (region_is_on_multiloop(pseudoloop.region, seed_pairs)
                    || region_is_on_multiloop(pseudoloop.region, added_pairs))
        })
        .count();
    let exterior_pseudoloop_count = topology
        .pseudoloops
        .len()
        .saturating_sub(nested_pseudoknot_count + multiloop_pseudoknot_count);
    let initiation = options.pseudoloop_initiation_kcal_mol * exterior_pseudoloop_count as f64;
    let multiloop_pseudoknot =
        options.multiloop_pseudoknot_kcal_mol * multiloop_pseudoknot_count as f64;
    let nested_pseudoknot = options.nested_pseudoknot_kcal_mol * nested_pseudoknot_count as f64;
    let band = options.band_kcal_mol * topology.pseudoknot_bands.len() as f64;
    let unpaired = options.pseudoloop_unpaired_kcal_mol * topology.pseudoloop_unpaired.len() as f64;
    let closed_subregion_count = topology
        .pseudoloops
        .iter()
        .map(|pseudoloop| pseudoloop.direct_closed_subregions.len())
        .sum::<usize>();
    let closed_subregion = options.closed_subregion_kcal_mol * closed_subregion_count as f64;
    let total = seed.total_kcal_mol
        + added.total_kcal_mol
        + stack_correction
        + internal_correction
        + multiloop_correction
        + initiation
        + multiloop_pseudoknot
        + nested_pseudoknot
        + band
        + unpaired
        + closed_subregion
        + constraint_kcal_mol;
    ConditionalDensity2Energy {
        seed_layer_kcal_mol: seed.total_kcal_mol,
        added_layer_kcal_mol: added.total_kcal_mol,
        spanning_stack_correction_kcal_mol: stack_correction,
        spanning_internal_correction_kcal_mol: internal_correction,
        spanning_multiloop_correction_kcal_mol: multiloop_correction,
        pseudoloop_initiation_kcal_mol: initiation,
        multiloop_pseudoknot_kcal_mol: multiloop_pseudoknot,
        nested_pseudoknot_kcal_mol: nested_pseudoknot,
        band_kcal_mol: band,
        pseudoloop_unpaired_kcal_mol: unpaired,
        closed_subregion_kcal_mol: closed_subregion,
        constraint_kcal_mol,
        decomposition_alignment_kcal_mol: 0.0,
        total_kcal_mol: total,
        crossing_component_count: topology.crossing_components.len(),
        exterior_pseudoloop_count,
        multiloop_pseudoknot_count,
        nested_pseudoknot_count,
        band_count: topology.band_count,
        pseudoknot_band_count: topology.pseudoknot_bands.len(),
        pseudoloop_unpaired_count: topology.pseudoloop_unpaired.len(),
        closed_subregion_count,
        maximum_band_density: topology.maximum_density,
        model: model_name,
    }
}

fn region_is_on_multiloop(region: (usize, usize), layer: &[(usize, usize)]) -> bool {
    let mut enclosing = layer
        .iter()
        .copied()
        .filter(|&(i, j)| i < region.0 && region.1 < j)
        .collect::<Vec<_>>();
    enclosing.sort_by_key(|&(i, j)| j - i);
    let Some(outer) = enclosing.first().copied() else {
        return false;
    };
    let (branches, _) = direct_children(0, outer, layer);
    if branches < 2 {
        return false;
    }
    let direct = layer
        .iter()
        .copied()
        .filter(|&(i, j)| outer.0 < i && j < outer.1)
        .filter(|&(i, j)| {
            !layer
                .iter()
                .any(|&(k, l)| outer.0 < k && k < i && j < l && l < outer.1)
        });
    !direct
        .into_iter()
        .any(|(i, j)| i <= region.0 && region.1 <= j)
}

fn spanning_loop_correction(
    length: usize,
    layer: &[(usize, usize)],
    other: &[(usize, usize)],
    breakdown: &EnergyBreakdown,
    options: &ConditionalDensity2Options,
    kind: &str,
) -> f64 {
    breakdown
        .loop_energies
        .iter()
        .filter(|entry| entry.kind == kind && entry.i > 0)
        .filter(|entry| {
            let pair = (entry.i - 1, entry.j - 1);
            other.iter().any(|&candidate| crosses(pair, candidate))
        })
        .map(|entry| match kind {
            "stack" => (options.spanning_stack_factor - 1.0) * entry.energy_kcal_mol,
            "internal" => (options.spanning_internal_factor - 1.0) * entry.energy_kcal_mol,
            "multiloop" => {
                let pair = (entry.i - 1, entry.j - 1);
                let (branches, unpaired) = direct_children(length, pair, layer);
                let replacement = options.spanning_multiloop_init_kcal_mol
                    + options.spanning_multiloop_branch_kcal_mol * branches as f64
                    + options.spanning_multiloop_unpaired_kcal_mol * unpaired as f64;
                replacement - entry.energy_kcal_mol
            }
            _ => 0.0,
        })
        .sum()
}

fn direct_children(
    _length: usize,
    outer: (usize, usize),
    pairs: &[(usize, usize)],
) -> (usize, usize) {
    let mut children = pairs
        .iter()
        .copied()
        .filter(|&(i, j)| outer.0 < i && j < outer.1)
        .filter(|&(i, j)| {
            !pairs
                .iter()
                .any(|&(k, l)| outer.0 < k && k < i && j < l && l < outer.1)
        })
        .collect::<Vec<_>>();
    children.sort_unstable();
    let occupied = children.iter().map(|&(i, j)| j - i + 1).sum::<usize>();
    let unpaired = (outer.1 - outer.0 - 1).saturating_sub(occupied);
    (children.len(), unpaired)
}

#[allow(clippy::too_many_arguments)]
fn finish_ensemble(
    sequence: String,
    seed_structure: String,
    seed_pairs: Vec<(usize, usize)>,
    states: Vec<State>,
    maximum_density: usize,
    rt: f64,
    gamma: f64,
    constraints: ConstraintSummary,
) -> Result<ConditionalDensity2Result, RnaError> {
    let maximum = states
        .iter()
        .map(|state| state.log_weight)
        .fold(f64::NEG_INFINITY, f64::max);
    let log_z = maximum
        + states
            .iter()
            .map(|state| (state.log_weight - maximum).exp())
            .sum::<f64>()
            .ln();
    let mut pair_mass = HashMap::<(usize, usize), f64>::new();
    for &pair in &seed_pairs {
        pair_mass.insert(pair, 1.0);
    }
    for state in &states {
        let probability = (state.log_weight - log_z).exp();
        for &pair in &state.added {
            *pair_mass.entry(pair).or_default() += probability;
        }
    }
    let mut pair_probabilities = pair_mass
        .into_iter()
        .map(|((i, j), probability)| PairProbability {
            i: i + 1,
            j: j + 1,
            probability: probability.clamp(0.0, 1.0),
        })
        .collect::<Vec<_>>();
    pair_probabilities.sort_by_key(|pair| (pair.i, pair.j));
    let mut unpaired_probabilities = vec![1.0; sequence.len()];
    for pair in &pair_probabilities {
        unpaired_probabilities[pair.i - 1] -= pair.probability;
        unpaired_probabilities[pair.j - 1] -= pair.probability;
    }
    for value in &mut unpaired_probabilities {
        *value = value.clamp(0.0, 1.0);
    }

    let variable_probabilities = pair_probabilities
        .iter()
        .filter(|pair| !seed_pairs.contains(&(pair.i - 1, pair.j - 1)))
        .map(|pair| ((pair.i - 1, pair.j - 1), pair.probability))
        .collect::<HashMap<_, _>>();
    let total_variable_probability = variable_probabilities.values().sum::<f64>();
    let mut centroid_index = 0usize;
    let mut centroid_distance = f64::INFINITY;
    let seed_occupied = seed_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .collect::<HashSet<_>>();
    let mut mea_index = 0usize;
    let mut mea_score = f64::NEG_INFINITY;
    for (index, state) in states.iter().enumerate() {
        let pair_sum = state
            .added
            .iter()
            .map(|pair| variable_probabilities.get(pair).copied().unwrap_or(0.0))
            .sum::<f64>();
        let distance = total_variable_probability + state.added.len() as f64 - 2.0 * pair_sum;
        if distance < centroid_distance {
            centroid_distance = distance;
            centroid_index = index;
        }
        let occupied = state
            .added
            .iter()
            .flat_map(|&(i, j)| [i, j])
            .chain(seed_occupied.iter().copied())
            .collect::<HashSet<_>>();
        let score = 2.0 * gamma * pair_sum
            + unpaired_probabilities
                .iter()
                .enumerate()
                .filter(|(position, _)| !occupied.contains(position))
                .map(|(_, probability)| probability)
                .sum::<f64>();
        if score > mea_score {
            mea_score = score;
            mea_index = index;
        }
    }
    let mfe = &states[0];
    Ok(ConditionalDensity2Result {
        mfe_structure: layered_structure(sequence.len(), &seed_pairs, &mfe.added),
        mfe_added_structure: planar_structure(sequence.len(), &mfe.added, '(', ')'),
        mfe_energy_kcal_mol: mfe.energy.total_kcal_mol,
        mfe_energy: mfe.energy.clone(),
        ensemble_free_energy_kcal_mol: -rt * log_z,
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        pair_probabilities,
        unpaired_probabilities,
        centroid_structure: layered_structure(
            sequence.len(),
            &seed_pairs,
            &states[centroid_index].added,
        ),
        centroid_distance,
        mea_structure: layered_structure(sequence.len(), &seed_pairs, &states[mea_index].added),
        mea_score,
        constraints,
        state_count: states.len(),
        state_count_exact: states.len().to_string(),
        state_space_complete: true,
        maximum_density,
        time_complexity: "exponential",
        space_complexity: "exponential",
        algorithm: "complete exhaustive conditional two-planar-layer density-2 oracle",
        model: mfe.energy.model,
        sequence,
        seed_structure,
    })
}

fn planar_structure(length: usize, pairs: &[(usize, usize)], open: char, close: char) -> String {
    let mut structure = vec!['.'; length];
    for &(i, j) in pairs {
        structure[i] = open;
        structure[j] = close;
    }
    structure.into_iter().collect()
}

fn layered_structure(length: usize, seed: &[(usize, usize)], added: &[(usize, usize)]) -> String {
    let mut structure = vec!['.'; length];
    for &(i, j) in seed {
        structure[i] = '(';
        structure[j] = ')';
    }
    for &(i, j) in added {
        structure[i] = '[';
        structure[j] = ']';
    }
    structure.into_iter().collect()
}

fn is_noncrossing(pairs: &[(usize, usize)]) -> bool {
    pairs.iter().enumerate().all(|(index, &left)| {
        pairs[index + 1..]
            .iter()
            .all(|&right| !crosses(left, right))
    })
}

fn crosses(left: (usize, usize), right: (usize, usize)) -> bool {
    (left.0 < right.0 && right.0 < left.1 && left.1 < right.1)
        || (right.0 < left.0 && left.0 < right.1 && right.1 < left.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition;

    #[test]
    fn ordinary_stacked_pairs_do_not_form_a_pseudoknot_band() {
        let topology = band_topology(10, &[(0, 9), (1, 8), (2, 7)]);
        assert_eq!(topology.band_count, 0);
        assert_eq!(topology.maximum_density, 0);
    }

    #[test]
    fn stacked_pseudoknotted_pairs_form_one_band() {
        let topology = band_topology(10, &[(0, 6), (1, 5), (3, 8)]);
        assert_eq!(topology.band_count, 2);
        assert_eq!(topology.crossing_components.len(), 1);
        assert_eq!(topology.pseudoloops.len(), 1);
        assert_eq!(topology.pseudoknot_bands.len(), 2);
        assert_eq!(topology.maximum_density, 2);
    }

    #[test]
    fn three_mutually_crossing_bands_have_density_three() {
        let topology = band_topology(10, &[(0, 6), (1, 7), (2, 8)]);
        assert_eq!(topology.band_count, 3);
        assert_eq!(topology.crossing_components.len(), 1);
        assert_eq!(topology.maximum_density, 3);
    }

    #[test]
    fn direct_closed_subregion_is_associated_with_its_pseudoloop() {
        let topology = band_topology(10, &[(0, 6), (3, 9), (1, 2)]);
        assert_eq!(topology.pseudoloops.len(), 1);
        assert_eq!(
            topology.pseudoloops[0].direct_closed_subregions,
            vec![(1, 2)]
        );
    }

    #[test]
    fn fixed_seed_pairs_have_probability_one() {
        let result = conditional_density2_ensemble(
            "GCGAAACGCU",
            "(........)",
            37.0,
            0,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        let seed = result
            .pair_probabilities
            .iter()
            .find(|pair| (pair.i, pair.j) == (1, 10))
            .unwrap();
        assert_eq!(seed.probability, 1.0);
        assert_eq!(result.unpaired_probabilities[0], 0.0);
        assert_eq!(result.unpaired_probabilities[9], 0.0);
        assert!(result.log_partition_function.is_finite());
        assert!(result.state_space_complete);
    }

    #[test]
    fn supplied_crossing_state_reports_published_penalty_components() {
        let result = evaluate_conditional_density2_structure(
            "GAGAACAACU",
            "(....)....",
            "..(.....).",
            37.0,
            0,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        assert_eq!(result.structure, "(.[..)..].");
        assert_eq!(result.energy.exterior_pseudoloop_count, 1);
        assert_eq!(result.energy.multiloop_pseudoknot_count, 0);
        assert_eq!(result.energy.nested_pseudoknot_count, 0);
        assert_eq!(result.energy.pseudoknot_band_count, 2);
        assert_eq!(result.energy.maximum_band_density, 2);
        assert!((result.energy.pseudoloop_initiation_kcal_mol + 1.38).abs() < 1e-12);
        assert!((result.energy.band_kcal_mol - 4.92).abs() < 1e-12);
    }

    #[test]
    fn empty_seed_reduces_state_energies_to_turner_layer_energy() {
        let result = conditional_density2_ensemble(
            "GGGAAACCC",
            ".........",
            37.0,
            3,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        assert_eq!(result.mfe_energy.crossing_component_count, 0);
        assert_eq!(result.mfe_energy.pseudoknot_band_count, 0);
        assert_eq!(result.mfe_energy.spanning_stack_correction_kcal_mol, 0.0);
        assert_eq!(
            result.mfe_energy.total_kcal_mol,
            result.mfe_energy.added_layer_kcal_mol
        );
    }

    #[test]
    fn empty_seed_matches_polynomial_planar_partition_and_marginals() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let expected = partition(sequence, 37.0, 3, &model).unwrap();
        let actual = conditional_density2_ensemble(
            sequence,
            ".........",
            37.0,
            3,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        assert!((actual.log_partition_function - expected.log_partition_function).abs() < 1e-10);
        let actual_pairs = actual
            .pair_probabilities
            .iter()
            .map(|pair| ((pair.i, pair.j), pair.probability))
            .collect::<HashMap<_, _>>();
        for pair in expected.pair_probabilities {
            let observed = actual_pairs.get(&(pair.i, pair.j)).copied().unwrap_or(0.0);
            assert!(
                (observed - pair.probability).abs() < 1e-10,
                "pair {}-{}: {observed} != {}",
                pair.i,
                pair.j,
                pair.probability
            );
        }
        for (observed, expected) in actual
            .unpaired_probabilities
            .iter()
            .zip(expected.unpaired_probabilities)
        {
            assert!((observed - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn short_fixed_seed_state_count_is_complete() {
        let result = conditional_density2_ensemble(
            "GCAUGC",
            "(....)",
            37.0,
            0,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        assert_eq!(result.state_count, 5);
        assert_eq!(result.state_count_exact, "5");
    }
}
