//! Probability-directed pseudoknot prediction.
//!
//! Pair marginals from the pseudoknot-free Turner ensemble are decoded with a
//! deterministic iterative mutual-maximum rule (the central ProbKnot idea).
//! Selected crossing pairs are serialized with extended dot-bracket levels.
//! A separate, explicitly parameterized diagnostic energy for the decoded
//! matching is reported; it is never presented as a complete DP09 loop
//! decomposition or a Turner nearest-neighbor free energy.

use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, parse_structure, RnaError};
use crate::ConstraintConfig;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct PseudoknotOptions {
    pub threshold: f64,
    pub iterations: usize,
    pub min_helix: usize,
    pub gamma: f64,
    pub initiation_kcal_mol: f64,
    pub crossing_kcal_mol: f64,
    pub unpaired_kcal_mol: f64,
    /// Favour helices whose individual pairs are supported by the planar
    /// partition ensemble. Units are kcal/mol per unit summed probability.
    pub evidence_weight_kcal_mol: f64,
    /// Optional user-requested restriction on the number of disjoint H-type
    /// components in one thermodynamic ensemble state. `None` is exhaustive.
    pub max_components: Option<usize>,
    /// Optional user-requested restriction on retained precursor helices,
    /// H-type cores, and combined states. `None` is exhaustive; the library
    /// never applies an implicit beam.
    pub max_ensemble_states: Option<usize>,
    /// Enumerate the complete canonical base-pair matching space, including
    /// arbitrary crossing topology.  This is an exact exponential algorithm
    /// intended for finite-system analysis and correctness oracles; it is
    /// never enabled implicitly.
    pub exact_arbitrary_ensemble: bool,
}

impl Default for PseudoknotOptions {
    fn default() -> Self {
        Self {
            threshold: 0.0,
            iterations: 1,
            min_helix: 3,
            gamma: 1.0,
            // HotKnots 2.0 DP09 exterior pseudoloop, per-band, and
            // pseudoloop-unpaired parameters at 37 C and 1 M NaCl.
            initiation_kcal_mol: -1.38,
            crossing_kcal_mol: 2.46,
            unpaired_kcal_mol: 0.06,
            evidence_weight_kcal_mol: 0.0,
            max_components: None,
            max_ensemble_states: None,
            exact_arbitrary_ensemble: false,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct PseudoknotPair {
    pub i: usize,
    pub j: usize,
    pub probability: f64,
    pub level: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct PseudoknotEnergyBreakdown {
    pub total_kcal_mol: f64,
    pub stack_kcal_mol: f64,
    pub terminal_kcal_mol: f64,
    pub initiation_kcal_mol: f64,
    pub crossing_kcal_mol: f64,
    pub unpaired_kcal_mol: f64,
    pub helix_count: usize,
    pub crossing_count: usize,
    pub enclosed_unpaired_count: usize,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct PseudoknotEvaluationResult {
    pub sequence: String,
    pub structure: String,
    pub pair_count: usize,
    pub crossing_count: usize,
    pub crossing_component_count: usize,
    pub energy: PseudoknotEnergyBreakdown,
}

#[derive(Clone, Debug, Serialize)]
pub struct PseudoknotResult {
    pub sequence: String,
    pub structure: String,
    pub pairs: Vec<PseudoknotPair>,
    pub pair_count: usize,
    pub crossing_count: usize,
    pub bracket_levels: usize,
    pub expected_accuracy_score: f64,
    pub decoded_diagnostic_energy: PseudoknotEnergyBreakdown,
    pub hybrid_structure: String,
    pub hybrid_pairs: Vec<PseudoknotPair>,
    pub hybrid_pair_count: usize,
    pub hybrid_crossing_count: usize,
    pub hybrid_bracket_levels: usize,
    pub hybrid_expected_accuracy_score: f64,
    pub matching_centroid_structure: String,
    pub matching_centroid_pairs: Vec<PseudoknotPair>,
    pub matching_centroid_pair_count: usize,
    pub matching_centroid_crossing_count: usize,
    pub matching_centroid_distance: f64,
    pub matching_mea_structure: String,
    pub matching_mea_pairs: Vec<PseudoknotPair>,
    pub matching_mea_pair_count: usize,
    pub matching_mea_crossing_count: usize,
    pub matching_mea_score: f64,
    pub thermodynamic_core_pairs: Vec<PairProbability>,
    pub thermodynamic_core_score_kcal_mol: Option<f64>,
    pub thermodynamic_component_count: usize,
    pub restricted_mfe_structure: String,
    pub restricted_mfe_energy_kcal_mol: f64,
    pub restricted_ensemble_free_energy_kcal_mol: f64,
    pub restricted_partition_function: f64,
    pub restricted_log_partition_function: f64,
    pub restricted_pair_probabilities: Vec<PairProbability>,
    pub restricted_unpaired_probabilities: Vec<f64>,
    pub restricted_centroid_structure: String,
    pub restricted_centroid_distance: f64,
    pub restricted_mea_structure: String,
    pub restricted_mea_score: f64,
    /// State count saturated at `usize::MAX` for compact consumers.
    pub restricted_state_count: usize,
    /// Exact decimal state count, including the open state.
    pub restricted_state_count_exact: String,
    pub restricted_state_space_complete: bool,
    pub restricted_ensemble_model: &'static str,
    /// Complete arbitrary-topology ensemble when explicitly requested.
    pub exact_arbitrary_ensemble: Option<ExactArbitraryEnsemble>,
    pub source_pair_probabilities: Vec<PairProbability>,
    pub threshold: f64,
    pub iterations: usize,
    pub min_helix: usize,
    pub decoder: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ExactArbitraryEnsemble {
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub centroid_structure: String,
    pub centroid_distance: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub state_count: usize,
    pub state_count_exact: String,
    pub state_space_complete: bool,
    pub model: &'static str,
}

#[allow(clippy::too_many_arguments)]
pub fn predict_pseudoknot(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &PseudoknotOptions,
) -> Result<PseudoknotResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    predict_pseudoknot_with_model(sequence, min_loop, &model, options)
}

pub fn predict_pseudoknot_with_model(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> Result<PseudoknotResult, RnaError> {
    validate_options(options)?;
    let analysis = crate::analyze_with_model(
        sequence.into(),
        min_loop,
        options.gamma,
        model.clone(),
        &ConstraintConfig::default(),
    )?;
    decode_pseudoknot_with_model(
        &analysis.sequence,
        &analysis.pair_probabilities,
        &analysis.unpaired_probabilities,
        min_loop,
        model,
        options,
    )
}

/// Evaluate any canonical extended-dot-bracket matching with the same
/// explicitly defined generalized pseudoknot diagnostic used by the
/// arbitrary-topology ensemble.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_pseudoknot_structure(
    sequence: &str,
    structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &PseudoknotOptions,
) -> Result<PseudoknotEvaluationResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    evaluate_pseudoknot_structure_with_model(sequence, structure, min_loop, &model, options)
}

pub fn evaluate_pseudoknot_structure_with_model(
    sequence: &str,
    structure: &str,
    min_loop: usize,
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> Result<PseudoknotEvaluationResult, RnaError> {
    validate_options(options)?;
    let parsed = parse_structure(sequence, structure)?;
    if !parsed.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "pseudoknot energy evaluation accepts one strand".into(),
        ));
    }
    if let Some(pair) = parsed.pairs.iter().find(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(format!(
            "pseudoknot energy evaluation does not parameterize the noncanonical pair {}-{}",
            pair.i, pair.j
        )));
    }
    if let Some(pair) = parsed.pairs.iter().find(|pair| pair.j - pair.i <= min_loop) {
        return Err(RnaError::InvalidOption(format!(
            "pseudoknot pair {}-{} violates min-loop={min_loop}",
            pair.i, pair.j
        )));
    }
    let pairs = parsed
        .pairs
        .iter()
        .map(|pair| PairProbability {
            i: pair.i,
            j: pair.j,
            probability: 0.0,
        })
        .collect::<Vec<_>>();
    let energy = decoded_diagnostic_energy(parsed.sequence.as_bytes(), &pairs, model, options);
    Ok(PseudoknotEvaluationResult {
        sequence: parsed.sequence,
        structure: parsed.structure,
        pair_count: pairs.len(),
        crossing_count: count_crossings(&pairs),
        crossing_component_count: crossing_component_count(&pairs),
        energy,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn decode_pseudoknot(
    sequence: &str,
    probabilities: &[PairProbability],
    unpaired_probabilities: &[f64],
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &PseudoknotOptions,
) -> Result<PseudoknotResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    decode_pseudoknot_with_model(
        sequence,
        probabilities,
        unpaired_probabilities,
        min_loop,
        &model,
        options,
    )
}

pub fn decode_pseudoknot_with_model(
    sequence: &str,
    probabilities: &[PairProbability],
    unpaired_probabilities: &[f64],
    min_loop: usize,
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> Result<PseudoknotResult, RnaError> {
    validate_options(options)?;
    let temperature_celsius = model.temperature_celsius();
    let sequence = normalize_sequence(sequence)?;
    let n = sequence.len();
    if unpaired_probabilities.len() != n {
        return Err(RnaError::LengthMismatch {
            sequence: n,
            structure: unpaired_probabilities.len(),
        });
    }
    if unpaired_probabilities
        .iter()
        .any(|probability| !probability.is_finite() || !(0.0..=1.0).contains(probability))
    {
        return Err(RnaError::InvalidOption(
            "pseudoknot unpaired probabilities must be finite values in [0, 1]".into(),
        ));
    }
    let mut probability_keys = HashSet::new();
    for entry in probabilities {
        if entry.i == 0
            || entry.i >= entry.j
            || entry.j > n
            || !entry.probability.is_finite()
            || !(0.0..=1.0).contains(&entry.probability)
        {
            return Err(RnaError::InvalidOption(
                "pseudoknot pair probabilities must be unique valid pairs with values in [0, 1]"
                    .into(),
            ));
        }
        if !model.can_pair(
            sequence.as_bytes()[entry.i - 1],
            sequence.as_bytes()[entry.j - 1],
        ) {
            return Err(RnaError::InvalidOption(format!(
                "pseudoknot source pair ({}, {}) is noncanonical",
                entry.i, entry.j
            )));
        }
        if !probability_keys.insert((entry.i, entry.j)) {
            return Err(RnaError::InvalidOption(format!(
                "duplicate pseudoknot source pair ({}, {})",
                entry.i, entry.j
            )));
        }
    }
    let mut candidates: Vec<_> = probabilities
        .iter()
        .filter(|entry| entry.probability > options.threshold)
        .cloned()
        .collect();
    candidates.sort_by(|a, b| {
        b.probability
            .total_cmp(&a.probability)
            .then_with(|| a.i.cmp(&b.i))
            .then_with(|| a.j.cmp(&b.j))
    });
    let mut available = vec![true; n];
    let mut selected = Vec::new();
    for _ in 0..options.iterations {
        let mut best: Vec<Option<(usize, f64)>> = vec![None; n];
        for entry in &candidates {
            let i = entry.i - 1;
            let j = entry.j - 1;
            if !available[i] || !available[j] {
                continue;
            }
            update_best(&mut best[i], j, entry.probability);
            update_best(&mut best[j], i, entry.probability);
        }
        let mut round = Vec::new();
        for entry in &candidates {
            let i = entry.i - 1;
            let j = entry.j - 1;
            if available[i]
                && available[j]
                && best[i].is_some_and(|(partner, _)| partner == j)
                && best[j].is_some_and(|(partner, _)| partner == i)
            {
                round.push(entry.clone());
            }
        }
        if round.is_empty() {
            break;
        }
        round.sort_by(|a, b| {
            b.probability
                .total_cmp(&a.probability)
                .then_with(|| a.i.cmp(&b.i))
                .then_with(|| a.j.cmp(&b.j))
        });
        for entry in round {
            let i = entry.i - 1;
            let j = entry.j - 1;
            if available[i] && available[j] {
                available[i] = false;
                available[j] = false;
                selected.push(entry);
            }
        }
    }
    selected = remove_short_helices(selected, options.min_helix);
    let cores = h_type_cores(sequence.as_bytes(), probabilities, model, options);
    let (restricted, restricted_ensemble_model) = if options.max_ensemble_states.is_some() {
        let states = combine_compatible_cores(&cores, options);
        (
            restricted_ensemble(n, &states, temperature_celsius, options.gamma)?,
            "caller-limited nucleotide-disjoint H-type component enumeration plus the open state",
        )
    } else {
        (
            interval_component_ensemble(
                n,
                &cores,
                temperature_celsius,
                options.gamma,
                options.max_components,
            )?,
            "exact span-disjoint independent H-type component ensemble plus the open state",
        )
    };
    let state = restricted.mfe_state.clone();
    let thermodynamic_core_pairs = state
        .as_ref()
        .map(|state| state.pairs.clone())
        .unwrap_or_default();
    let thermodynamic_core_score_kcal_mol = state.as_ref().map(|state| state.score_kcal_mol);
    let thermodynamic_component_count = state.as_ref().map_or(0, |state| state.component_count);
    selected.sort_by_key(|entry| (entry.i, entry.j));
    let (structure, pairs, crossing_count, bracket_levels) = serialize_prediction(n, &selected)?;
    let expected_accuracy_score = expected_accuracy(
        &selected,
        unpaired_probabilities,
        options.gamma,
        sequence.len(),
    );
    let energy = decoded_diagnostic_energy(sequence.as_bytes(), &selected, model, options);

    let exact_arbitrary_ensemble = options
        .exact_arbitrary_ensemble
        .then(|| {
            exact_arbitrary_matching_ensemble(
                sequence.as_bytes(),
                probabilities,
                temperature_celsius,
                min_loop,
                model,
                options,
            )
        })
        .transpose()?
        .map(ExactArbitraryEnsemble::from);

    let mut hybrid_selected = selected.clone();
    let decoder_seed = cores
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| {
            left.score_kcal_mol
                .total_cmp(&right.score_kcal_mol)
                .then_with(|| left.pairs.len().cmp(&right.pairs.len()).reverse())
        })
        .filter(|(_, core)| core.score_kcal_mol < 0.0)
        .map(|(index, core)| PseudoknotState {
            pairs: core.pairs.clone(),
            score_kcal_mol: core.score_kcal_mol,
            component_count: 1,
            last_core: index,
        });
    if let Some(state) = decoder_seed {
        // Seed the probability decoder with one jointly favourable H-type
        // core. The exact multi-component MFE remains available separately;
        // forcing all of its independent components into this heuristic
        // decoder would suppress otherwise well-supported planar marginals.
        let mut occupied = vec![false; n];
        let mut hybrid = Vec::with_capacity(state.pairs.len() + hybrid_selected.len());
        for entry in state.pairs {
            occupied[entry.i - 1] = true;
            occupied[entry.j - 1] = true;
            hybrid.push(entry);
        }
        for entry in hybrid_selected {
            if !occupied[entry.i - 1] && !occupied[entry.j - 1] {
                occupied[entry.i - 1] = true;
                occupied[entry.j - 1] = true;
                hybrid.push(entry);
            }
        }
        hybrid_selected = hybrid;
    }
    hybrid_selected.sort_by_key(|entry| (entry.i, entry.j));
    let (hybrid_structure, hybrid_pairs, hybrid_crossing_count, hybrid_bracket_levels) =
        serialize_prediction(n, &hybrid_selected)?;
    let hybrid_expected_accuracy_score = expected_accuracy(
        &hybrid_selected,
        unpaired_probabilities,
        options.gamma,
        sequence.len(),
    );

    let matching_centroid_selected =
        maximum_weight_pair_matching(sequence.as_bytes(), probabilities, |pair| {
            2.0 * pair.probability - 1.0
        });
    let (matching_centroid_structure, matching_centroid_pairs, matching_centroid_crossing_count, _) =
        serialize_prediction(n, &matching_centroid_selected)?;
    let total_pair_probability = probabilities
        .iter()
        .map(|pair| pair.probability)
        .sum::<f64>();
    let matching_centroid_distance = total_pair_probability
        - matching_centroid_selected
            .iter()
            .map(|pair| 2.0 * pair.probability - 1.0)
            .sum::<f64>();

    let matching_mea_selected =
        maximum_weight_pair_matching(sequence.as_bytes(), probabilities, |pair| {
            2.0 * options.gamma * pair.probability
                - unpaired_probabilities[pair.i - 1]
                - unpaired_probabilities[pair.j - 1]
        });
    let (matching_mea_structure, matching_mea_pairs, matching_mea_crossing_count, _) =
        serialize_prediction(n, &matching_mea_selected)?;
    let matching_mea_score = expected_accuracy(
        &matching_mea_selected,
        unpaired_probabilities,
        options.gamma,
        sequence.len(),
    );
    Ok(PseudoknotResult {
        sequence,
        structure,
        pair_count: pairs.len(),
        crossing_count,
        bracket_levels,
        pairs,
        expected_accuracy_score,
        decoded_diagnostic_energy: energy,
        hybrid_structure,
        hybrid_pair_count: hybrid_pairs.len(),
        hybrid_crossing_count,
        hybrid_bracket_levels,
        hybrid_pairs,
        hybrid_expected_accuracy_score,
        matching_centroid_structure,
        matching_centroid_pair_count: matching_centroid_pairs.len(),
        matching_centroid_crossing_count,
        matching_centroid_pairs,
        matching_centroid_distance,
        matching_mea_structure,
        matching_mea_pair_count: matching_mea_pairs.len(),
        matching_mea_crossing_count,
        matching_mea_pairs,
        matching_mea_score,
        thermodynamic_core_pairs,
        thermodynamic_core_score_kcal_mol,
        thermodynamic_component_count,
        restricted_mfe_structure: restricted.mfe_structure,
        restricted_mfe_energy_kcal_mol: restricted.mfe_energy_kcal_mol,
        restricted_ensemble_free_energy_kcal_mol: restricted.ensemble_free_energy_kcal_mol,
        restricted_partition_function: restricted.partition_function,
        restricted_log_partition_function: restricted.log_partition_function,
        restricted_pair_probabilities: restricted.pair_probabilities,
        restricted_unpaired_probabilities: restricted.unpaired_probabilities,
        restricted_centroid_structure: restricted.centroid_structure,
        restricted_centroid_distance: restricted.centroid_distance,
        restricted_mea_structure: restricted.mea_structure,
        restricted_mea_score: restricted.mea_score,
        restricted_state_count: restricted.state_count,
        restricted_state_count_exact: restricted.state_count_exact,
        restricted_state_space_complete: options.max_components.is_none()
            && options.max_ensemble_states.is_none(),
        restricted_ensemble_model,
        exact_arbitrary_ensemble,
        source_pair_probabilities: probabilities.to_vec(),
        threshold: options.threshold,
        iterations: options.iterations,
        min_helix: options.min_helix,
        decoder: "ProbKnot mutual-maximum probability decoder with one-nucleotide-bulge-aware helix filtering",
    })
}

fn serialize_prediction(
    length: usize,
    selected: &[PairProbability],
) -> Result<(String, Vec<PseudoknotPair>, usize, usize), RnaError> {
    let levels = color_crossing_graph(selected)?;
    let structure = extended_dot_bracket(length, selected, &levels)?;
    let crossing_count = count_crossings(selected);
    let pairs = selected
        .iter()
        .zip(levels)
        .map(|(entry, level)| PseudoknotPair {
            i: entry.i,
            j: entry.j,
            probability: entry.probability,
            level,
        })
        .collect::<Vec<_>>();
    let bracket_levels = pairs.iter().map(|pair| pair.level + 1).max().unwrap_or(0);
    Ok((structure, pairs, crossing_count, bracket_levels))
}

fn expected_accuracy(
    selected: &[PairProbability],
    unpaired_probabilities: &[f64],
    gamma: f64,
    length: usize,
) -> f64 {
    let mut paired = vec![false; length];
    let pair_score = selected
        .iter()
        .map(|entry| {
            paired[entry.i - 1] = true;
            paired[entry.j - 1] = true;
            2.0 * gamma * entry.probability
        })
        .sum::<f64>();
    pair_score
        + unpaired_probabilities
            .iter()
            .enumerate()
            .filter(|(index, _)| !paired[*index])
            .map(|(_, probability)| probability)
            .sum::<f64>()
}

/// Exact maximum-weight matching over arbitrary crossing topologies.
///
/// Every canonical RNA pair joins a purine (A/G) to a pyrimidine (C/U), so
/// the candidate graph is bipartite. Padding both sides with zero-weight dummy
/// vertices turns optional pairing into a square assignment problem. The
/// Hungarian primal-dual algorithm then obtains the global optimum in O(n^3).
fn maximum_weight_pair_matching<F>(
    bases: &[u8],
    probabilities: &[PairProbability],
    gain: F,
) -> Vec<PairProbability>
where
    F: Fn(&PairProbability) -> f64,
{
    let purines = bases
        .iter()
        .enumerate()
        .filter_map(|(index, &base)| matches!(base, b'A' | b'G').then_some(index))
        .collect::<Vec<_>>();
    let pyrimidines = bases
        .iter()
        .enumerate()
        .filter_map(|(index, &base)| matches!(base, b'C' | b'U' | b'T').then_some(index))
        .collect::<Vec<_>>();
    let size = purines.len() + pyrimidines.len();
    if size == 0 {
        return Vec::new();
    }
    let purine_slot = purines
        .iter()
        .enumerate()
        .map(|(slot, &position)| (position, slot))
        .collect::<std::collections::HashMap<_, _>>();
    let pyrimidine_slot = pyrimidines
        .iter()
        .enumerate()
        .map(|(slot, &position)| (position, slot))
        .collect::<std::collections::HashMap<_, _>>();
    let mut weights = vec![vec![0.0; size]; size];
    let mut edges = std::collections::HashMap::<(usize, usize), PairProbability>::new();
    for pair in probabilities {
        let first = pair.i - 1;
        let second = pair.j - 1;
        let slots = purine_slot
            .get(&first)
            .zip(pyrimidine_slot.get(&second))
            .or_else(|| purine_slot.get(&second).zip(pyrimidine_slot.get(&first)));
        let Some((&left, &right)) = slots else {
            continue;
        };
        let weight = gain(pair);
        if weight > 0.0 {
            weights[left][right] = weight;
            edges.insert((left, right), pair.clone());
        }
    }
    let assignment = hungarian_maximize(&weights);
    let mut selected = assignment
        .into_iter()
        .enumerate()
        .filter_map(|(left, right)| {
            (left < purines.len() && right < pyrimidines.len())
                .then(|| edges.get(&(left, right)).cloned())
                .flatten()
        })
        .collect::<Vec<_>>();
    selected.sort_by_key(|pair| (pair.i, pair.j));
    selected
}

fn hungarian_maximize(weights: &[Vec<f64>]) -> Vec<usize> {
    let size = weights.len();
    debug_assert!(weights.iter().all(|row| row.len() == size));
    let mut row_potential = vec![0.0; size + 1];
    let mut column_potential = vec![0.0; size + 1];
    let mut matched_row = vec![0usize; size + 1];
    let mut predecessor = vec![0usize; size + 1];
    for row in 1..=size {
        matched_row[0] = row;
        let mut column = 0usize;
        let mut minimum = vec![f64::INFINITY; size + 1];
        let mut used = vec![false; size + 1];
        loop {
            used[column] = true;
            let current_row = matched_row[column];
            let mut delta = f64::INFINITY;
            let mut next_column = 0usize;
            for candidate_column in 1..=size {
                if used[candidate_column] {
                    continue;
                }
                let reduced = -weights[current_row - 1][candidate_column - 1]
                    - row_potential[current_row]
                    - column_potential[candidate_column];
                if reduced < minimum[candidate_column] {
                    minimum[candidate_column] = reduced;
                    predecessor[candidate_column] = column;
                }
                if minimum[candidate_column] < delta {
                    delta = minimum[candidate_column];
                    next_column = candidate_column;
                }
            }
            for candidate_column in 0..=size {
                if used[candidate_column] {
                    row_potential[matched_row[candidate_column]] += delta;
                    column_potential[candidate_column] -= delta;
                } else {
                    minimum[candidate_column] -= delta;
                }
            }
            column = next_column;
            if matched_row[column] == 0 {
                break;
            }
        }
        loop {
            let previous_column = predecessor[column];
            matched_row[column] = matched_row[previous_column];
            column = previous_column;
            if column == 0 {
                break;
            }
        }
    }
    let mut assignment = vec![0usize; size];
    for column in 1..=size {
        assignment[matched_row[column] - 1] = column - 1;
    }
    assignment
}

fn validate_options(options: &PseudoknotOptions) -> Result<(), RnaError> {
    if !options.threshold.is_finite() || !(0.0..=1.0).contains(&options.threshold) {
        return Err(RnaError::InvalidOption(
            "pseudoknot threshold must be between zero and one".into(),
        ));
    }
    if options.iterations == 0 || options.min_helix == 0 {
        return Err(RnaError::InvalidOption(
            "pseudoknot iterations and min_helix must be positive".into(),
        ));
    }
    if options.max_components == Some(0) || options.max_ensemble_states == Some(0) {
        return Err(RnaError::InvalidOption(
            "pseudoknot max-components and max-ensemble-states must be positive".into(),
        ));
    }
    for (name, value) in [
        ("gamma", options.gamma),
        ("initiation", options.initiation_kcal_mol),
        ("crossing", options.crossing_kcal_mol),
        ("unpaired", options.unpaired_kcal_mol),
        ("evidence weight", options.evidence_weight_kcal_mol),
    ] {
        if !value.is_finite()
            || (name == "gamma" && value <= 0.0)
            || (name == "evidence weight" && value < 0.0)
        {
            return Err(RnaError::InvalidOption(format!(
                "pseudoknot {name} parameter is invalid"
            )));
        }
    }
    Ok(())
}

fn update_best(slot: &mut Option<(usize, f64)>, partner: usize, probability: f64) {
    if slot.is_none_or(|(old_partner, old_probability)| {
        probability > old_probability || (probability == old_probability && partner < old_partner)
    }) {
        *slot = Some((partner, probability));
    }
}

fn remove_short_helices(selected: Vec<PairProbability>, min_helix: usize) -> Vec<PairProbability> {
    if min_helix <= 1 {
        return selected;
    }
    let mut ordered = selected;
    ordered.sort_by_key(|entry| (entry.i, std::cmp::Reverse(entry.j)));
    let mut outer_depth = vec![1usize; ordered.len()];
    for inner in 0..ordered.len() {
        for outer in 0..inner {
            if helix_neighbors(&ordered[outer], &ordered[inner]) {
                outer_depth[inner] = outer_depth[inner].max(outer_depth[outer] + 1);
            }
        }
    }
    let mut inner_depth = vec![1usize; ordered.len()];
    for outer in (0..ordered.len()).rev() {
        for inner in outer + 1..ordered.len() {
            if helix_neighbors(&ordered[outer], &ordered[inner]) {
                inner_depth[outer] = inner_depth[outer].max(inner_depth[inner] + 1);
            }
        }
    }
    ordered
        .into_iter()
        .zip(outer_depth.into_iter().zip(inner_depth))
        .filter_map(|(entry, (outside, inside))| (outside + inside > min_helix).then_some(entry))
        .collect()
}

/// ProbKnot treats one-nucleotide bulges as continuous helical stacking.
fn helix_neighbors(outer: &PairProbability, inner: &PairProbability) -> bool {
    if outer.i >= inner.i || inner.j >= outer.j {
        return false;
    }
    let left = inner.i - outer.i;
    let right = outer.j - inner.j;
    matches!((left, right), (1, 1) | (1, 2) | (2, 1))
}

#[derive(Clone, Debug)]
struct HelixCandidate {
    pairs: Vec<(usize, usize)>,
    energy_kcal_mol: f64,
    stack_kcal_mol: f64,
    terminal_kcal_mol: f64,
    probability_sum: f64,
}

#[derive(Clone, Debug)]
struct HTypeCore {
    pairs: Vec<PairProbability>,
    score_kcal_mol: f64,
}

#[derive(Clone, Debug)]
struct PseudoknotState {
    pairs: Vec<PairProbability>,
    score_kcal_mol: f64,
    component_count: usize,
    last_core: usize,
}

/// Enumerate density-2 H-type pseudoknot cores directly from the sequence.
///
/// This is an independent implementation of the physical stem ordering
/// `A5 < B5 < A3 < B3`. Each stem is a contiguous antiparallel helix scored
/// with the same nearest-neighbour stack and terminal terms as the planar
/// model. The three connecting regions use explicit initiation, topology, and
/// unpaired penalties from [`PseudoknotOptions`]. Keeping this stage separate
/// from the marginal decoder matters: two crossing stems are mutually
/// exclusive in a planar partition function even when their joint
/// pseudoknotted state is favourable.
fn h_type_cores(
    bases: &[u8],
    probabilities: &[PairProbability],
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> Vec<HTypeCore> {
    let probability = probabilities
        .iter()
        .map(|entry| ((entry.i - 1, entry.j - 1), entry.probability))
        .collect::<std::collections::HashMap<_, _>>();
    let n = bases.len();
    let mut helices = Vec::new();
    for i in 0..n {
        for j in i + 3..n {
            if !model.can_pair(bases[i], bases[j]) {
                continue;
            }
            let mut maximum = 1usize;
            while i + maximum < j.saturating_sub(maximum)
                && model.can_pair(bases[i + maximum], bases[j - maximum])
            {
                maximum += 1;
            }
            if maximum < 2 {
                continue;
            }
            let mut stack = 0.0;
            let mut probability_sum = probability.get(&(i, j)).copied().unwrap_or(0.0);
            for length in 2..=maximum {
                let inner = length - 1;
                stack += model.stack_energy(
                    bases[i + inner - 1],
                    bases[j - inner + 1],
                    bases[i + inner],
                    bases[j - inner],
                );
                probability_sum += probability
                    .get(&(i + inner, j - inner))
                    .copied()
                    .unwrap_or(0.0);
                let terminal = model.terminal_pair_energy(bases[i], bases[j])
                    + model.terminal_pair_energy(bases[j - inner], bases[i + inner]);
                let pairs = (0..length).map(|offset| (i + offset, j - offset)).collect();
                helices.push(HelixCandidate {
                    pairs,
                    energy_kcal_mol: stack + terminal,
                    stack_kcal_mol: stack,
                    terminal_kcal_mol: terminal,
                    probability_sum,
                });
            }
        }
    }
    // Retain the complete helix set. Ordering is total and deterministic so
    // native and WASM results remain bitwise reproducible.
    helices.sort_by(|a, b| {
        helix_rank(a, options.evidence_weight_kcal_mol)
            .total_cmp(&helix_rank(b, options.evidence_weight_kcal_mol))
            .then_with(|| b.pairs.len().cmp(&a.pairs.len()))
            .then_with(|| a.pairs.cmp(&b.pairs))
    });
    if let Some(limit) = options.max_ensemble_states {
        helices.truncate(limit);
    }
    let mut helices_by_outer_left = vec![Vec::new(); n];
    for helix in &helices {
        helices_by_outer_left[helix.pairs[0].0].push(helix);
    }
    let mut cores = Vec::new();
    for a in &helices {
        let &(a_left_end, a_right_start) = a
            .pairs
            .last()
            .expect("helix candidates contain at least two pairs");
        for candidates in helices_by_outer_left
            .iter()
            .take(a_right_start)
            .skip(a_left_end + 1)
        {
            for &b in candidates {
                let Some(unpaired) = h_type_unpaired_ordered(a, b) else {
                    continue;
                };
                let probability_bonus =
                    options.evidence_weight_kcal_mol * (a.probability_sum + b.probability_sum);
                // DP09 scales only stacks that span a band.  Terminal pair
                // terms are ordinary Turner contributions and must not be
                // multiplied by the 0.89 spanning-stack factor.
                let score = 0.89 * (a.stack_kcal_mol + b.stack_kcal_mol)
                    + a.terminal_kcal_mol
                    + b.terminal_kcal_mol
                    + options.initiation_kcal_mol
                    + 2.0 * options.crossing_kcal_mol
                    + options.unpaired_kcal_mol * unpaired as f64
                    - probability_bonus;
                let mut pairs = a
                    .pairs
                    .iter()
                    .chain(&b.pairs)
                    .map(|&(i, j)| PairProbability {
                        i: i + 1,
                        j: j + 1,
                        probability: probability.get(&(i, j)).copied().unwrap_or(0.0),
                    })
                    .collect::<Vec<_>>();
                pairs.sort_by_key(|entry| (entry.i, entry.j));
                cores.push(HTypeCore {
                    pairs,
                    score_kcal_mol: score,
                });
            }
        }
    }
    cores.sort_by(|a, b| {
        a.score_kcal_mol
            .total_cmp(&b.score_kcal_mol)
            .then_with(|| b.pairs.len().cmp(&a.pairs.len()))
            .then_with(|| {
                a.pairs
                    .iter()
                    .map(|pair| (pair.i, pair.j))
                    .cmp(b.pairs.iter().map(|pair| (pair.i, pair.j)))
            })
    });
    let mut seen = HashSet::new();
    cores.retain(|core| {
        seen.insert(
            core.pairs
                .iter()
                .map(|pair| (pair.i, pair.j))
                .collect::<Vec<_>>(),
        )
    });
    if let Some(limit) = options.max_ensemble_states {
        cores.truncate(limit);
    }
    cores
}

fn helix_rank(helix: &HelixCandidate, evidence_weight_kcal_mol: f64) -> f64 {
    helix.energy_kcal_mol - evidence_weight_kcal_mol * helix.probability_sum
}

fn h_type_unpaired_ordered(first: &HelixCandidate, second: &HelixCandidate) -> Option<usize> {
    let &(_, a_right_end) = first.pairs.first()?;
    let &(a_left_end, a_right_start) = first.pairs.last()?;
    let &(b_left_start, _) = second.pairs.first()?;
    let &(b_left_end, b_right_start) = second.pairs.last()?;
    if a_left_end < b_left_start && b_left_end < a_right_start && a_right_end < b_right_start {
        Some(
            (b_left_start - a_left_end - 1)
                + (a_right_start - b_left_end - 1)
                + (b_right_start - a_right_end - 1),
        )
    } else {
        None
    }
}

/// Form a finite thermodynamic ensemble from mutually nucleotide-disjoint
/// H-type components.  Enumerating components instead of treating every core
/// as a globally exclusive state permits multiple pseudoknots and arbitrary
/// crossing-graph bracket depth while retaining exact Boltzmann sums over the
/// complete enumerated component set (unless the caller requests a limit).
fn combine_compatible_cores(
    cores: &[HTypeCore],
    options: &PseudoknotOptions,
) -> Vec<PseudoknotState> {
    let mut all = cores
        .iter()
        .enumerate()
        .map(|(index, core)| PseudoknotState {
            pairs: core.pairs.clone(),
            score_kcal_mol: core.score_kcal_mol,
            component_count: 1,
            last_core: index,
        })
        .collect::<Vec<_>>();
    let mut frontier = cores
        .iter()
        .enumerate()
        .map(|(index, core)| PseudoknotState {
            pairs: core.pairs.clone(),
            score_kcal_mol: core.score_kcal_mol,
            component_count: 1,
            last_core: index,
        })
        .collect::<Vec<_>>();
    sort_and_deduplicate_states(&mut frontier, options.max_ensemble_states);

    let mut component_count = 2usize;
    while options
        .max_components
        .is_none_or(|maximum| component_count <= maximum)
    {
        let mut next = Vec::new();
        for state in &frontier {
            let occupied = state
                .pairs
                .iter()
                .flat_map(|pair| [pair.i, pair.j])
                .collect::<HashSet<_>>();
            for (index, core) in cores.iter().enumerate().skip(state.last_core + 1) {
                if core
                    .pairs
                    .iter()
                    .any(|pair| occupied.contains(&pair.i) || occupied.contains(&pair.j))
                {
                    continue;
                }
                let mut pairs = state.pairs.clone();
                pairs.extend(core.pairs.iter().cloned());
                pairs.sort_by_key(|pair| (pair.i, pair.j));
                next.push(PseudoknotState {
                    pairs,
                    score_kcal_mol: state.score_kcal_mol + core.score_kcal_mol,
                    component_count,
                    last_core: index,
                });
            }
        }
        if next.is_empty() {
            break;
        }
        sort_and_deduplicate_states(&mut next, options.max_ensemble_states);
        all.extend(next.iter().cloned());
        frontier = next;
        component_count += 1;
    }
    sort_and_deduplicate_states(&mut all, options.max_ensemble_states);
    all
}

fn sort_and_deduplicate_states(states: &mut Vec<PseudoknotState>, limit: Option<usize>) {
    states.sort_by(|a, b| {
        a.score_kcal_mol
            .total_cmp(&b.score_kcal_mol)
            .then_with(|| b.pairs.len().cmp(&a.pairs.len()))
            .then_with(|| {
                a.pairs
                    .iter()
                    .map(|pair| (pair.i, pair.j))
                    .cmp(b.pairs.iter().map(|pair| (pair.i, pair.j)))
            })
    });
    let mut seen = HashSet::new();
    states.retain(|state| {
        seen.insert(
            state
                .pairs
                .iter()
                .map(|pair| (pair.i, pair.j))
                .collect::<Vec<_>>(),
        )
    });
    if let Some(limit) = limit {
        states.truncate(limit);
    }
}

struct RestrictedEnsemble {
    mfe_state: Option<PseudoknotState>,
    mfe_structure: String,
    mfe_energy_kcal_mol: f64,
    ensemble_free_energy_kcal_mol: f64,
    partition_function: f64,
    log_partition_function: f64,
    pair_probabilities: Vec<PairProbability>,
    unpaired_probabilities: Vec<f64>,
    centroid_structure: String,
    centroid_distance: f64,
    mea_structure: String,
    mea_score: f64,
    state_count: usize,
    state_count_exact: String,
}

impl From<RestrictedEnsemble> for ExactArbitraryEnsemble {
    fn from(value: RestrictedEnsemble) -> Self {
        Self {
            mfe_structure: value.mfe_structure,
            mfe_energy_kcal_mol: value.mfe_energy_kcal_mol,
            ensemble_free_energy_kcal_mol: value.ensemble_free_energy_kcal_mol,
            partition_function: value.partition_function,
            log_partition_function: value.log_partition_function,
            pair_probabilities: value.pair_probabilities,
            unpaired_probabilities: value.unpaired_probabilities,
            centroid_structure: value.centroid_structure,
            centroid_distance: value.centroid_distance,
            mea_structure: value.mea_structure,
            mea_score: value.mea_score,
            state_count: value.state_count,
            state_count_exact: value.state_count_exact,
            state_space_complete: true,
            model:
                "complete canonical matching ensemble with Ribon generalized DP09 diagnostic energy",
        }
    }
}

/// Exhaustively enumerate every canonical matching, without a crossing or
/// topology restriction, and evaluate the resulting finite Boltzmann sum.
///
/// This deliberately remains exponential: exact arbitrary-pseudoknot
/// thermodynamics is NP-hard, and silently pruning the state space would make
/// the returned partition function ill-defined.  The polynomial H-type
/// interval ensemble is the separate default path.
fn exact_arbitrary_matching_ensemble(
    bases: &[u8],
    source_probabilities: &[PairProbability],
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> Result<RestrictedEnsemble, RnaError> {
    let evidence = source_probabilities
        .iter()
        .map(|pair| ((pair.i, pair.j), pair.probability))
        .collect::<std::collections::HashMap<_, _>>();
    let mut states = Vec::new();
    let mut occupied = vec![false; bases.len()];
    let mut pairs = Vec::new();
    enumerate_arbitrary_matchings(
        bases,
        min_loop,
        model,
        &evidence,
        options,
        0,
        &mut occupied,
        &mut pairs,
        &mut states,
    );
    sort_and_deduplicate_states(&mut states, None);
    restricted_ensemble(bases.len(), &states, temperature_celsius, options.gamma)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_arbitrary_matchings(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    evidence: &std::collections::HashMap<(usize, usize), f64>,
    options: &PseudoknotOptions,
    mut index: usize,
    occupied: &mut [bool],
    pairs: &mut Vec<PairProbability>,
    states: &mut Vec<PseudoknotState>,
) {
    while index < bases.len() && occupied[index] {
        index += 1;
    }
    if index == bases.len() {
        if !pairs.is_empty() {
            let breakdown = decoded_diagnostic_energy(bases, pairs, model, options);
            let probability_sum = pairs.iter().map(|pair| pair.probability).sum::<f64>();
            states.push(PseudoknotState {
                pairs: pairs.clone(),
                score_kcal_mol: breakdown.total_kcal_mol
                    - options.evidence_weight_kcal_mol * probability_sum,
                component_count: crossing_component_count(pairs),
                last_core: 0,
            });
        }
        return;
    }

    // The smallest unprocessed base is unpaired.
    enumerate_arbitrary_matchings(
        bases,
        min_loop,
        model,
        evidence,
        options,
        index + 1,
        occupied,
        pairs,
        states,
    );

    let first_partner = index.saturating_add(min_loop).saturating_add(1);
    for partner in first_partner..bases.len() {
        if occupied[partner] || !model.can_pair(bases[index], bases[partner]) {
            continue;
        }
        occupied[index] = true;
        occupied[partner] = true;
        pairs.push(PairProbability {
            i: index + 1,
            j: partner + 1,
            probability: evidence
                .get(&(index + 1, partner + 1))
                .copied()
                .unwrap_or(0.0),
        });
        enumerate_arbitrary_matchings(
            bases,
            min_loop,
            model,
            evidence,
            options,
            index + 1,
            occupied,
            pairs,
            states,
        );
        pairs.pop();
        occupied[index] = false;
        occupied[partner] = false;
    }
}

fn crossing_component_count(pairs: &[PairProbability]) -> usize {
    let crossing = pairs
        .iter()
        .enumerate()
        .filter(|(index, pair)| {
            pairs
                .iter()
                .enumerate()
                .any(|(other_index, other)| index != &other_index && crosses(pair, other))
        })
        .map(|(index, _)| index)
        .collect::<HashSet<_>>();
    if crossing.is_empty() {
        return 0;
    }
    let mut seen = HashSet::new();
    let mut components = 0usize;
    for &start in &crossing {
        if !seen.insert(start) {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        while let Some(index) = stack.pop() {
            for &other in &crossing {
                if !seen.contains(&other) && crosses(&pairs[index], &pairs[other]) {
                    seen.insert(other);
                    stack.push(other);
                }
            }
        }
    }
    components
}

fn restricted_ensemble(
    length: usize,
    states: &[PseudoknotState],
    temperature_celsius: f64,
    gamma: f64,
) -> Result<RestrictedEnsemble, RnaError> {
    const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let log_weights = states
        .iter()
        .map(|state| -state.score_kcal_mol / rt)
        .collect::<Vec<_>>();
    let maximum = log_weights.iter().copied().fold(0.0, f64::max);
    let scaled_sum = (-maximum).exp()
        + log_weights
            .iter()
            .map(|weight| (weight - maximum).exp())
            .sum::<f64>();
    let log_partition_function = maximum + scaled_sum.ln();

    let mut pair_mass = std::collections::HashMap::<(usize, usize), f64>::new();
    for (state, log_weight) in states.iter().zip(&log_weights) {
        let probability = (*log_weight - log_partition_function).exp();
        for pair in &state.pairs {
            *pair_mass.entry((pair.i, pair.j)).or_default() += probability;
        }
    }
    let mut pair_probabilities = pair_mass
        .into_iter()
        .map(|((i, j), probability)| PairProbability {
            i,
            j,
            probability: probability.clamp(0.0, 1.0),
        })
        .collect::<Vec<_>>();
    pair_probabilities.sort_by_key(|pair| (pair.i, pair.j));
    let mut unpaired_probabilities = vec![1.0; length];
    for pair in &pair_probabilities {
        unpaired_probabilities[pair.i - 1] -= pair.probability;
        unpaired_probabilities[pair.j - 1] -= pair.probability;
    }
    for probability in &mut unpaired_probabilities {
        *probability = probability.clamp(0.0, 1.0);
    }

    let open_structure = ".".repeat(length);
    let (mfe_structure, mfe_energy_kcal_mol) = states
        .first()
        .filter(|state| state.score_kcal_mol < 0.0)
        .map(|state| Ok((state_structure(length, state)?, state.score_kcal_mol)))
        .transpose()?
        .unwrap_or_else(|| (open_structure.clone(), 0.0));

    let total_pair_probability = pair_probabilities
        .iter()
        .map(|pair| pair.probability)
        .sum::<f64>();
    let probability = pair_probabilities
        .iter()
        .map(|pair| ((pair.i, pair.j), pair.probability))
        .collect::<std::collections::HashMap<_, _>>();
    let mut centroid_structure = open_structure.clone();
    let mut centroid_distance = total_pair_probability;
    let mut centroid_state = None;
    let mut mea_structure = open_structure.clone();
    let mut mea_score = unpaired_probabilities.iter().sum::<f64>();
    let mut mea_state = None;
    for state in states {
        let pair_sum = state
            .pairs
            .iter()
            .map(|pair| probability.get(&(pair.i, pair.j)).copied().unwrap_or(0.0))
            .sum::<f64>();
        let distance = total_pair_probability + state.pairs.len() as f64 - 2.0 * pair_sum;
        let occupied = state
            .pairs
            .iter()
            .flat_map(|pair| [pair.i - 1, pair.j - 1])
            .collect::<HashSet<_>>();
        let score = 2.0 * gamma * pair_sum
            + unpaired_probabilities
                .iter()
                .enumerate()
                .filter(|(index, _)| !occupied.contains(index))
                .map(|(_, value)| value)
                .sum::<f64>();
        if distance < centroid_distance {
            centroid_distance = distance;
            centroid_state = Some(state);
        }
        if score > mea_score {
            mea_score = score;
            mea_state = Some(state);
        }
    }
    if let Some(state) = centroid_state {
        centroid_structure = state_structure(length, state)?;
    }
    if let Some(state) = mea_state {
        mea_structure = state_structure(length, state)?;
    }

    Ok(RestrictedEnsemble {
        mfe_state: states
            .first()
            .filter(|state| state.score_kcal_mol < 0.0)
            .cloned(),
        mfe_structure,
        mfe_energy_kcal_mol,
        ensemble_free_energy_kcal_mol: -rt * log_partition_function,
        partition_function: if log_partition_function < f64::MAX.ln() {
            log_partition_function.exp()
        } else {
            f64::MAX
        },
        log_partition_function,
        pair_probabilities,
        unpaired_probabilities,
        centroid_structure,
        centroid_distance,
        mea_structure,
        mea_score,
        state_count: states.len() + 1,
        state_count_exact: (states.len() + 1).to_string(),
    })
}

/// Exact Boltzmann ensemble of independent H-type components.
///
/// Independence is defined topologically: component spans may not overlap.
/// This prevents two nominal “components” from forming one higher-order,
/// interleaved pseudoknot. Once cores are ordered by their right endpoint,
/// the complete ensemble is a weighted interval family. The forward/reverse
/// dynamic program below therefore obtains the partition function and every
/// core marginal without materialising an exponential number of states.
fn interval_component_ensemble(
    length: usize,
    cores: &[HTypeCore],
    temperature_celsius: f64,
    gamma: f64,
    max_components: Option<usize>,
) -> Result<RestrictedEnsemble, RnaError> {
    const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut ordered = cores.iter().collect::<Vec<_>>();
    ordered.sort_by(|a, b| {
        core_span(a)
            .1
            .cmp(&core_span(b).1)
            .then_with(|| core_span(a).0.cmp(&core_span(b).0))
            .then_with(|| a.score_kcal_mol.total_cmp(&b.score_kcal_mol))
            .then_with(|| {
                a.pairs
                    .iter()
                    .map(|pair| (pair.i, pair.j))
                    .cmp(b.pairs.iter().map(|pair| (pair.i, pair.j)))
            })
    });
    let ends = ordered
        .iter()
        .map(|core| core_span(core).1)
        .collect::<Vec<_>>();
    let predecessors = ordered
        .iter()
        .enumerate()
        .map(|(index, core)| {
            let start = core_span(core).0;
            ends[..index].partition_point(|&end| end < start)
        })
        .collect::<Vec<_>>();
    let mut maximum_components = vec![0usize; ordered.len() + 1];
    for index in 1..=ordered.len() {
        maximum_components[index] =
            maximum_components[index - 1].max(maximum_components[predecessors[index - 1]] + 1);
    }
    let component_limit = max_components
        .unwrap_or(maximum_components[ordered.len()])
        .min(maximum_components[ordered.len()]);

    let mut forward = vec![vec![f64::NEG_INFINITY; component_limit + 1]; ordered.len() + 1];
    forward[0][0] = 0.0;
    let log_weights = ordered
        .iter()
        .map(|core| -core.score_kcal_mol / rt)
        .collect::<Vec<_>>();
    for index in 1..=ordered.len() {
        let predecessor = predecessors[index - 1];
        for count in 0..=component_limit {
            let skip = forward[index - 1][count];
            let take = if count > 0 {
                forward[predecessor][count - 1] + log_weights[index - 1]
            } else {
                f64::NEG_INFINITY
            };
            forward[index][count] = log_add_exp(skip, take);
        }
    }
    let log_partition_function = forward[ordered.len()]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, log_add_exp);

    // Reverse-mode differentiation of the interval recurrence gives each
    // take-edge marginal. Summing those marginals for a base pair is exact,
    // because overlapping spans (and hence repeated-pair cores) are mutually
    // exclusive in this ensemble.
    let mut outside = vec![vec![f64::NEG_INFINITY; component_limit + 1]; ordered.len() + 1];
    for count in 0..=component_limit {
        if forward[ordered.len()][count].is_finite() {
            outside[ordered.len()][count] = 0.0;
        }
    }
    let mut core_log_mass = vec![f64::NEG_INFINITY; ordered.len()];
    for index in (1..=ordered.len()).rev() {
        let predecessor = predecessors[index - 1];
        for count in 0..=component_limit {
            let suffix = outside[index][count];
            if !suffix.is_finite() {
                continue;
            }
            outside[index - 1][count] = log_add_exp(outside[index - 1][count], suffix);
            if count > 0 && forward[predecessor][count - 1].is_finite() {
                outside[predecessor][count - 1] = log_add_exp(
                    outside[predecessor][count - 1],
                    suffix + log_weights[index - 1],
                );
                core_log_mass[index - 1] = log_add_exp(
                    core_log_mass[index - 1],
                    forward[predecessor][count - 1] + log_weights[index - 1] + suffix,
                );
            }
        }
    }

    let mut pair_mass = std::collections::HashMap::<(usize, usize), f64>::new();
    for (core, log_mass) in ordered.iter().zip(core_log_mass) {
        let probability = (log_mass - log_partition_function).exp();
        for pair in &core.pairs {
            *pair_mass.entry((pair.i, pair.j)).or_default() += probability;
        }
    }
    let mut pair_probabilities = pair_mass
        .into_iter()
        .map(|((i, j), probability)| PairProbability {
            i,
            j,
            probability: probability.clamp(0.0, 1.0),
        })
        .collect::<Vec<_>>();
    pair_probabilities.sort_by_key(|pair| (pair.i, pair.j));
    let mut unpaired_probabilities = vec![1.0; length];
    for pair in &pair_probabilities {
        unpaired_probabilities[pair.i - 1] -= pair.probability;
        unpaired_probabilities[pair.j - 1] -= pair.probability;
    }
    for probability in &mut unpaired_probabilities {
        *probability = probability.clamp(0.0, 1.0);
    }

    let mfe_values = ordered
        .iter()
        .map(|core| core.score_kcal_mol)
        .collect::<Vec<_>>();
    let (mfe_energy_kcal_mol, mfe_indices) =
        select_intervals(&predecessors, &mfe_values, component_limit, false);
    let mfe_state = state_from_interval_indices(&ordered, &mfe_indices);
    let open_structure = ".".repeat(length);
    let mfe_structure = mfe_state
        .as_ref()
        .map(|state| state_structure(length, state))
        .transpose()?
        .unwrap_or_else(|| open_structure.clone());

    let probability = pair_probabilities
        .iter()
        .map(|pair| ((pair.i, pair.j), pair.probability))
        .collect::<std::collections::HashMap<_, _>>();
    let total_pair_probability = pair_probabilities
        .iter()
        .map(|pair| pair.probability)
        .sum::<f64>();
    let centroid_values = ordered
        .iter()
        .map(|core| {
            core.pairs.len() as f64
                - 2.0
                    * core
                        .pairs
                        .iter()
                        .map(|pair| probability.get(&(pair.i, pair.j)).copied().unwrap_or(0.0))
                        .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let (centroid_increment, centroid_indices) =
        select_intervals(&predecessors, &centroid_values, component_limit, false);
    let centroid_state = state_from_interval_indices(&ordered, &centroid_indices);
    let centroid_structure = centroid_state
        .as_ref()
        .map(|state| state_structure(length, state))
        .transpose()?
        .unwrap_or_else(|| open_structure.clone());

    let mea_values = ordered
        .iter()
        .map(|core| {
            let pair_reward = 2.0
                * gamma
                * core
                    .pairs
                    .iter()
                    .map(|pair| probability.get(&(pair.i, pair.j)).copied().unwrap_or(0.0))
                    .sum::<f64>();
            let unpaired_loss = core
                .pairs
                .iter()
                .flat_map(|pair| [pair.i - 1, pair.j - 1])
                .map(|position| unpaired_probabilities[position])
                .sum::<f64>();
            pair_reward - unpaired_loss
        })
        .collect::<Vec<_>>();
    let (mea_increment, mea_indices) =
        select_intervals(&predecessors, &mea_values, component_limit, true);
    let mea_state = state_from_interval_indices(&ordered, &mea_indices);
    let mea_structure = mea_state
        .as_ref()
        .map(|state| state_structure(length, state))
        .transpose()?
        .unwrap_or_else(|| open_structure.clone());

    let mut counts = vec![vec!["0".to_string(); component_limit + 1]; ordered.len() + 1];
    counts[0][0] = "1".to_string();
    for index in 1..=ordered.len() {
        let predecessor = predecessors[index - 1];
        for count in 0..=component_limit {
            counts[index][count] = if count == 0 {
                counts[index - 1][count].clone()
            } else {
                add_decimal(&counts[index - 1][count], &counts[predecessor][count - 1])
            };
        }
    }
    let state_count_exact = counts[ordered.len()]
        .iter()
        .fold("0".to_string(), |sum, value| add_decimal(&sum, value));
    let state_count = state_count_exact.parse().unwrap_or(usize::MAX);

    Ok(RestrictedEnsemble {
        mfe_state,
        mfe_structure,
        mfe_energy_kcal_mol,
        ensemble_free_energy_kcal_mol: -rt * log_partition_function,
        partition_function: if log_partition_function < f64::MAX.ln() {
            log_partition_function.exp()
        } else {
            f64::MAX
        },
        log_partition_function,
        pair_probabilities,
        unpaired_probabilities: unpaired_probabilities.clone(),
        centroid_structure,
        centroid_distance: total_pair_probability + centroid_increment,
        mea_structure,
        mea_score: unpaired_probabilities.iter().sum::<f64>() + mea_increment,
        state_count,
        state_count_exact,
    })
}

fn core_span(core: &HTypeCore) -> (usize, usize) {
    (
        core.pairs.iter().map(|pair| pair.i).min().unwrap_or(0),
        core.pairs.iter().map(|pair| pair.j).max().unwrap_or(0),
    )
}

fn log_add_exp(left: f64, right: f64) -> f64 {
    if left == f64::NEG_INFINITY {
        return right;
    }
    if right == f64::NEG_INFINITY {
        return left;
    }
    let maximum = left.max(right);
    maximum + ((left - maximum).exp() + (right - maximum).exp()).ln()
}

/// Select a minimum- or maximum-weight span-disjoint set with an optional
/// cardinality limit. Ties retain the skip edge for deterministic output.
fn select_intervals(
    predecessors: &[usize],
    values: &[f64],
    component_limit: usize,
    maximize: bool,
) -> (f64, Vec<usize>) {
    let unreachable = if maximize {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    };
    let mut dynamic = vec![vec![unreachable; component_limit + 1]; values.len() + 1];
    let mut take = vec![vec![false; component_limit + 1]; values.len() + 1];
    dynamic[0][0] = 0.0;
    for index in 1..=values.len() {
        let predecessor = predecessors[index - 1];
        for count in 0..=component_limit {
            dynamic[index][count] = dynamic[index - 1][count];
            if count == 0 || !dynamic[predecessor][count - 1].is_finite() {
                continue;
            }
            let candidate = dynamic[predecessor][count - 1] + values[index - 1];
            let better = if maximize {
                candidate > dynamic[index][count]
            } else {
                candidate < dynamic[index][count]
            };
            if better {
                dynamic[index][count] = candidate;
                take[index][count] = true;
            }
        }
    }
    let mut best_count = 0usize;
    for count in 1..=component_limit {
        let better = if maximize {
            dynamic[values.len()][count] > dynamic[values.len()][best_count]
        } else {
            dynamic[values.len()][count] < dynamic[values.len()][best_count]
        };
        if better {
            best_count = count;
        }
    }
    let best = dynamic[values.len()][best_count];
    let mut indices = Vec::with_capacity(best_count);
    let mut index = values.len();
    let mut count = best_count;
    while index > 0 && count > 0 {
        if take[index][count] {
            indices.push(index - 1);
            index = predecessors[index - 1];
            count -= 1;
        } else {
            index -= 1;
        }
    }
    indices.reverse();
    (best, indices)
}

fn state_from_interval_indices(
    ordered: &[&HTypeCore],
    indices: &[usize],
) -> Option<PseudoknotState> {
    if indices.is_empty() {
        return None;
    }
    let mut pairs = indices
        .iter()
        .flat_map(|&index| ordered[index].pairs.iter().cloned())
        .collect::<Vec<_>>();
    pairs.sort_by_key(|pair| (pair.i, pair.j));
    Some(PseudoknotState {
        pairs,
        score_kcal_mol: indices
            .iter()
            .map(|&index| ordered[index].score_kcal_mol)
            .sum(),
        component_count: indices.len(),
        last_core: *indices.last().expect("nonempty interval selection"),
    })
}

fn add_decimal(left: &str, right: &str) -> String {
    let mut result = Vec::with_capacity(left.len().max(right.len()) + 1);
    let mut left = left.bytes().rev();
    let mut right = right.bytes().rev();
    let mut carry = 0u8;
    loop {
        let a = left.next().map(|digit| digit - b'0');
        let b = right.next().map(|digit| digit - b'0');
        if a.is_none() && b.is_none() && carry == 0 {
            break;
        }
        let sum = a.unwrap_or(0) + b.unwrap_or(0) + carry;
        result.push(b'0' + sum % 10);
        carry = sum / 10;
    }
    result.reverse();
    String::from_utf8(result).expect("decimal addition emits ASCII digits")
}

fn state_structure(length: usize, state: &PseudoknotState) -> Result<String, RnaError> {
    let levels = color_crossing_graph(&state.pairs)?;
    extended_dot_bracket(length, &state.pairs, &levels)
}

fn crosses(a: &PairProbability, b: &PairProbability) -> bool {
    (a.i < b.i && b.i < a.j && a.j < b.j) || (b.i < a.i && a.i < b.j && b.j < a.j)
}

fn count_crossings(pairs: &[PairProbability]) -> usize {
    let mut count = 0;
    for (index, a) in pairs.iter().enumerate() {
        count += pairs[index + 1..].iter().filter(|b| crosses(a, b)).count();
    }
    count
}

fn color_crossing_graph(pairs: &[PairProbability]) -> Result<Vec<usize>, RnaError> {
    if pairs.is_empty() {
        return Ok(Vec::new());
    }
    let n = pairs.len();
    let mut adjacency = vec![vec![false; n]; n];
    let mut degrees = vec![0usize; n];
    for i in 0..n {
        for j in i + 1..n {
            if crosses(&pairs[i], &pairs[j]) {
                adjacency[i][j] = true;
                adjacency[j][i] = true;
                degrees[i] += 1;
                degrees[j] += 1;
            }
        }
    }

    // A valid <=30-color greedy result needs no optimization. Only a greedy
    // overflow can create a false format failure, so that uncommon boundary
    // case is resolved by exact 30-colorability search below.
    let mut greedy = vec![usize::MAX; n];
    for vertex in 0..n {
        let mut forbidden = vec![false; n];
        for neighbor in 0..vertex {
            if adjacency[vertex][neighbor] {
                forbidden[greedy[neighbor]] = true;
            }
        }
        greedy[vertex] = forbidden
            .iter()
            .position(|&value| !value)
            .expect("n colors always suffice");
    }
    if greedy.iter().copied().max().unwrap_or(0) < 30 {
        return Ok(greedy);
    }

    // DSATUR now answers the exact question imposed by the one-symbol format:
    // whether any 30-coloring exists. There is no search or iteration cap.
    let mut colors = vec![None; n];
    if color_dsatur(&adjacency, &degrees, 30, &mut colors, 0) {
        return Ok(colors
            .into_iter()
            .map(|color| color.expect("all vertices colored"))
            .collect());
    }
    Err(RnaError::InvalidOption(
        "pseudoknot crossing graph is not 30-colorable; the pair table is exact, but the standard one-symbol extended dot-bracket alphabet has only 30 bracket levels"
            .into(),
    ))
}

fn color_dsatur(
    adjacency: &[Vec<bool>],
    degrees: &[usize],
    color_count: usize,
    colors: &mut [Option<usize>],
    colored: usize,
) -> bool {
    if colored == colors.len() {
        return true;
    }
    let vertex = (0..colors.len())
        .filter(|&index| colors[index].is_none())
        .max_by_key(|&index| {
            let mut used = vec![false; color_count];
            for (neighbor, &connected) in adjacency[index].iter().enumerate() {
                if connected {
                    if let Some(color) = colors[neighbor] {
                        used[color] = true;
                    }
                }
            }
            (
                used.into_iter().filter(|&value| value).count(),
                degrees[index],
                std::cmp::Reverse(index),
            )
        })
        .expect("an uncolored vertex exists");
    let mut forbidden = vec![false; color_count];
    for (neighbor, &connected) in adjacency[vertex].iter().enumerate() {
        if connected {
            if let Some(color) = colors[neighbor] {
                forbidden[color] = true;
            }
        }
    }
    let used = colors
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(0, |maximum| maximum + 1);
    let available = (used + 1).min(color_count);
    for (color, &is_forbidden) in forbidden.iter().take(available).enumerate() {
        if is_forbidden {
            continue;
        }
        colors[vertex] = Some(color);
        if color_dsatur(adjacency, degrees, color_count, colors, colored + 1) {
            return true;
        }
        colors[vertex] = None;
    }
    false
}

fn brackets(level: usize) -> Option<(char, char)> {
    match level {
        0 => Some(('(', ')')),
        1 => Some(('[', ']')),
        2 => Some(('{', '}')),
        3 => Some(('<', '>')),
        4..=29 => {
            let upper = char::from_u32('A' as u32 + (level - 4) as u32)?;
            Some((upper, upper.to_ascii_lowercase()))
        }
        _ => None,
    }
}

fn extended_dot_bracket(
    length: usize,
    pairs: &[PairProbability],
    levels: &[usize],
) -> Result<String, RnaError> {
    let mut result = vec!['.'; length];
    for (pair, &level) in pairs.iter().zip(levels) {
        if pair.i == 0 || pair.i >= pair.j || pair.j > length {
            return Err(RnaError::InvalidOption(
                "pseudoknot pair lies outside the sequence".into(),
            ));
        }
        if result[pair.i - 1] != '.' || result[pair.j - 1] != '.' {
            return Err(RnaError::MultiplePartners { position: pair.i });
        }
        let (open, close) = brackets(level).ok_or_else(|| {
            RnaError::InvalidOption("pseudoknot bracket level is unavailable".into())
        })?;
        result[pair.i - 1] = open;
        result[pair.j - 1] = close;
    }
    Ok(result.into_iter().collect())
}

fn decoded_diagnostic_energy(
    bases: &[u8],
    pairs: &[PairProbability],
    model: &EnergyModel,
    options: &PseudoknotOptions,
) -> PseudoknotEnergyBreakdown {
    let lookup: HashSet<(usize, usize)> = pairs
        .iter()
        .map(|entry| (entry.i - 1, entry.j - 1))
        .collect();
    let crossing_pairs = pairs
        .iter()
        .filter(|pair| pairs.iter().any(|other| crosses(pair, other)))
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<HashSet<_>>();
    let mut stack = 0.0;
    let mut terminal = 0.0;
    let mut helices = 0usize;
    for &(i, j) in &lookup {
        if i == 0 || !lookup.contains(&(i - 1, j + 1)) {
            helices += 1;
            terminal += model.terminal_pair_energy(bases[i], bases[j]);
        }
        if lookup.contains(&(i + 1, j.saturating_sub(1))) {
            let energy = model.stack_energy(bases[i], bases[j], bases[i + 1], bases[j - 1]);
            stack += if crossing_pairs.contains(&(i, j)) {
                0.89 * energy
            } else {
                energy
            };
        } else {
            terminal += model.terminal_pair_energy(bases[j], bases[i]);
        }
    }
    let crossing_count = count_crossings(pairs);
    let enclosed_positions: HashSet<usize> = pairs
        .iter()
        .flat_map(|entry| entry.i..entry.j - 1)
        .collect();
    let paired_positions: HashSet<usize> = pairs
        .iter()
        .flat_map(|entry| [entry.i - 1, entry.j - 1])
        .collect();
    let enclosed_unpaired = enclosed_positions
        .iter()
        .filter(|position| !paired_positions.contains(position))
        .count();
    let initiation = options.initiation_kcal_mol * crossing_component_count(pairs) as f64;
    let crossing = options.crossing_kcal_mol
        * crossing_pairs
            .iter()
            .filter(|&&(i, j)| i == 0 || !crossing_pairs.contains(&(i - 1, j + 1)))
            .count() as f64;
    let unpaired = options.unpaired_kcal_mol * enclosed_unpaired as f64;
    PseudoknotEnergyBreakdown {
        total_kcal_mol: stack + terminal + initiation + crossing + unpaired,
        stack_kcal_mol: stack,
        terminal_kcal_mol: terminal,
        initiation_kcal_mol: initiation,
        crossing_kcal_mol: crossing,
        unpaired_kcal_mol: unpaired,
        helix_count: helices,
        crossing_count,
        enclosed_unpaired_count: enclosed_unpaired,
        model: "Ribon generalized pseudoknot diagnostic v1: DP09 Ps/Pb/Pup terms, 0.89 spanning stacks, and RNAstructure 6.6 RNA terminal pairs",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_structure;

    #[test]
    fn mutual_maximum_decoder_emits_valid_crossing_brackets() {
        let probabilities = vec![
            PairProbability {
                i: 1,
                j: 6,
                probability: 0.9,
            },
            PairProbability {
                i: 3,
                j: 8,
                probability: 0.8,
            },
            PairProbability {
                i: 2,
                j: 5,
                probability: 0.7,
            },
            PairProbability {
                i: 4,
                j: 7,
                probability: 0.6,
            },
        ];
        let options = PseudoknotOptions {
            min_helix: 1,
            ..PseudoknotOptions::default()
        };
        let result = decode_pseudoknot(
            "GCGCGCGC",
            &probabilities,
            &[0.1; 8],
            37.0,
            3,
            0,
            1.021,
            &options,
        )
        .unwrap();
        assert!(result.crossing_count > 0);
        let parsed = parse_structure(&result.sequence, &result.structure).unwrap();
        assert_eq!(parsed.pairs.len(), result.pair_count);
        assert!(result.bracket_levels >= 2);
        assert!(result.restricted_state_count > 1);
        assert!(result.restricted_log_partition_function.is_finite());
        for position in 1..=result.sequence.len() {
            let paired = result
                .restricted_pair_probabilities
                .iter()
                .filter(|pair| pair.i == position || pair.j == position)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!(
                (paired + result.restricted_unpaired_probabilities[position - 1] - 1.0).abs()
                    < 1.0e-10
            );
        }
        parse_structure(&result.sequence, &result.restricted_mfe_structure).unwrap();
        parse_structure(&result.sequence, &result.restricted_centroid_structure).unwrap();
        parse_structure(&result.sequence, &result.restricted_mea_structure).unwrap();
    }

    #[test]
    fn deterministic_ties_choose_the_lower_partner() {
        let probabilities = vec![
            PairProbability {
                i: 1,
                j: 5,
                probability: 0.8,
            },
            PairProbability {
                i: 1,
                j: 6,
                probability: 0.8,
            },
            PairProbability {
                i: 2,
                j: 6,
                probability: 0.7,
            },
        ];
        let options = PseudoknotOptions {
            min_helix: 1,
            ..PseudoknotOptions::default()
        };
        let result = decode_pseudoknot(
            "GGAACC",
            &probabilities,
            &[0.5; 6],
            37.0,
            3,
            0,
            1.021,
            &options,
        )
        .unwrap();
        assert!(result.pairs.iter().any(|pair| pair.i == 1 && pair.j == 5));
        assert!(!result.pairs.iter().any(|pair| pair.i == 1 && pair.j == 6));
    }

    #[test]
    fn minimum_helix_filter_keeps_every_pair_in_a_valid_helix() {
        let helix = vec![
            PairProbability {
                i: 2,
                j: 8,
                probability: 0.9,
            },
            PairProbability {
                i: 3,
                j: 7,
                probability: 0.9,
            },
            PairProbability {
                i: 4,
                j: 6,
                probability: 0.9,
            },
        ];
        assert_eq!(remove_short_helices(helix.clone(), 3).len(), 3);
        assert!(remove_short_helices(helix, 4).is_empty());

        let one_base_bulge = vec![
            PairProbability {
                i: 1,
                j: 9,
                probability: 0.9,
            },
            PairProbability {
                i: 3,
                j: 8,
                probability: 0.9,
            },
            PairProbability {
                i: 4,
                j: 7,
                probability: 0.9,
            },
        ];
        assert_eq!(remove_short_helices(one_base_bulge, 3).len(), 3);
    }

    #[test]
    fn malformed_probability_inputs_are_rejected() {
        let invalid = vec![PairProbability {
            i: 1,
            j: 5,
            probability: 1.1,
        }];
        assert!(decode_pseudoknot(
            "GAAAC",
            &invalid,
            &[0.5; 5],
            37.0,
            3,
            0,
            1.021,
            &PseudoknotOptions::default(),
        )
        .is_err());
    }

    #[test]
    fn exact_arbitrary_ensemble_enumerates_the_complete_matching_space() {
        let options = PseudoknotOptions {
            exact_arbitrary_ensemble: true,
            ..PseudoknotOptions::default()
        };
        // With min_loop=3, GAAAC has exactly the open matching and (1,5).
        let result = predict_pseudoknot("GAAAC", 37.0, 3, 0, 1.021, &options).unwrap();
        let exact = result
            .exact_arbitrary_ensemble
            .expect("the explicitly requested ensemble must be returned");
        assert_eq!(exact.state_count, 2);
        assert_eq!(exact.state_count_exact, "2");
        assert!(exact.state_space_complete);
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let pair = PairProbability {
            i: 1,
            j: 5,
            probability: result
                .source_pair_probabilities
                .iter()
                .find(|pair| (pair.i, pair.j) == (1, 5))
                .map_or(0.0, |pair| pair.probability),
        };
        let pair_energy =
            decoded_diagnostic_energy(b"GAAAC", std::slice::from_ref(&pair), &model, &options)
                .total_kcal_mol;
        let rt = 0.001_987_17 * (37.0 + 273.15);
        let expected_log_z = log_add_exp(0.0, -pair_energy / rt);
        assert!((exact.log_partition_function - expected_log_z).abs() < 1.0e-12);
        let probability = exact
            .pair_probabilities
            .iter()
            .find(|pair| (pair.i, pair.j) == (1, 5))
            .unwrap()
            .probability;
        assert!((probability - (-pair_energy / rt - expected_log_z).exp()).abs() < 1.0e-12);
        for position in 0..5 {
            let paired = exact
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == position + 1 || pair.j == position + 1)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!((paired + exact.unpaired_probabilities[position] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn exact_arbitrary_ensemble_includes_crossing_perfect_matchings() {
        // GGCC induces K2,2: empty, four one-edge, and two perfect
        // matchings. One perfect matching is (1,3)(2,4), which crosses.
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let ensemble = exact_arbitrary_matching_ensemble(
            b"GGCC",
            &[],
            37.0,
            0,
            &model,
            &PseudoknotOptions::default(),
        )
        .unwrap();
        assert_eq!(ensemble.state_count, 7);
        assert_eq!(ensemble.state_count_exact, "7");
        assert!(ensemble.log_partition_function.is_finite());
    }

    #[test]
    fn supplied_crossing_structure_has_a_finite_reproducible_energy() {
        let first = evaluate_pseudoknot_structure(
            "GGGAAACCCGGGAAACCC",
            "([......)......]..",
            37.0,
            3,
            2,
            1.021,
            &PseudoknotOptions::default(),
        )
        .unwrap();
        let second = evaluate_pseudoknot_structure(
            "GGGAAACCCGGGAAACCC",
            "([......)......]..",
            37.0,
            3,
            2,
            1.021,
            &PseudoknotOptions::default(),
        )
        .unwrap();
        assert_eq!(first.crossing_count, 1);
        assert_eq!(first.crossing_component_count, 1);
        assert_eq!(first.energy.total_kcal_mol, second.energy.total_kcal_mol);
        assert!(first.energy.total_kcal_mol.is_finite());
    }

    #[test]
    fn compatible_components_coexist_in_the_boltzmann_ensemble() {
        let pair = |i, j| PairProbability {
            i,
            j,
            probability: 0.8,
        };
        let cores = vec![
            HTypeCore {
                pairs: vec![pair(1, 6), pair(3, 8)],
                score_kcal_mol: -1.0,
            },
            HTypeCore {
                pairs: vec![pair(9, 14), pair(11, 16)],
                score_kcal_mol: -2.0,
            },
            HTypeCore {
                pairs: vec![pair(1, 10), pair(4, 12)],
                score_kcal_mol: -4.0,
            },
        ];
        let options = PseudoknotOptions {
            max_components: Some(3),
            max_ensemble_states: Some(64),
            ..PseudoknotOptions::default()
        };
        let states = combine_compatible_cores(&cores, &options);
        let combined = states
            .iter()
            .find(|state| {
                state.component_count == 2
                    && state.pairs.iter().any(|pair| (pair.i, pair.j) == (1, 6))
                    && state.pairs.iter().any(|pair| (pair.i, pair.j) == (9, 14))
            })
            .expect("the two disjoint pseudoknot components must coexist");
        assert_eq!(combined.score_kcal_mol, -3.0);
        let ensemble = restricted_ensemble(16, &states, 37.0, 1.0).unwrap();
        assert!(ensemble.state_count > cores.len() + 1);
        for position in 0..16 {
            let paired = ensemble
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == position + 1 || pair.j == position + 1)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!((paired + ensemble.unpaired_probabilities[position] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn interval_dynamic_program_matches_complete_small_enumeration() {
        let pair = |i, j| PairProbability {
            i,
            j,
            probability: 0.5,
        };
        // The third core conflicts with both shorter-span cores. Thus the
        // complete state set is {}, {A}, {B}, {C}, and {A,B}.
        let cores = vec![
            HTypeCore {
                pairs: vec![pair(1, 6), pair(3, 8)],
                score_kcal_mol: -1.25,
            },
            HTypeCore {
                pairs: vec![pair(9, 14), pair(11, 16)],
                score_kcal_mol: -0.75,
            },
            HTypeCore {
                pairs: vec![pair(1, 10), pair(4, 16)],
                score_kcal_mol: -1.6,
            },
        ];
        let options = PseudoknotOptions {
            max_components: Some(3),
            max_ensemble_states: Some(64),
            ..PseudoknotOptions::default()
        };
        let states = combine_compatible_cores(&cores, &options);
        let enumerated = restricted_ensemble(16, &states, 37.0, 1.0).unwrap();
        let dynamic = interval_component_ensemble(16, &cores, 37.0, 1.0, None).unwrap();
        assert_eq!(dynamic.state_count_exact, "5");
        assert_eq!(dynamic.state_count, enumerated.state_count);
        assert_eq!(dynamic.mfe_structure, enumerated.mfe_structure);
        assert!((dynamic.mfe_energy_kcal_mol - enumerated.mfe_energy_kcal_mol).abs() < 1e-12);
        assert!((dynamic.log_partition_function - enumerated.log_partition_function).abs() < 1e-12);
        assert_eq!(
            dynamic.pair_probabilities.len(),
            enumerated.pair_probabilities.len()
        );
        for (actual, expected) in dynamic
            .pair_probabilities
            .iter()
            .zip(&enumerated.pair_probabilities)
        {
            assert_eq!((actual.i, actual.j), (expected.i, expected.j));
            assert!((actual.probability - expected.probability).abs() < 1e-12);
        }
        assert!((dynamic.centroid_distance - enumerated.centroid_distance).abs() < 1e-12);
        assert!((dynamic.mea_score - enumerated.mea_score).abs() < 1e-12);
    }

    #[test]
    fn interval_state_count_remains_exact_past_machine_integer_range() {
        let pair = |i, j| PairProbability {
            i,
            j,
            probability: 0.5,
        };
        // Sixty-four pairwise span-disjoint choices have 2^64 states.
        let cores = (0..64)
            .map(|index| {
                let offset = 8 * index;
                HTypeCore {
                    pairs: vec![pair(offset + 1, offset + 6), pair(offset + 3, offset + 8)],
                    score_kcal_mol: -0.1,
                }
            })
            .collect::<Vec<_>>();
        let ensemble = interval_component_ensemble(512, &cores, 37.0, 1.0, None).unwrap();
        assert_eq!(ensemble.state_count_exact, "18446744073709551616");
        assert_eq!(ensemble.state_count, usize::MAX);
    }

    #[test]
    fn arbitrary_topology_decoder_is_the_global_bipartite_matching() {
        let pair = |i, j| PairProbability {
            i,
            j,
            probability: 0.5,
        };
        let candidates = vec![pair(1, 3), pair(1, 4), pair(2, 3), pair(2, 4)];
        let selected =
            maximum_weight_pair_matching(b"GGCC", &candidates, |entry| match (entry.i, entry.j) {
                (1, 3) => 10.0,
                (1, 4) => 8.0,
                (2, 3) => 9.0,
                (2, 4) => 1.0,
                _ => unreachable!(),
            });
        assert_eq!(
            selected
                .iter()
                .map(|entry| (entry.i, entry.j))
                .collect::<Vec<_>>(),
            vec![(1, 4), (2, 3)]
        );

        let complete = (1..=3)
            .flat_map(|i| (4..=6).map(move |j| pair(i, j)))
            .collect::<Vec<_>>();
        let weight = |entry: &PairProbability| {
            ((entry.i * 17 + entry.j * 29 + entry.i * entry.j * 7) % 31) as f64 / 10.0
        };
        let actual = maximum_weight_pair_matching(b"GGGCCC", &complete, weight)
            .iter()
            .map(weight)
            .sum::<f64>();
        let mut expected = 0.0f64;
        for mask in 0usize..1usize << complete.len() {
            let mut occupied = [false; 6];
            let mut score = 0.0;
            let mut valid = true;
            for (index, entry) in complete.iter().enumerate() {
                if mask & (1 << index) == 0 {
                    continue;
                }
                if occupied[entry.i - 1] || occupied[entry.j - 1] {
                    valid = false;
                    break;
                }
                occupied[entry.i - 1] = true;
                occupied[entry.j - 1] = true;
                score += weight(entry);
            }
            if valid {
                expected = expected.max(score);
            }
        }
        assert!((actual - expected).abs() < 1e-12);
    }
}
