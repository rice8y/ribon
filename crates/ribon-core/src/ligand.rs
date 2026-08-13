//! Exact ligand-coupled secondary-structure ensemble.
//!
//! Every admitted planar RNA structure is combined with every independent set
//! of compatible ligand sites. Overlapping site intervals cannot be occupied
//! simultaneously. Binding free energies are standard-state values corrected
//! by `-RT ln(c / 1 M)`. No structure, site, or microstate cap is applied.

use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::decode::{decode_centroid_with_constraints, decode_mea_with_constraints};
use crate::energy::EnergyModel;
use crate::exact_enumeration::for_each_noncrossing_structure;
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::{Deserialize, Serialize};

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct LigandMotif {
    pub id: String,
    pub start: usize,
    pub sequence: String,
    pub structure: String,
    pub standard_binding_energy_kcal_mol: f64,
    pub concentration_molar: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LigandOccupancy {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub sequence_pattern: String,
    pub structure_pattern: String,
    pub standard_binding_energy_kcal_mol: f64,
    pub concentration_molar: f64,
    pub effective_binding_energy_kcal_mol: f64,
    pub sequence_matches: bool,
    pub occupancy_probability: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LigandEnsembleResult {
    pub sequence: String,
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub mfe_bound_motif_ids: Vec<String>,
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub centroid_structure: String,
    pub centroid_score: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub mea_gamma: f64,
    pub motifs: Vec<LigandOccupancy>,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub constraints: ConstraintSummary,
    pub state_space_complete: bool,
    pub model: &'static str,
    pub algorithm: &'static str,
    pub time_complexity: &'static str,
    pub space_complexity: &'static str,
}

struct CompiledMotif {
    input: LigandMotif,
    start: usize,
    end: usize,
    partner: Vec<Option<usize>>,
    sequence_matches: bool,
    effective_energy: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn ligand_ensemble_exact(
    sequence: &str,
    motifs: &[LigandMotif],
    min_loop: usize,
    gamma: f64,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<LigandEnsembleResult, RnaError> {
    if motifs.is_empty() {
        return Err(RnaError::InvalidOption(
            "ligand ensemble requires at least one motif".into(),
        ));
    }
    if !gamma.is_finite() || gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "ligand MEA gamma must be finite and positive".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    let rt = GAS_CONSTANT_KCAL * (model.temperature_celsius() + 273.15);
    let compiled = compile_motifs(&sequence, motifs, rt)?;
    let n = sequence.len();
    let mut log_z = NEG_INF;
    let mut pair_log_mass = vec![vec![NEG_INF; n]; n];
    let mut motif_log_mass = vec![NEG_INF; compiled.len()];
    let mut mfe_energy = f64::INFINITY;
    let mut mfe_structure = String::new();
    let mut mfe_bound = Vec::new();

    for_each_noncrossing_structure(sequence.as_bytes(), min_loop, model, constraints, |pairs| {
        let structure = pairs_to_dot_bracket(n, pairs);
        let structural_energy = model
            .evaluate_with_constraints(&sequence, &structure, min_loop, constraints)?
            .total_kcal_mol;
        let mut partner = vec![None; n];
        for &(i, j) in pairs {
            partner[i] = Some(j);
            partner[j] = Some(i);
        }
        let compatible = compiled
            .iter()
            .enumerate()
            .filter_map(|(index, motif)| {
                (motif.sequence_matches && motif_matches_structure(motif, &partner))
                    .then_some(index)
            })
            .collect::<Vec<_>>();
        enumerate_binding_microstates(
            &compatible,
            &compiled,
            0,
            &mut Vec::new(),
            structural_energy,
            &mut |bound, energy| {
                let log_weight = -energy / rt;
                log_z = log_add(log_z, log_weight);
                for &(i, j) in pairs {
                    pair_log_mass[i][j] = log_add(pair_log_mass[i][j], log_weight);
                }
                for &motif in bound {
                    motif_log_mass[motif] = log_add(motif_log_mass[motif], log_weight);
                }
                let ids = bound
                    .iter()
                    .map(|&index| compiled[index].input.id.clone())
                    .collect::<Vec<_>>();
                if energy + 1.0e-12 < mfe_energy
                    || ((energy - mfe_energy).abs() <= 1.0e-12
                        && (structure.as_str(), ids.as_slice())
                            < (mfe_structure.as_str(), mfe_bound.as_slice()))
                {
                    mfe_energy = energy;
                    mfe_structure.clone_from(&structure);
                    mfe_bound = ids;
                }
                Ok(())
            },
        )
    })?;
    if log_z == NEG_INF || !mfe_energy.is_finite() {
        return Err(RnaError::InvalidOption(
            "constraints admit no ligand-coupled structure".into(),
        ));
    }

    let mut pair_probabilities = Vec::new();
    let mut paired_mass = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            if pair_log_mass[i][j] == NEG_INF {
                continue;
            }
            let probability = (pair_log_mass[i][j] - log_z).exp().clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                pair_probabilities.push(PairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                });
                paired_mass[i] += probability;
                paired_mass[j] += probability;
            }
        }
    }
    let unpaired_probabilities = paired_mass
        .into_iter()
        .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
        .collect::<Vec<_>>();
    let (centroid_structure, centroid_score) = decode_centroid_with_constraints(
        sequence.as_bytes(),
        min_loop,
        &pair_probabilities,
        constraints,
        model,
    )?;
    let (mea_structure, mea_score) = decode_mea_with_constraints(
        sequence.as_bytes(),
        min_loop,
        gamma,
        &pair_probabilities,
        &unpaired_probabilities,
        constraints,
        model,
    )?;
    let motifs = compiled
        .into_iter()
        .enumerate()
        .map(|(index, motif)| LigandOccupancy {
            id: motif.input.id,
            start: motif.start + 1,
            end: motif.end,
            sequence_pattern: motif.input.sequence,
            structure_pattern: motif.input.structure,
            standard_binding_energy_kcal_mol: motif.input.standard_binding_energy_kcal_mol,
            concentration_molar: motif.input.concentration_molar,
            effective_binding_energy_kcal_mol: motif.effective_energy,
            sequence_matches: motif.sequence_matches,
            occupancy_probability: if motif_log_mass[index] == NEG_INF {
                0.0
            } else {
                (motif_log_mass[index] - log_z).exp().clamp(0.0, 1.0)
            },
        })
        .collect();

    Ok(LigandEnsembleResult {
        sequence,
        mfe_structure,
        mfe_energy_kcal_mol: mfe_energy,
        mfe_bound_motif_ids: mfe_bound,
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        pair_probabilities,
        unpaired_probabilities,
        centroid_structure,
        centroid_score,
        mea_structure,
        mea_score,
        mea_gamma: gamma,
        motifs,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        constraints: constraints.summary(),
        state_space_complete: true,
        model: model.model_name(),
        algorithm: "complete planar RNA structures times complete non-overlapping ligand binding microstates",
        time_complexity: "exponential in RNA structures and compatible overlapping ligand sites",
        space_complexity: "quadratic probability tables plus enumeration stacks",
    })
}

fn compile_motifs(
    sequence: &str,
    motifs: &[LigandMotif],
    rt: f64,
) -> Result<Vec<CompiledMotif>, RnaError> {
    let mut ids = std::collections::HashSet::new();
    motifs
        .iter()
        .cloned()
        .map(|input| {
            if input.id.is_empty() || !ids.insert(input.id.clone()) {
                return Err(RnaError::InvalidOption(
                    "ligand motif ids must be non-empty and unique".into(),
                ));
            }
            if input.start == 0 {
                return Err(RnaError::InvalidOption(
                    "ligand motif start is one-based and must be positive".into(),
                ));
            }
            if !input.standard_binding_energy_kcal_mol.is_finite()
                || !input.concentration_molar.is_finite()
                || input.concentration_molar <= 0.0
            {
                return Err(RnaError::InvalidOption(
                    "ligand energy must be finite and concentration must be finite and positive"
                        .into(),
                ));
            }
            let pattern = normalize_sequence(&input.sequence)?;
            let (structure, partner) = parse_local_structure(&input.structure)?;
            if pattern.len() != structure.len() {
                return Err(RnaError::LengthMismatch {
                    sequence: pattern.len(),
                    structure: structure.len(),
                });
            }
            let start = input.start - 1;
            let end = start
                .checked_add(pattern.len())
                .ok_or_else(|| RnaError::InvalidOption("ligand motif interval overflow".into()))?;
            if end > sequence.len() {
                return Err(RnaError::InvalidOption(format!(
                    "ligand motif {:?} interval exceeds sequence length {}",
                    input.id,
                    sequence.len()
                )));
            }
            let sequence_matches = pattern
                .bytes()
                .zip(sequence.as_bytes()[start..end].iter().copied())
                .all(|(symbol, base)| iupac_matches(symbol, base));
            let effective_energy =
                input.standard_binding_energy_kcal_mol - rt * input.concentration_molar.ln();
            Ok(CompiledMotif {
                input: LigandMotif {
                    sequence: pattern,
                    structure,
                    ..input
                },
                start,
                end,
                partner,
                sequence_matches,
                effective_energy,
            })
        })
        .collect()
}

fn parse_local_structure(structure: &str) -> Result<(String, Vec<Option<usize>>), RnaError> {
    let symbols = structure
        .chars()
        .filter(|symbol| !symbol.is_whitespace())
        .collect::<Vec<_>>();
    if symbols.is_empty() {
        return Err(RnaError::InvalidOption(
            "ligand structure pattern is empty".into(),
        ));
    }
    let mut partner = vec![None; symbols.len()];
    let mut stack = Vec::new();
    for (position, symbol) in symbols.iter().copied().enumerate() {
        match symbol {
            '.' => {}
            '(' => stack.push(position),
            ')' => {
                let opening = stack.pop().ok_or(RnaError::UnmatchedClosing {
                    position: position + 1,
                    symbol,
                })?;
                partner[opening] = Some(position);
                partner[position] = Some(opening);
            }
            _ => {
                return Err(RnaError::InvalidStructure {
                    position: position + 1,
                    symbol,
                })
            }
        }
    }
    if let Some(opening) = stack.pop() {
        return Err(RnaError::UnmatchedOpening {
            position: opening + 1,
            symbol: '(',
        });
    }
    Ok((symbols.into_iter().collect(), partner))
}

fn motif_matches_structure(motif: &CompiledMotif, global: &[Option<usize>]) -> bool {
    motif.partner.iter().enumerate().all(|(offset, partner)| {
        let position = motif.start + offset;
        match partner {
            Some(mate) => global[position] == Some(motif.start + mate),
            None => global[position].is_none(),
        }
    })
}

#[allow(clippy::too_many_arguments)]
fn enumerate_binding_microstates(
    compatible: &[usize],
    motifs: &[CompiledMotif],
    position: usize,
    bound: &mut Vec<usize>,
    energy: f64,
    visit: &mut impl FnMut(&[usize], f64) -> Result<(), RnaError>,
) -> Result<(), RnaError> {
    if position == compatible.len() {
        return visit(bound, energy);
    }
    enumerate_binding_microstates(compatible, motifs, position + 1, bound, energy, visit)?;
    let candidate = compatible[position];
    let site = &motifs[candidate];
    if bound.iter().all(|&occupied| {
        let other = &motifs[occupied];
        site.end <= other.start || other.end <= site.start
    }) {
        bound.push(candidate);
        enumerate_binding_microstates(
            compatible,
            motifs,
            position + 1,
            bound,
            energy + site.effective_energy,
            visit,
        )?;
        bound.pop();
    }
    Ok(())
}

fn iupac_matches(symbol: u8, base: u8) -> bool {
    match symbol {
        b'A' => base == b'A',
        b'C' => base == b'C',
        b'G' => base == b'G',
        b'U' => base == b'U',
        b'R' => matches!(base, b'A' | b'G'),
        b'Y' => matches!(base, b'C' | b'U'),
        b'S' => matches!(base, b'C' | b'G'),
        b'W' => matches!(base, b'A' | b'U'),
        b'K' => matches!(base, b'G' | b'U'),
        b'M' => matches!(base, b'A' | b'C'),
        b'B' => matches!(base, b'C' | b'G' | b'U'),
        b'D' => matches!(base, b'A' | b'G' | b'U'),
        b'H' => matches!(base, b'A' | b'C' | b'U'),
        b'V' => matches!(base, b'A' | b'C' | b'G'),
        b'N' => matches!(base, b'A' | b'C' | b'G' | b'U'),
        _ => false,
    }
}

fn log_add(left: f64, right: f64) -> f64 {
    if left == NEG_INF {
        return right;
    }
    if right == NEG_INF {
        return left;
    }
    let maximum = left.max(right);
    maximum + ((left - maximum).exp() + (right - maximum).exp()).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::ConstraintConfig;
    use crate::partition::partition_with_constraints;

    fn motif(id: &str) -> LigandMotif {
        LigandMotif {
            id: id.into(),
            start: 1,
            sequence: "GGGAAACCC".into(),
            structure: "(((...)))".into(),
            standard_binding_energy_kcal_mol: -3.0,
            concentration_molar: 1.0,
        }
    }

    #[test]
    fn one_site_partition_matches_an_independent_two_macrostate_identity() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let constraints =
            ConstraintModel::compile(sequence.len(), &ConstraintConfig::default()).unwrap();
        let ordinary = partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        let result =
            ligand_ensemble_exact(sequence, &[motif("aptamer")], 3, 1.0, &model, &constraints)
                .unwrap();

        let target_config = ConstraintConfig {
            force_unpaired: vec![4, 5, 6],
            force_pairs: vec![
                crate::constraints::PairConstraint { i: 1, j: 9 },
                crate::constraints::PairConstraint { i: 2, j: 8 },
                crate::constraints::PairConstraint { i: 3, j: 7 },
            ],
            ..ConstraintConfig::default()
        };
        let target_constraints = ConstraintModel::compile(sequence.len(), &target_config).unwrap();
        let target =
            partition_with_constraints(sequence, 37.0, 3, &model, &target_constraints).unwrap();
        let rt = GAS_CONSTANT_KCAL * 310.15;
        let expected_log_z = log_add(
            ordinary.log_partition_function,
            target.log_partition_function + 3.0 / rt,
        );
        assert!((result.log_partition_function - expected_log_z).abs() < 1.0e-10);
        let expected_occupancy = (target.log_partition_function + 3.0 / rt - expected_log_z).exp();
        assert!((result.motifs[0].occupancy_probability - expected_occupancy).abs() < 1.0e-10);
    }

    #[test]
    fn overlapping_sites_are_mutually_exclusive_microstates() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let result = ligand_ensemble_exact(
            sequence,
            &[motif("one"), motif("two")],
            3,
            1.0,
            &model,
            &constraints,
        )
        .unwrap();
        assert!(
            (result.motifs[0].occupancy_probability - result.motifs[1].occupancy_probability).abs()
                < 1.0e-12
        );
        assert!(
            result
                .motifs
                .iter()
                .map(|motif| motif.occupancy_probability)
                .sum::<f64>()
                < 1.0
        );
    }
}
