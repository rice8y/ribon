//! Exact inverse folding by complete sequence-space search.
//!
//! Every sequence admitted by an IUPAC template, target canonical-pair
//! requirements, and GC interval is evaluated. Ranking uses the exact target
//! structure probability obtained as a constrained-partition ratio. The
//! requested output count truncates only the returned ranking, never search.

use crate::constraints::{ConstraintConfig, ConstraintModel, PairConstraint};
use crate::energy::{EnergyModel, NucleicAcid};
use crate::fold::fold_mfe_with_constraints;
use crate::partition::partition_with_constraints;
use crate::structure::{normalize_sequence, RnaError};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct InverseDesignOptions {
    pub minimum_gc_fraction: f64,
    pub maximum_gc_fraction: f64,
    pub return_count: usize,
}

impl Default for InverseDesignOptions {
    fn default() -> Self {
        Self {
            minimum_gc_fraction: 0.0,
            maximum_gc_fraction: 1.0,
            return_count: 10,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct DesignCandidate {
    pub rank: usize,
    pub sequence: String,
    pub target_structure: String,
    pub target_probability: f64,
    pub log_target_probability: f64,
    pub target_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub target_is_mfe: bool,
    pub gc_fraction: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct InverseDesignResult {
    pub target_structure: String,
    pub template: String,
    pub length: usize,
    pub candidate_sequence_count: usize,
    pub evaluated_sequence_count: usize,
    pub returned: usize,
    pub candidates: Vec<DesignCandidate>,
    pub objective: &'static str,
    pub search_complete: bool,
    pub output_truncated: bool,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub model: &'static str,
    pub algorithm: &'static str,
    pub time_complexity: &'static str,
    pub space_complexity: &'static str,
}

#[derive(Clone)]
struct RankedSequence {
    sequence: String,
    log_probability: f64,
    gc_fraction: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn inverse_fold_exact(
    target_structure: &str,
    template: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraint_config: &ConstraintConfig,
    options: &InverseDesignOptions,
) -> Result<InverseDesignResult, RnaError> {
    validate_options(options)?;
    let ParsedPlanarTarget {
        structure: target,
        pairs,
        partner,
    } = parse_planar_target(target_structure)?;
    let template = normalize_sequence(template)?;
    if template.len() != target.len() {
        return Err(RnaError::LengthMismatch {
            sequence: template.len(),
            structure: target.len(),
        });
    }
    for &(i, j) in &pairs {
        if j - i <= min_loop {
            return Err(RnaError::InvalidOption(format!(
                "target pair {}-{} violates min_loop={min_loop}",
                i + 1,
                j + 1
            )));
        }
    }
    let allowed = template
        .bytes()
        .map(allowed_bases)
        .collect::<Result<Vec<_>, _>>()?;
    let base_constraints = ConstraintModel::compile(target.len(), constraint_config)?;
    let target_config = exact_target_constraints(constraint_config, &pairs, &partner);
    let target_constraints = ConstraintModel::compile(target.len(), &target_config)?;

    let mut ranked = Vec::new();
    let mut sequence = vec![b'A'; target.len()];
    let mut candidate_sequence_count = 0usize;
    let mut evaluated_sequence_count = 0usize;
    enumerate_sequences(
        0,
        &allowed,
        &partner,
        model,
        &mut sequence,
        &mut |candidate| {
            candidate_sequence_count =
                candidate_sequence_count.checked_add(1).ok_or_else(|| {
                    RnaError::Numerical("inverse-design sequence count overflow".into())
                })?;
            let gc_count = candidate
                .iter()
                .filter(|&&base| matches!(base, b'G' | b'C'))
                .count();
            let gc_fraction = gc_count as f64 / candidate.len() as f64;
            if gc_fraction + 1.0e-15 < options.minimum_gc_fraction
                || gc_fraction - 1.0e-15 > options.maximum_gc_fraction
            {
                return Ok(());
            }
            let candidate = std::str::from_utf8(candidate)
                .map_err(|_| RnaError::Numerical("generated a non-ASCII sequence".into()))?;
            if model
                .evaluate_with_constraints(candidate, &target, min_loop, &target_constraints)
                .is_err()
            {
                return Ok(());
            }
            let ensemble = partition_with_constraints(
                candidate,
                model.temperature_celsius(),
                min_loop,
                model,
                &base_constraints,
            )?;
            let target_ensemble = partition_with_constraints(
                candidate,
                model.temperature_celsius(),
                min_loop,
                model,
                &target_constraints,
            )?;
            let log_probability =
                target_ensemble.log_partition_function - ensemble.log_partition_function;
            if !log_probability.is_finite() || log_probability > 1.0e-9 {
                return Err(RnaError::Numerical(format!(
                    "invalid target log probability {log_probability}"
                )));
            }
            evaluated_sequence_count =
                evaluated_sequence_count.checked_add(1).ok_or_else(|| {
                    RnaError::Numerical("inverse-design evaluated count overflow".into())
                })?;
            ranked.push(RankedSequence {
                sequence: candidate.into(),
                log_probability,
                gc_fraction,
            });
            Ok(())
        },
    )?;
    if ranked.is_empty() {
        return Err(RnaError::InvalidOption(
            "no sequence satisfies the template, target, GC interval, and constraints".into(),
        ));
    }
    ranked.sort_by(|left, right| {
        right
            .log_probability
            .total_cmp(&left.log_probability)
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    let output_truncated = ranked.len() > options.return_count;
    ranked.truncate(options.return_count);

    let mut candidates = Vec::with_capacity(ranked.len());
    for (rank, entry) in ranked.into_iter().enumerate() {
        let target_energy = model
            .evaluate_with_constraints(&entry.sequence, &target, min_loop, &target_constraints)?
            .total_kcal_mol;
        let ensemble = partition_with_constraints(
            &entry.sequence,
            model.temperature_celsius(),
            min_loop,
            model,
            &base_constraints,
        )?;
        let mfe = fold_mfe_with_constraints(&entry.sequence, min_loop, model, &base_constraints)?;
        let display_sequence = display_sequence(&entry.sequence, model.nucleic_acid());
        candidates.push(DesignCandidate {
            rank: rank + 1,
            sequence: display_sequence,
            target_structure: target.clone(),
            target_probability: entry.log_probability.exp().clamp(0.0, 1.0),
            log_target_probability: entry.log_probability,
            target_energy_kcal_mol: target_energy,
            ensemble_free_energy_kcal_mol: ensemble.ensemble_free_energy_kcal_mol,
            mfe_structure: mfe.structure.clone(),
            mfe_energy_kcal_mol: mfe.energy_kcal_mol,
            target_is_mfe: (target_energy - mfe.energy_kcal_mol).abs() <= 1.0e-9,
            gc_fraction: entry.gc_fraction,
        });
    }

    Ok(InverseDesignResult {
        target_structure: target,
        template: display_sequence(&template, model.nucleic_acid()),
        length: template.len(),
        candidate_sequence_count,
        evaluated_sequence_count,
        returned: candidates.len(),
        candidates,
        objective: "maximum exact target-structure ensemble probability",
        search_complete: true,
        output_truncated,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        model: model.model_name(),
        algorithm: "complete IUPAC sequence enumeration with canonical-pair pruning and exact constrained-partition scoring",
        time_complexity: "exponential in unconstrained template positions",
        space_complexity: "returned ranking plus one thermodynamic DP",
    })
}

fn validate_options(options: &InverseDesignOptions) -> Result<(), RnaError> {
    if !options.minimum_gc_fraction.is_finite()
        || !options.maximum_gc_fraction.is_finite()
        || options.minimum_gc_fraction < 0.0
        || options.maximum_gc_fraction > 1.0
        || options.minimum_gc_fraction > options.maximum_gc_fraction
    {
        return Err(RnaError::InvalidOption(
            "inverse-design GC fractions must satisfy 0 <= minimum <= maximum <= 1".into(),
        ));
    }
    if options.return_count == 0 {
        return Err(RnaError::InvalidOption(
            "inverse-design return_count must be positive".into(),
        ));
    }
    Ok(())
}

fn exact_target_constraints(
    base: &ConstraintConfig,
    pairs: &[(usize, usize)],
    partner: &[Option<usize>],
) -> ConstraintConfig {
    let mut target = base.clone();
    for &(i, j) in pairs {
        target
            .force_pairs
            .push(PairConstraint { i: i + 1, j: j + 1 });
    }
    for (position, paired) in partner.iter().enumerate() {
        if paired.is_none() {
            target.force_unpaired.push(position + 1);
        }
    }
    target
}

struct ParsedPlanarTarget {
    structure: String,
    pairs: Vec<(usize, usize)>,
    partner: Vec<Option<usize>>,
}

fn parse_planar_target(structure: &str) -> Result<ParsedPlanarTarget, RnaError> {
    let symbols = structure
        .chars()
        .filter(|symbol| !symbol.is_whitespace())
        .collect::<Vec<_>>();
    if symbols.is_empty() {
        return Err(RnaError::InvalidOption(
            "inverse-design target structure is empty".into(),
        ));
    }
    let mut stack = Vec::new();
    let mut pairs = Vec::new();
    let mut partner = vec![None; symbols.len()];
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
                pairs.push((opening, position));
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
    pairs.sort_unstable();
    Ok(ParsedPlanarTarget {
        structure: symbols.into_iter().collect(),
        pairs,
        partner,
    })
}

fn allowed_bases(symbol: u8) -> Result<Vec<u8>, RnaError> {
    let bases = match symbol {
        b'A' => b"A".as_slice(),
        b'C' => b"C".as_slice(),
        b'G' => b"G".as_slice(),
        b'U' => b"U".as_slice(),
        b'T' => b"T".as_slice(),
        b'R' => b"AG".as_slice(),
        b'Y' => b"CU".as_slice(),
        b'S' => b"CG".as_slice(),
        b'W' => b"AU".as_slice(),
        b'K' => b"GU".as_slice(),
        b'M' => b"AC".as_slice(),
        b'B' => b"CGU".as_slice(),
        b'D' => b"AGU".as_slice(),
        b'H' => b"ACU".as_slice(),
        b'V' => b"ACG".as_slice(),
        b'N' => b"ACGU".as_slice(),
        _ => {
            return Err(RnaError::InvalidSequence {
                position: 0,
                symbol: symbol as char,
            })
        }
    };
    Ok(bases.to_vec())
}

fn enumerate_sequences(
    position: usize,
    allowed: &[Vec<u8>],
    partner: &[Option<usize>],
    model: &EnergyModel,
    sequence: &mut [u8],
    visit: &mut impl FnMut(&[u8]) -> Result<(), RnaError>,
) -> Result<(), RnaError> {
    if position == sequence.len() {
        return visit(sequence);
    }
    if let Some(mate) = partner[position] {
        if mate < position {
            return enumerate_sequences(position + 1, allowed, partner, model, sequence, visit);
        }
        for &left in &allowed[position] {
            for &right in &allowed[mate] {
                if model.can_pair(left, right) {
                    sequence[position] = left;
                    sequence[mate] = right;
                    enumerate_sequences(position + 1, allowed, partner, model, sequence, visit)?;
                }
            }
        }
    } else {
        for &base in &allowed[position] {
            sequence[position] = base;
            enumerate_sequences(position + 1, allowed, partner, model, sequence, visit)?;
        }
    }
    Ok(())
}

fn display_sequence(sequence: &str, family: NucleicAcid) -> String {
    if family == NucleicAcid::Dna {
        sequence.replace('U', "T")
    } else {
        sequence.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_design_enumerates_every_target_compatible_sequence() {
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let result = inverse_fold_exact(
            "(...)",
            "NNNNN",
            3,
            &model,
            &ConstraintConfig::default(),
            &InverseDesignOptions {
                return_count: 4,
                ..InverseDesignOptions::default()
            },
        )
        .unwrap();
        // Six oriented canonical pairs and 4^3 unpaired assignments.
        assert_eq!(result.candidate_sequence_count, 6 * 4usize.pow(3));
        assert_eq!(
            result.evaluated_sequence_count,
            result.candidate_sequence_count
        );
        assert!(result.search_complete);
        assert!(result.output_truncated);
        assert_eq!(result.returned, 4);
        for candidate in &result.candidates {
            assert!(candidate.target_probability > 0.0);
            assert!(candidate.target_probability <= 1.0);
        }
        assert!(result.candidates.windows(2).all(|entries| {
            entries[0].log_target_probability + 1.0e-14 >= entries[1].log_target_probability
        }));
    }

    #[test]
    fn gc_filter_changes_only_the_explicit_admissible_space() {
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let result = inverse_fold_exact(
            "(...)",
            "NNNNN",
            3,
            &model,
            &ConstraintConfig::default(),
            &InverseDesignOptions {
                minimum_gc_fraction: 1.0,
                maximum_gc_fraction: 1.0,
                return_count: 2,
            },
        )
        .unwrap();
        // Two GC/CG pair orientations and 2^3 G/C unpaired assignments.
        assert_eq!(result.evaluated_sequence_count, 2 * 2usize.pow(3));
        assert!(result
            .candidates
            .iter()
            .all(|candidate| candidate.gc_fraction == 1.0));
    }

    #[test]
    fn dna_template_accepts_and_preserves_thymine() {
        let model = EnergyModel::with_parameter_family(37.0, 0, 1.021, NucleicAcid::Dna).unwrap();
        let result = inverse_fold_exact(
            "(...)",
            "TNNNA",
            3,
            &model,
            &ConstraintConfig::default(),
            &InverseDesignOptions {
                return_count: 3,
                ..InverseDesignOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.candidate_sequence_count, 4usize.pow(3));
        assert!(result.candidates.iter().all(|candidate| {
            candidate.sequence.starts_with('T')
                && candidate.sequence.ends_with('A')
                && !candidate.sequence.contains('U')
        }));
    }
}
