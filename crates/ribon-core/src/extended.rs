//! Extended alphabets and G-quadruplex analysis.

#[cfg(test)]
use crate::analyze_with_options;
use crate::constraints::{ConstraintConfig, ConstraintModel, PositionEnergy};
use crate::energy::EnergyModel;
use crate::gquad_integrated;
use crate::modified_parameters::{self, ModifiedBaseKind};
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, RnaError};
use crate::{analyze_with_compiled_constraints, AnalysisResult};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const T_MEASURE_KELVIN: f64 = 310.15;

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct ModifiedBase {
    pub position: usize,
    pub symbol: String,
    pub canonical_base: char,
    /// Selects a calibrated sparse nearest-neighbor parameter set.  `None`
    /// retains the generic explicit pseudo-energy API.
    #[serde(default)]
    pub kind: Option<ModifiedBaseKind>,
    #[serde(default)]
    pub paired_energy_kcal_mol: f64,
    #[serde(default)]
    pub unpaired_energy_kcal_mol: f64,
    #[serde(default)]
    pub stack_energy_kcal_mol: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModifiedFoldResult {
    pub canonical_sequence: String,
    /// Canonical alphabet actually used by the folding grammar.  This differs
    /// only for inosine, whose supported I-C/I-U pairs use G as the explicitly
    /// defined canonical folding surrogate.
    pub folding_sequence: String,
    pub display_symbols: Vec<String>,
    pub modifications: Vec<ModifiedBase>,
    pub parameter_usage: Vec<ModifiedParameterUse>,
    pub analysis: AnalysisResult,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModifiedParameterUse {
    pub kind: ModifiedBaseKind,
    pub name: &'static str,
    pub source_url: &'static str,
    pub calibration: &'static str,
    pub matched_stacks: usize,
    pub canonical_reference_stacks: usize,
    pub matched_helix_ends: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct ParameterCounts {
    matched_stacks: usize,
    canonical_reference_stacks: usize,
    matched_helix_ends: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct GQuadruplex {
    pub from: usize,
    pub to: usize,
    pub stack_size: usize,
    pub linker_lengths: [usize; 3],
    pub energy_kcal_mol: f64,
    pub probability: f64,
    pub guanine_positions: Vec<usize>,
    pub structure: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct GQuadResult {
    pub sequence: String,
    pub structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub secondary_structure_mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub log_partition_function: f64,
    pub gquad_probability: f64,
    pub expected_gquad_count: f64,
    pub gquad_elements: Vec<GQuadProbability>,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub gquad_position_probabilities: Vec<f64>,
    pub candidates: Vec<GQuadruplex>,
    pub temperature_celsius: f64,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct GQuadProbability {
    pub from: usize,
    pub to: usize,
    pub probability: f64,
}

#[allow(clippy::too_many_arguments)]
pub fn fold_modified(
    sequence: &str,
    modifications: &[ModifiedBase],
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
) -> Result<ModifiedFoldResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    fold_modified_with_model(sequence, modifications, min_loop, gamma, model)
}

pub fn fold_modified_with_model(
    sequence: &str,
    modifications: &[ModifiedBase],
    min_loop: usize,
    gamma: f64,
    energy_model: EnergyModel,
) -> Result<ModifiedFoldResult, RnaError> {
    if energy_model.nucleic_acid() != crate::energy::NucleicAcid::Rna {
        return Err(RnaError::InvalidOption(
            "modified-base calibrations are defined for RNA parameter families".into(),
        ));
    }
    let temperature_celsius = energy_model.temperature_celsius();
    let custom_profile = energy_model.parameter_profile_name().is_some();
    let canonical = normalize_sequence(sequence)?;
    let n = canonical.len();
    let canonical_bases = canonical.as_bytes();
    let mut folding_bases = canonical_bases.to_vec();
    let mut display_symbols = canonical.chars().map(String::from).collect::<Vec<_>>();
    let mut config = ConstraintConfig::default();
    let mut kinds = vec![None; n];
    let mut occupied = HashSet::new();
    for modification in modifications {
        if modification.position == 0 || modification.position > n {
            return Err(RnaError::InvalidOption(format!(
                "modified base position {} is outside 1-{n}",
                modification.position
            )));
        }
        if modification.symbol.trim().is_empty() {
            return Err(RnaError::InvalidOption(
                "modified base symbol must not be empty".into(),
            ));
        }
        if !occupied.insert(modification.position) {
            return Err(RnaError::InvalidOption(format!(
                "more than one modified base was supplied at position {}",
                modification.position
            )));
        }
        let base = modification.canonical_base.to_ascii_uppercase();
        let base = if base == 'T' { 'U' } else { base };
        if !matches!(base, 'A' | 'C' | 'G' | 'U') {
            return Err(RnaError::InvalidOption(format!(
                "modified base {} must map to A, C, G, or U",
                modification.symbol
            )));
        }
        if canonical.as_bytes()[modification.position - 1] != base as u8 {
            return Err(RnaError::InvalidOption(format!(
                "modified base {} maps to {base}, but canonical sequence position {} is {}",
                modification.symbol,
                modification.position,
                canonical.as_bytes()[modification.position - 1] as char
            )));
        }
        if let Some(kind) = modification.kind {
            if kind.precursor() != base as u8 {
                return Err(RnaError::InvalidOption(format!(
                    "{} requires precursor {}, not {base}",
                    kind.label(),
                    kind.precursor() as char
                )));
            }
            kinds[modification.position - 1] = Some(kind);
            folding_bases[modification.position - 1] = kind.folding_base();
        }
        for (name, value) in [
            ("paired", modification.paired_energy_kcal_mol),
            ("unpaired", modification.unpaired_energy_kcal_mol),
            ("stack", modification.stack_energy_kcal_mol),
        ] {
            if !value.is_finite() {
                return Err(RnaError::InvalidOption(format!(
                    "modified base {name} energy must be finite"
                )));
            }
        }
        display_symbols[modification.position - 1] = modification.symbol.clone();
        config.soft.paired.push(PositionEnergy {
            position: modification.position,
            energy_kcal_mol: modification.paired_energy_kcal_mol,
        });
        config.soft.unpaired.push(PositionEnergy {
            position: modification.position,
            energy_kcal_mol: modification.unpaired_energy_kcal_mol,
        });
        config.soft.stack.push(PositionEnergy {
            position: modification.position,
            energy_kcal_mol: modification.stack_energy_kcal_mol,
        });
    }
    let folding_sequence = String::from_utf8(folding_bases.clone())
        .expect("normalized canonical bases and folding substitutions are ASCII");
    let reference_model = energy_model.with_salt_molar(1.021)?;
    let mut constraints = ConstraintModel::compile(n, &config)?;
    let mut counts = HashMap::<ModifiedBaseKind, ParameterCounts>::new();

    // Apply sparse measured stack parameters as delta energies relative to the
    // model's explicitly defined canonical reference state. This preserves
    // the salt term and makes the same correction visible to MFE and PF.
    for i in 0..n {
        let k = i + 1;
        if k >= n {
            break;
        }
        for j in (k + 1)..n {
            let l = j - 1;
            if k >= l
                || !energy_model.can_pair(folding_bases[i], folding_bases[j])
                || !energy_model.can_pair(folding_bases[k], folding_bases[l])
            {
                continue;
            }
            let positions = [i, k, j, l];
            let canonical_key = [
                folding_bases[i],
                folding_bases[k],
                folding_bases[j],
                folding_bases[l],
            ];
            let canonical_reference = reference_model.internal_energy(&folding_bases, i, j, k, l);
            let present = positions
                .iter()
                .enumerate()
                .filter_map(|(offset, &position)| kinds[position].map(|kind| (kind, offset)))
                .collect::<Vec<_>>();
            if present.is_empty() || !canonical_reference.is_finite() {
                continue;
            }

            let unique = present
                .iter()
                .map(|&(kind, _)| kind)
                .collect::<HashSet<_>>();
            let mut correction = 0.0;
            for kind in unique {
                let offsets = present
                    .iter()
                    .filter_map(|&(candidate, offset)| (candidate == kind).then_some(offset))
                    .collect::<Vec<_>>();
                let stats = counts.entry(kind).or_default();
                if kind == ModifiedBaseKind::Dihydrouridine {
                    correction += 1.5;
                    stats.matched_stacks += 1;
                    continue;
                }
                if let Some(measured) = modified_parameters::stack_energy(
                    kind,
                    canonical_key,
                    &offsets,
                    temperature_celsius,
                ) {
                    correction += measured - canonical_reference;
                    stats.matched_stacks += 1;
                } else if offsets.len() > 1 {
                    let mut matched = 0;
                    for offset in offsets {
                        if let Some(measured) = modified_parameters::stack_energy(
                            kind,
                            canonical_key,
                            &[offset],
                            temperature_celsius,
                        ) {
                            correction += measured - canonical_reference;
                            matched += 1;
                        }
                    }
                    if matched == 0 {
                        stats.canonical_reference_stacks += 1;
                    } else {
                        stats.matched_stacks += 1;
                    }
                } else {
                    stats.canonical_reference_stacks += 1;
                }
            }
            constraints.add_context_stack_energy(i, j, k, l, correction);
        }
    }

    // Published terminal terms replace the canonical terminal AU/GU-style
    // contribution whenever the modified pair occurs at a helix end.
    for i in 0..n {
        for j in (i + 1)..n {
            if !energy_model.can_pair(folding_bases[i], folding_bases[j]) {
                continue;
            }
            let mut correction = 0.0;
            for (position, other) in [(i, folding_bases[j]), (j, folding_bases[i])] {
                let Some(kind) = kinds[position] else {
                    continue;
                };
                if let Some(measured) =
                    modified_parameters::terminal_energy(kind, other, temperature_celsius)
                {
                    correction += measured
                        - reference_model.terminal_pair_energy(folding_bases[i], folding_bases[j]);
                    counts.entry(kind).or_default().matched_helix_ends += 1;
                }
            }
            constraints.add_terminal_pair_energy(i, j, correction);
        }
    }

    let analysis = analyze_with_compiled_constraints(
        folding_sequence.clone(),
        temperature_celsius,
        min_loop,
        gamma,
        energy_model,
        constraints,
    )?;
    let mut parameter_usage = counts
        .into_iter()
        .map(|(kind, counts)| ModifiedParameterUse {
            kind,
            name: kind.label(),
            source_url: kind.source_url(),
            calibration: kind.calibration(),
            matched_stacks: counts.matched_stacks,
            canonical_reference_stacks: counts.canonical_reference_stacks,
            matched_helix_ends: counts.matched_helix_ends,
        })
        .collect::<Vec<_>>();
    parameter_usage.sort_by_key(|entry| entry.name);
    Ok(ModifiedFoldResult {
        canonical_sequence: canonical,
        folding_sequence,
        display_symbols,
        modifications: modifications.to_vec(),
        parameter_usage,
        analysis,
        model: if custom_profile {
            "custom normalized RNA model with experimental sparse modified-base nearest-neighbor corrections"
        } else {
            "Turner 2004 with experimental sparse modified-base nearest-neighbor corrections"
        },
    })
}

#[allow(clippy::too_many_arguments)]
pub fn fold_gquad(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
) -> Result<GQuadResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    fold_gquad_with_model(sequence, min_loop, gamma, model)
}

pub fn fold_gquad_with_model(
    sequence: &str,
    min_loop: usize,
    gamma: f64,
    model: EnergyModel,
) -> Result<GQuadResult, RnaError> {
    if model.nucleic_acid() != crate::energy::NucleicAcid::Rna {
        return Err(RnaError::InvalidOption(
            "G-quadruplex calibration is defined for RNA parameter families".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    let custom_profile = model.parameter_profile_name().is_some();
    let temperature_celsius = model.temperature_celsius();
    let secondary = crate::analyze_with_model(
        sequence.clone(),
        min_loop,
        gamma,
        model.clone(),
        &ConstraintConfig::default(),
    )?;
    let mut candidates = enumerate_gquads(&sequence, temperature_celsius);
    let integrated =
        gquad_integrated::predict_with_model(&sequence, min_loop, &model, &candidates)?;
    for (candidate, &probability) in candidates
        .iter_mut()
        .zip(&integrated.candidate_probabilities)
    {
        candidate.probability = probability;
    }
    candidates.sort_by(|a, b| {
        a.energy_kcal_mol
            .total_cmp(&b.energy_kcal_mol)
            .then_with(|| a.from.cmp(&b.from))
            .then_with(|| a.to.cmp(&b.to))
    });
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let log_total = integrated.log_partition_function;
    let gquad_probability = if secondary.log_partition_function >= log_total {
        0.0
    } else {
        (-(secondary.log_partition_function - log_total).exp_m1()).clamp(0.0, 1.0)
    };
    let expected_gquad_count = integrated
        .interval_probabilities
        .iter()
        .map(|entry| entry.2)
        .sum();
    Ok(GQuadResult {
        sequence,
        structure: integrated.structure,
        mfe_energy_kcal_mol: integrated.mfe_energy,
        secondary_structure_mfe_energy_kcal_mol: secondary.mfe_energy_kcal_mol,
        ensemble_free_energy_kcal_mol: -rt * log_total,
        log_partition_function: log_total,
        gquad_probability,
        expected_gquad_count,
        gquad_elements: integrated
            .interval_probabilities
            .into_iter()
            .map(|(from, to, probability)| GQuadProbability {
                from,
                to,
                probability,
            })
            .collect(),
        pair_probabilities: integrated.pair_probabilities,
        unpaired_probabilities: integrated.unpaired_probabilities,
        gquad_position_probabilities: integrated.gquad_position_probabilities,
        candidates,
        temperature_celsius,
        model: if custom_profile {
            "integrated custom RNA/G-quadruplex exterior and generalized-multiloop grammar"
        } else {
            "integrated Turner/G-quadruplex exterior and generalized-multiloop grammar"
        },
    })
}

fn enumerate_gquads(sequence: &str, temperature_celsius: f64) -> Vec<GQuadruplex> {
    let bases = sequence.as_bytes();
    let n = bases.len();
    let mut result = Vec::new();
    for start in 0..n {
        for stack_size in 2..=7 {
            for linker_1 in 1..=15 {
                for linker_2 in 1..=15 {
                    for linker_3 in 1..=15 {
                        let second = start + stack_size + linker_1;
                        let third = second + stack_size + linker_2;
                        let fourth = third + stack_size + linker_3;
                        let end = fourth + stack_size;
                        if end > n {
                            continue;
                        }
                        if [start, second, third, fourth].iter().all(|&position| {
                            bases[position..position + stack_size]
                                .iter()
                                .all(|&b| b == b'G')
                        }) {
                            let positions = [start, second, third, fourth]
                                .into_iter()
                                .flat_map(|position| position..position + stack_size)
                                .map(|position| position + 1)
                                .collect::<Vec<_>>();
                            let mut structure = vec!['.'; n];
                            for &position in &positions {
                                structure[position - 1] = '+';
                            }
                            result.push(GQuadruplex {
                                from: start + 1,
                                to: end,
                                stack_size,
                                linker_lengths: [linker_1, linker_2, linker_3],
                                energy_kcal_mol: gquad_energy(
                                    stack_size,
                                    linker_1 + linker_2 + linker_3,
                                    temperature_celsius,
                                ),
                                probability: 0.0,
                                guanine_positions: positions,
                                structure: structure.into_iter().collect(),
                            });
                        }
                    }
                }
            }
        }
    }
    result
}

fn gquad_energy(stack_size: usize, total_linker: usize, temperature_celsius: f64) -> f64 {
    let ratio = (temperature_celsius + 273.15) / T_MEASURE_KELVIN;
    // ViennaRNA retains the rescaled alpha/beta terms as doubles and casts
    // only their final linear/logarithmic combination to centi-kcal.
    let alpha = -11_934.0 - (-11_934.0 + 1_800.0) * ratio;
    let beta = 0.0 - (0.0 - 1_200.0) * ratio;
    (alpha * (stack_size - 1) as f64 + beta * ((total_linker - 2) as f64).ln()).trunc() / 100.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gquad_enumeration_matches_vienna_parameter_formula_at_37c() {
        let result = fold_gquad("GGAGGAGGAGG", 37.0, 3, 1.0, 0, 1.021).unwrap();
        let candidate = result
            .candidates
            .iter()
            .find(|candidate| candidate.stack_size == 2 && candidate.linker_lengths == [1, 1, 1])
            .unwrap();
        let expected = (-1800.0 + 1200.0 * (1.0f64).ln()) / 100.0;
        assert!((candidate.energy_kcal_mol - expected).abs() < 1.0e-12);
        assert_eq!(candidate.structure, "++.++.++.++");
    }

    #[test]
    fn integrated_gquad_grammar_reduces_to_ordinary_analysis_without_candidates() {
        let sequence = "GCGAAAGC";
        let result = fold_gquad(sequence, 37.0, 3, 1.0, 2, 1.021).unwrap();
        let ordinary = analyze_with_options(
            sequence,
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap();
        assert!(result.candidates.is_empty());
        assert_eq!(result.structure, ordinary.mfe_structure);
        assert!((result.mfe_energy_kcal_mol - ordinary.mfe_energy_kcal_mol).abs() < 1.0e-12);
        assert!((result.log_partition_function - ordinary.log_partition_function).abs() < 1.0e-12);
        assert_eq!(
            result.pair_probabilities.len(),
            ordinary.pair_probabilities.len()
        );
        for (integrated, canonical) in result
            .pair_probabilities
            .iter()
            .zip(&ordinary.pair_probabilities)
        {
            assert_eq!((integrated.i, integrated.j), (canonical.i, canonical.j));
            assert!((integrated.probability - canonical.probability).abs() < 1.0e-12);
        }
    }

    #[test]
    fn multiple_gquadruplexes_coexist_and_probabilities_are_normalized() {
        let sequence = "GGAGGAGGAGGCCCCCGGAGGAGGAGG";
        let result = fold_gquad(sequence, 37.0, 3, 1.0, 0, 1.021).unwrap();
        assert_eq!(
            result
                .structure
                .chars()
                .filter(|&symbol| symbol == '+')
                .count(),
            16
        );
        assert!(result.mfe_energy_kcal_mol <= -36.0 + 1.0e-12);
        assert!(result.expected_gquad_count > 1.5);
        let candidate_sum: f64 = result
            .candidates
            .iter()
            .map(|candidate| candidate.probability)
            .sum();
        assert!((candidate_sum - result.expected_gquad_count).abs() < 1.0e-10);
        for index in 0..sequence.len() {
            let paired_mass: f64 = result
                .pair_probabilities
                .iter()
                .filter_map(|pair| {
                    (pair.i == index + 1 || pair.j == index + 1).then_some(pair.probability)
                })
                .sum();
            let total = paired_mass
                + result.gquad_position_probabilities[index]
                + result.unpaired_probabilities[index];
            assert!(
                (total - 1.0).abs() < 1.0e-8,
                "position {}: {total}",
                index + 1
            );
        }
    }

    #[test]
    fn integrated_ensemble_contains_more_than_the_old_exclusive_state_space() {
        let sequence = "GCGGAGGAGGAGGCGC";
        let result = fold_gquad(sequence, 37.0, 3, 1.0, 0, 1.021).unwrap();
        let ordinary = analyze_with_options(
            sequence,
            37.0,
            3,
            1.0,
            0,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap();
        let rt = GAS_CONSTANT_KCAL * (37.0 + 273.15);
        let exclusive = result
            .candidates
            .iter()
            .fold(ordinary.log_partition_function, |total, candidate| {
                log_sum(total, -candidate.energy_kcal_mol / rt)
            });
        assert!(result.log_partition_function > exclusive + 1.0e-8);
    }

    fn log_sum(a: f64, b: f64) -> f64 {
        let high = a.max(b);
        high + (a.min(b) - high).exp().ln_1p()
    }

    #[test]
    fn zero_delta_modifications_equal_the_canonical_analysis() {
        let modified = fold_modified(
            "GGGAAACCC",
            &[ModifiedBase {
                position: 4,
                symbol: "m6A".into(),
                canonical_base: 'A',
                kind: None,
                paired_energy_kcal_mol: 0.0,
                unpaired_energy_kcal_mol: 0.0,
                stack_energy_kcal_mol: 0.0,
            }],
            37.0,
            3,
            1.0,
            2,
            1.021,
        )
        .unwrap();
        let canonical = analyze_with_options(
            "GGGAAACCC",
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap();
        assert_eq!(modified.analysis.mfe_structure, canonical.mfe_structure);
        assert!(
            (modified.analysis.mfe_energy_kcal_mol - canonical.mfe_energy_kcal_mol).abs() < 1.0e-12
        );
        assert_eq!(modified.display_symbols[3], "m6A");
    }

    #[test]
    fn modified_profile_is_exactly_the_corresponding_soft_constraint() {
        let modification = ModifiedBase {
            position: 3,
            symbol: "mG".into(),
            canonical_base: 'G',
            kind: None,
            paired_energy_kcal_mol: 0.35,
            unpaired_energy_kcal_mol: -0.20,
            stack_energy_kcal_mol: -0.45,
        };
        let modified = fold_modified(
            "GGGAAACCC",
            std::slice::from_ref(&modification),
            25.0,
            3,
            1.0,
            0,
            1.021,
        )
        .unwrap();
        let mut config = ConstraintConfig::default();
        config.soft.paired.push(PositionEnergy {
            position: 3,
            energy_kcal_mol: 0.35,
        });
        config.soft.unpaired.push(PositionEnergy {
            position: 3,
            energy_kcal_mol: -0.20,
        });
        config.soft.stack.push(PositionEnergy {
            position: 3,
            energy_kcal_mol: -0.45,
        });
        let direct = analyze_with_options("GGGAAACCC", 25.0, 3, 1.0, 0, 1.021, &config).unwrap();
        assert_eq!(modified.analysis.mfe_structure, direct.mfe_structure);
        assert_eq!(
            modified.analysis.mfe_energy_kcal_mol,
            direct.mfe_energy_kcal_mol
        );
        assert_eq!(
            modified.analysis.log_partition_function,
            direct.log_partition_function
        );
        assert_eq!(
            modified.analysis.pair_probabilities.len(),
            direct.pair_probabilities.len()
        );
        for (left, right) in modified
            .analysis
            .pair_probabilities
            .iter()
            .zip(&direct.pair_probabilities)
        {
            assert_eq!((left.i, left.j), (right.i, right.j));
            assert_eq!(left.probability, right.probability);
        }
    }

    #[test]
    fn calibrated_m6a_stack_changes_mfe_and_partition_with_provenance() {
        let canonical = analyze_with_options(
            "ACGAAACGU",
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap();
        let modified = fold_modified(
            "ACGAAACGU",
            &[ModifiedBase {
                position: 1,
                symbol: "m6A".into(),
                canonical_base: 'A',
                kind: Some(ModifiedBaseKind::M6a),
                paired_energy_kcal_mol: 0.0,
                unpaired_energy_kcal_mol: 0.0,
                stack_energy_kcal_mol: 0.0,
            }],
            37.0,
            3,
            1.0,
            2,
            1.021,
        )
        .unwrap();
        assert_eq!(modified.folding_sequence, "ACGAAACGU");
        assert_eq!(modified.parameter_usage.len(), 1);
        assert!(modified.parameter_usage[0].matched_stacks > 0);
        assert!(modified.parameter_usage[0].source_url.contains("doi.org"));
        assert_ne!(
            modified.analysis.log_partition_function,
            canonical.log_partition_function
        );
    }

    #[test]
    fn inosine_uses_g_surrogate_and_experimental_i_c_stack() {
        let modified = fold_modified(
            "AGGAAACCC",
            &[ModifiedBase {
                position: 1,
                symbol: "I".into(),
                canonical_base: 'A',
                kind: Some(ModifiedBaseKind::Inosine),
                paired_energy_kcal_mol: 0.0,
                unpaired_energy_kcal_mol: 0.0,
                stack_energy_kcal_mol: 0.0,
            }],
            37.0,
            3,
            1.0,
            2,
            1.021,
        )
        .unwrap();
        assert_eq!(modified.canonical_sequence, "AGGAAACCC");
        assert_eq!(modified.folding_sequence, "GGGAAACCC");
        assert!(modified.parameter_usage[0].matched_stacks > 0);
        assert!(modified
            .analysis
            .pair_probabilities
            .iter()
            .any(|pair| pair.i == 1 && pair.j == 9));
    }
}
