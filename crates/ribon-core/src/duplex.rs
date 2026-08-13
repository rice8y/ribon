//! Connected intermolecular RNA duplex prediction.
//!
//! This is the RNAduplex model: all pairs are between the two antiparallel
//! strands and form one connected interaction through stacks, bulges, and
//! internal loops. Intramolecular base pairs are intentionally excluded.

use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, RnaError};
use serde::Serialize;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const INF: f64 = 1.0e100;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Debug, Serialize)]
pub struct DuplexResult {
    pub sequence_a: String,
    pub sequence_b: String,
    pub structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub bound_partition_function: f64,
    pub log_bound_partition_function: f64,
    pub association_ensemble_free_energy_kcal_mol: f64,
    pub standard_state_bound_probability: f64,
    pub standard_state_ensemble_free_energy_kcal_mol: f64,
    pub conditional_pair_probabilities: Vec<PairProbability>,
    pub standard_state_pair_probabilities: Vec<PairProbability>,
    pub temperature_celsius: f64,
    pub salt_molar: f64,
    pub model: &'static str,
}

/// Predict a connected intermolecular duplex under Turner internal-loop rules.
pub fn duplex(
    sequence_a: &str,
    sequence_b: &str,
    temperature_celsius: f64,
    salt_molar: f64,
) -> Result<DuplexResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, 0, salt_molar)?;
    duplex_with_model(sequence_a, sequence_b, &model)
}

/// Predict a connected duplex with an explicitly selected RNA or DNA table
/// family. This connected grammar has no free exterior dangle choice.
pub fn duplex_with_model(
    sequence_a: &str,
    sequence_b: &str,
    model: &EnergyModel,
) -> Result<DuplexResult, RnaError> {
    let a = normalize_sequence(sequence_a)?;
    let b = normalize_sequence(sequence_b)?;
    let temperature_celsius = model.temperature_celsius();
    let salt_molar = model.salt_molar();
    let na = a.len();
    let nb = b.len();
    let mut combined = Vec::with_capacity(na + nb);
    combined.extend_from_slice(a.as_bytes());
    combined.extend_from_slice(b.as_bytes());

    let mut mfe = vec![vec![INF; nb]; na];
    let mut child = vec![vec![None; nb]; na];
    let mut log_inside = vec![vec![NEG_INF; nb]; na];
    for i in (0..na).rev() {
        for j in 0..nb {
            if !model.can_pair(a.as_bytes()[i], b.as_bytes()[j]) {
                continue;
            }
            // Inner interaction end, expressed in the reversed strand order
            // used by ViennaRNA's second exterior-stem contribution.
            let terminal = model.duplex_end_energy(
                b.as_bytes()[j],
                a.as_bytes()[i],
                j.checked_sub(1).map(|index| b.as_bytes()[index]),
                a.as_bytes().get(i + 1).copied(),
            );
            let terminal_pf = model.duplex_end_boltzmann_energy(
                b.as_bytes()[j],
                a.as_bytes()[i],
                j.checked_sub(1).map(|index| b.as_bytes()[index]),
                a.as_bytes().get(i + 1).copied(),
            );
            mfe[i][j] = terminal;
            let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
            log_inside[i][j] = -terminal_pf / rt;
            let geometric_limit = (na - i - 1).saturating_add(j).saturating_sub(2);
            let internal_limit = model.internal_loop_limit(geometric_limit);
            for delta_a in 1..=(internal_limit + 1) {
                let k = i + delta_a;
                if k >= na {
                    break;
                }
                let remaining = internal_limit + 2 - delta_a;
                for delta_b in 1..=remaining {
                    let Some(l) = j.checked_sub(delta_b) else {
                        continue;
                    };
                    if mfe[k][l] >= INF / 2.0 {
                        continue;
                    }
                    let transition = model.internal_energy(&combined, i, na + j, k, na + l);
                    let candidate = transition + mfe[k][l];
                    if candidate < mfe[i][j] {
                        mfe[i][j] = candidate;
                        child[i][j] = Some((k, l));
                    }
                    let transition_pf =
                        model.internal_boltzmann_energy(&combined, i, na + j, k, na + l);
                    log_inside[i][j] =
                        log_add(log_inside[i][j], log_inside[k][l] - transition_pf / rt);
                }
            }
        }
    }

    let initiation = model.duplex_initialization_energy();
    let initiation_pf = model.duplex_initialization_boltzmann_energy();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut best = INF;
    let mut root = None;
    let mut log_bound = NEG_INF;
    let mut outside = vec![vec![NEG_INF; nb]; na];
    for i in 0..na {
        for j in 0..nb {
            if mfe[i][j] >= INF / 2.0 {
                continue;
            }
            let outer = model.duplex_end_energy(
                a.as_bytes()[i],
                b.as_bytes()[j],
                i.checked_sub(1).map(|index| a.as_bytes()[index]),
                b.as_bytes().get(j + 1).copied(),
            );
            let candidate = initiation + outer + mfe[i][j];
            if candidate < best {
                best = candidate;
                root = Some((i, j));
            }
            let outer_pf = model.duplex_end_boltzmann_energy(
                a.as_bytes()[i],
                b.as_bytes()[j],
                i.checked_sub(1).map(|index| a.as_bytes()[index]),
                b.as_bytes().get(j + 1).copied(),
            );
            let root_transition = -(initiation_pf + outer_pf) / rt;
            log_bound = log_add(log_bound, root_transition + log_inside[i][j]);
            outside[i][j] = log_add(outside[i][j], root_transition);
        }
    }
    let Some((mut i, mut j)) = root else {
        return Err(RnaError::InvalidOption(
            "the two strands have no canonical intermolecular base pair".into(),
        ));
    };

    let mut mfe_pairs = Vec::new();
    loop {
        mfe_pairs.push((i, j));
        let Some((k, l)) = child[i][j] else {
            break;
        };
        i = k;
        j = l;
    }
    let mut left = vec!['.'; na];
    let mut right = vec!['.'; nb];
    for &(i, j) in &mfe_pairs {
        left[i] = '(';
        right[j] = ')';
    }
    let structure = format!(
        "{}&{}",
        left.into_iter().collect::<String>(),
        right.into_iter().collect::<String>()
    );

    // Outside propagation follows outer pairs to their enclosed pair.
    for i in 0..na {
        for j in (0..nb).rev() {
            if outside[i][j] == NEG_INF || log_inside[i][j] == NEG_INF {
                continue;
            }
            let geometric_limit = (na - i - 1).saturating_add(j).saturating_sub(2);
            let internal_limit = model.internal_loop_limit(geometric_limit);
            for delta_a in 1..=(internal_limit + 1) {
                let k = i + delta_a;
                if k >= na {
                    break;
                }
                let remaining = internal_limit + 2 - delta_a;
                for delta_b in 1..=remaining {
                    let Some(l) = j.checked_sub(delta_b) else {
                        continue;
                    };
                    if log_inside[k][l] == NEG_INF {
                        continue;
                    }
                    let transition =
                        model.internal_boltzmann_energy(&combined, i, na + j, k, na + l);
                    outside[k][l] = log_add(outside[k][l], outside[i][j] - transition / rt);
                }
            }
        }
    }

    let log_total = log_add(0.0, log_bound);
    let bound_probability = (log_bound - log_total).exp();
    let mut conditional = Vec::new();
    let mut standard = Vec::new();
    for i in 0..na {
        for j in 0..nb {
            if outside[i][j] == NEG_INF || log_inside[i][j] == NEG_INF {
                continue;
            }
            let mass = outside[i][j] + log_inside[i][j];
            let conditional_probability = (mass - log_bound).exp().clamp(0.0, 1.0);
            let standard_probability = (mass - log_total).exp().clamp(0.0, 1.0);
            // Global one-based index for strand B is offset by strand A.
            conditional.push(PairProbability {
                i: i + 1,
                j: na + j + 1,
                probability: conditional_probability,
            });
            standard.push(PairProbability {
                i: i + 1,
                j: na + j + 1,
                probability: standard_probability,
            });
        }
    }

    Ok(DuplexResult {
        sequence_a: a,
        sequence_b: b,
        structure,
        mfe_energy_kcal_mol: best,
        bound_partition_function: if log_bound < f64::MAX.ln() {
            log_bound.exp()
        } else {
            f64::MAX
        },
        log_bound_partition_function: log_bound,
        association_ensemble_free_energy_kcal_mol: -rt * log_bound,
        standard_state_bound_probability: bound_probability,
        standard_state_ensemble_free_energy_kcal_mol: -rt * log_total,
        conditional_pair_probabilities: conditional,
        standard_state_pair_probabilities: standard,
        temperature_celsius,
        salt_molar,
        model: if model.parameter_profile_name().is_some() {
            "connected intermolecular custom nearest-neighbor duplex, dangles=0"
        } else {
            match model.nucleic_acid() {
                crate::energy::NucleicAcid::Rna => {
                    "connected intermolecular Turner 2004 RNA duplex, dangles=0"
                }
                crate::energy::NucleicAcid::Dna => {
                    "connected intermolecular Mathews 2004 DNA duplex, dangles=0"
                }
            }
        },
    })
}

fn log_add(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    let high = a.max(b);
    high + (a.min(b) - high).exp().ln_1p()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complementary_strands_form_antiparallel_connected_duplex() {
        let result = duplex("GGG", "CCC", 37.0, 1.021).unwrap();
        assert_eq!(result.structure, "(((&)))");
        assert!(result.mfe_energy_kcal_mol < 0.0);
        assert!((0.0..=1.0).contains(&result.standard_state_bound_probability));
        assert!(result
            .conditional_pair_probabilities
            .iter()
            .all(|pair| pair.probability <= 1.0));
    }

    #[test]
    fn rejects_strands_without_complementarity() {
        assert!(duplex("AAAA", "AAAA", 37.0, 1.021).is_err());
    }
}
