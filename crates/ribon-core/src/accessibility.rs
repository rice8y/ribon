//! Exact equilibrium opening probabilities for selected contiguous regions.
//!
//! A region is accessible when every nucleotide in it is unpaired.  Its
//! probability is evaluated as `Z(force-unpaired region) / Z`, using the same
//! Turner model, hard/soft constraints, and log-domain partition function as
//! the rest of Ribon.  Single-nucleotide queries reuse the already exact
//! unpaired marginals from the baseline partition.

use crate::constraints::{ConstraintConfig, ConstraintModel};
use crate::energy::EnergyModel;
use crate::partition::partition_with_constraints;
use crate::structure::{normalize_sequence, RnaError};
use serde::{Deserialize, Serialize};

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct AccessibilityWindow {
    /// One-based inclusive start.
    pub from: usize,
    /// One-based inclusive end.
    pub to: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct AccessibilityEntry {
    pub from: usize,
    pub to: usize,
    pub length: usize,
    pub probability_unpaired: f64,
    /// `None` denotes an impossible opening (infinite opening energy).
    pub opening_energy_kcal_mol: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct AccessibilityResult {
    pub sequence: String,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub baseline_log_partition_function: f64,
    pub windows: Vec<AccessibilityEntry>,
    pub method: &'static str,
}

/// Compute exact joint unpaired probabilities for explicitly selected windows.
pub fn accessibility(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    config: &ConstraintConfig,
    windows: &[AccessibilityWindow],
) -> Result<AccessibilityResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    accessibility_with_model(sequence, min_loop, &model, config, windows)
}

/// Compute accessibility with an explicitly selected thermodynamic family.
pub fn accessibility_with_model(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    config: &ConstraintConfig,
    windows: &[AccessibilityWindow],
) -> Result<AccessibilityResult, RnaError> {
    let sequence = normalize_sequence(sequence)?;
    let n = sequence.len();
    for window in windows {
        if window.from == 0 || window.from > window.to || window.to > n {
            return Err(RnaError::InvalidOption(format!(
                "accessibility window {}-{} is outside 1-{n}",
                window.from, window.to
            )));
        }
    }

    let temperature_celsius = model.temperature_celsius();
    let dangles = model.dangles();
    let salt_molar = model.salt_molar();
    let baseline_constraints = ConstraintModel::compile(n, config)?;
    let baseline = partition_with_constraints(
        &sequence,
        temperature_celsius,
        min_loop,
        model,
        &baseline_constraints,
    )?;
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut entries = Vec::with_capacity(windows.len());
    for window in windows {
        let length = window.to - window.from + 1;
        let log_probability = if length == 1 {
            let probability = baseline.unpaired_probabilities[window.from - 1];
            if probability > 0.0 {
                probability.ln()
            } else {
                f64::NEG_INFINITY
            }
        } else {
            let mut conditioned = config.clone();
            conditioned.force_unpaired.extend(window.from..=window.to);
            conditioned.force_unpaired.sort_unstable();
            conditioned.force_unpaired.dedup();
            match ConstraintModel::compile(n, &conditioned) {
                Ok(constraints) => {
                    let partition = partition_with_constraints(
                        &sequence,
                        temperature_celsius,
                        min_loop,
                        model,
                        &constraints,
                    )?;
                    let difference =
                        partition.log_partition_function - baseline.log_partition_function;
                    if difference > 1.0e-10 {
                        return Err(RnaError::Numerical(format!(
                            "conditioned accessibility partition exceeds baseline by {difference:e}"
                        )));
                    }
                    difference.min(0.0)
                }
                // The baseline configuration has already compiled. Therefore
                // an error after adding only force-unpaired positions denotes
                // an impossible opening, not a numerical PF failure.
                Err(_) => f64::NEG_INFINITY,
            }
        };
        let probability = if log_probability.is_finite() {
            log_probability.exp().clamp(0.0, 1.0)
        } else {
            0.0
        };
        entries.push(AccessibilityEntry {
            from: window.from,
            to: window.to,
            length,
            probability_unpaired: probability,
            opening_energy_kcal_mol: log_probability.is_finite().then_some(-rt * log_probability),
        });
    }

    Ok(AccessibilityResult {
        sequence,
        temperature_celsius,
        dangles,
        salt_molar,
        baseline_log_partition_function: baseline.log_partition_function,
        windows: entries,
        method: "exact constrained-partition ratio",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::PairConstraint;

    #[test]
    fn single_base_queries_equal_partition_marginals() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::default();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let partition =
            partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        let windows: Vec<_> = (1..=sequence.len())
            .map(|position| AccessibilityWindow {
                from: position,
                to: position,
            })
            .collect();
        let result = accessibility(
            sequence,
            37.0,
            3,
            2,
            1.021,
            &ConstraintConfig::default(),
            &windows,
        )
        .unwrap();
        for (entry, expected) in result.windows.iter().zip(partition.unpaired_probabilities) {
            assert!((entry.probability_unpaired - expected).abs() < 1.0e-12);
        }
    }

    #[test]
    fn joint_opening_probability_obeys_marginals_and_hard_constraints() {
        let windows = [
            AccessibilityWindow { from: 1, to: 2 },
            AccessibilityWindow { from: 1, to: 1 },
            AccessibilityWindow { from: 2, to: 2 },
        ];
        let result = accessibility(
            "GGGAAACCC",
            37.0,
            3,
            2,
            1.021,
            &ConstraintConfig::default(),
            &windows,
        )
        .unwrap();
        assert!(result.windows[0].probability_unpaired <= result.windows[1].probability_unpaired);
        assert!(result.windows[0].probability_unpaired <= result.windows[2].probability_unpaired);

        let constrained = ConstraintConfig {
            force_pairs: vec![PairConstraint { i: 1, j: 9 }],
            ..ConstraintConfig::default()
        };
        let impossible = accessibility(
            "GGGAAACCC",
            37.0,
            3,
            2,
            1.021,
            &constrained,
            &[AccessibilityWindow { from: 1, to: 2 }],
        )
        .unwrap();
        assert_eq!(impossible.windows[0].probability_unpaired, 0.0);
        assert_eq!(impossible.windows[0].opening_energy_kcal_mol, None);
    }
}
