//! Observables derived exactly from McCaskill base-pair probabilities.

use crate::partition::PartitionResult;
use crate::structure::{parse_structure, RnaError};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct EnsembleSummary {
    /// Expected number of base pairs in a Boltzmann-distributed structure.
    pub expected_base_pairs: f64,
    /// Expected base-pair distance between two independent ensemble samples.
    pub mean_base_pair_distance: f64,
    /// Shannon entropy in bits for every nucleotide's unpaired/partner state.
    pub positional_entropy_bits: Vec<f64>,
    pub mean_positional_entropy_bits: f64,
    pub max_positional_entropy_bits: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct EnsembleDefectResult {
    pub target_structure: String,
    /// Nucleotide ensemble defect divided by sequence length.
    pub normalized_ensemble_defect: f64,
    pub expected_correct_nucleotides: f64,
}

pub fn summarize(partition: &PartitionResult) -> EnsembleSummary {
    let n = partition.unpaired_probabilities.len();
    let expected_base_pairs = partition
        .pair_probabilities
        .iter()
        .map(|pair| pair.probability)
        .sum();
    let mean_base_pair_distance = 2.0
        * partition
            .pair_probabilities
            .iter()
            .map(|pair| pair.probability * (1.0 - pair.probability))
            .sum::<f64>();
    let mut positional_entropy_bits = vec![0.0; n];
    for (position, &probability) in partition.unpaired_probabilities.iter().enumerate() {
        positional_entropy_bits[position] -= entropy_term(probability);
    }
    for pair in &partition.pair_probabilities {
        let term = entropy_term(pair.probability);
        positional_entropy_bits[pair.i - 1] -= term;
        positional_entropy_bits[pair.j - 1] -= term;
    }
    let mean_positional_entropy_bits = if n == 0 {
        0.0
    } else {
        positional_entropy_bits.iter().sum::<f64>() / n as f64
    };
    let max_positional_entropy_bits = positional_entropy_bits.iter().copied().fold(0.0, f64::max);
    EnsembleSummary {
        expected_base_pairs,
        mean_base_pair_distance,
        positional_entropy_bits,
        mean_positional_entropy_bits,
        max_positional_entropy_bits,
    }
}

pub fn ensemble_defect(
    sequence: &str,
    target_structure: &str,
    partition: &PartitionResult,
) -> Result<EnsembleDefectResult, RnaError> {
    let parsed = parse_structure(sequence, target_structure)?;
    if parsed.length != partition.unpaired_probabilities.len() {
        return Err(RnaError::InvalidOption(
            "partition result length does not match target structure".into(),
        ));
    }
    let mut correct = 0.0;
    for position in 0..parsed.length {
        if let Some(partner) = parsed.partner[position] {
            correct += pair_probability(partition, position + 1, partner + 1);
        } else {
            correct += partition.unpaired_probabilities[position];
        }
    }
    let defect = 1.0 - correct / parsed.length as f64;
    Ok(EnsembleDefectResult {
        target_structure: parsed.structure,
        normalized_ensemble_defect: defect.clamp(0.0, 1.0),
        expected_correct_nucleotides: correct,
    })
}

fn pair_probability(partition: &PartitionResult, i: usize, j: usize) -> f64 {
    let (i, j) = (i.min(j), i.max(j));
    partition
        .pair_probabilities
        .iter()
        .find(|pair| pair.i == i && pair.j == j)
        .map(|pair| pair.probability)
        .unwrap_or(0.0)
}

fn entropy_term(probability: f64) -> f64 {
    if probability > 0.0 {
        probability * probability.log2()
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{partition, EnergyModel};

    #[test]
    fn summary_and_defect_match_probability_definitions() {
        let sequence = "GGGAAACCC";
        let result = partition(sequence, 37.0, 3, &EnergyModel::default()).unwrap();
        let summary = summarize(&result);
        assert!(summary.expected_base_pairs > 0.0);
        assert!(summary.mean_base_pair_distance >= 0.0);
        assert_eq!(summary.positional_entropy_bits.len(), sequence.len());
        let defect = ensemble_defect(sequence, "(((...)))", &result).unwrap();
        assert!((0.0..=1.0).contains(&defect.normalized_ensemble_defect));
        assert!(defect.expected_correct_nucleotides > 0.0);
    }
}
