//! Exact k-best suboptimal structures for the even-dangle Turner grammar.

use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::energy::EnergyModel;
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::Serialize;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;

#[derive(Clone, Debug, Serialize)]
pub struct SuboptimalStructure {
    pub rank: usize,
    pub structure: String,
    pub energy_kcal_mol: f64,
    pub delta_energy_kcal_mol: f64,
    pub relative_boltzmann_weight: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct SuboptimalResult {
    pub sequence: String,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub energy_band_kcal_mol: f64,
    pub requested_limit: usize,
    pub truncated: bool,
    pub structures: Vec<SuboptimalStructure>,
    pub constraints: ConstraintSummary,
    pub method: &'static str,
}

#[derive(Clone, Debug)]
struct Candidate {
    energy: f64,
    pairs: Vec<(usize, usize)>,
}

type Table = Vec<Vec<Vec<Candidate>>>;

/// Return the exact lowest-energy structures retained by a k-best dynamic
/// program, filtered to `energy_band_kcal_mol` above the MFE.
pub fn suboptimal_structures(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    energy_band_kcal_mol: f64,
    limit: usize,
) -> Result<SuboptimalResult, RnaError> {
    if !energy_band_kcal_mol.is_finite() || energy_band_kcal_mol < 0.0 {
        return Err(RnaError::InvalidOption(
            "suboptimal energy band must be finite and non-negative".into(),
        ));
    }
    if limit == 0 {
        return Err(RnaError::InvalidOption(
            "suboptimal structure limit must be positive".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    if model.dangles() % 2 == 1 {
        return suboptimal_odd_exact(
            &sequence,
            min_loop,
            model,
            constraints,
            energy_band_kcal_mol,
            limit,
        );
    }
    let bases = sequence.as_bytes();
    let n = bases.len();
    let keep = limit.saturating_add(1);
    let mut f: Table = vec![vec![Vec::new(); n]; n];
    let mut v: Table = vec![vec![Vec::new(); n]; n];
    let mut v_secure: Table = vec![vec![Vec::new(); n]; n];
    let mut m1: Table = vec![vec![Vec::new(); n]; n];
    let mut m2: Table = vec![vec![Vec::new(); n]; n];
    for (i, row) in f.iter_mut().enumerate() {
        if constraints.allows_unpaired(i) {
            row[i].push(Candidate {
                energy: constraints.unpaired_energy(i),
                pairs: Vec::new(),
            });
        }
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            if span > min_loop && constraints.allows_pair(bases, i, j, model) {
                let pair_energy =
                    model.pair_energy(bases[i], bases[j]) + constraints.pair_energy(i, j);
                let mut choices = Vec::new();
                if let Some(soft) = constraints.unpaired_range_energy(i + 1, j - 1) {
                    let energy = pair_energy + model.hairpin_energy(bases, i, j) + soft;
                    if energy.is_finite() {
                        choices.push(Candidate {
                            energy,
                            pairs: vec![(i, j)],
                        });
                    }
                }
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    for right in 0..=(internal_limit - left) {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        let stacked = left == 0 && right == 0;
                        let children = if constraints.no_lonely_pairs() && !stacked {
                            &v_secure[k][l]
                        } else {
                            &v[k][l]
                        };
                        if k >= l || l - k <= min_loop || children.is_empty() {
                            continue;
                        }
                        let Some(left_soft) = constraints.unpaired_range_energy(i + 1, k - 1)
                        else {
                            continue;
                        };
                        let Some(right_soft) = constraints.unpaired_range_energy(l + 1, j - 1)
                        else {
                            continue;
                        };
                        let transition = pair_energy
                            + model.internal_energy(bases, i, j, k, l)
                            + left_soft
                            + right_soft
                            + if stacked {
                                constraints.stack_energy(i, j, k, l)
                            } else {
                                0.0
                            };
                        for child in children {
                            choices.push(wrap_pair(child, i, j, transition));
                        }
                    }
                }
                if i + 1 < j {
                    let transition = pair_energy
                        + model.multiloop_closing()
                        + model.multiloop_closing_stem_energy(bases, i, j);
                    for child in &m2[i + 1][j - 1] {
                        choices.push(wrap_pair(child, i, j, transition));
                    }
                }
                v[i][j] = best(choices, keep);

                if i + 1 < j {
                    let transition = pair_energy
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1);
                    v_secure[i][j] = best(
                        v[i + 1][j - 1]
                            .iter()
                            .map(|child| wrap_pair(child, i, j, transition))
                            .collect(),
                        keep,
                    );
                }
            }

            let mut exterior = Vec::new();
            if constraints.allows_unpaired(j) {
                exterior.extend(f[i][j - 1].iter().map(|child| Candidate {
                    energy: child.energy + constraints.unpaired_energy(j),
                    pairs: child.pairs.clone(),
                }));
            }
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branches = if constraints.no_lonely_pairs() {
                        &v_secure[k][j]
                    } else {
                        &v[k][j]
                    };
                    let transition = model.exterior_stem_energy(bases, k, j);
                    if k == i {
                        exterior.extend(branches.iter().map(|branch| Candidate {
                            energy: branch.energy + transition,
                            pairs: branch.pairs.clone(),
                        }));
                    } else {
                        combine(&mut exterior, &f[i][k - 1], branches, transition);
                    }
                }
            }
            f[i][j] = best(exterior, keep);

            let mut one = Vec::new();
            let mut two = Vec::new();
            if constraints.allows_unpaired(j) {
                let transition = model.multiloop_unpaired() + constraints.unpaired_energy(j);
                one.extend(m1[i][j - 1].iter().map(|child| Candidate {
                    energy: child.energy + transition,
                    pairs: child.pairs.clone(),
                }));
                two.extend(m2[i][j - 1].iter().map(|child| Candidate {
                    energy: child.energy + transition,
                    pairs: child.pairs.clone(),
                }));
            }
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branches = if constraints.no_lonely_pairs() {
                        &v_secure[k][j]
                    } else {
                        &v[k][j]
                    };
                    let branch_transition = model.multiloop_stem_energy(bases, k, j);
                    let leading_soft = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    };
                    if let Some(soft) = leading_soft {
                        let transition =
                            branch_transition + (k - i) as f64 * model.multiloop_unpaired() + soft;
                        one.extend(branches.iter().map(|branch| Candidate {
                            energy: branch.energy + transition,
                            pairs: branch.pairs.clone(),
                        }));
                    }
                    if k > i {
                        combine(&mut one, &m1[i][k - 1], branches, branch_transition);
                        combine(&mut two, &m1[i][k - 1], branches, branch_transition);
                    }
                }
            }
            m1[i][j] = best(one, keep);
            m2[i][j] = best(two, keep);
        }
    }

    let root = &f[0][n - 1];
    let Some(mfe) = root.first().map(|entry| entry.energy) else {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    };
    let rt = GAS_CONSTANT_KCAL * (model.temperature_celsius() + 273.15);
    let retained: Vec<_> = root
        .iter()
        .take(limit)
        .filter(|entry| entry.energy <= mfe + energy_band_kcal_mol + 1.0e-9)
        .collect();
    let structures = retained
        .iter()
        .enumerate()
        .map(|(rank, entry)| {
            let delta = entry.energy - mfe;
            SuboptimalStructure {
                rank: rank + 1,
                structure: pairs_to_dot_bracket(n, &entry.pairs),
                energy_kcal_mol: entry.energy,
                delta_energy_kcal_mol: delta,
                relative_boltzmann_weight: (-delta / rt).exp(),
            }
        })
        .collect();
    let truncated = root.len() > limit
        && root
            .get(limit)
            .is_some_and(|entry| entry.energy <= mfe + energy_band_kcal_mol + 1.0e-9);

    Ok(SuboptimalResult {
        sequence,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        energy_band_kcal_mol,
        requested_limit: limit,
        truncated,
        structures,
        constraints: constraints.summary(),
        method: "exact k-best Turner dynamic program",
    })
}

fn suboptimal_odd_exact(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    energy_band_kcal_mol: f64,
    limit: usize,
) -> Result<SuboptimalResult, RnaError> {
    let n = sequence.len();
    let mut candidates = Vec::new();
    crate::exact_enumeration::for_each_noncrossing_structure(
        sequence.as_bytes(),
        min_loop,
        model,
        constraints,
        |pairs| {
            let mut partner = vec![None; n];
            for &(i, j) in pairs {
                partner[i] = Some(j);
                partner[j] = Some(i);
            }
            if constraints
                .validate_structure(sequence.as_bytes(), &partner, min_loop, model)
                .is_err()
            {
                return Ok(());
            }
            let structure = pairs_to_dot_bracket(n, pairs);
            let energy = model
                .evaluate_with_constraints(sequence, &structure, min_loop, constraints)?
                .total_kcal_mol;
            candidates.push((energy, structure));
            Ok(())
        },
    )?;
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    let Some(mfe) = candidates.first().map(|candidate| candidate.0) else {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    };
    let in_band = candidates
        .into_iter()
        .take_while(|(energy, _)| *energy <= mfe + energy_band_kcal_mol + 1.0e-12)
        .collect::<Vec<_>>();
    let truncated = in_band.len() > limit;
    let rt = GAS_CONSTANT_KCAL * (model.temperature_celsius() + 273.15);
    let structures = in_band
        .into_iter()
        .take(limit)
        .enumerate()
        .map(|(index, (energy, structure))| SuboptimalStructure {
            rank: index + 1,
            structure,
            energy_kcal_mol: energy,
            delta_energy_kcal_mol: energy - mfe,
            relative_boltzmann_weight: (-(energy - mfe) / rt).exp(),
        })
        .collect();
    Ok(SuboptimalResult {
        sequence: sequence.into(),
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        energy_band_kcal_mol,
        requested_limit: limit,
        truncated,
        structures,
        constraints: constraints.summary(),
        method: "exact exhaustive fixed-structure odd-dangle energy ordering",
    })
}

fn wrap_pair(child: &Candidate, i: usize, j: usize, transition_energy: f64) -> Candidate {
    let mut pairs = Vec::with_capacity(child.pairs.len() + 1);
    pairs.push((i, j));
    pairs.extend_from_slice(&child.pairs);
    Candidate {
        energy: child.energy + transition_energy,
        pairs,
    }
}

fn combine(
    target: &mut Vec<Candidate>,
    left: &[Candidate],
    right: &[Candidate],
    transition_energy: f64,
) {
    for first in left {
        for second in right {
            let mut pairs = Vec::with_capacity(first.pairs.len() + second.pairs.len());
            pairs.extend_from_slice(&first.pairs);
            pairs.extend_from_slice(&second.pairs);
            target.push(Candidate {
                energy: first.energy + second.energy + transition_energy,
                pairs,
            });
        }
    }
}

fn best(mut candidates: Vec<Candidate>, keep: usize) -> Vec<Candidate> {
    candidates.retain(|candidate| candidate.energy.is_finite());
    candidates.sort_by(|a, b| {
        a.energy
            .total_cmp(&b.energy)
            .then_with(|| a.pairs.cmp(&b.pairs))
    });
    candidates.dedup_by(|a, b| a.pairs == b.pairs);
    candidates.truncate(keep);
    candidates
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::parse_structure;

    fn enumerate(
        sequence: &str,
        model: &EnergyModel,
        i: usize,
        j: usize,
    ) -> Vec<Vec<(usize, usize)>> {
        if i >= j {
            return vec![Vec::new()];
        }
        let mut result = enumerate(sequence, model, i, j - 1);
        for k in i..j.saturating_sub(3) {
            if !model.can_pair(sequence.as_bytes()[k], sequence.as_bytes()[j]) {
                continue;
            }
            for mut inside in enumerate(sequence, model, k + 1, j - 1) {
                inside.push((k, j));
                for mut left in if k == i {
                    vec![Vec::new()]
                } else {
                    enumerate(sequence, model, i, k - 1)
                } {
                    left.extend_from_slice(&inside);
                    left.sort_unstable();
                    result.push(left);
                }
            }
        }
        result.sort();
        result.dedup();
        result
    }

    #[test]
    fn k_best_matches_exhaustive_structure_energies() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(0.0 + 37.0, 0).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let result = suboptimal_structures(sequence, 3, &model, &constraints, 100.0, 200).unwrap();
        let mut expected: Vec<_> = enumerate(sequence, &model, 0, sequence.len() - 1)
            .into_iter()
            .map(|pairs| {
                let structure = pairs_to_dot_bracket(sequence.len(), &pairs);
                let parsed = parse_structure(sequence, &structure).unwrap();
                let energy = model
                    .evaluate(sequence, &parsed.structure)
                    .unwrap()
                    .total_kcal_mol;
                (energy, structure)
            })
            .collect();
        expected.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        expected.dedup_by(|a, b| a.1 == b.1);
        assert_eq!(result.structures.len(), expected.len());
        for (actual, (energy, structure)) in result.structures.iter().zip(expected) {
            assert_eq!(actual.structure, structure);
            assert!((actual.energy_kcal_mol - energy).abs() < 1.0e-9);
        }
    }

    #[test]
    fn energy_band_and_limit_are_reported() {
        let model = EnergyModel::with_dangles(37.0, 2).unwrap();
        let constraints = ConstraintModel::unconstrained(9);
        let result = suboptimal_structures("GGGAAACCC", 3, &model, &constraints, 0.0, 1).unwrap();
        assert_eq!(result.structures.len(), 1);
        assert!(result.structures[0].delta_energy_kcal_mol.abs() < 1.0e-12);
    }

    #[test]
    fn odd_dangle_suboptimals_are_ordered_by_fixed_structure_energy() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 3).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let result = suboptimal_structures(sequence, 3, &model, &constraints, 20.0, 32).unwrap();
        assert!(!result.structures.is_empty());
        assert!(result
            .structures
            .windows(2)
            .all(|window| window[0].energy_kcal_mol <= window[1].energy_kcal_mol));
        for entry in result.structures {
            let evaluated = model.evaluate(sequence, &entry.structure).unwrap();
            assert!((entry.energy_kcal_mol - evaluated.total_kcal_mol).abs() < 1e-12);
        }
    }
}
