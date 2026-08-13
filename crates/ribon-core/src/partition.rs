use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::energy::EnergyModel;
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::Serialize;

// Molar gas constant in kcal mol^-1 K^-1.
const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Debug, Serialize)]
pub struct PairProbability {
    /// One-based nucleotide index.
    pub i: usize,
    /// One-based nucleotide index.
    pub j: usize,
    pub probability: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct PartitionResult {
    /// The ordinary partition function. Values larger than `f64::MAX` are
    /// saturated; use `log_partition_function` for numerically stable work.
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub model: &'static str,
    pub salt_molar: f64,
    pub constraints: ConstraintSummary,
}

fn log_empty_or(table: &[Vec<f64>], i: usize, j: usize) -> f64 {
    if i > j {
        0.0
    } else {
        table[i][j]
    }
}

fn log_add(a: f64, b: f64) -> f64 {
    if a == NEG_INF {
        return b;
    }
    if b == NEG_INF {
        return a;
    }
    let high = a.max(b);
    let low = a.min(b);
    high + (low - high).exp().ln_1p()
}

fn update(target: &mut f64, value: f64) {
    *target = log_add(*target, value);
}

/// McCaskill-style partition function over the same Turner loop decomposition
/// used by MFE and supplied-structure evaluation.
///
/// The grammar is unambiguous: exterior and multiloop segments are decomposed
/// by their rightmost branch. `m1` contains one-or-more branches and `m2`
/// contains two-or-more branches. All calculations are performed in log space.
pub fn partition(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<PartitionResult, RnaError> {
    let length = normalize_sequence(sequence)?.len();
    let constraints = ConstraintModel::unconstrained(length);
    partition_with_constraints(sequence, temperature_celsius, min_loop, model, &constraints)
}

pub fn partition_with_constraints(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<PartitionResult, RnaError> {
    debug_assert!(model.supports_partition());
    let sequence = normalize_sequence(sequence)?;
    let bases = sequence.as_bytes();
    let n = bases.len();
    if (model.temperature_celsius() - temperature_celsius).abs() > 1.0e-9 {
        return Err(RnaError::InvalidOption(
            "partition temperature and energy-model temperature differ".into(),
        ));
    }
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    if !rt.is_finite() || rt <= 0.0 {
        return Err(RnaError::InvalidOption(
            "temperature must be finite and above absolute zero".into(),
        ));
    }
    if model.dangles() % 2 == 1 {
        return partition_odd_exact(&sequence, min_loop, model, constraints, rt);
    }

    // q: arbitrary exterior structure, qb: structure closed by i-j,
    // m1/m2: multiloop segments with >=1 / >=2 branches.
    let mut q = vec![vec![NEG_INF; n]; n];
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut qb_secure = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    // ViennaRNA applies noLP as a stack-capable pair filter in PF recursions,
    // not as the secure-pair state used by its MFE grammar.
    let strict_no_lp = false;
    for (i, row) in q.iter_mut().enumerate() {
        if constraints.allows_unpaired(i) {
            row[i] = -constraints.unpaired_energy(i) / rt;
        }
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;

            if span > min_loop
                && constraints.allows_pair_for_partition(bases, i, j, min_loop, model)
            {
                let pair_soft = constraints.pair_energy(i, j);
                let hairpin = model.hairpin_boltzmann_energy(bases, i, j);
                if let Some(unpaired_soft) = constraints.unpaired_range_energy(i + 1, j - 1) {
                    if hairpin.is_finite() {
                        qb[i][j] = -(hairpin + pair_soft + unpaired_soft) / rt;
                    }
                }

                // All geometrically possible internal loops, including
                // stacks and bulges, unless the caller requests a bound.
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    let remaining = internal_limit - left;
                    for right in 0..=remaining {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        let stacked = left == 0 && right == 0;
                        let child = if strict_no_lp && !stacked {
                            qb_secure[k][l]
                        } else {
                            qb[k][l]
                        };
                        if k >= l || l - k <= min_loop || child == NEG_INF {
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
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l)
                            + pair_soft
                            + left_soft
                            + right_soft
                            + constraints.stack_energy(i, j, k, l);
                        if energy.is_finite() {
                            let candidate = child - energy / rt;
                            update(&mut qb[i][j], candidate);
                        }
                    }
                }

                if span >= 2 && m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
                        + constraints.terminal_pair_energy(i, j)
                        + pair_soft;
                    update(&mut qb[i][j], m2[i + 1][j - 1] - energy / rt);
                }
                if i + 1 < j && qb[i + 1][j - 1] != NEG_INF {
                    let energy = model.internal_boltzmann_energy(bases, i, j, i + 1, j - 1)
                        + pair_soft
                        + constraints.stack_energy(i, j, i + 1, j - 1);
                    qb_secure[i][j] = qb[i + 1][j - 1] - energy / rt;
                }
            }

            // Multiloop segments. A rightmost branch k-j has its own stem
            // energy; qb supplies everything enclosed by that branch.
            let log_unpaired_ml = if constraints.allows_unpaired(j) {
                -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
            } else {
                NEG_INF
            };
            if m1[i][j - 1] != NEG_INF && log_unpaired_ml != NEG_INF {
                m1[i][j] = m1[i][j - 1] + log_unpaired_ml;
            }
            if m2[i][j - 1] != NEG_INF && log_unpaired_ml != NEG_INF {
                m2[i][j] = m2[i][j - 1] + log_unpaired_ml;
            }
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branch_qb = if strict_no_lp {
                        qb_secure[k][j]
                    } else {
                        qb[k][j]
                    };
                    if branch_qb == NEG_INF {
                        continue;
                    }
                    let branch = branch_qb
                        - (model.multiloop_stem_boltzmann_energy(bases, k, j)
                            + constraints.terminal_pair_energy(k, j))
                            / rt;
                    let leading_soft = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    };
                    if let Some(soft) = leading_soft {
                        let unpaired =
                            -((k - i) as f64 * model.multiloop_unpaired_boltzmann() + soft) / rt;
                        update(&mut m1[i][j], unpaired + branch);
                    }
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        update(&mut m1[i][j], additional);
                        update(&mut m2[i][j], additional);
                    }
                }
            }

            // Exterior rightmost decomposition.
            if constraints.allows_unpaired(j) && q[i][j - 1] != NEG_INF {
                q[i][j] = q[i][j - 1] - constraints.unpaired_energy(j) / rt;
            }
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    let branch_qb = if strict_no_lp {
                        qb_secure[k][j]
                    } else {
                        qb_row[j]
                    };
                    if branch_qb == NEG_INF {
                        continue;
                    }
                    let left = if k == i {
                        0.0
                    } else {
                        log_empty_or(&q, i, k - 1)
                    };
                    let stem = -(model.exterior_stem_boltzmann_energy(bases, k, j)
                        + constraints.terminal_pair_energy(k, j))
                        / rt;
                    update(&mut q[i][j], left + stem + branch_qb);
                }
            }
        }
    }

    let log_z = q[0][n - 1];
    if log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    }

    // Reverse-mode differentiation of the four inside tables. Each outside
    // table stores dZ(root)/d(state) in log space.
    let mut oq = vec![vec![NEG_INF; n]; n];
    let mut oqb = vec![vec![NEG_INF; n]; n];
    let mut oqb_secure = vec![vec![NEG_INF; n]; n];
    let mut om1 = vec![vec![NEG_INF; n]; n];
    let mut om2 = vec![vec![NEG_INF; n]; n];
    oq[0][n - 1] = 0.0;

    for span in (1..n).rev() {
        for i in 0..(n - span) {
            let j = i + span;

            // q(i,j)
            let parent = oq[i][j];
            if parent != NEG_INF {
                if constraints.allows_unpaired(j) && q[i][j - 1] != NEG_INF {
                    update(
                        &mut oq[i][j - 1],
                        parent - constraints.unpaired_energy(j) / rt,
                    );
                }
                if j > min_loop {
                    for k in i..=(j - min_loop - 1) {
                        let branch_qb = if strict_no_lp {
                            qb_secure[k][j]
                        } else {
                            qb[k][j]
                        };
                        if branch_qb == NEG_INF {
                            continue;
                        }
                        let left = if k == i {
                            0.0
                        } else {
                            log_empty_or(&q, i, k - 1)
                        };
                        let stem = -(model.exterior_stem_boltzmann_energy(bases, k, j)
                            + constraints.terminal_pair_energy(k, j))
                            / rt;
                        if strict_no_lp {
                            update(&mut oqb_secure[k][j], parent + left + stem);
                        } else {
                            update(&mut oqb[k][j], parent + left + stem);
                        }
                        if k > i {
                            update(&mut oq[i][k - 1], parent + stem + branch_qb);
                        }
                    }
                }
            }

            // m2(i,j)
            let parent = om2[i][j];
            if parent != NEG_INF {
                let log_unpaired_ml = if constraints.allows_unpaired(j) {
                    -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
                } else {
                    NEG_INF
                };
                if m2[i][j - 1] != NEG_INF && log_unpaired_ml != NEG_INF {
                    update(&mut om2[i][j - 1], parent + log_unpaired_ml);
                }
                if j > min_loop {
                    for k in (i + 1)..=(j - min_loop - 1) {
                        let branch_qb = if strict_no_lp {
                            qb_secure[k][j]
                        } else {
                            qb[k][j]
                        };
                        if m1[i][k - 1] == NEG_INF || branch_qb == NEG_INF {
                            continue;
                        }
                        let stem = -(model.multiloop_stem_boltzmann_energy(bases, k, j)
                            + constraints.terminal_pair_energy(k, j))
                            / rt;
                        update(&mut om1[i][k - 1], parent + stem + branch_qb);
                        if strict_no_lp {
                            update(&mut oqb_secure[k][j], parent + m1[i][k - 1] + stem);
                        } else {
                            update(&mut oqb[k][j], parent + m1[i][k - 1] + stem);
                        }
                    }
                }
            }

            // m1(i,j)
            let parent = om1[i][j];
            if parent != NEG_INF {
                let log_unpaired_ml = if constraints.allows_unpaired(j) {
                    -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
                } else {
                    NEG_INF
                };
                if m1[i][j - 1] != NEG_INF && log_unpaired_ml != NEG_INF {
                    update(&mut om1[i][j - 1], parent + log_unpaired_ml);
                }
                if j > min_loop {
                    for k in i..=(j - min_loop - 1) {
                        let branch_qb = if strict_no_lp {
                            qb_secure[k][j]
                        } else {
                            qb[k][j]
                        };
                        if branch_qb == NEG_INF {
                            continue;
                        }
                        let stem = -(model.multiloop_stem_boltzmann_energy(bases, k, j)
                            + constraints.terminal_pair_energy(k, j))
                            / rt;
                        let leading_soft = if k == i {
                            Some(0.0)
                        } else {
                            constraints.unpaired_range_energy(i, k - 1)
                        };
                        if let Some(soft) = leading_soft {
                            let unpaired = -((k - i) as f64 * model.multiloop_unpaired_boltzmann()
                                + soft)
                                / rt;
                            if strict_no_lp {
                                update(&mut oqb_secure[k][j], parent + unpaired + stem);
                            } else {
                                update(&mut oqb[k][j], parent + unpaired + stem);
                            }
                        }
                        if k > i && m1[i][k - 1] != NEG_INF {
                            update(&mut om1[i][k - 1], parent + stem + branch_qb);
                            if strict_no_lp {
                                update(&mut oqb_secure[k][j], parent + m1[i][k - 1] + stem);
                            } else {
                                update(&mut oqb[k][j], parent + m1[i][k - 1] + stem);
                            }
                        }
                    }
                }
            }

            // qb_secure(i,j): the root pair is guaranteed non-lonely by its
            // immediately enclosed stack; the child may rely on that outer
            // stack and therefore uses the ordinary qb state.
            let secure_parent = oqb_secure[i][j];
            if strict_no_lp
                && secure_parent != NEG_INF
                && qb_secure[i][j] != NEG_INF
                && i + 1 < j
                && qb[i + 1][j - 1] != NEG_INF
            {
                let energy = model.internal_boltzmann_energy(bases, i, j, i + 1, j - 1)
                    + constraints.pair_energy(i, j)
                    + constraints.stack_energy(i, j, i + 1, j - 1);
                update(&mut oqb[i + 1][j - 1], secure_parent - energy / rt);
            }

            // qb(i,j)
            let parent = oqb[i][j];
            if parent != NEG_INF && qb[i][j] != NEG_INF {
                let pair_soft = constraints.pair_energy(i, j);
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    let remaining = internal_limit - left;
                    for right in 0..=remaining {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        let stacked = left == 0 && right == 0;
                        let child = if strict_no_lp && !stacked {
                            qb_secure[k][l]
                        } else {
                            qb[k][l]
                        };
                        if k >= l || l - k <= min_loop || child == NEG_INF {
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
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l)
                            + pair_soft
                            + left_soft
                            + right_soft
                            + constraints.stack_energy(i, j, k, l);
                        if energy.is_finite() {
                            if strict_no_lp && !stacked {
                                update(&mut oqb_secure[k][l], parent - energy / rt);
                            } else {
                                update(&mut oqb[k][l], parent - energy / rt);
                            }
                        }
                    }
                }
                if span >= 2 && m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
                        + constraints.terminal_pair_energy(i, j)
                        + pair_soft;
                    update(&mut om2[i + 1][j - 1], parent - energy / rt);
                }
            }
        }
    }

    let mut pair_probabilities = Vec::new();
    let mut paired_mass = vec![0.0f64; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let ordinary = if oqb[i][j] == NEG_INF || qb[i][j] == NEG_INF {
                NEG_INF
            } else {
                oqb[i][j] + qb[i][j]
            };
            let secure = if oqb_secure[i][j] == NEG_INF || qb_secure[i][j] == NEG_INF {
                NEG_INF
            } else {
                oqb_secure[i][j] + qb_secure[i][j]
            };
            let pair_log_weight = log_add(ordinary, secure);
            if pair_log_weight == NEG_INF {
                continue;
            }
            let probability = (pair_log_weight - log_z).exp().clamp(0.0, 1.0);
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
        .map(|probability| (1.0 - probability).clamp(0.0, 1.0))
        .collect();
    let partition_function = if log_z < f64::MAX.ln() {
        log_z.exp()
    } else {
        f64::MAX
    };

    Ok(PartitionResult {
        partition_function,
        log_partition_function: log_z,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        pair_probabilities,
        unpaired_probabilities,
        model: model.ensemble_model_name(),
        salt_molar: model.salt_molar(),
        constraints: constraints.summary(),
    })
}

fn partition_odd_exact(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    rt: f64,
) -> Result<PartitionResult, RnaError> {
    let n = sequence.len();
    let mut log_z = NEG_INF;
    let mut pair_log_mass = vec![vec![NEG_INF; n]; n];
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
            let log_weight = -energy / rt;
            log_z = log_add(log_z, log_weight);
            for &(i, j) in pairs {
                pair_log_mass[i][j] = log_add(pair_log_mass[i][j], log_weight);
            }
            Ok(())
        },
    )?;
    if log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
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
        .collect();
    Ok(PartitionResult {
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        pair_probabilities,
        unpaired_probabilities,
        model: if model.dangles() == 1 {
            "exact fixed-structure Turner ensemble, exclusive single dangles"
        } else {
            "exact fixed-structure Turner ensemble, single dangles and coaxial stacking"
        },
        salt_molar: model.salt_molar(),
        constraints: constraints.summary(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::ConstraintConfig;
    use crate::structure::pairs_to_dot_bracket;

    #[test]
    fn odd_dangle_partition_uses_the_exact_fixed_structure_ensemble() {
        for dangles in [1, 3] {
            let model = EnergyModel::with_dangles(37.0, dangles).unwrap();
            let odd = partition("GGGAAACCC", 37.0, 3, &model).unwrap();
            assert!(odd.log_partition_function.is_finite());
            assert!(odd.model.contains("exact fixed-structure"));
            for (index, &unpaired) in odd.unpaired_probabilities.iter().enumerate() {
                let paired = odd
                    .pair_probabilities
                    .iter()
                    .filter(|pair| pair.i == index + 1 || pair.j == index + 1)
                    .map(|pair| pair.probability)
                    .sum::<f64>();
                assert!((paired + unpaired - 1.0).abs() < 1.0e-10);
            }
        }
    }

    fn enumerate_matchings(
        bases: &[u8],
        model: &EnergyModel,
        min_loop: usize,
        i: usize,
        j: usize,
    ) -> Vec<Vec<(usize, usize)>> {
        if i >= j {
            return vec![Vec::new()];
        }
        let mut results = enumerate_matchings(bases, model, min_loop, i, j - 1);
        if j > min_loop {
            for k in i..=(j - min_loop - 1) {
                if !model.can_pair(bases[k], bases[j]) {
                    continue;
                }
                let left = if k == i {
                    vec![Vec::new()]
                } else {
                    enumerate_matchings(bases, model, min_loop, i, k - 1)
                };
                let inner = enumerate_matchings(bases, model, min_loop, k + 1, j - 1);
                for prefix in &left {
                    for enclosed in &inner {
                        let mut pairs = prefix.clone();
                        pairs.extend(enclosed.iter().copied());
                        pairs.push((k, j));
                        results.push(pairs);
                    }
                }
            }
        }
        results
    }

    #[test]
    fn probabilities_are_normalized_per_base() {
        let model = EnergyModel::default();
        let result = partition("GGGAAACCC", 37.0, 3, &model).unwrap();
        assert!(result.partition_function >= 1.0);
        for (index, &unpaired) in result.unpaired_probabilities.iter().enumerate() {
            let paired: f64 = result
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == index + 1 || pair.j == index + 1)
                .map(|pair| pair.probability)
                .sum();
            assert!((paired + unpaired - 1.0).abs() < 1.0e-8);
        }
    }

    #[test]
    fn zero_dangle_salt_grammar_matches_complete_structure_enumeration() {
        let sequence = "GGAAACGAAACGAAACC";
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let rt = GAS_CONSTANT_KCAL * (37.0 + 273.15);
        for salt_molar in [0.05, 0.2, 0.5, 2.0] {
            let model = EnergyModel::with_dangles_and_salt(37.0, 0, salt_molar).unwrap();
            let dynamic =
                partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
            let enumerated = partition_odd_exact(sequence, 3, &model, &constraints, rt).unwrap();
            assert!(
                (dynamic.log_partition_function - enumerated.log_partition_function).abs()
                    < 1.0e-10,
                "salt={salt_molar}: dynamic={} exhaustive={}",
                dynamic.log_partition_function,
                enumerated.log_partition_function
            );
            for pair in &dynamic.pair_probabilities {
                let expected = enumerated
                    .pair_probabilities
                    .iter()
                    .find(|candidate| candidate.i == pair.i && candidate.j == pair.j)
                    .map_or(0.0, |candidate| candidate.probability);
                assert!(
                    (pair.probability - expected).abs() < 1.0e-10,
                    "salt={salt_molar}, pair=({}, {}): dynamic={} exhaustive={expected}",
                    pair.i,
                    pair.j,
                    pair.probability
                );
            }
        }
    }

    #[test]
    fn log_domain_remains_finite_for_long_gc_rich_sequences() {
        let sequence = "GGGGCCCC".repeat(35);
        let model = EnergyModel::default();
        let result = partition(&sequence, 37.0, 3, &model).unwrap();
        assert!(result.log_partition_function.is_finite());
        assert!(result.ensemble_free_energy_kcal_mol.is_finite());
        assert!(result
            .pair_probabilities
            .iter()
            .all(|pair| pair.probability.is_finite()));
    }

    #[test]
    fn inside_outside_matches_exhaustive_turner_ensemble() {
        let sequence = "GGGAAACCC";
        // dangles=0 has identical MFE and PF loop weights. ViennaRNA's
        // dangles=2 PF path deliberately smooths dangling-end Boltzmann
        // factors, so it is checked against ViennaRNA's own corpus instead.
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let temperature = 37.0;
        let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
        let structures = enumerate_matchings(sequence.as_bytes(), &model, 3, 0, sequence.len() - 1);
        let mut z = 0.0;
        let mut mass = vec![vec![0.0; sequence.len()]; sequence.len()];
        for pairs in structures {
            let structure = pairs_to_dot_bracket(sequence.len(), &pairs);
            let energy = model.evaluate(sequence, &structure).unwrap().total_kcal_mol;
            let weight = (-energy / rt).exp();
            z += weight;
            for (i, j) in pairs {
                mass[i][j] += weight;
            }
        }

        let result = partition(sequence, temperature, 3, &model).unwrap();
        assert!((result.log_partition_function - z.ln()).abs() < 1.0e-10);
        for pair in result.pair_probabilities {
            assert!((pair.probability - mass[pair.i - 1][pair.j - 1] / z).abs() < 1.0e-10);
        }
    }

    #[test]
    fn no_lonely_pair_partition_matches_vienna_stack_capable_filter() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 2).unwrap();
        let constraints = ConstraintModel::compile(
            sequence.len(),
            &ConstraintConfig {
                no_lonely_pairs: true,
                ..ConstraintConfig::default()
            },
        )
        .unwrap();
        let result = partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        assert!((result.ensemble_free_energy_kcal_mol + 1.618_238_26).abs() < 1.0e-7);
        assert!(!result
            .pair_probabilities
            .iter()
            .any(|pair| pair.i == 1 && pair.j == 7));
    }

    #[test]
    fn no_lp_stack_capability_uses_vienna_raw_pair_matrix_with_no_gu() {
        let sequence = "AUUCAAUAAACAUCUAAAAAGGAUUCUACGACCUACUAUGC";
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let constraints = ConstraintModel::compile(
            sequence.len(),
            &ConstraintConfig {
                no_lonely_pairs: true,
                no_gu: true,
                ..ConstraintConfig::default()
            },
        )
        .unwrap();
        let result = partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        assert!((result.ensemble_free_energy_kcal_mol + 0.404_761_117).abs() < 1.0e-6);
    }
}
