//! Exact stochastic backtracking through the McCaskill inside grammar.

use crate::constraints::ConstraintModel;
use crate::energy::EnergyModel;
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::Serialize;
use std::collections::HashSet;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Debug, Serialize)]
pub struct BoltzmannSample {
    pub structure: String,
    pub effective_energy_kcal_mol: f64,
    pub probability: f64,
    pub log_probability: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct SamplingResult {
    pub sequence: String,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub seed: u64,
    pub requested: usize,
    pub returned: usize,
    pub unique: bool,
    pub log_partition_function: f64,
    pub samples: Vec<BoltzmannSample>,
}

struct Inside {
    q: Vec<Vec<f64>>,
    qb: Vec<Vec<f64>>,
    m1: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
    rt: f64,
}

#[derive(Clone, Copy)]
enum QChoice {
    Unpaired,
    Pair(usize),
}

#[derive(Clone, Copy)]
enum QbChoice {
    Hairpin,
    Internal(usize, usize),
    Multiloop,
}

#[derive(Clone, Copy)]
enum MChoice {
    Unpaired,
    LeadingBranch(usize),
    SplitBranch(usize),
}

#[allow(clippy::too_many_arguments)]
pub fn sample_boltzmann(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    count: usize,
    seed: u64,
    unique: bool,
) -> Result<SamplingResult, RnaError> {
    if count == 0 {
        return Err(RnaError::InvalidOption(
            "sample count must be positive".into(),
        ));
    }
    if !model.supports_partition() {
        return Err(RnaError::InvalidOption(
            "Boltzmann sampling requires dangles=0, 1, 2, or 3".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    if model.dangles() % 2 == 1 {
        return sample_boltzmann_odd_exact(
            &sequence,
            temperature_celsius,
            min_loop,
            model,
            constraints,
            count,
            seed,
            unique,
        );
    }
    let bases = sequence.as_bytes();
    let inside = compute_inside(bases, temperature_celsius, min_loop, model, constraints)?;
    let log_z = inside.q[0][bases.len() - 1];
    let mut rng = SplitMix64::new(seed);
    let mut samples = Vec::with_capacity(count);
    let mut seen = HashSet::new();
    // `count` always means the exact number of independent Boltzmann draws.
    // With `unique=true`, duplicate outcomes are removed after drawing; this
    // avoids an undocumented rejection-sampling attempt cap and cannot bias
    // the distribution of the retained structures.
    for _ in 0..count {
        let mut pairs = Vec::new();
        let mut log_probability = 0.0;
        sample_q(
            0,
            bases.len() - 1,
            bases,
            min_loop,
            model,
            constraints,
            &inside,
            &mut rng,
            &mut pairs,
            &mut log_probability,
        );
        pairs.sort_unstable();
        let structure = pairs_to_dot_bracket(bases.len(), &pairs);
        if unique && !seen.insert(structure.clone()) {
            continue;
        }
        let log_weight = log_z + log_probability;
        samples.push(BoltzmannSample {
            structure,
            effective_energy_kcal_mol: -inside.rt * log_weight,
            probability: log_probability.exp(),
            log_probability,
        });
    }
    Ok(SamplingResult {
        sequence,
        temperature_celsius,
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        seed,
        requested: count,
        returned: samples.len(),
        unique,
        log_partition_function: log_z,
        samples,
    })
}

#[allow(clippy::too_many_arguments)]
fn sample_boltzmann_odd_exact(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    count: usize,
    seed: u64,
    unique: bool,
) -> Result<SamplingResult, RnaError> {
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let n = sequence.len();
    let mut states = Vec::<(String, f64, f64)>::new();
    let mut log_z = NEG_INF;
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
            states.push((structure, energy, log_weight));
            Ok(())
        },
    )?;
    if states.is_empty() || log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    }
    let mut cumulative = Vec::with_capacity(states.len());
    let mut total = 0.0;
    for (_, _, log_weight) in &states {
        total += (*log_weight - log_z).exp();
        cumulative.push(total);
    }
    let mut rng = SplitMix64::new(seed);
    let mut samples = Vec::with_capacity(count);
    let mut seen = HashSet::new();
    for _ in 0..count {
        let draw = rng.next_f64() * total;
        let index = cumulative.partition_point(|&boundary| boundary < draw);
        let (structure, energy, log_weight) = &states[index.min(states.len() - 1)];
        if unique && !seen.insert(structure.clone()) {
            continue;
        }
        let log_probability = *log_weight - log_z;
        samples.push(BoltzmannSample {
            structure: structure.clone(),
            effective_energy_kcal_mol: *energy,
            probability: log_probability.exp(),
            log_probability,
        });
    }
    Ok(SamplingResult {
        sequence: sequence.into(),
        temperature_celsius,
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        seed,
        requested: count,
        returned: samples.len(),
        unique,
        log_partition_function: log_z,
        samples,
    })
}

fn compute_inside(
    bases: &[u8],
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<Inside, RnaError> {
    let n = bases.len();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    if n == 0 || !rt.is_finite() || rt <= 0.0 {
        return Err(RnaError::InvalidOption("invalid sampling input".into()));
    }
    let mut q = vec![vec![NEG_INF; n]; n];
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    for (i, row) in q.iter_mut().enumerate() {
        if constraints.allows_unpaired(i) {
            row[i] = -constraints.unpaired_energy(i) / rt;
        }
    }
    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if span > min_loop
                && constraints.allows_pair_for_partition(bases, i, j, min_loop, model)
            {
                let pair_soft = constraints.pair_energy(i, j);
                let hairpin = model.hairpin_boltzmann_energy(bases, i, j);
                if let Some(unpaired) = constraints.unpaired_range_energy(i + 1, j - 1) {
                    if hairpin.is_finite() {
                        qb[i][j] = -(hairpin + pair_soft + unpaired) / rt;
                    }
                }
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    for right in 0..=internal_limit - left {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        if k >= l || l - k <= min_loop || qb[k][l] == NEG_INF {
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
                            qb[i][j] = log_add(qb[i][j], qb[k][l] - energy / rt);
                        }
                    }
                }
                if span >= 2 && m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
                        + pair_soft;
                    qb[i][j] = log_add(qb[i][j], m2[i + 1][j - 1] - energy / rt);
                }
            }
            let unpaired_ml = if constraints.allows_unpaired(j) {
                -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
            } else {
                NEG_INF
            };
            if m1[i][j - 1] != NEG_INF && unpaired_ml != NEG_INF {
                m1[i][j] = m1[i][j - 1] + unpaired_ml;
            }
            if m2[i][j - 1] != NEG_INF && unpaired_ml != NEG_INF {
                m2[i][j] = m2[i][j - 1] + unpaired_ml;
            }
            if j > min_loop {
                for k in i..=j - min_loop - 1 {
                    if qb[k][j] == NEG_INF {
                        continue;
                    }
                    let branch = qb[k][j] - model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                    let leading = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    };
                    if let Some(soft) = leading {
                        let unpaired =
                            -((k - i) as f64 * model.multiloop_unpaired_boltzmann() + soft) / rt;
                        m1[i][j] = log_add(m1[i][j], unpaired + branch);
                    }
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let value = m1[i][k - 1] + branch;
                        m1[i][j] = log_add(m1[i][j], value);
                        m2[i][j] = log_add(m2[i][j], value);
                    }
                }
            }
            if constraints.allows_unpaired(j) && q[i][j - 1] != NEG_INF {
                q[i][j] = q[i][j - 1] - constraints.unpaired_energy(j) / rt;
            }
            if j > min_loop {
                for k in i..=j - min_loop - 1 {
                    if qb[k][j] == NEG_INF {
                        continue;
                    }
                    let left = if k == i { 0.0 } else { q[i][k - 1] };
                    let stem = -model.exterior_stem_boltzmann_energy(bases, k, j) / rt;
                    q[i][j] = log_add(q[i][j], left + stem + qb[k][j]);
                }
            }
        }
    }
    if q[0][n - 1] == NEG_INF {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    }
    Ok(Inside { q, qb, m1, m2, rt })
}

#[allow(clippy::too_many_arguments)]
fn sample_q(
    i: usize,
    j: usize,
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    inside: &Inside,
    rng: &mut SplitMix64,
    pairs: &mut Vec<(usize, usize)>,
    log_probability: &mut f64,
) {
    if i >= j {
        return;
    }
    let mut candidates = Vec::new();
    if constraints.allows_unpaired(j) && inside.q[i][j - 1] != NEG_INF {
        candidates.push((
            inside.q[i][j - 1] - constraints.unpaired_energy(j) / inside.rt,
            QChoice::Unpaired,
        ));
    }
    if j > min_loop {
        for k in i..=j - min_loop - 1 {
            if inside.qb[k][j] == NEG_INF {
                continue;
            }
            let left = if k == i { 0.0 } else { inside.q[i][k - 1] };
            let stem = -model.exterior_stem_boltzmann_energy(bases, k, j) / inside.rt;
            candidates.push((left + stem + inside.qb[k][j], QChoice::Pair(k)));
        }
    }
    let (weight, choice) = choose(&candidates, rng);
    *log_probability += weight - inside.q[i][j];
    match choice {
        QChoice::Unpaired => sample_q(
            i,
            j - 1,
            bases,
            min_loop,
            model,
            constraints,
            inside,
            rng,
            pairs,
            log_probability,
        ),
        QChoice::Pair(k) => {
            if k > i {
                sample_q(
                    i,
                    k - 1,
                    bases,
                    min_loop,
                    model,
                    constraints,
                    inside,
                    rng,
                    pairs,
                    log_probability,
                );
            }
            sample_qb(
                k,
                j,
                bases,
                min_loop,
                model,
                constraints,
                inside,
                rng,
                pairs,
                log_probability,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn sample_qb(
    i: usize,
    j: usize,
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    inside: &Inside,
    rng: &mut SplitMix64,
    pairs: &mut Vec<(usize, usize)>,
    log_probability: &mut f64,
) {
    pairs.push((i, j));
    let pair_soft = constraints.pair_energy(i, j);
    let mut candidates = Vec::new();
    if let Some(unpaired) = constraints.unpaired_range_energy(i + 1, j - 1) {
        let hairpin = model.hairpin_boltzmann_energy(bases, i, j);
        if hairpin.is_finite() {
            candidates.push((
                -(hairpin + pair_soft + unpaired) / inside.rt,
                QbChoice::Hairpin,
            ));
        }
    }
    let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
    for left in 0..=internal_limit {
        let k = i + 1 + left;
        if k >= j {
            break;
        }
        for right in 0..=internal_limit - left {
            let Some(l) = j.checked_sub(right + 1) else {
                continue;
            };
            if k >= l || l - k <= min_loop || inside.qb[k][l] == NEG_INF {
                continue;
            }
            let Some(left_soft) = constraints.unpaired_range_energy(i + 1, k - 1) else {
                continue;
            };
            let Some(right_soft) = constraints.unpaired_range_energy(l + 1, j - 1) else {
                continue;
            };
            let energy = model.internal_boltzmann_energy(bases, i, j, k, l)
                + pair_soft
                + left_soft
                + right_soft
                + constraints.stack_energy(i, j, k, l);
            if energy.is_finite() {
                candidates.push((
                    inside.qb[k][l] - energy / inside.rt,
                    QbChoice::Internal(k, l),
                ));
            }
        }
    }
    if i + 1 < j && inside.m2[i + 1][j - 1] != NEG_INF {
        let energy = model.multiloop_closing_boltzmann()
            + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
            + pair_soft;
        candidates.push((
            inside.m2[i + 1][j - 1] - energy / inside.rt,
            QbChoice::Multiloop,
        ));
    }
    let (weight, choice) = choose(&candidates, rng);
    *log_probability += weight - inside.qb[i][j];
    match choice {
        QbChoice::Hairpin => {}
        QbChoice::Internal(k, l) => sample_qb(
            k,
            l,
            bases,
            min_loop,
            model,
            constraints,
            inside,
            rng,
            pairs,
            log_probability,
        ),
        QbChoice::Multiloop => sample_m(
            i + 1,
            j - 1,
            true,
            bases,
            min_loop,
            model,
            constraints,
            inside,
            rng,
            pairs,
            log_probability,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn sample_m(
    i: usize,
    j: usize,
    two_or_more: bool,
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    inside: &Inside,
    rng: &mut SplitMix64,
    pairs: &mut Vec<(usize, usize)>,
    log_probability: &mut f64,
) {
    let table = if two_or_more { &inside.m2 } else { &inside.m1 };
    let parent = table[i][j];
    let mut candidates = Vec::new();
    let unpaired = if constraints.allows_unpaired(j) {
        -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / inside.rt
    } else {
        NEG_INF
    };
    if j > i && table[i][j - 1] != NEG_INF && unpaired != NEG_INF {
        candidates.push((table[i][j - 1] + unpaired, MChoice::Unpaired));
    }
    if j > min_loop {
        for k in i..=j - min_loop - 1 {
            if inside.qb[k][j] == NEG_INF {
                continue;
            }
            let branch =
                inside.qb[k][j] - model.multiloop_stem_boltzmann_energy(bases, k, j) / inside.rt;
            if !two_or_more {
                let leading = if k == i {
                    Some(0.0)
                } else {
                    constraints.unpaired_range_energy(i, k - 1)
                };
                if let Some(soft) = leading {
                    let value = branch
                        - ((k - i) as f64 * model.multiloop_unpaired_boltzmann() + soft)
                            / inside.rt;
                    candidates.push((value, MChoice::LeadingBranch(k)));
                }
            }
            if k > i && inside.m1[i][k - 1] != NEG_INF {
                candidates.push((inside.m1[i][k - 1] + branch, MChoice::SplitBranch(k)));
            }
        }
    }
    let (weight, choice) = choose(&candidates, rng);
    *log_probability += weight - parent;
    match choice {
        MChoice::Unpaired => sample_m(
            i,
            j - 1,
            two_or_more,
            bases,
            min_loop,
            model,
            constraints,
            inside,
            rng,
            pairs,
            log_probability,
        ),
        MChoice::LeadingBranch(k) => sample_qb(
            k,
            j,
            bases,
            min_loop,
            model,
            constraints,
            inside,
            rng,
            pairs,
            log_probability,
        ),
        MChoice::SplitBranch(k) => {
            sample_m(
                i,
                k - 1,
                false,
                bases,
                min_loop,
                model,
                constraints,
                inside,
                rng,
                pairs,
                log_probability,
            );
            sample_qb(
                k,
                j,
                bases,
                min_loop,
                model,
                constraints,
                inside,
                rng,
                pairs,
                log_probability,
            );
        }
    }
}

fn choose<T: Copy>(candidates: &[(f64, T)], rng: &mut SplitMix64) -> (f64, T) {
    debug_assert!(!candidates.is_empty());
    let maximum = candidates
        .iter()
        .map(|entry| entry.0)
        .fold(NEG_INF, f64::max);
    let total = candidates
        .iter()
        .map(|entry| (entry.0 - maximum).exp())
        .sum::<f64>();
    let mut target = rng.next_f64() * total;
    for &(weight, choice) in candidates {
        target -= (weight - maximum).exp();
        if target <= 0.0 {
            return (weight, choice);
        }
    }
    *candidates.last().expect("non-empty candidates")
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

struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D049BB133111EB);
        value ^ (value >> 31)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition;

    #[test]
    fn sampling_is_reproducible_and_probabilities_are_valid() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::default();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let first =
            sample_boltzmann(sequence, 37.0, 3, &model, &constraints, 64, 42, false).unwrap();
        let second =
            sample_boltzmann(sequence, 37.0, 3, &model, &constraints, 64, 42, false).unwrap();
        assert_eq!(first.samples[0].structure, second.samples[0].structure);
        assert!(first
            .samples
            .iter()
            .all(|sample| sample.log_probability <= 1.0e-12));
    }

    #[test]
    fn odd_dangle_sampling_uses_the_exact_partition() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 3).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let partition =
            partition::partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        let sampled =
            sample_boltzmann(sequence, 37.0, 3, &model, &constraints, 32, 7, false).unwrap();
        assert!((sampled.log_partition_function - partition.log_partition_function).abs() < 1e-12);
        assert_eq!(sampled.returned, 32);
    }

    #[test]
    fn sampled_pair_frequencies_track_inside_outside_probabilities() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let partition = partition(sequence, 37.0, 3, &model).unwrap();
        let sampled =
            sample_boltzmann(sequence, 37.0, 3, &model, &constraints, 10_000, 7, false).unwrap();
        for pair in partition
            .pair_probabilities
            .iter()
            .filter(|pair| pair.probability > 0.02)
        {
            let frequency = sampled
                .samples
                .iter()
                .filter(|sample| {
                    let parsed = crate::parse_structure(sequence, &sample.structure).unwrap();
                    parsed.partner[pair.i - 1] == Some(pair.j - 1)
                })
                .count() as f64
                / sampled.samples.len() as f64;
            assert!((frequency - pair.probability).abs() < 0.02);
        }
        let rt = GAS_CONSTANT_KCAL * (37.0 + 273.15);
        for sample in sampled.samples.iter().take(100) {
            let energy = model
                .evaluate(sequence, &sample.structure)
                .unwrap()
                .total_kcal_mol;
            let expected = -energy / rt - partition.log_partition_function;
            assert!((sample.log_probability - expected).abs() < 1.0e-10);
        }
    }
}
