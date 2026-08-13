//! Sliding-window local folding in the RNAplfold sense.
//!
//! Every reported marginal is the arithmetic mean over all complete windows
//! that contain the corresponding pair or unpaired interval.  The window
//! ensembles use the same log-domain Turner partition implementation as the
//! global API.  Joint unpaired probabilities are exact constrained-partition
//! ratios inside each window.

use crate::energy::EnergyModel;
use crate::structure::{normalize_sequence, RnaError};
use crate::{partition_with_constraints, ConstraintConfig, ConstraintModel};
use serde::Serialize;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const NEG_INF: f64 = f64::NEG_INFINITY;

struct LocalInside {
    q_open: Vec<Vec<f64>>,
    q_closed: Vec<Vec<f64>>,
    qb: Vec<Vec<f64>>,
    m1: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
}

struct LocalMarginals {
    pair_sums: Vec<Vec<f64>>,
    oq_open: Vec<Vec<f64>>,
    oq_closed: Vec<Vec<f64>>,
    oqb: Vec<Vec<f64>>,
    om1: Vec<Vec<f64>>,
    om2: Vec<Vec<f64>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalPairProbability {
    pub i: usize,
    pub j: usize,
    pub probability: f64,
    pub containing_windows: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalAccessibility {
    pub from: usize,
    pub to: usize,
    pub length: usize,
    pub probability_unpaired: f64,
    pub opening_energy_kcal_mol: Option<f64>,
    pub containing_windows: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalWindowSummary {
    pub from: usize,
    pub to: usize,
    pub log_partition_function: f64,
    pub ensemble_free_energy_kcal_mol: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalFoldResult {
    pub sequence: String,
    pub window_size: usize,
    pub max_pair_span: usize,
    pub max_unpaired: usize,
    pub pair_probabilities: Vec<LocalPairProbability>,
    pub accessibility: Vec<LocalAccessibility>,
    pub windows: Vec<LocalWindowSummary>,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub method: &'static str,
}

#[allow(clippy::too_many_arguments, unreachable_code, unused_variables)]
pub fn local_fold(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    window_size: usize,
    max_pair_span: usize,
    max_unpaired: usize,
) -> Result<LocalFoldResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    local_fold_with_model(
        sequence,
        min_loop,
        &model,
        window_size,
        max_pair_span,
        max_unpaired,
    )
}

/// Local/windowed analysis with an explicitly selected thermodynamic family.
pub fn local_fold_with_model(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    window_size: usize,
    max_pair_span: usize,
    max_unpaired: usize,
) -> Result<LocalFoldResult, RnaError> {
    let sequence = normalize_sequence(sequence)?;
    let n = sequence.len();
    if window_size == 0 || max_pair_span == 0 {
        return Err(RnaError::InvalidOption(
            "window_size and max_pair_span must be positive".into(),
        ));
    }
    let width = window_size.min(n);
    let span = max_pair_span.min(width);
    let starts: Vec<usize> = if n <= width {
        vec![0]
    } else {
        (0..=n - width).collect()
    };
    let temperature_celsius = model.temperature_celsius();
    let dangles = model.dangles();
    let salt_molar = model.salt_molar();
    if dangles % 2 == 1 {
        return local_fold_odd_exact(
            sequence,
            temperature_celsius,
            min_loop,
            salt_molar,
            width,
            span,
            max_unpaired,
            &starts,
            model,
        );
    }
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let inside = fill_banded_inside(sequence.as_bytes(), min_loop, width, span, rt, model);
    let marginals = accumulate_window_pair_marginals(
        sequence.as_bytes(),
        min_loop,
        width,
        span,
        rt,
        model,
        &inside,
    );
    let pair_sums = &marginals.pair_sums;
    let access_sums = accumulate_joint_unpaired_marginals(
        sequence.as_bytes(),
        min_loop,
        width,
        max_unpaired,
        rt,
        model,
        &inside,
        &marginals,
    );
    let mut window_summaries = Vec::with_capacity(starts.len());

    for &start in &starts {
        let end = start + width;
        let log_z = inside.q_closed[start][end - 1];
        window_summaries.push(LocalWindowSummary {
            from: start + 1,
            to: end,
            log_partition_function: log_z,
            ensemble_free_energy_kcal_mol: -rt * log_z,
        });
    }

    let mut pair_probabilities = Vec::new();
    for (i, row) in pair_sums.iter().enumerate() {
        for (j, &sum) in row.iter().enumerate().skip(i + 1) {
            let count = if j - i < span {
                containing_window_count(n, width, i, j)
            } else {
                0
            };
            if count > 0 {
                let probability = (sum / count as f64).clamp(0.0, 1.0);
                if probability <= 1.0e-12 {
                    continue;
                }
                pair_probabilities.push(LocalPairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                    containing_windows: count,
                });
            }
        }
    }
    let mut accessibility = Vec::new();
    for (from, row) in access_sums.iter().enumerate() {
        for (length, &sum) in row
            .iter()
            .enumerate()
            .take(max_unpaired.min(n - from) + 1)
            .skip(1)
        {
            let count = containing_window_count(n, width, from, from + length - 1);
            if count == 0 {
                continue;
            }
            let probability = (sum / count as f64).clamp(0.0, 1.0);
            accessibility.push(LocalAccessibility {
                from: from + 1,
                to: from + length,
                length,
                probability_unpaired: probability,
                opening_energy_kcal_mol: (probability > 0.0).then_some(-rt * probability.ln()),
                containing_windows: count,
            });
        }
    }

    Ok(LocalFoldResult {
        sequence,
        window_size: width,
        max_pair_span: span,
        max_unpaired,
        pair_probabilities,
        accessibility,
        windows: window_summaries,
        temperature_celsius,
        dangles,
        salt_molar,
        method: "exact banded sliding-window Turner inside/outside ensembles; arithmetic mean over containing windows",
    })
}

#[allow(clippy::too_many_arguments)]
fn local_fold_odd_exact(
    sequence: String,
    temperature_celsius: f64,
    min_loop: usize,
    salt_molar: f64,
    width: usize,
    span: usize,
    max_unpaired: usize,
    starts: &[usize],
    model: &EnergyModel,
) -> Result<LocalFoldResult, RnaError> {
    let n = sequence.len();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut pair_sums = vec![vec![0.0; n]; n];
    let mut access_sums = vec![vec![0.0; max_unpaired + 1]; n];
    let mut windows = Vec::with_capacity(starts.len());
    for &start in starts {
        let end = start + width;
        let subsequence = &sequence[start..end];
        let base_config = ConstraintConfig {
            max_span: Some(span),
            ..ConstraintConfig::default()
        };
        let base_constraints = ConstraintModel::compile(width, &base_config)?;
        let ensemble = partition_with_constraints(
            subsequence,
            temperature_celsius,
            min_loop,
            model,
            &base_constraints,
        )?;
        windows.push(LocalWindowSummary {
            from: start + 1,
            to: end,
            log_partition_function: ensemble.log_partition_function,
            ensemble_free_energy_kcal_mol: ensemble.ensemble_free_energy_kcal_mol,
        });
        for pair in &ensemble.pair_probabilities {
            pair_sums[start + pair.i - 1][start + pair.j - 1] += pair.probability;
        }
        for from in 0..width {
            let global_from = start + from;
            for (length, sum) in access_sums[global_from]
                .iter_mut()
                .enumerate()
                .take(max_unpaired.min(width - from) + 1)
                .skip(1)
            {
                let probability = if length == 1 {
                    ensemble.unpaired_probabilities[from]
                } else {
                    let mut conditioned = base_config.clone();
                    conditioned.force_unpaired.extend(from + 1..=from + length);
                    let conditioned = ConstraintModel::compile(width, &conditioned)?;
                    let result = partition_with_constraints(
                        subsequence,
                        temperature_celsius,
                        min_loop,
                        model,
                        &conditioned,
                    )?;
                    (result.log_partition_function - ensemble.log_partition_function)
                        .exp()
                        .clamp(0.0, 1.0)
                };
                *sum += probability;
            }
        }
    }
    let mut pair_probabilities = Vec::new();
    for (i, row) in pair_sums.iter().enumerate() {
        for (j, &sum) in row.iter().enumerate().skip(i + 1) {
            let count = if j - i < span {
                containing_window_count(n, width, i, j)
            } else {
                0
            };
            if count == 0 {
                continue;
            }
            let probability = (sum / count as f64).clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                pair_probabilities.push(LocalPairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                    containing_windows: count,
                });
            }
        }
    }
    let mut accessibility = Vec::new();
    for (from, row) in access_sums.iter().enumerate() {
        for (length, &sum) in row
            .iter()
            .enumerate()
            .take(max_unpaired.min(n - from) + 1)
            .skip(1)
        {
            let count = containing_window_count(n, width, from, from + length - 1);
            if count == 0 {
                continue;
            }
            let probability = (sum / count as f64).clamp(0.0, 1.0);
            accessibility.push(LocalAccessibility {
                from: from + 1,
                to: from + length,
                length,
                probability_unpaired: probability,
                opening_energy_kcal_mol: (probability > 0.0).then_some(-rt * probability.ln()),
                containing_windows: count,
            });
        }
    }
    Ok(LocalFoldResult {
        sequence,
        window_size: width,
        max_pair_span: span,
        max_unpaired,
        pair_probabilities,
        accessibility,
        windows,
        temperature_celsius,
        dangles: model.dangles(),
        salt_molar,
        method: "exact fixed-structure odd-dangle sliding-window ensembles; arithmetic mean over containing windows",
    })
}

fn fill_banded_inside(
    bases: &[u8],
    min_loop: usize,
    width: usize,
    max_pair_span: usize,
    rt: f64,
    model: &EnergyModel,
) -> LocalInside {
    let n = bases.len();
    let mut q_open = vec![vec![NEG_INF; n]; n];
    let mut q_closed = vec![vec![NEG_INF; n]; n];
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    for i in 0..n {
        q_open[i][i] = 0.0;
        q_closed[i][i] = 0.0;
    }
    for distance in 1..width {
        for i in 0..n - distance {
            let j = i + distance;
            if distance > min_loop && distance < max_pair_span && model.can_pair(bases[i], bases[j])
            {
                let hairpin = model.hairpin_boltzmann_energy(bases, i, j);
                if hairpin.is_finite() {
                    qb[i][j] = -hairpin / rt;
                }
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
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
                        let child = qb[k][l];
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        if energy.is_finite() {
                            log_update(&mut qb[i][j], child - energy / rt);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                    log_update(&mut qb[i][j], m2[i + 1][j - 1] - energy / rt);
                }
            }

            let log_unpaired_ml = -model.multiloop_unpaired_boltzmann() / rt;
            if m1[i][j - 1] != NEG_INF {
                m1[i][j] = m1[i][j - 1] + log_unpaired_ml;
            }
            if m2[i][j - 1] != NEG_INF {
                m2[i][j] = m2[i][j - 1] + log_unpaired_ml;
            }
            if j > min_loop {
                for k in i..=j - min_loop - 1 {
                    if qb[k][j] == NEG_INF {
                        continue;
                    }
                    let branch = qb[k][j] - model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                    log_update(
                        &mut m1[i][j],
                        branch - (k - i) as f64 * model.multiloop_unpaired_boltzmann() / rt,
                    );
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        log_update(&mut m1[i][j], additional);
                        log_update(&mut m2[i][j], additional);
                    }
                }
            }

            q_open[i][j] = q_open[i][j - 1];
            q_closed[i][j] = q_open[i][j - 1];
            if j > min_loop {
                for k in i..=j - min_loop - 1 {
                    if qb[k][j] == NEG_INF {
                        continue;
                    }
                    let left = if k == i { 0.0 } else { q_open[i][k - 1] };
                    let five_prime = (k > i).then(|| bases[k - 1]);
                    let closed_stem =
                        -model.exterior_stem_boltzmann_selected(bases, k, j, five_prime, None) / rt;
                    log_update(&mut q_closed[i][j], left + closed_stem + qb[k][j]);
                    if j + 1 < n {
                        let open_stem = -model.exterior_stem_boltzmann_selected(
                            bases,
                            k,
                            j,
                            five_prime,
                            Some(bases[j + 1]),
                        ) / rt;
                        log_update(&mut q_open[i][j], left + open_stem + qb[k][j]);
                    }
                }
            }
        }
    }
    LocalInside {
        q_open,
        q_closed,
        qb,
        m1,
        m2,
    }
}

#[allow(clippy::too_many_arguments)]
fn accumulate_window_pair_marginals(
    bases: &[u8],
    min_loop: usize,
    width: usize,
    max_pair_span: usize,
    rt: f64,
    model: &EnergyModel,
    inside: &LocalInside,
) -> LocalMarginals {
    let n = bases.len();
    let mut oq_open = vec![vec![NEG_INF; n]; n];
    let mut oq_closed = vec![vec![NEG_INF; n]; n];
    let mut oqb = vec![vec![NEG_INF; n]; n];
    let mut om1 = vec![vec![NEG_INF; n]; n];
    let mut om2 = vec![vec![NEG_INF; n]; n];
    let starts = if n <= width { 0..1 } else { 0..n - width + 1 };
    for start in starts {
        let end = start + width - 1;
        log_update(&mut oq_closed[start][end], -inside.q_closed[start][end]);
    }

    for distance in (1..width).rev() {
        for i in 0..n - distance {
            let j = i + distance;
            let closed_parent = oq_closed[i][j];
            if closed_parent != NEG_INF {
                log_update(&mut oq_open[i][j - 1], closed_parent);
                if j > min_loop {
                    for k in i..=j - min_loop - 1 {
                        if inside.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let left = if k == i { 0.0 } else { inside.q_open[i][k - 1] };
                        let five_prime = (k > i).then(|| bases[k - 1]);
                        let stem = -model
                            .exterior_stem_boltzmann_selected(bases, k, j, five_prime, None)
                            / rt;
                        log_update(&mut oqb[k][j], closed_parent + left + stem);
                        if k > i {
                            log_update(
                                &mut oq_open[i][k - 1],
                                closed_parent + stem + inside.qb[k][j],
                            );
                        }
                    }
                }
            }

            let open_parent = oq_open[i][j];
            if open_parent != NEG_INF {
                log_update(&mut oq_open[i][j - 1], open_parent);
                if j > min_loop && j + 1 < n {
                    for k in i..=j - min_loop - 1 {
                        if inside.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let left = if k == i { 0.0 } else { inside.q_open[i][k - 1] };
                        let five_prime = (k > i).then(|| bases[k - 1]);
                        let stem = -model.exterior_stem_boltzmann_selected(
                            bases,
                            k,
                            j,
                            five_prime,
                            Some(bases[j + 1]),
                        ) / rt;
                        log_update(&mut oqb[k][j], open_parent + left + stem);
                        if k > i {
                            log_update(
                                &mut oq_open[i][k - 1],
                                open_parent + stem + inside.qb[k][j],
                            );
                        }
                    }
                }
            }

            let parent2 = om2[i][j];
            if parent2 != NEG_INF {
                let log_unpaired = -model.multiloop_unpaired_boltzmann() / rt;
                if inside.m2[i][j - 1] != NEG_INF {
                    log_update(&mut om2[i][j - 1], parent2 + log_unpaired);
                }
                if j > min_loop {
                    for k in i + 1..=j - min_loop - 1 {
                        if inside.m1[i][k - 1] == NEG_INF || inside.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        log_update(&mut om1[i][k - 1], parent2 + stem + inside.qb[k][j]);
                        log_update(&mut oqb[k][j], parent2 + inside.m1[i][k - 1] + stem);
                    }
                }
            }

            let parent1 = om1[i][j];
            if parent1 != NEG_INF {
                let log_unpaired = -model.multiloop_unpaired_boltzmann() / rt;
                if inside.m1[i][j - 1] != NEG_INF {
                    log_update(&mut om1[i][j - 1], parent1 + log_unpaired);
                }
                if j > min_loop {
                    for k in i..=j - min_loop - 1 {
                        if inside.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        log_update(
                            &mut oqb[k][j],
                            parent1 + stem
                                - (k - i) as f64 * model.multiloop_unpaired_boltzmann() / rt,
                        );
                        if k > i && inside.m1[i][k - 1] != NEG_INF {
                            log_update(&mut om1[i][k - 1], parent1 + stem + inside.qb[k][j]);
                            log_update(&mut oqb[k][j], parent1 + inside.m1[i][k - 1] + stem);
                        }
                    }
                }
            }

            let pair_parent = oqb[i][j];
            if pair_parent != NEG_INF && inside.qb[i][j] != NEG_INF {
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
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
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        if energy.is_finite() {
                            log_update(&mut oqb[k][l], pair_parent - energy / rt);
                        }
                    }
                }
                if i + 1 < j && inside.m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                    log_update(&mut om2[i + 1][j - 1], pair_parent - energy / rt);
                }
            }
        }
    }

    let mut pair_sums = vec![vec![0.0; n]; n];
    for i in 0..n {
        let maximum_j = (i + max_pair_span.saturating_sub(1)).min(n - 1);
        for j in i + 1..=maximum_j {
            if oqb[i][j] != NEG_INF && inside.qb[i][j] != NEG_INF {
                pair_sums[i][j] = (oqb[i][j] + inside.qb[i][j]).exp();
            }
        }
    }
    LocalMarginals {
        pair_sums,
        oq_open,
        oq_closed,
        oqb,
        om1,
        om2,
    }
}

#[allow(clippy::too_many_arguments)]
fn accumulate_joint_unpaired_marginals(
    bases: &[u8],
    min_loop: usize,
    width: usize,
    max_unpaired: usize,
    rt: f64,
    model: &EnergyModel,
    inside: &LocalInside,
    outside: &LocalMarginals,
) -> Vec<Vec<f64>> {
    let n = bases.len();
    let mut differences = vec![vec![0.0; n + 1]; max_unpaired + 1];
    if max_unpaired == 0 {
        return vec![vec![0.0; 1]; n];
    }

    // Exterior-loop runs. q(i,b) -> q(i,a-1) emits a..b as unpaired;
    // q(i,i) is the one-base all-unpaired base case.
    for i in 0..n {
        let maximum_b = (i + width - 1).min(n - 1);
        let diagonal_outside = log_add(outside.oq_open[i][i], outside.oq_closed[i][i]);
        if diagonal_outside != NEG_INF {
            add_unpaired_range(&mut differences, i, i, diagonal_outside.exp());
        }
        for b in i + 1..=maximum_b {
            let parent = log_add(outside.oq_open[i][b], outside.oq_closed[i][b]);
            if parent == NEG_INF {
                continue;
            }
            for length in 1..=max_unpaired.min(b - i + 1) {
                let a = b + 1 - length;
                let child = if a == i {
                    inside.q_open[i][i]
                } else {
                    inside.q_open[i][a - 1]
                };
                if child != NEG_INF {
                    add_unpaired_interval(&mut differences, a, length, (parent + child).exp());
                }
            }
        }
    }

    let log_unpaired_ml = -model.multiloop_unpaired_boltzmann() / rt;
    for distance in 1..width {
        for i in 0..n - distance {
            let j = i + distance;
            if outside.oqb[i][j] != NEG_INF && inside.qb[i][j] != NEG_INF {
                let hairpin = model.hairpin_boltzmann_energy(bases, i, j);
                if hairpin.is_finite() && j > i + 1 {
                    add_unpaired_range(
                        &mut differences,
                        i + 1,
                        j - 1,
                        (outside.oqb[i][j] - hairpin / rt).exp(),
                    );
                }

                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
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
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        if !energy.is_finite() {
                            continue;
                        }
                        let mass = (outside.oqb[i][j] + inside.qb[k][l] - energy / rt).exp();
                        if i + 1 < k {
                            add_unpaired_range(&mut differences, i + 1, k - 1, mass);
                        }
                        if l + 1 < j {
                            add_unpaired_range(&mut differences, l + 1, j - 1, mass);
                        }
                    }
                }
            }

            // Consecutive unpaired suffixes in multiloop segment states.
            for (parent, table) in [
                (outside.om1[i][j], &inside.m1),
                (outside.om2[i][j], &inside.m2),
            ] {
                if parent == NEG_INF {
                    continue;
                }
                for length in 1..=max_unpaired.min(j - i) {
                    let a = j + 1 - length;
                    let child = table[i][a - 1];
                    if child != NEG_INF {
                        add_unpaired_interval(
                            &mut differences,
                            a,
                            length,
                            (parent + child + length as f64 * log_unpaired_ml).exp(),
                        );
                    }
                }
            }

            // Leading unpaired run before the first branch in m1(i,j).
            let parent = outside.om1[i][j];
            if parent != NEG_INF && j > min_loop {
                for k in i..=j - min_loop - 1 {
                    if inside.qb[k][j] == NEG_INF {
                        continue;
                    }
                    let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                    let mass =
                        (parent + inside.qb[k][j] + stem + (k - i) as f64 * log_unpaired_ml).exp();
                    if i < k {
                        add_unpaired_range(&mut differences, i, k - 1, mass);
                    }
                }
            }
        }
    }

    let mut sums = vec![vec![0.0; max_unpaired + 1]; n];
    for length in 1..=max_unpaired {
        let mut running = 0.0;
        for (start, row) in sums.iter_mut().enumerate() {
            running += differences[length][start];
            if start + length <= n {
                row[length] = running.max(0.0);
            }
        }
    }
    sums
}

fn add_unpaired_range(differences: &mut [Vec<f64>], left: usize, right: usize, mass: f64) {
    if left > right || !mass.is_finite() || mass <= 0.0 {
        return;
    }
    let maximum = (right - left + 1).min(differences.len() - 1);
    for (length, row) in differences.iter_mut().enumerate().take(maximum + 1).skip(1) {
        let last_start = right + 1 - length;
        row[left] += mass;
        row[last_start + 1] -= mass;
    }
}

fn add_unpaired_interval(differences: &mut [Vec<f64>], start: usize, length: usize, mass: f64) {
    if length == 0 || length >= differences.len() || !mass.is_finite() || mass <= 0.0 {
        return;
    }
    differences[length][start] += mass;
    differences[length][start + 1] -= mass;
}

fn containing_window_count(n: usize, width: usize, from: usize, to: usize) -> usize {
    if n <= width {
        return 1;
    }
    let first = to.saturating_add(1).saturating_sub(width);
    let last = from.min(n - width);
    last.saturating_sub(first) + usize::from(first <= last)
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

fn log_update(target: &mut f64, value: f64) {
    *target = log_add(*target, value);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_full_window_is_identical_to_global_partition() {
        let sequence = "GGGAAACCC";
        let local = local_fold(sequence, 37.0, 3, 2, 1.021, 99, 99, 2).unwrap();
        let model = EnergyModel::default();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let global = partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        for expected in global.pair_probabilities {
            let actual = local
                .pair_probabilities
                .iter()
                .find(|entry| entry.i == expected.i && entry.j == expected.j)
                .unwrap();
            assert!((actual.probability - expected.probability).abs() < 1.0e-12);
            assert_eq!(actual.containing_windows, 1);
        }
        for entry in local.accessibility.iter().filter(|entry| entry.length == 1) {
            assert!(
                (entry.probability_unpaired - global.unpaired_probabilities[entry.from - 1]).abs()
                    < 1.0e-12
            );
        }
    }

    #[test]
    fn odd_dangle_full_window_is_identical_to_exact_global_partition() {
        let sequence = "GGGAAACCC";
        let local = local_fold(
            sequence,
            37.0,
            3,
            3,
            1.021,
            sequence.len(),
            sequence.len(),
            3,
        )
        .unwrap();
        let model = EnergyModel::with_dangles(37.0, 3).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let global = partition_with_constraints(sequence, 37.0, 3, &model, &constraints).unwrap();
        assert!(
            (local.windows[0].log_partition_function - global.log_partition_function).abs() < 1e-12
        );
        for pair in global.pair_probabilities {
            let observed = local
                .pair_probabilities
                .iter()
                .find(|entry| entry.i == pair.i && entry.j == pair.j)
                .unwrap();
            assert!((observed.probability - pair.probability).abs() < 1e-12);
        }
    }

    #[test]
    fn denominators_count_only_containing_windows() {
        let result = local_fold("GGGAAACCC", 37.0, 3, 0, 1.021, 6, 5, 1).unwrap();
        let pair = result
            .pair_probabilities
            .iter()
            .find(|entry| entry.i == 3 && entry.j == 7)
            .unwrap();
        assert_eq!(pair.containing_windows, 2);
        let left = result
            .accessibility
            .iter()
            .find(|entry| entry.from == 1 && entry.length == 1)
            .unwrap();
        let middle = result
            .accessibility
            .iter()
            .find(|entry| entry.from == 4 && entry.length == 1)
            .unwrap();
        assert_eq!(left.containing_windows, 1);
        assert_eq!(middle.containing_windows, 4);
    }

    #[test]
    fn banded_inside_outside_and_pu_match_independent_window_partitions() {
        let sequence = "GGGAAACCCGGGAA";
        let width = 10;
        let span = 8;
        let max_unpaired = 3;
        let result = local_fold(sequence, 37.0, 3, 2, 1.021, width, span, max_unpaired).unwrap();
        let model = EnergyModel::default();
        let starts = 0..=sequence.len() - width;
        let baseline_config = ConstraintConfig {
            max_span: Some(span),
            ..ConstraintConfig::default()
        };

        for actual in &result.pair_probabilities {
            let mut sum = 0.0;
            let mut count = 0usize;
            for start in starts.clone() {
                if start < actual.i && actual.j - 1 < start + width {
                    let constraints = ConstraintModel::compile(width, &baseline_config).unwrap();
                    let window = partition_with_constraints(
                        &sequence[start..start + width],
                        37.0,
                        3,
                        &model,
                        &constraints,
                    )
                    .unwrap();
                    sum += window
                        .pair_probabilities
                        .iter()
                        .find(|pair| pair.i == actual.i - start && pair.j == actual.j - start)
                        .map_or(0.0, |pair| pair.probability);
                    count += 1;
                }
            }
            assert_eq!(actual.containing_windows, count);
            assert!(
                (actual.probability - sum / count as f64).abs() < 1.0e-10,
                "pair {}-{}: {} != {}",
                actual.i,
                actual.j,
                actual.probability,
                sum / count as f64
            );
        }

        for actual in &result.accessibility {
            let mut sum = 0.0;
            let mut count = 0usize;
            for start in starts.clone() {
                if start < actual.from && actual.to - 1 < start + width {
                    let subsequence = &sequence[start..start + width];
                    let baseline_constraints =
                        ConstraintModel::compile(width, &baseline_config).unwrap();
                    let baseline = partition_with_constraints(
                        subsequence,
                        37.0,
                        3,
                        &model,
                        &baseline_constraints,
                    )
                    .unwrap();
                    let mut conditioned_config = baseline_config.clone();
                    conditioned_config.force_unpaired =
                        (actual.from - start..=actual.to - start).collect();
                    let conditioned_constraints =
                        ConstraintModel::compile(width, &conditioned_config).unwrap();
                    let conditioned = partition_with_constraints(
                        subsequence,
                        37.0,
                        3,
                        &model,
                        &conditioned_constraints,
                    )
                    .unwrap();
                    sum += (conditioned.log_partition_function - baseline.log_partition_function)
                        .exp();
                    count += 1;
                }
            }
            assert_eq!(actual.containing_windows, count);
            assert!(
                (actual.probability_unpaired - sum / count as f64).abs() < 1.0e-9,
                "{}-{}: {} != {}",
                actual.from,
                actual.to,
                actual.probability_unpaired,
                sum / count as f64
            );
        }
    }
}
