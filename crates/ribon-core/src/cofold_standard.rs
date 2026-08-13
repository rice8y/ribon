//! Exact cut-point secondary-structure grammar for two interacting strands.
//!
//! State index 0 contains no intermolecular pair and state index 1 contains at
//! least one pair crossing `cut`. A crossing pair may close another crossing
//! pair through an internal loop, a multiloop containing a crossing branch, or
//! the inner intermolecular exterior face. The latter contains arbitrary
//! intramolecular exterior structures on both strand suffixes. This is the
//! missing state that distinguishes full cofolding from a bare RNAduplex.

use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{pairs_to_dot_bracket, RnaError};

const INF: f64 = 1.0e100;
const NEG_INF: f64 = f64::NEG_INFINITY;
const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;

#[derive(Clone, Debug)]
pub(crate) struct StandardCofold {
    pub structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub bound_mfe_energy_kcal_mol: Option<f64>,
    pub log_partition_function: f64,
    pub log_unbound_partition_function: f64,
    pub log_bound_partition_function: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
}

#[derive(Clone, Copy, Default)]
enum QChoice {
    #[default]
    Invalid,
    Unpaired,
    Pair {
        left: usize,
        left_state: usize,
        branch_state: usize,
    },
}

#[derive(Clone, Copy, Default)]
enum VChoice {
    #[default]
    Invalid,
    Hairpin,
    IntermolecularExterior,
    Internal(usize, usize),
    Multiloop,
}

#[derive(Clone, Copy, Default)]
enum MChoice {
    #[default]
    Invalid,
    Unpaired,
    First {
        left: usize,
        branch_state: usize,
    },
    Add {
        left: usize,
        prior_state: usize,
        branch_state: usize,
    },
}

#[derive(Clone, Copy, Default)]
enum M2Choice {
    #[default]
    Invalid,
    Unpaired,
    Add {
        left: usize,
        prior_state: usize,
        branch_state: usize,
    },
}

pub(crate) fn cofold_standard(
    sequence_a: &str,
    sequence_b: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<StandardCofold, RnaError> {
    if model.dangles() % 2 == 1 {
        return cofold_odd_exact(sequence_a, sequence_b, temperature_celsius, min_loop, model);
    }
    let cut = sequence_a.len();
    let sequence = format!("{sequence_a}{sequence_b}");
    let bases = sequence.as_bytes();
    let (mfe_structure, mfe_energy, bound_mfe) = mfe(bases, cut, min_loop, model)?;
    let ensemble = partition(bases, cut, temperature_celsius, min_loop, model)?;
    Ok(StandardCofold {
        structure: insert_cut(&mfe_structure, cut),
        mfe_energy_kcal_mol: mfe_energy,
        bound_mfe_energy_kcal_mol: bound_mfe,
        log_partition_function: ensemble.log_z,
        log_unbound_partition_function: ensemble.log_unbound,
        log_bound_partition_function: ensemble.log_bound,
        pair_probabilities: ensemble.pairs,
        unpaired_probabilities: ensemble.unpaired,
    })
}

fn cofold_odd_exact(
    sequence_a: &str,
    sequence_b: &str,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<StandardCofold, RnaError> {
    let cut = sequence_a.len();
    let sequence = format!("{sequence_a}{sequence_b}");
    let bases = sequence.as_bytes();
    let n = bases.len();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut log_z = NEG_INF;
    let mut log_unbound = NEG_INF;
    let mut log_bound = NEG_INF;
    let mut pair_mass = vec![vec![NEG_INF; n]; n];
    let mut best_energy = INF;
    let mut best_bound = INF;
    let mut best_structure = ".".repeat(n);
    enumerate_cofold_structures(bases, cut, min_loop, model, |pairs| {
        let flat = pairs_to_dot_bracket(n, pairs);
        let structure = insert_cut(&flat, cut);
        let energy = model
            .evaluate_cofold(sequence_a, sequence_b, &structure)?
            .total_kcal_mol;
        let log_weight = -energy / rt;
        let bound = pairs.iter().any(|&(i, j)| crosses_cut(i, j, cut));
        log_z = log_add(log_z, log_weight);
        if bound {
            log_bound = log_add(log_bound, log_weight);
            best_bound = best_bound.min(energy);
        } else {
            log_unbound = log_add(log_unbound, log_weight);
        }
        if energy < best_energy {
            best_energy = energy;
            best_structure = flat;
        }
        for &(i, j) in pairs {
            pair_mass[i][j] = log_add(pair_mass[i][j], log_weight);
        }
        Ok(())
    })?;
    let mut probabilities = Vec::new();
    let mut paired = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            if pair_mass[i][j] == NEG_INF {
                continue;
            }
            let probability = (pair_mass[i][j] - log_z).exp().clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                probabilities.push(PairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                });
                paired[i] += probability;
                paired[j] += probability;
            }
        }
    }
    Ok(StandardCofold {
        structure: insert_cut(&best_structure, cut),
        mfe_energy_kcal_mol: best_energy,
        bound_mfe_energy_kcal_mol: (best_bound < INF / 2.0).then_some(best_bound),
        log_partition_function: log_z,
        log_unbound_partition_function: log_unbound,
        log_bound_partition_function: log_bound,
        pair_probabilities: probabilities,
        unpaired_probabilities: paired
            .into_iter()
            .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
            .collect(),
    })
}

fn enumerate_cofold_structures(
    bases: &[u8],
    cut: usize,
    min_loop: usize,
    model: &EnergyModel,
    mut visit: impl FnMut(&[(usize, usize)]) -> Result<(), RnaError>,
) -> Result<(), RnaError> {
    fn recurse(
        bases: &[u8],
        cut: usize,
        min_loop: usize,
        model: &EnergyModel,
        intervals: &mut Vec<(usize, usize)>,
        pairs: &mut Vec<(usize, usize)>,
        visit: &mut impl FnMut(&[(usize, usize)]) -> Result<(), RnaError>,
    ) -> Result<(), RnaError> {
        let Some((i, j)) = intervals.pop() else {
            return visit(pairs);
        };
        if i > j {
            recurse(bases, cut, min_loop, model, intervals, pairs, visit)?;
            intervals.push((i, j));
            return Ok(());
        }
        if i < j {
            intervals.push((i + 1, j));
        }
        recurse(bases, cut, min_loop, model, intervals, pairs, visit)?;
        if i < j {
            intervals.pop();
        }
        for k in i + 1..=j {
            if !pair_allowed(bases, i, k, cut, min_loop, model) {
                continue;
            }
            pairs.push((i, k));
            let mut pushed = 0;
            if k < j {
                intervals.push((k + 1, j));
                pushed += 1;
            }
            if i + 1 < k {
                intervals.push((i + 1, k - 1));
                pushed += 1;
            }
            recurse(bases, cut, min_loop, model, intervals, pairs, visit)?;
            for _ in 0..pushed {
                intervals.pop();
            }
            pairs.pop();
        }
        intervals.push((i, j));
        Ok(())
    }
    let mut intervals = vec![(0, bases.len() - 1)];
    recurse(
        bases,
        cut,
        min_loop,
        model,
        &mut intervals,
        &mut Vec::new(),
        &mut visit,
    )
}

fn crosses_cut(i: usize, j: usize, cut: usize) -> bool {
    i < cut && cut <= j
}

fn pair_allowed(
    bases: &[u8],
    i: usize,
    j: usize,
    cut: usize,
    min_loop: usize,
    model: &EnergyModel,
) -> bool {
    i < j && model.can_pair(bases[i], bases[j]) && (crosses_cut(i, j, cut) || j - i > min_loop)
}

fn interval(table: &[Vec<[f64; 2]>], i: usize, j: usize, state: usize, empty: f64) -> f64 {
    if i > j {
        if state == 0 {
            empty
        } else {
            NEG_INF
        }
    } else {
        table[i][j][state]
    }
}

fn min_interval(table: &[Vec<[f64; 2]>], i: usize, j: usize, state: usize) -> f64 {
    if i > j {
        if state == 0 {
            0.0
        } else {
            INF
        }
    } else {
        table[i][j][state]
    }
}

fn mfe(
    bases: &[u8],
    cut: usize,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<(String, f64, Option<f64>), RnaError> {
    let n = bases.len();
    let mut q = vec![vec![[INF; 2]; n]; n];
    let mut qb = vec![vec![[INF; 2]; n]; n];
    let mut m1 = vec![vec![[INF; 2]; n]; n];
    let mut m2 = vec![vec![[INF; 2]; n]; n];
    let mut qc = vec![vec![[QChoice::Invalid; 2]; n]; n];
    let mut vc = vec![vec![[VChoice::Invalid; 2]; n]; n];
    let mut m1c = vec![vec![[MChoice::Invalid; 2]; n]; n];
    let mut m2c = vec![vec![[M2Choice::Invalid; 2]; n]; n];
    for (i, row) in q.iter_mut().enumerate() {
        row[i][0] = 0.0;
        qc[i][i][0] = QChoice::Unpaired;
    }

    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if pair_allowed(bases, i, j, cut, min_loop, model) {
                let state = usize::from(crosses_cut(i, j, cut));
                if state == 0 {
                    qb[i][j][0] = model.hairpin_energy(bases, i, j);
                    vc[i][j][0] = VChoice::Hairpin;
                } else {
                    let left = min_interval(&q, i + 1, cut - 1, 0);
                    let right = min_interval(&q, cut, j - 1, 0);
                    let end = model.oriented_exterior_stem_energy(
                        bases[j],
                        bases[i],
                        j.checked_sub(1).filter(|&p| p >= cut).map(|p| bases[p]),
                        (i + 1 < cut).then_some(bases[i + 1]),
                    );
                    qb[i][j][1] = left + right + end;
                    vc[i][j][1] = VChoice::IntermolecularExterior;
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
                        if k >= l
                            || usize::from(crosses_cut(k, l, cut)) != state
                            || !pair_allowed(bases, k, l, cut, min_loop, model)
                            || qb[k][l][state] >= INF / 2.0
                        {
                            continue;
                        }
                        let energy = model.internal_energy(bases, i, j, k, l) + qb[k][l][state];
                        if energy < qb[i][j][state] {
                            qb[i][j][state] = energy;
                            vc[i][j][state] = VChoice::Internal(k, l);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1][state] < INF / 2.0 {
                    let energy = model.multiloop_closing()
                        + model.multiloop_closing_stem_energy(bases, i, j)
                        + m2[i + 1][j - 1][state];
                    if energy < qb[i][j][state] {
                        qb[i][j][state] = energy;
                        vc[i][j][state] = VChoice::Multiloop;
                    }
                }
            }

            for state in 0..2 {
                if q[i][j - 1][state] < q[i][j][state] {
                    q[i][j][state] = q[i][j - 1][state];
                    qc[i][j][state] = QChoice::Unpaired;
                }
                if m1[i][j - 1][state] < INF / 2.0 {
                    m1[i][j][state] = m1[i][j - 1][state] + model.multiloop_unpaired();
                    m1c[i][j][state] = MChoice::Unpaired;
                }
                if m2[i][j - 1][state] < INF / 2.0 {
                    m2[i][j][state] = m2[i][j - 1][state] + model.multiloop_unpaired();
                    m2c[i][j][state] = M2Choice::Unpaired;
                }
            }

            for k in i..j {
                if !pair_allowed(bases, k, j, cut, min_loop, model) {
                    continue;
                }
                let branch_state = usize::from(crosses_cut(k, j, cut));
                let branch = qb[k][j][branch_state];
                if branch >= INF / 2.0 {
                    continue;
                }
                let exterior = branch
                    + model.cofold_exterior_stem_energy(bases, k, j, cut)
                    + if branch_state == 1 {
                        model.duplex_initialization_energy()
                    } else {
                        0.0
                    };
                for left_state in 0..2 {
                    let left = min_interval(&q, i, k.saturating_sub(1), left_state);
                    if k == i && left_state == 1 || left >= INF / 2.0 {
                        continue;
                    }
                    let state = left_state | branch_state;
                    let candidate = left + exterior;
                    if candidate < q[i][j][state] {
                        q[i][j][state] = candidate;
                        qc[i][j][state] = QChoice::Pair {
                            left: k,
                            left_state,
                            branch_state,
                        };
                    }
                }

                let ml_branch = branch + model.multiloop_stem_energy(bases, k, j);
                let leading = (k - i) as f64 * model.multiloop_unpaired() + ml_branch;
                if leading < m1[i][j][branch_state] {
                    m1[i][j][branch_state] = leading;
                    m1c[i][j][branch_state] = MChoice::First {
                        left: k,
                        branch_state,
                    };
                }
                if k > i {
                    for prior_state in 0..2 {
                        let prior = m1[i][k - 1][prior_state];
                        if prior >= INF / 2.0 {
                            continue;
                        }
                        let state = prior_state | branch_state;
                        let candidate = prior + ml_branch;
                        if candidate < m1[i][j][state] {
                            m1[i][j][state] = candidate;
                            m1c[i][j][state] = MChoice::Add {
                                left: k,
                                prior_state,
                                branch_state,
                            };
                        }
                        if candidate < m2[i][j][state] {
                            m2[i][j][state] = candidate;
                            m2c[i][j][state] = M2Choice::Add {
                                left: k,
                                prior_state,
                                branch_state,
                            };
                        }
                    }
                }
            }
        }
    }
    let selected_state = usize::from(q[0][n - 1][1] < q[0][n - 1][0]);
    let mut pairs = Vec::new();
    Trace {
        cut,
        q: &qc,
        v: &vc,
        m1: &m1c,
        m2: &m2c,
    }
    .q(0, n - 1, selected_state, &mut pairs);
    pairs.sort_unstable();
    let structure = pairs_to_dot_bracket(n, &pairs);
    let bound = (q[0][n - 1][1] < INF / 2.0).then_some(q[0][n - 1][1]);
    Ok((structure, q[0][n - 1][selected_state], bound))
}

struct Trace<'a> {
    cut: usize,
    q: &'a [Vec<[QChoice; 2]>],
    v: &'a [Vec<[VChoice; 2]>],
    m1: &'a [Vec<[MChoice; 2]>],
    m2: &'a [Vec<[M2Choice; 2]>],
}

impl Trace<'_> {
    fn q(&self, i: usize, j: usize, state: usize, pairs: &mut Vec<(usize, usize)>) {
        if i > j || j >= self.q.len() || i == j {
            return;
        }
        match self.q[i][j][state] {
            QChoice::Unpaired => self.q(i, j - 1, state, pairs),
            QChoice::Pair {
                left,
                left_state,
                branch_state,
            } => {
                if left > i {
                    self.q(i, left - 1, left_state, pairs);
                }
                self.v(left, j, branch_state, pairs);
            }
            QChoice::Invalid => {}
        }
    }

    fn v(&self, i: usize, j: usize, state: usize, pairs: &mut Vec<(usize, usize)>) {
        pairs.push((i, j));
        match self.v[i][j][state] {
            VChoice::Hairpin | VChoice::Invalid => {}
            VChoice::IntermolecularExterior => {
                if i + 1 < self.cut {
                    self.q(i + 1, self.cut - 1, 0, pairs);
                }
                if self.cut < j {
                    self.q(self.cut, j - 1, 0, pairs);
                }
            }
            VChoice::Internal(k, l) => self.v(k, l, state, pairs),
            VChoice::Multiloop => self.m2(i + 1, j - 1, state, pairs),
        }
    }

    fn m1(&self, i: usize, j: usize, state: usize, pairs: &mut Vec<(usize, usize)>) {
        if i > j || j >= self.m1.len() {
            return;
        }
        match self.m1[i][j][state] {
            MChoice::Unpaired => self.m1(i, j - 1, state, pairs),
            MChoice::First { left, branch_state } => self.v(left, j, branch_state, pairs),
            MChoice::Add {
                left,
                prior_state,
                branch_state,
            } => {
                self.m1(i, left - 1, prior_state, pairs);
                self.v(left, j, branch_state, pairs);
            }
            MChoice::Invalid => {}
        }
    }

    fn m2(&self, i: usize, j: usize, state: usize, pairs: &mut Vec<(usize, usize)>) {
        if i > j || j >= self.m2.len() {
            return;
        }
        match self.m2[i][j][state] {
            M2Choice::Unpaired => self.m2(i, j - 1, state, pairs),
            M2Choice::Add {
                left,
                prior_state,
                branch_state,
            } => {
                self.m1(i, left - 1, prior_state, pairs);
                self.v(left, j, branch_state, pairs);
            }
            M2Choice::Invalid => {}
        }
    }
}

struct Ensemble {
    log_z: f64,
    log_unbound: f64,
    log_bound: f64,
    pairs: Vec<PairProbability>,
    unpaired: Vec<f64>,
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

fn update(value: &mut f64, candidate: f64) {
    *value = log_add(*value, candidate);
}

fn partition(
    bases: &[u8],
    cut: usize,
    temperature_celsius: f64,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<Ensemble, RnaError> {
    let n = bases.len();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut q = vec![vec![[NEG_INF; 2]; n]; n];
    let mut qb = vec![vec![[NEG_INF; 2]; n]; n];
    let mut m1 = vec![vec![[NEG_INF; 2]; n]; n];
    let mut m2 = vec![vec![[NEG_INF; 2]; n]; n];
    for (i, row) in q.iter_mut().enumerate() {
        row[i][0] = 0.0;
    }
    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if pair_allowed(bases, i, j, cut, min_loop, model) {
                let state = usize::from(crosses_cut(i, j, cut));
                if state == 0 {
                    qb[i][j][0] = -model.hairpin_boltzmann_energy(bases, i, j) / rt;
                } else {
                    let left = interval(&q, i + 1, cut - 1, 0, 0.0);
                    let right = interval(&q, cut, j - 1, 0, 0.0);
                    let end = model.oriented_exterior_stem_boltzmann_energy(
                        bases[j],
                        bases[i],
                        j.checked_sub(1).filter(|&p| p >= cut).map(|p| bases[p]),
                        (i + 1 < cut).then_some(bases[i + 1]),
                    );
                    qb[i][j][1] = left + right - end / rt;
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
                        if k >= l
                            || usize::from(crosses_cut(k, l, cut)) != state
                            || !pair_allowed(bases, k, l, cut, min_loop, model)
                            || qb[k][l][state] == NEG_INF
                        {
                            continue;
                        }
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        let candidate = qb[k][l][state] - energy / rt;
                        update(&mut qb[i][j][state], candidate);
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1][state] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                    update(&mut qb[i][j][state], m2[i + 1][j - 1][state] - energy / rt);
                }
            }
            let ml_unpaired = -model.multiloop_unpaired_boltzmann() / rt;
            for state in 0..2 {
                if q[i][j - 1][state] != NEG_INF {
                    q[i][j][state] = q[i][j - 1][state];
                }
                if m1[i][j - 1][state] != NEG_INF {
                    m1[i][j][state] = m1[i][j - 1][state] + ml_unpaired;
                }
                if m2[i][j - 1][state] != NEG_INF {
                    m2[i][j][state] = m2[i][j - 1][state] + ml_unpaired;
                }
            }
            for k in i..j {
                if !pair_allowed(bases, k, j, cut, min_loop, model) {
                    continue;
                }
                let bs = usize::from(crosses_cut(k, j, cut));
                if qb[k][j][bs] == NEG_INF {
                    continue;
                }
                let exterior = qb[k][j][bs]
                    - (model.cofold_exterior_stem_boltzmann_energy(bases, k, j, cut)
                        + if bs == 1 {
                            model.duplex_initialization_boltzmann_energy()
                        } else {
                            0.0
                        })
                        / rt;
                for ls in 0..2 {
                    let left = interval(&q, i, k.saturating_sub(1), ls, 0.0);
                    if k == i && ls == 1 || left == NEG_INF {
                        continue;
                    }
                    update(&mut q[i][j][ls | bs], left + exterior);
                }
                let branch = qb[k][j][bs] - model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                update(
                    &mut m1[i][j][bs],
                    -((k - i) as f64 * model.multiloop_unpaired_boltzmann()) / rt + branch,
                );
                if k > i {
                    for ps in 0..2 {
                        if m1[i][k - 1][ps] == NEG_INF {
                            continue;
                        }
                        let candidate = m1[i][k - 1][ps] + branch;
                        update(&mut m1[i][j][ps | bs], candidate);
                        update(&mut m2[i][j][ps | bs], candidate);
                    }
                }
            }
        }
    }
    let log_z = log_add(q[0][n - 1][0], q[0][n - 1][1]);
    let log_bound = q[0][n - 1][1];

    let mut oq = vec![vec![[NEG_INF; 2]; n]; n];
    let mut oqb = vec![vec![[NEG_INF; 2]; n]; n];
    let mut om1 = vec![vec![[NEG_INF; 2]; n]; n];
    let mut om2 = vec![vec![[NEG_INF; 2]; n]; n];
    oq[0][n - 1] = [0.0, 0.0];
    for span in (1..n).rev() {
        for i in 0..n - span {
            let j = i + span;
            let ml_unpaired = -model.multiloop_unpaired_boltzmann() / rt;
            for state in 0..2 {
                let parent = oq[i][j][state];
                if parent != NEG_INF {
                    if q[i][j - 1][state] != NEG_INF {
                        update(&mut oq[i][j - 1][state], parent);
                    }
                    for k in i..j {
                        if !pair_allowed(bases, k, j, cut, min_loop, model) {
                            continue;
                        }
                        let bs = usize::from(crosses_cut(k, j, cut));
                        if qb[k][j][bs] == NEG_INF {
                            continue;
                        }
                        let stem = -(model.cofold_exterior_stem_boltzmann_energy(bases, k, j, cut)
                            + if bs == 1 {
                                model.duplex_initialization_boltzmann_energy()
                            } else {
                                0.0
                            })
                            / rt;
                        for ls in [0usize, 1] {
                            if ls | bs != state {
                                continue;
                            }
                            let left = interval(&q, i, k.saturating_sub(1), ls, 0.0);
                            if k == i && ls == 1 || left == NEG_INF {
                                continue;
                            }
                            update(&mut oqb[k][j][bs], parent + left + stem);
                            if k > i {
                                update(&mut oq[i][k - 1][ls], parent + qb[k][j][bs] + stem);
                            }
                        }
                    }
                }

                let parent = om2[i][j][state];
                if parent != NEG_INF {
                    if m2[i][j - 1][state] != NEG_INF {
                        update(&mut om2[i][j - 1][state], parent + ml_unpaired);
                    }
                    for k in i + 1..j {
                        if !pair_allowed(bases, k, j, cut, min_loop, model) {
                            continue;
                        }
                        let bs = usize::from(crosses_cut(k, j, cut));
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        for ps in 0..2 {
                            if ps | bs != state
                                || m1[i][k - 1][ps] == NEG_INF
                                || qb[k][j][bs] == NEG_INF
                            {
                                continue;
                            }
                            update(&mut om1[i][k - 1][ps], parent + qb[k][j][bs] + stem);
                            update(&mut oqb[k][j][bs], parent + m1[i][k - 1][ps] + stem);
                        }
                    }
                }

                let parent = om1[i][j][state];
                if parent != NEG_INF {
                    if m1[i][j - 1][state] != NEG_INF {
                        update(&mut om1[i][j - 1][state], parent + ml_unpaired);
                    }
                    for k in i..j {
                        if !pair_allowed(bases, k, j, cut, min_loop, model) {
                            continue;
                        }
                        let bs = usize::from(crosses_cut(k, j, cut));
                        if qb[k][j][bs] == NEG_INF {
                            continue;
                        }
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        if bs == state {
                            let leading =
                                -((k - i) as f64 * model.multiloop_unpaired_boltzmann()) / rt;
                            update(&mut oqb[k][j][bs], parent + leading + stem);
                        }
                        if k > i {
                            for ps in 0..2 {
                                if ps | bs != state || m1[i][k - 1][ps] == NEG_INF {
                                    continue;
                                }
                                update(&mut om1[i][k - 1][ps], parent + qb[k][j][bs] + stem);
                                update(&mut oqb[k][j][bs], parent + m1[i][k - 1][ps] + stem);
                            }
                        }
                    }
                }

                let parent = oqb[i][j][state];
                if parent == NEG_INF || qb[i][j][state] == NEG_INF {
                    continue;
                }
                if state == 1 {
                    let left = interval(&q, i + 1, cut - 1, 0, 0.0);
                    let right = interval(&q, cut, j - 1, 0, 0.0);
                    let end = model.oriented_exterior_stem_boltzmann_energy(
                        bases[j],
                        bases[i],
                        j.checked_sub(1).filter(|&p| p >= cut).map(|p| bases[p]),
                        (i + 1 < cut).then_some(bases[i + 1]),
                    );
                    let root = parent - end / rt;
                    if i + 1 < cut {
                        update(&mut oq[i + 1][cut - 1][0], root + right);
                    }
                    if cut < j {
                        update(&mut oq[cut][j - 1][0], root + left);
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
                        if k >= l
                            || usize::from(crosses_cut(k, l, cut)) != state
                            || !pair_allowed(bases, k, l, cut, min_loop, model)
                            || qb[k][l][state] == NEG_INF
                        {
                            continue;
                        }
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        update(&mut oqb[k][l][state], parent - energy / rt);
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1][state] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                    update(&mut om2[i + 1][j - 1][state], parent - energy / rt);
                }
            }
        }
    }

    let mut pairs = Vec::new();
    let mut paired = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            let state = usize::from(crosses_cut(i, j, cut));
            if oqb[i][j][state] == NEG_INF || qb[i][j][state] == NEG_INF {
                continue;
            }
            let probability = (oqb[i][j][state] + qb[i][j][state] - log_z)
                .exp()
                .clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                pairs.push(PairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                });
                paired[i] += probability;
                paired[j] += probability;
            }
        }
    }
    let unpaired = paired
        .into_iter()
        .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
        .collect();
    Ok(Ensemble {
        log_z,
        log_unbound: q[0][n - 1][0],
        log_bound,
        pairs,
        unpaired,
    })
}

fn insert_cut(structure: &str, cut: usize) -> String {
    format!("{}&{}", &structure[..cut], &structure[cut..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_grammar_reduces_to_duplex_for_three_complementary_bases() {
        let model = EnergyModel::with_dangles(37.0, 0).unwrap();
        let result = cofold_standard("GGG", "CCC", 37.0, 3, &model).unwrap();
        let duplex = crate::duplex::duplex("GGG", "CCC", 37.0, 1.021).unwrap();
        assert_eq!(result.structure, "(((&)))");
        assert!(
            (result.bound_mfe_energy_kcal_mol.unwrap() - duplex.mfe_energy_kcal_mol).abs()
                < 1.0e-10
        );
        let mut exhaustive = NEG_INF;
        let rt = GAS_CONSTANT_KCAL * (37.0 + 273.15);
        enumerate_cofold_structures(b"GGGCCC", 3, 3, &model, |pairs| {
            if pairs.iter().any(|&(i, j)| crosses_cut(i, j, 3)) {
                let flat = pairs_to_dot_bracket(6, pairs);
                let structure = insert_cut(&flat, 3);
                let energy = model
                    .evaluate_cofold("GGG", "CCC", &structure)?
                    .total_kcal_mol;
                exhaustive = log_add(exhaustive, -energy / rt);
            }
            Ok(())
        })
        .unwrap();
        assert!((result.log_bound_partition_function - exhaustive).abs() < 1.0e-10);
    }

    #[test]
    fn bound_state_contains_intramolecular_branches_and_normalized_marginals() {
        let model = EnergyModel::with_dangles(37.0, 2).unwrap();
        let result = cofold_standard("GGGAAACCCG", "CGGGAAACCC", 37.0, 3, &model).unwrap();
        assert!(result.log_bound_partition_function.is_finite());
        for i in 0..result.unpaired_probabilities.len() {
            let paired = result
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == i + 1 || pair.j == i + 1)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!((paired + result.unpaired_probabilities[i] - 1.0).abs() < 1.0e-8);
        }
    }

    #[test]
    fn unbound_state_is_exactly_the_product_of_the_two_monomer_ensembles() {
        let a = "GGGAAACCC";
        let b = "GCGAAACGC";
        for dangles in [0, 2] {
            let model = EnergyModel::with_dangles(37.0, dangles).unwrap();
            let result = cofold_standard(a, b, 37.0, 3, &model).unwrap();
            let left = crate::partition::partition(a, 37.0, 3, &model).unwrap();
            let right = crate::partition::partition(b, 37.0, 3, &model).unwrap();
            assert!(
                (result.log_unbound_partition_function
                    - left.log_partition_function
                    - right.log_partition_function)
                    .abs()
                    < 1.0e-10,
                "dangles={dangles}",
            );
        }
    }

    #[test]
    fn odd_dangle_full_cofold_is_exhaustive_and_self_consistent() {
        for dangles in [1, 3] {
            let model = EnergyModel::with_dangles(37.0, dangles).unwrap();
            let result = cofold_standard("GGG", "CCC", 37.0, 3, &model).unwrap();
            let evaluated = model
                .evaluate_cofold("GGG", "CCC", &result.structure)
                .unwrap();
            assert!(
                (evaluated.total_kcal_mol - result.mfe_energy_kcal_mol).abs() < 1.0e-10,
                "dangles={dangles}",
            );
            assert!(result.log_bound_partition_function.is_finite());
            for (index, unpaired) in result.unpaired_probabilities.iter().enumerate() {
                let paired = result
                    .pair_probabilities
                    .iter()
                    .filter(|pair| pair.i == index + 1 || pair.j == index + 1)
                    .map(|pair| pair.probability)
                    .sum::<f64>();
                assert!((paired + unpaired - 1.0).abs() < 1.0e-10);
            }
        }
    }
}
