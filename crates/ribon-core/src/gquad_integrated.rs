//! Integrated G-quadruplex/secondary-structure dynamic programming.
//!
//! The grammar follows the public recurrences described by Lorenz et al.:
//! every interval-spanning G-quadruplex is a terminal substructure that may
//! occur in the exterior loop or as a multiloop branch.  A single G-quadruplex
//! enclosed by a canonical pair is treated as a generalized multiloop.  This
//! makes multiple and nested occurrences possible without enumerating whole
//! structures.  MFE and partition-function calculations share the grammar.

use crate::energy::EnergyModel;
use crate::extended::GQuadruplex;
use crate::partition::PairProbability;
use crate::structure::RnaError;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const INF: f64 = 1.0e100;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Debug)]
pub(crate) struct IntegratedGquadResult {
    pub structure: String,
    pub mfe_energy: f64,
    pub log_partition_function: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub gquad_position_probabilities: Vec<f64>,
    pub interval_probabilities: Vec<(usize, usize, f64)>,
    pub candidate_probabilities: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
enum ExteriorChoice {
    Unpaired,
    Pair(usize),
    Gquad(usize),
    Invalid,
}

#[derive(Clone, Copy, Debug)]
enum PairChoice {
    Hairpin,
    Internal(usize, usize),
    Multiloop,
    SingleGquad,
    Invalid,
}

#[derive(Clone, Copy, Debug)]
enum SegmentChoice {
    Unpaired,
    FirstPair(usize),
    FirstGquad(usize),
    AddPair(usize),
    AddGquad(usize),
    Invalid,
}

#[derive(Clone, Copy, Debug)]
enum SingleGquadChoice {
    Unpaired,
    First(usize),
    Invalid,
}

fn interval(table: &[Vec<f64>], i: usize, j: usize) -> f64 {
    if i > j {
        0.0
    } else {
        table[i][j]
    }
}

fn log_interval(table: &[Vec<f64>], i: usize, j: usize) -> f64 {
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
    high + (a.min(b) - high).exp().ln_1p()
}

fn update(target: &mut f64, value: f64) {
    *target = log_add(*target, value);
}

pub(crate) fn predict_with_model(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    candidates: &[GQuadruplex],
) -> Result<IntegratedGquadResult, RnaError> {
    let bases = sequence.as_bytes();
    let n = bases.len();
    let temperature_celsius = model.temperature_celsius();
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);

    let mut g_mfe = vec![vec![INF; n]; n];
    let mut g_choice = vec![vec![None; n]; n];
    let mut g_log = vec![vec![NEG_INF; n]; n];
    for (candidate_index, candidate) in candidates.iter().enumerate() {
        let i = candidate.from - 1;
        let j = candidate.to - 1;
        if candidate.energy_kcal_mol < g_mfe[i][j] {
            g_mfe[i][j] = candidate.energy_kcal_mol;
            g_choice[i][j] = Some(candidate_index);
        }
        update(&mut g_log[i][j], -candidate.energy_kcal_mol / rt);
    }

    let (mfe_energy, structure) = mfe(bases, min_loop, model, &g_mfe, &g_choice, candidates)?;
    let ensemble = partition(bases, min_loop, model, rt, &g_log, candidates)?;
    Ok(IntegratedGquadResult {
        structure,
        mfe_energy,
        log_partition_function: ensemble.log_z,
        pair_probabilities: ensemble.pair_probabilities,
        unpaired_probabilities: ensemble.unpaired_probabilities,
        gquad_position_probabilities: ensemble.gquad_position_probabilities,
        interval_probabilities: ensemble.interval_probabilities,
        candidate_probabilities: ensemble.candidate_probabilities,
    })
}

fn mfe(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    g: &[Vec<f64>],
    g_choice: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
) -> Result<(f64, String), RnaError> {
    let n = bases.len();
    let mut q = vec![vec![0.0; n]; n];
    let mut qb = vec![vec![INF; n]; n];
    let mut m1 = vec![vec![INF; n]; n];
    let mut m2 = vec![vec![INF; n]; n];
    // Exactly one G-quadruplex branch, plus any flanking unpaired bases.
    let mut sg = vec![vec![INF; n]; n];
    let mut q_choice = vec![vec![ExteriorChoice::Invalid; n]; n];
    let mut qb_choice = vec![vec![PairChoice::Invalid; n]; n];
    let mut m1_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    let mut m2_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    let mut sg_choice = vec![vec![SingleGquadChoice::Invalid; n]; n];
    for (i, row) in q_choice.iter_mut().enumerate() {
        row[i] = ExteriorChoice::Unpaired;
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            if span > min_loop && model.can_pair(bases[i], bases[j]) {
                let mut best = model.hairpin_energy(bases, i, j);
                let mut choice = if best.is_finite() {
                    PairChoice::Hairpin
                } else {
                    PairChoice::Invalid
                };
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
                    if k >= j {
                        break;
                    }
                    for right in 0..=(internal_limit - left) {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        if k >= l || l - k <= min_loop || qb[k][l] >= INF / 2.0 {
                            continue;
                        }
                        let energy = model.internal_energy(bases, i, j, k, l);
                        if energy.is_finite() && qb[k][l] + energy < best {
                            best = qb[k][l] + energy;
                            choice = PairChoice::Internal(k, l);
                        }
                    }
                }
                let closing =
                    model.multiloop_closing() + model.multiloop_closing_stem_energy(bases, i, j);
                if m2[i + 1][j - 1] + closing < best {
                    best = m2[i + 1][j - 1] + closing;
                    choice = PairChoice::Multiloop;
                }
                if sg[i + 1][j - 1] + closing < best {
                    best = sg[i + 1][j - 1] + closing;
                    choice = PairChoice::SingleGquad;
                }
                qb[i][j] = best;
                qb_choice[i][j] = choice;
            }

            let mut best_m1 = if m1[i][j - 1] < INF / 2.0 {
                m1[i][j - 1] + model.multiloop_unpaired()
            } else {
                INF
            };
            let mut best_m1_choice = if best_m1 < INF / 2.0 {
                SegmentChoice::Unpaired
            } else {
                SegmentChoice::Invalid
            };
            let mut best_m2 = if m2[i][j - 1] < INF / 2.0 {
                m2[i][j - 1] + model.multiloop_unpaired()
            } else {
                INF
            };
            let mut best_m2_choice = if best_m2 < INF / 2.0 {
                SegmentChoice::Unpaired
            } else {
                SegmentChoice::Invalid
            };
            let mut best_sg = if sg[i][j - 1] < INF / 2.0 {
                sg[i][j - 1] + model.multiloop_unpaired()
            } else {
                INF
            };
            let mut best_sg_choice = if best_sg < INF / 2.0 {
                SingleGquadChoice::Unpaired
            } else {
                SingleGquadChoice::Invalid
            };

            for k in i..=j {
                if j - k > min_loop && qb[k][j] < INF / 2.0 {
                    let branch = qb[k][j] + model.multiloop_stem_energy(bases, k, j);
                    let first = (k - i) as f64 * model.multiloop_unpaired() + branch;
                    if first < best_m1 {
                        best_m1 = first;
                        best_m1_choice = SegmentChoice::FirstPair(k);
                    }
                    if k > i && m1[i][k - 1] + branch < best_m1 {
                        best_m1 = m1[i][k - 1] + branch;
                        best_m1_choice = SegmentChoice::AddPair(k);
                    }
                    if k > i && m1[i][k - 1] + branch < best_m2 {
                        best_m2 = m1[i][k - 1] + branch;
                        best_m2_choice = SegmentChoice::AddPair(k);
                    }
                }
                if g[k][j] < INF / 2.0 {
                    let branch = g[k][j] + model.multiloop_branch();
                    let first = (k - i) as f64 * model.multiloop_unpaired() + branch;
                    if first < best_m1 {
                        best_m1 = first;
                        best_m1_choice = SegmentChoice::FirstGquad(k);
                    }
                    if first < best_sg {
                        best_sg = first;
                        best_sg_choice = SingleGquadChoice::First(k);
                    }
                    if k > i && m1[i][k - 1] + branch < best_m1 {
                        best_m1 = m1[i][k - 1] + branch;
                        best_m1_choice = SegmentChoice::AddGquad(k);
                    }
                    if k > i && m1[i][k - 1] + branch < best_m2 {
                        best_m2 = m1[i][k - 1] + branch;
                        best_m2_choice = SegmentChoice::AddGquad(k);
                    }
                }
            }
            m1[i][j] = best_m1;
            m2[i][j] = best_m2;
            sg[i][j] = best_sg;
            m1_choice[i][j] = best_m1_choice;
            m2_choice[i][j] = best_m2_choice;
            sg_choice[i][j] = best_sg_choice;

            let mut best_q = q[i][j - 1];
            let mut best_q_choice = ExteriorChoice::Unpaired;
            for k in i..=j {
                if j - k > min_loop && qb[k][j] < INF / 2.0 {
                    let candidate = interval(&q, i, k.saturating_sub(1))
                        + qb[k][j]
                        + model.exterior_stem_energy(bases, k, j);
                    if candidate < best_q {
                        best_q = candidate;
                        best_q_choice = ExteriorChoice::Pair(k);
                    }
                }
                if g[k][j] < INF / 2.0 {
                    let candidate = interval(&q, i, k.saturating_sub(1)) + g[k][j];
                    if candidate < best_q {
                        best_q = candidate;
                        best_q_choice = ExteriorChoice::Gquad(k);
                    }
                }
            }
            q[i][j] = best_q;
            q_choice[i][j] = best_q_choice;
        }
    }

    let mut symbols = vec!['.'; n];
    trace_q(
        0,
        n - 1,
        &q_choice,
        &qb_choice,
        &m1_choice,
        &m2_choice,
        &sg_choice,
        g_choice,
        candidates,
        &mut symbols,
    );
    Ok((q[0][n - 1], symbols.into_iter().collect()))
}

#[allow(clippy::too_many_arguments)]
fn trace_q(
    i: usize,
    j: usize,
    q_choice: &[Vec<ExteriorChoice>],
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    sg_choice: &[Vec<SingleGquadChoice>],
    g_choice: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
    symbols: &mut [char],
) {
    if i >= j {
        return;
    }
    match q_choice[i][j] {
        ExteriorChoice::Unpaired => trace_q(
            i,
            j - 1,
            q_choice,
            qb_choice,
            m1_choice,
            m2_choice,
            sg_choice,
            g_choice,
            candidates,
            symbols,
        ),
        ExteriorChoice::Pair(k) => {
            if k > i {
                trace_q(
                    i,
                    k - 1,
                    q_choice,
                    qb_choice,
                    m1_choice,
                    m2_choice,
                    sg_choice,
                    g_choice,
                    candidates,
                    symbols,
                );
            }
            trace_pair(
                k, j, qb_choice, m1_choice, m2_choice, sg_choice, g_choice, candidates, symbols,
            );
        }
        ExteriorChoice::Gquad(k) => {
            if k > i {
                trace_q(
                    i,
                    k - 1,
                    q_choice,
                    qb_choice,
                    m1_choice,
                    m2_choice,
                    sg_choice,
                    g_choice,
                    candidates,
                    symbols,
                );
            }
            trace_gquad(k, j, g_choice, candidates, symbols);
        }
        ExteriorChoice::Invalid => {}
    }
}

#[allow(clippy::too_many_arguments)]
fn trace_pair(
    i: usize,
    j: usize,
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    sg_choice: &[Vec<SingleGquadChoice>],
    g_choice: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
    symbols: &mut [char],
) {
    symbols[i] = '(';
    symbols[j] = ')';
    match qb_choice[i][j] {
        PairChoice::Hairpin | PairChoice::Invalid => {}
        PairChoice::Internal(k, l) => trace_pair(
            k, l, qb_choice, m1_choice, m2_choice, sg_choice, g_choice, candidates, symbols,
        ),
        PairChoice::Multiloop => trace_segment(
            i + 1,
            j - 1,
            true,
            qb_choice,
            m1_choice,
            m2_choice,
            sg_choice,
            g_choice,
            candidates,
            symbols,
        ),
        PairChoice::SingleGquad => {
            trace_single_gquad(i + 1, j - 1, sg_choice, g_choice, candidates, symbols)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn trace_segment(
    i: usize,
    j: usize,
    second: bool,
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    sg_choice: &[Vec<SingleGquadChoice>],
    g_choice: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
    symbols: &mut [char],
) {
    if i > j {
        return;
    }
    let choice = if second {
        m2_choice[i][j]
    } else {
        m1_choice[i][j]
    };
    match choice {
        SegmentChoice::Unpaired => trace_segment(
            i,
            j - 1,
            second,
            qb_choice,
            m1_choice,
            m2_choice,
            sg_choice,
            g_choice,
            candidates,
            symbols,
        ),
        SegmentChoice::FirstPair(k) => trace_pair(
            k, j, qb_choice, m1_choice, m2_choice, sg_choice, g_choice, candidates, symbols,
        ),
        SegmentChoice::FirstGquad(k) => trace_gquad(k, j, g_choice, candidates, symbols),
        SegmentChoice::AddPair(k) => {
            trace_segment(
                i,
                k - 1,
                false,
                qb_choice,
                m1_choice,
                m2_choice,
                sg_choice,
                g_choice,
                candidates,
                symbols,
            );
            trace_pair(
                k, j, qb_choice, m1_choice, m2_choice, sg_choice, g_choice, candidates, symbols,
            );
        }
        SegmentChoice::AddGquad(k) => {
            trace_segment(
                i,
                k - 1,
                false,
                qb_choice,
                m1_choice,
                m2_choice,
                sg_choice,
                g_choice,
                candidates,
                symbols,
            );
            trace_gquad(k, j, g_choice, candidates, symbols);
        }
        SegmentChoice::Invalid => {}
    }
}

fn trace_single_gquad(
    i: usize,
    j: usize,
    choices: &[Vec<SingleGquadChoice>],
    g_choice: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
    symbols: &mut [char],
) {
    if i > j {
        return;
    }
    match choices[i][j] {
        SingleGquadChoice::Unpaired => {
            trace_single_gquad(i, j - 1, choices, g_choice, candidates, symbols)
        }
        SingleGquadChoice::First(k) => trace_gquad(k, j, g_choice, candidates, symbols),
        SingleGquadChoice::Invalid => {}
    }
}

fn trace_gquad(
    i: usize,
    j: usize,
    choices: &[Vec<Option<usize>>],
    candidates: &[GQuadruplex],
    symbols: &mut [char],
) {
    if let Some(index) = choices[i][j] {
        for &position in &candidates[index].guanine_positions {
            symbols[position - 1] = '+';
        }
    }
}

struct EnsembleResult {
    log_z: f64,
    pair_probabilities: Vec<PairProbability>,
    unpaired_probabilities: Vec<f64>,
    gquad_position_probabilities: Vec<f64>,
    interval_probabilities: Vec<(usize, usize, f64)>,
    candidate_probabilities: Vec<f64>,
}

fn partition(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    rt: f64,
    g: &[Vec<f64>],
    candidates: &[GQuadruplex],
) -> Result<EnsembleResult, RnaError> {
    let n = bases.len();
    let mut q = vec![vec![NEG_INF; n]; n];
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    let mut sg = vec![vec![NEG_INF; n]; n];
    for (i, row) in q.iter_mut().enumerate() {
        row[i] = 0.0;
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            if span > min_loop && model.can_pair(bases[i], bases[j]) {
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
                    for right in 0..=(internal_limit - left) {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        if k >= l || l - k <= min_loop || qb[k][l] == NEG_INF {
                            continue;
                        }
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        if energy.is_finite() {
                            let child = qb[k][l];
                            update(&mut qb[i][j], child - energy / rt);
                        }
                    }
                }
                let closing = model.multiloop_closing_boltzmann()
                    + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                if m2[i + 1][j - 1] != NEG_INF {
                    update(&mut qb[i][j], m2[i + 1][j - 1] - closing / rt);
                }
                if sg[i + 1][j - 1] != NEG_INF {
                    update(&mut qb[i][j], sg[i + 1][j - 1] - closing / rt);
                }
            }

            let unpaired_ml = -model.multiloop_unpaired_boltzmann() / rt;
            if m1[i][j - 1] != NEG_INF {
                m1[i][j] = m1[i][j - 1] + unpaired_ml;
            }
            if m2[i][j - 1] != NEG_INF {
                m2[i][j] = m2[i][j - 1] + unpaired_ml;
            }
            if sg[i][j - 1] != NEG_INF {
                sg[i][j] = sg[i][j - 1] + unpaired_ml;
            }
            for k in i..=j {
                if j - k > min_loop && qb[k][j] != NEG_INF {
                    let branch = qb[k][j] - model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                    let first = branch + (k - i) as f64 * unpaired_ml;
                    update(&mut m1[i][j], first);
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        update(&mut m1[i][j], additional);
                        update(&mut m2[i][j], additional);
                    }
                }
                if g[k][j] != NEG_INF {
                    let branch = g[k][j] - model.multiloop_branch_boltzmann() / rt;
                    let first = branch + (k - i) as f64 * unpaired_ml;
                    update(&mut m1[i][j], first);
                    update(&mut sg[i][j], first);
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        update(&mut m1[i][j], additional);
                        update(&mut m2[i][j], additional);
                    }
                }
            }

            q[i][j] = q[i][j - 1];
            for k in i..=j {
                let left = if k == i {
                    0.0
                } else {
                    log_interval(&q, i, k - 1)
                };
                if j - k > min_loop && qb[k][j] != NEG_INF {
                    let branch = qb[k][j] - model.exterior_stem_boltzmann_energy(bases, k, j) / rt;
                    update(&mut q[i][j], left + branch);
                }
                if g[k][j] != NEG_INF {
                    update(&mut q[i][j], left + g[k][j]);
                }
            }
        }
    }
    let log_z = q[0][n - 1];
    if log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "integrated G-quadruplex grammar admits no structure".into(),
        ));
    }

    let mut oq = vec![vec![NEG_INF; n]; n];
    let mut oqb = vec![vec![NEG_INF; n]; n];
    let mut om1 = vec![vec![NEG_INF; n]; n];
    let mut om2 = vec![vec![NEG_INF; n]; n];
    let mut osg = vec![vec![NEG_INF; n]; n];
    let mut og = vec![vec![NEG_INF; n]; n];
    oq[0][n - 1] = 0.0;

    for span in (1..n).rev() {
        for i in 0..(n - span) {
            let j = i + span;
            let q_parent = oq[i][j];
            if q_parent != NEG_INF {
                update(&mut oq[i][j - 1], q_parent);
                for k in i..=j {
                    let left = if k == i {
                        0.0
                    } else {
                        log_interval(&q, i, k - 1)
                    };
                    if j - k > min_loop && qb[k][j] != NEG_INF {
                        let stem = -model.exterior_stem_boltzmann_energy(bases, k, j) / rt;
                        update(&mut oqb[k][j], q_parent + left + stem);
                        if k > i {
                            update(&mut oq[i][k - 1], q_parent + stem + qb[k][j]);
                        }
                    }
                    if g[k][j] != NEG_INF {
                        update(&mut og[k][j], q_parent + left);
                        if k > i {
                            update(&mut oq[i][k - 1], q_parent + g[k][j]);
                        }
                    }
                }
            }

            let unpaired_ml = -model.multiloop_unpaired_boltzmann() / rt;
            let m2_parent = om2[i][j];
            if m2_parent != NEG_INF {
                if m2[i][j - 1] != NEG_INF {
                    update(&mut om2[i][j - 1], m2_parent + unpaired_ml);
                }
                for k in (i + 1)..=j {
                    if m1[i][k - 1] == NEG_INF {
                        continue;
                    }
                    if j - k > min_loop && qb[k][j] != NEG_INF {
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        update(&mut om1[i][k - 1], m2_parent + stem + qb[k][j]);
                        update(&mut oqb[k][j], m2_parent + m1[i][k - 1] + stem);
                    }
                    if g[k][j] != NEG_INF {
                        let stem = -model.multiloop_branch_boltzmann() / rt;
                        update(&mut om1[i][k - 1], m2_parent + stem + g[k][j]);
                        update(&mut og[k][j], m2_parent + m1[i][k - 1] + stem);
                    }
                }
            }

            let m1_parent = om1[i][j];
            if m1_parent != NEG_INF {
                if m1[i][j - 1] != NEG_INF {
                    update(&mut om1[i][j - 1], m1_parent + unpaired_ml);
                }
                for k in i..=j {
                    if j - k > min_loop && qb[k][j] != NEG_INF {
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        let leading = (k - i) as f64 * unpaired_ml;
                        update(&mut oqb[k][j], m1_parent + leading + stem);
                        if k > i && m1[i][k - 1] != NEG_INF {
                            update(&mut om1[i][k - 1], m1_parent + stem + qb[k][j]);
                            update(&mut oqb[k][j], m1_parent + m1[i][k - 1] + stem);
                        }
                    }
                    if g[k][j] != NEG_INF {
                        let stem = -model.multiloop_branch_boltzmann() / rt;
                        let leading = (k - i) as f64 * unpaired_ml;
                        update(&mut og[k][j], m1_parent + leading + stem);
                        if k > i && m1[i][k - 1] != NEG_INF {
                            update(&mut om1[i][k - 1], m1_parent + stem + g[k][j]);
                            update(&mut og[k][j], m1_parent + m1[i][k - 1] + stem);
                        }
                    }
                }
            }

            let sg_parent = osg[i][j];
            if sg_parent != NEG_INF {
                if sg[i][j - 1] != NEG_INF {
                    update(&mut osg[i][j - 1], sg_parent + unpaired_ml);
                }
                for k in i..=j {
                    if g[k][j] == NEG_INF {
                        continue;
                    }
                    let stem = -model.multiloop_branch_boltzmann() / rt;
                    let leading = (k - i) as f64 * unpaired_ml;
                    update(&mut og[k][j], sg_parent + leading + stem);
                }
            }

            let qb_parent = oqb[i][j];
            if qb_parent != NEG_INF && qb[i][j] != NEG_INF {
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
                    if k >= j {
                        break;
                    }
                    for right in 0..=(internal_limit - left) {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        if k >= l || l - k <= min_loop || qb[k][l] == NEG_INF {
                            continue;
                        }
                        let energy = model.internal_boltzmann_energy(bases, i, j, k, l);
                        if energy.is_finite() {
                            update(&mut oqb[k][l], qb_parent - energy / rt);
                        }
                    }
                }
                let closing = model.multiloop_closing_boltzmann()
                    + model.multiloop_closing_stem_boltzmann_energy(bases, i, j);
                if m2[i + 1][j - 1] != NEG_INF {
                    update(&mut om2[i + 1][j - 1], qb_parent - closing / rt);
                }
                if sg[i + 1][j - 1] != NEG_INF {
                    update(&mut osg[i + 1][j - 1], qb_parent - closing / rt);
                }
            }
        }
    }

    let mut pair_probabilities = Vec::new();
    let mut occupied = vec![0.0; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if oqb[i][j] == NEG_INF || qb[i][j] == NEG_INF {
                continue;
            }
            let probability = (oqb[i][j] + qb[i][j] - log_z).exp().clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                pair_probabilities.push(PairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability,
                });
                occupied[i] += probability;
                occupied[j] += probability;
            }
        }
    }

    let mut interval_probabilities = Vec::new();
    for i in 0..n {
        for j in i..n {
            if og[i][j] == NEG_INF || g[i][j] == NEG_INF {
                continue;
            }
            let probability = (og[i][j] + g[i][j] - log_z).exp().clamp(0.0, 1.0);
            if probability > 1.0e-12 {
                interval_probabilities.push((i + 1, j + 1, probability));
            }
        }
    }

    let mut candidate_probabilities = vec![0.0; candidates.len()];
    let mut gquad_position_probabilities = vec![0.0; n];
    for (index, candidate) in candidates.iter().enumerate() {
        let i = candidate.from - 1;
        let j = candidate.to - 1;
        if og[i][j] == NEG_INF {
            continue;
        }
        let probability = (og[i][j] - candidate.energy_kcal_mol / rt - log_z)
            .exp()
            .clamp(0.0, 1.0);
        candidate_probabilities[index] = probability;
        for &position in &candidate.guanine_positions {
            gquad_position_probabilities[position - 1] += probability;
        }
    }
    for (mass, &gquad_mass) in occupied.iter_mut().zip(&gquad_position_probabilities) {
        *mass += gquad_mass;
    }
    let unpaired_probabilities = occupied
        .into_iter()
        .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
        .collect();

    Ok(EnsembleResult {
        log_z,
        pair_probabilities,
        unpaired_probabilities,
        gquad_position_probabilities,
        interval_probabilities,
        candidate_probabilities,
    })
}
