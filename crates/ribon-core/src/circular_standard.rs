//! Circular RNA MFE and partition functions.
//!
//! The ordinary paired and multibranch interval states are closed by a
//! circular root that is either fully unpaired, an exterior hairpin, an
//! exterior internal loop, or a multiloop with at least three branches.

use crate::constraints::{ConstraintConfig, ConstraintModel, ConstraintSummary};
use crate::decode::{decode_centroid_with_constraints, decode_mea_with_constraints};
use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::Serialize;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const INF: f64 = 1.0e100;
const NEG_INF: f64 = f64::NEG_INFINITY;
const CIRC_ALPHA0: f64 = 4.385;

#[derive(Clone, Debug, Serialize)]
pub struct CircularResult {
    pub sequence: String,
    pub length: usize,
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub centroid_structure: String,
    pub centroid_score: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub mea_gamma: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub constraints: ConstraintSummary,
    pub model: &'static str,
}

#[derive(Clone, Copy)]
enum VChoice {
    Hairpin,
    Stack,
    Internal(usize, usize),
    Multiloop,
    Invalid,
}

#[derive(Clone, Copy)]
enum MChoice {
    Invalid,
    Unpaired,
    FirstPair(usize),
    AddPair(usize),
}

#[derive(Clone, Copy)]
enum RootChoice {
    Open,
    Hairpin(usize, usize),
    Internal(usize, usize, usize, usize),
    Multiloop(usize, usize),
}

struct MfeWorkspace {
    v: Vec<Vec<f64>>,
    v_secure: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
    v_choice: Vec<Vec<VChoice>>,
    m1_choice: Vec<Vec<MChoice>>,
    m2_choice: Vec<Vec<MChoice>>,
}

struct PfWorkspace {
    qb: Vec<Vec<f64>>,
    m1: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
}

#[allow(clippy::too_many_arguments)]
pub fn circular_fold(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
) -> Result<CircularResult, RnaError> {
    circular_fold_with_constraints(
        sequence,
        temperature_celsius,
        min_loop,
        gamma,
        dangles,
        salt_molar,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn circular_fold_with_constraints(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
    constraint_config: &ConstraintConfig,
) -> Result<CircularResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    circular_fold_with_model(sequence, min_loop, gamma, &model, constraint_config)
}

/// Circular ensemble with an explicitly selected RNA or DNA parameter family.
pub fn circular_fold_with_model(
    sequence: &str,
    min_loop: usize,
    gamma: f64,
    model: &EnergyModel,
    constraint_config: &ConstraintConfig,
) -> Result<CircularResult, RnaError> {
    if !gamma.is_finite() || gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "circular MEA gamma must be finite and positive".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    let n = sequence.len();
    let temperature_celsius = model.temperature_celsius();
    let dangles = model.dangles();
    let salt_molar = model.salt_molar();
    let constraints = ConstraintModel::compile(n, constraint_config)?;
    let mfe_workspace = fill_mfe(sequence.as_bytes(), min_loop, model, &constraints);
    let (mfe_energy, root) = circular_mfe_root(
        sequence.as_bytes(),
        min_loop,
        temperature_celsius,
        model,
        &mfe_workspace,
        &constraints,
    );
    let mut mfe_pairs = Vec::new();
    traceback_root(
        root,
        &mfe_workspace,
        constraints.no_lonely_pairs(),
        &mut mfe_pairs,
    );
    mfe_pairs.sort_unstable();
    let mfe_structure = pairs_to_dot_bracket(n, &mfe_pairs);

    let pf_workspace = fill_partition(
        sequence.as_bytes(),
        min_loop,
        temperature_celsius,
        model,
        &constraints,
    );
    let (log_z, pair_probabilities, unpaired_probabilities) = circular_probabilities(
        sequence.as_bytes(),
        min_loop,
        temperature_celsius,
        model,
        &pf_workspace,
        &constraints,
    )?;
    let (mfe_energy, mfe_structure, log_z, pair_probabilities, unpaired_probabilities) =
        if dangles % 2 == 1 {
            let exact = circular_odd_exact(
                &sequence,
                min_loop,
                temperature_celsius,
                model,
                &constraints,
            )?;
            (
                exact.mfe_energy_kcal_mol,
                exact.mfe_structure,
                exact.log_partition_function,
                exact.pair_probabilities,
                exact.unpaired_probabilities,
            )
        } else {
            (
                mfe_energy,
                mfe_structure,
                log_z,
                pair_probabilities,
                unpaired_probabilities,
            )
        };
    let (centroid_structure, centroid_score) = decode_centroid_with_constraints(
        sequence.as_bytes(),
        min_loop,
        &pair_probabilities,
        &constraints,
        model,
    )?;
    let (mea_structure, mea_score) = decode_mea_with_constraints(
        sequence.as_bytes(),
        min_loop,
        gamma,
        &pair_probabilities,
        &unpaired_probabilities,
        &constraints,
        model,
    )?;
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    Ok(CircularResult {
        sequence,
        length: n,
        mfe_structure,
        mfe_energy_kcal_mol: mfe_energy,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        centroid_structure,
        centroid_score,
        mea_structure,
        mea_score,
        mea_gamma: gamma,
        pair_probabilities,
        unpaired_probabilities,
        temperature_celsius,
        dangles,
        salt_molar,
        constraints: constraints.summary(),
        model: if model.parameter_profile_name().is_some() {
            if dangles % 2 == 1 {
                "exact constrained circular custom ensemble with fixed-structure single-dangle/coaxial cycle evaluation"
            } else {
                "constrained circular custom root grammar with circular entropy penalty"
            }
        } else if model.nucleic_acid() == crate::energy::NucleicAcid::Dna {
            if dangles % 2 == 1 {
                "exact constrained circular RNAstructure 6.6 DNA ensemble with fixed-structure single-dangle/coaxial cycle evaluation"
            } else {
                "constrained circular RNAstructure 6.6 DNA root grammar with circular entropy penalty"
            }
        } else if dangles % 2 == 1 {
            "exact constrained circular RNAstructure 6.6 RNA ensemble with fixed-structure single-dangle/coaxial cycle evaluation"
        } else {
            "constrained circular RNAstructure 6.6 RNA root grammar with circular entropy penalty"
        },
    })
}

struct CircularExactEnsemble {
    mfe_structure: String,
    mfe_energy_kcal_mol: f64,
    log_partition_function: f64,
    pair_probabilities: Vec<PairProbability>,
    unpaired_probabilities: Vec<f64>,
}

fn circular_odd_exact(
    sequence: &str,
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<CircularExactEnsemble, RnaError> {
    let n = sequence.len();
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    let mut log_z = NEG_INF;
    let mut pair_log_mass = vec![vec![NEG_INF; n]; n];
    let mut mfe_energy = INF;
    let mut mfe_structure = String::new();
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
            let energy = evaluate_circular_structure(
                sequence,
                &structure,
                min_loop,
                temperature,
                model,
                constraints,
            )?;
            if !energy.is_finite() {
                return Ok(());
            }
            if energy < mfe_energy
                || (energy == mfe_energy && (mfe_structure.is_empty() || structure < mfe_structure))
            {
                mfe_energy = energy;
                mfe_structure.clone_from(&structure);
            }
            let log_weight = -energy / rt;
            log_z = log_add(log_z, log_weight);
            for &(i, j) in pairs {
                pair_log_mass[i][j] = log_add(pair_log_mass[i][j], log_weight);
            }
            Ok(())
        },
    )?;
    if log_z == NEG_INF || mfe_structure.is_empty() {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid circular secondary structure".into(),
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
    Ok(CircularExactEnsemble {
        mfe_structure,
        mfe_energy_kcal_mol: mfe_energy,
        log_partition_function: log_z,
        pair_probabilities,
        unpaired_probabilities: paired_mass
            .into_iter()
            .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
            .collect(),
    })
}

fn evaluate_circular_structure(
    sequence: &str,
    structure: &str,
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<f64, RnaError> {
    let parsed = crate::structure::parse_structure(sequence, structure)?;
    let linear = model.evaluate_with_constraints(sequence, structure, min_loop, constraints)?;
    let mut top_level = Vec::new();
    for pair in &parsed.pairs {
        let i = pair.i - 1;
        let j = pair.j - 1;
        if !parsed
            .pairs
            .iter()
            .any(|outer| outer.i - 1 < i && j < outer.j - 1)
        {
            top_level.push((i, j));
        }
    }
    top_level.sort_unstable();
    let bases = sequence.as_bytes();
    let root = match top_level.as_slice() {
        [] => circular_exterior_penalty(sequence.len(), temperature),
        &[(i, j)] => circular_hairpin_energy(bases, i, j, model, false),
        &[(i, j), (p, q)] => circular_internal_energy(bases, i, j, p, q, model, false),
        stems => {
            let occupied = stems.iter().map(|&(i, j)| j - i + 1).sum::<usize>();
            model.multiloop_closing()
                + model.multiloop_unpaired() * (sequence.len() - occupied) as f64
                + model.circular_multiloop_stems_energy(bases, &parsed.partner, stems)
        }
    };
    Ok(linear.total_kcal_mol - linear.exterior_kcal_mol + root)
}

fn fill_mfe(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> MfeWorkspace {
    let n = bases.len();
    let mut v = vec![vec![INF; n]; n];
    let mut v_secure = vec![vec![INF; n]; n];
    let mut m1 = vec![vec![INF; n]; n];
    let mut m2 = vec![vec![INF; n]; n];
    let mut v_choice = vec![vec![VChoice::Invalid; n]; n];
    let mut m1_choice = vec![vec![MChoice::Invalid; n]; n];
    let mut m2_choice = vec![vec![MChoice::Invalid; n]; n];
    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if span > min_loop && constraints.allows_pair(bases, i, j, model) {
                let pair_soft = constraints.pair_energy(i, j);
                let mut best = constraints
                    .unpaired_range_energy(i + 1, j - 1)
                    .map(|soft| model.hairpin_energy(bases, i, j) + pair_soft + soft)
                    .unwrap_or(INF);
                let mut choice = if best < INF / 2.0 {
                    VChoice::Hairpin
                } else {
                    VChoice::Invalid
                };
                if i + 1 < j && v[i + 1][j - 1] < INF / 2.0 {
                    let energy = pair_soft
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1)
                        + v[i + 1][j - 1];
                    if energy < best {
                        best = energy;
                        choice = VChoice::Stack;
                    }
                    v_secure[i][j] = energy;
                }
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
                    if k >= j {
                        break;
                    }
                    for right in 0..=internal_limit - left {
                        if left == 0 && right == 0 {
                            continue;
                        }
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        let child = if constraints.no_lonely_pairs() {
                            v_secure[k][l]
                        } else {
                            v[k][l]
                        };
                        if k >= l || l - k <= min_loop || child >= INF / 2.0 {
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
                        let energy = pair_soft
                            + model.internal_energy(bases, i, j, k, l)
                            + constraints.stack_energy(i, j, k, l)
                            + left_soft
                            + right_soft
                            + child;
                        if energy < best {
                            best = energy;
                            choice = VChoice::Internal(k, l);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] < INF / 2.0 {
                    let energy = model.multiloop_closing()
                        + model.multiloop_closing_stem_energy(bases, i, j)
                        + pair_soft
                        + m2[i + 1][j - 1];
                    if energy < best {
                        best = energy;
                        choice = VChoice::Multiloop;
                    }
                }
                v[i][j] = best;
                v_choice[i][j] = choice;
            }

            let mut best_m1 = if m1[i][j - 1] < INF / 2.0 && constraints.allows_unpaired(j) {
                m1[i][j - 1] + model.multiloop_unpaired() + constraints.unpaired_energy(j)
            } else {
                INF
            };
            let mut choice_m1 = if best_m1 < INF / 2.0 {
                MChoice::Unpaired
            } else {
                MChoice::Invalid
            };
            let mut best_m2 = if m2[i][j - 1] < INF / 2.0 && constraints.allows_unpaired(j) {
                m2[i][j - 1] + model.multiloop_unpaired() + constraints.unpaired_energy(j)
            } else {
                INF
            };
            let mut choice_m2 = if best_m2 < INF / 2.0 {
                MChoice::Unpaired
            } else {
                MChoice::Invalid
            };
            if j > min_loop {
                for k in i..=j - min_loop - 1 {
                    let branch_v = if constraints.no_lonely_pairs() {
                        v_secure[k][j]
                    } else {
                        v[k][j]
                    };
                    if branch_v >= INF / 2.0 {
                        continue;
                    }
                    let branch = branch_v + model.multiloop_stem_energy(bases, k, j);
                    if let Some(soft) = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    } {
                        let first = (k - i) as f64 * model.multiloop_unpaired() + soft + branch;
                        if first < best_m1 {
                            best_m1 = first;
                            choice_m1 = MChoice::FirstPair(k);
                        }
                    }
                    if k > i && m1[i][k - 1] < INF / 2.0 {
                        let additional = m1[i][k - 1] + branch;
                        if additional < best_m1 {
                            best_m1 = additional;
                            choice_m1 = MChoice::AddPair(k);
                        }
                        if additional < best_m2 {
                            best_m2 = additional;
                            choice_m2 = MChoice::AddPair(k);
                        }
                    }
                }
            }
            m1[i][j] = best_m1;
            m2[i][j] = best_m2;
            m1_choice[i][j] = choice_m1;
            m2_choice[i][j] = choice_m2;
        }
    }
    MfeWorkspace {
        v,
        v_secure,
        m2,
        v_choice,
        m1_choice,
        m2_choice,
    }
}

fn circular_mfe_root(
    bases: &[u8],
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    workspace: &MfeWorkspace,
    constraints: &ConstraintModel,
) -> (f64, RootChoice) {
    let n = bases.len();
    let mut best = unpaired_segments_energy(constraints, &[(0, n)])
        .map(|soft| circular_exterior_penalty(n, temperature) + soft)
        .unwrap_or(INF);
    let mut root = RootChoice::Open;
    let internal_limit = model.internal_loop_limit(n.saturating_sub(4));
    for i in 0..n - 1 {
        for j in i + min_loop + 1..n {
            let first_branch = if constraints.no_lonely_pairs() {
                workspace.v_secure[i][j]
            } else {
                workspace.v[i][j]
            };
            if first_branch >= INF / 2.0 {
                continue;
            }
            let outside = n - j + i - 1;
            if outside >= min_loop {
                if let Some(soft) = unpaired_segments_energy(constraints, &[(j + 1, n), (0, i)]) {
                    let energy =
                        first_branch + circular_hairpin_energy(bases, i, j, model, false) + soft;
                    if energy < best {
                        best = energy;
                        root = RootChoice::Hairpin(i, j);
                    }
                }
            }
            for p in j + 1..n {
                let u1 = p - j - 1;
                if u1 + i > internal_limit {
                    break;
                }
                let qmin = p + 1;
                for q in (qmin..n).rev() {
                    let u2 = i + n - q - 1;
                    let second_branch = if constraints.no_lonely_pairs() {
                        workspace.v_secure[p][q]
                    } else {
                        workspace.v[p][q]
                    };
                    if u1 + u2 > internal_limit || second_branch >= INF / 2.0 {
                        continue;
                    }
                    if let Some(soft) =
                        unpaired_segments_energy(constraints, &[(j + 1, p), (q + 1, n), (0, i)])
                    {
                        let loop_energy = circular_internal_energy(bases, i, j, p, q, model, false);
                        let energy = first_branch + second_branch + loop_energy + soft;
                        if energy < best {
                            best = energy;
                            root = RootChoice::Internal(i, j, p, q);
                        }
                    }
                }
            }
        }
    }
    for end in 0..n.saturating_sub(1) {
        let Some((start, first)) =
            first_branch_mfe(bases, end, min_loop, model, workspace, constraints)
        else {
            continue;
        };
        if end + 1 < n && workspace.m2[end + 1][n - 1] < INF / 2.0 {
            let energy = first + workspace.m2[end + 1][n - 1] + model.multiloop_closing();
            if energy < best {
                best = energy;
                root = RootChoice::Multiloop(start, end);
            }
        }
    }
    (best, root)
}

fn first_branch_mfe(
    bases: &[u8],
    end: usize,
    min_loop: usize,
    model: &EnergyModel,
    workspace: &MfeWorkspace,
    constraints: &ConstraintModel,
) -> Option<(usize, f64)> {
    if end <= min_loop {
        return None;
    }
    let mut best = INF;
    let mut start = 0;
    for i in 0..=end - min_loop - 1 {
        let branch = if constraints.no_lonely_pairs() {
            workspace.v_secure[i][end]
        } else {
            workspace.v[i][end]
        };
        if branch < INF / 2.0 {
            let Some(soft) = unpaired_segments_energy(constraints, &[(0, i)]) else {
                continue;
            };
            let energy = i as f64 * model.multiloop_unpaired()
                + soft
                + model.multiloop_stem_energy(bases, i, end)
                + branch;
            if energy < best {
                best = energy;
                start = i;
            }
        }
    }
    (best < INF / 2.0).then_some((start, best))
}

fn traceback_root(
    root: RootChoice,
    workspace: &MfeWorkspace,
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    match root {
        RootChoice::Open => {}
        RootChoice::Hairpin(i, j) => traceback_v_branch(i, j, no_lonely_pairs, workspace, pairs),
        RootChoice::Internal(i, j, p, q) => {
            traceback_v_branch(i, j, no_lonely_pairs, workspace, pairs);
            traceback_v_branch(p, q, no_lonely_pairs, workspace, pairs);
        }
        RootChoice::Multiloop(i, j) => {
            traceback_v_branch(i, j, no_lonely_pairs, workspace, pairs);
            traceback_m(
                j + 1,
                workspace.v.len() - 1,
                false,
                no_lonely_pairs,
                workspace,
                pairs,
            );
        }
    }
}

fn traceback_v_branch(
    i: usize,
    j: usize,
    secure: bool,
    workspace: &MfeWorkspace,
    pairs: &mut Vec<(usize, usize)>,
) {
    if secure {
        pairs.push((i, j));
        traceback_v(i + 1, j - 1, true, workspace, pairs);
    } else {
        traceback_v(i, j, false, workspace, pairs);
    }
}

fn traceback_v(
    i: usize,
    j: usize,
    no_lonely_pairs: bool,
    workspace: &MfeWorkspace,
    pairs: &mut Vec<(usize, usize)>,
) {
    pairs.push((i, j));
    match workspace.v_choice[i][j] {
        VChoice::Hairpin | VChoice::Invalid => {}
        VChoice::Stack => traceback_v(i + 1, j - 1, no_lonely_pairs, workspace, pairs),
        VChoice::Internal(k, l) => traceback_v_branch(k, l, no_lonely_pairs, workspace, pairs),
        VChoice::Multiloop => traceback_m(i + 1, j - 1, false, no_lonely_pairs, workspace, pairs),
    }
}

fn traceback_m(
    i: usize,
    j: usize,
    one: bool,
    no_lonely_pairs: bool,
    workspace: &MfeWorkspace,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= workspace.v.len() || j >= workspace.v.len() || i > j {
        return;
    }
    let choice = if one {
        workspace.m1_choice[i][j]
    } else {
        workspace.m2_choice[i][j]
    };
    match choice {
        MChoice::Invalid => {}
        MChoice::Unpaired => traceback_m(i, j - 1, one, no_lonely_pairs, workspace, pairs),
        MChoice::FirstPair(k) => traceback_v_branch(k, j, no_lonely_pairs, workspace, pairs),
        MChoice::AddPair(k) => {
            traceback_m(i, k - 1, true, no_lonely_pairs, workspace, pairs);
            traceback_v_branch(k, j, no_lonely_pairs, workspace, pairs);
        }
    }
}

fn fill_partition(
    bases: &[u8],
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> PfWorkspace {
    let n = bases.len();
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    for span in 1..n {
        for i in 0..n - span {
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
                        let candidate = qb[k][l] - energy / rt;
                        update(&mut qb[i][j], candidate);
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
                        + pair_soft;
                    update(&mut qb[i][j], m2[i + 1][j - 1] - energy / rt);
                }
            }
            let log_unpaired = if constraints.allows_unpaired(j) {
                -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
            } else {
                NEG_INF
            };
            if m1[i][j - 1] != NEG_INF && log_unpaired != NEG_INF {
                m1[i][j] = m1[i][j - 1] + log_unpaired;
            }
            if m2[i][j - 1] != NEG_INF && log_unpaired != NEG_INF {
                m2[i][j] = m2[i][j - 1] + log_unpaired;
            }
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    if qb_row[j] == NEG_INF {
                        continue;
                    }
                    let branch =
                        qb_row[j] - model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                    if let Some(soft) = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    } {
                        update(
                            &mut m1[i][j],
                            branch
                                - ((k - i) as f64 * model.multiloop_unpaired_boltzmann() + soft)
                                    / rt,
                        );
                    }
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        update(&mut m1[i][j], additional);
                        update(&mut m2[i][j], additional);
                    }
                }
            }
        }
    }
    PfWorkspace { qb, m1, m2 }
}

fn circular_probabilities(
    bases: &[u8],
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    workspace: &PfWorkspace,
    constraints: &ConstraintModel,
) -> Result<(f64, Vec<PairProbability>, Vec<f64>), RnaError> {
    let n = bases.len();
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    let mut log_z = unpaired_segments_energy(constraints, &[(0, n)])
        .map(|soft| -(circular_exterior_penalty(n, temperature) + soft) / rt)
        .unwrap_or(NEG_INF);
    let mut oqb = vec![vec![NEG_INF; n]; n];
    let mut om1 = vec![vec![NEG_INF; n]; n];
    let mut om2 = vec![vec![NEG_INF; n]; n];
    let internal_limit = model.internal_loop_limit(n.saturating_sub(4));

    for i in 0..n - 1 {
        for j in i + min_loop + 1..n {
            if workspace.qb[i][j] == NEG_INF {
                continue;
            }
            let outside = n - j + i - 1;
            if outside >= min_loop {
                if let Some(soft) = unpaired_segments_energy(constraints, &[(j + 1, n), (0, i)]) {
                    let energy = circular_hairpin_energy(bases, i, j, model, true) + soft;
                    let root = workspace.qb[i][j] - energy / rt;
                    log_z = log_add(log_z, root);
                    update(&mut oqb[i][j], -energy / rt);
                }
            }
            for p in j + 1..n {
                let u1 = p - j - 1;
                if u1 + i > internal_limit {
                    break;
                }
                for q in ((p + 1)..n).rev() {
                    let u2 = i + n - q - 1;
                    if u1 + u2 > internal_limit || workspace.qb[p][q] == NEG_INF {
                        continue;
                    }
                    let Some(soft) =
                        unpaired_segments_energy(constraints, &[(j + 1, p), (q + 1, n), (0, i)])
                    else {
                        continue;
                    };
                    let energy = circular_internal_energy(bases, i, j, p, q, model, true) + soft;
                    let transition = -energy / rt;
                    let root = workspace.qb[i][j] + workspace.qb[p][q] + transition;
                    log_z = log_add(log_z, root);
                    update(&mut oqb[i][j], workspace.qb[p][q] + transition);
                    update(&mut oqb[p][q], workspace.qb[i][j] + transition);
                }
            }
        }
    }
    let closing = -model.multiloop_closing_boltzmann() / rt;
    for end in 0..n.saturating_sub(1) {
        let first = first_branch_log(
            bases,
            end,
            min_loop,
            temperature,
            model,
            workspace,
            constraints,
        );
        if first == NEG_INF || end + 1 >= n || workspace.m2[end + 1][n - 1] == NEG_INF {
            continue;
        }
        let root = first + workspace.m2[end + 1][n - 1] + closing;
        log_z = log_add(log_z, root);
        seed_first_branch_outside(
            bases,
            end,
            min_loop,
            temperature,
            model,
            workspace,
            constraints,
            workspace.m2[end + 1][n - 1] + closing,
            &mut oqb,
        );
        update(&mut om2[end + 1][n - 1], first + closing);
    }

    // Reverse differentiation of multiloop and paired interval states.
    for span in (1..n).rev() {
        for i in 0..n - span {
            let j = i + span;
            let log_unpaired = if constraints.allows_unpaired(j) {
                -(model.multiloop_unpaired_boltzmann() + constraints.unpaired_energy(j)) / rt
            } else {
                NEG_INF
            };
            let parent2 = om2[i][j];
            if parent2 != NEG_INF {
                if workspace.m2[i][j - 1] != NEG_INF && log_unpaired != NEG_INF {
                    update(&mut om2[i][j - 1], parent2 + log_unpaired);
                }
                if j > min_loop {
                    for (k, oqb_row) in oqb.iter_mut().enumerate().take(j - min_loop).skip(i + 1) {
                        if workspace.m1[i][k - 1] == NEG_INF || workspace.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        update(&mut om1[i][k - 1], parent2 + stem + workspace.qb[k][j]);
                        update(&mut oqb_row[j], parent2 + workspace.m1[i][k - 1] + stem);
                    }
                }
            }
            let parent1 = om1[i][j];
            if parent1 != NEG_INF {
                if workspace.m1[i][j - 1] != NEG_INF && log_unpaired != NEG_INF {
                    update(&mut om1[i][j - 1], parent1 + log_unpaired);
                }
                if j > min_loop {
                    for (k, oqb_row) in oqb.iter_mut().enumerate().take(j - min_loop).skip(i) {
                        if workspace.qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -model.multiloop_stem_boltzmann_energy(bases, k, j) / rt;
                        if let Some(soft) = if k == i {
                            Some(0.0)
                        } else {
                            constraints.unpaired_range_energy(i, k - 1)
                        } {
                            update(
                                &mut oqb_row[j],
                                parent1 + stem
                                    - ((k - i) as f64 * model.multiloop_unpaired_boltzmann()
                                        + soft)
                                        / rt,
                            );
                        }
                        if k > i && workspace.m1[i][k - 1] != NEG_INF {
                            update(&mut om1[i][k - 1], parent1 + stem + workspace.qb[k][j]);
                            update(&mut oqb_row[j], parent1 + workspace.m1[i][k - 1] + stem);
                        }
                    }
                }
            }
            let parent = oqb[i][j];
            if parent != NEG_INF && workspace.qb[i][j] != NEG_INF {
                let pair_soft = constraints.pair_energy(i, j);
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
                        if k >= l || l - k <= min_loop || workspace.qb[k][l] == NEG_INF {
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
                        update(&mut oqb[k][l], parent - energy / rt);
                    }
                }
                if i + 1 < j && workspace.m2[i + 1][j - 1] != NEG_INF {
                    let energy = model.multiloop_closing_boltzmann()
                        + model.multiloop_closing_stem_boltzmann_energy(bases, i, j)
                        + pair_soft;
                    update(&mut om2[i + 1][j - 1], parent - energy / rt);
                }
            }
        }
    }

    if log_z == NEG_INF {
        return Err(RnaError::Numerical(
            "circular partition function has zero weight".into(),
        ));
    }
    let mut pair_probabilities = Vec::new();
    let mut paired_mass = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            if oqb[i][j] == NEG_INF || workspace.qb[i][j] == NEG_INF {
                continue;
            }
            let probability = (oqb[i][j] + workspace.qb[i][j] - log_z)
                .exp()
                .clamp(0.0, 1.0);
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
    let unpaired = paired_mass
        .into_iter()
        .map(|value| (1.0 - value).clamp(0.0, 1.0))
        .collect();
    Ok((log_z, pair_probabilities, unpaired))
}

fn first_branch_log(
    bases: &[u8],
    end: usize,
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    workspace: &PfWorkspace,
    constraints: &ConstraintModel,
) -> f64 {
    if end <= min_loop {
        return NEG_INF;
    }
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    let mut value = NEG_INF;
    for i in 0..=end - min_loop - 1 {
        if workspace.qb[i][end] != NEG_INF {
            let Some(soft) = unpaired_segments_energy(constraints, &[(0, i)]) else {
                continue;
            };
            let branch = workspace.qb[i][end]
                - model.multiloop_stem_boltzmann_energy(bases, i, end) / rt
                - (i as f64 * model.multiloop_unpaired_boltzmann() + soft) / rt;
            value = log_add(value, branch);
        }
    }
    value
}

#[allow(clippy::too_many_arguments)]
fn seed_first_branch_outside(
    bases: &[u8],
    end: usize,
    min_loop: usize,
    temperature: f64,
    model: &EnergyModel,
    workspace: &PfWorkspace,
    constraints: &ConstraintModel,
    outside: f64,
    oqb: &mut [Vec<f64>],
) {
    if end <= min_loop {
        return;
    }
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    for (i, row) in oqb.iter_mut().enumerate().take(end - min_loop) {
        if workspace.qb[i][end] != NEG_INF {
            let Some(soft) = unpaired_segments_energy(constraints, &[(0, i)]) else {
                continue;
            };
            let transition = -model.multiloop_stem_boltzmann_energy(bases, i, end) / rt
                - (i as f64 * model.multiloop_unpaired_boltzmann() + soft) / rt;
            update(&mut row[end], outside + transition);
        }
    }
}

/// Sum pseudo-energies for disjoint half-open unpaired segments and reject a
/// decomposition as soon as it would leave a forced-paired nucleotide open.
fn unpaired_segments_energy(
    constraints: &ConstraintModel,
    segments: &[(usize, usize)],
) -> Option<f64> {
    let mut total = 0.0;
    for &(start, end) in segments {
        if start < end {
            total += constraints.unpaired_range_energy(start, end - 1)?;
        }
    }
    Some(total)
}

fn circular_exterior_penalty(length: usize, temperature: f64) -> f64 {
    let rt = GAS_CONSTANT_KCAL * (temperature + 273.15);
    ((100.0 * rt * (CIRC_ALPHA0 + 1.5 * (length as f64).ln()) + 0.5) as i64) as f64 / 100.0
}

fn circular_hairpin_energy(
    bases: &[u8],
    i: usize,
    j: usize,
    model: &EnergyModel,
    boltzmann: bool,
) -> f64 {
    let mut loop_bases = Vec::with_capacity(bases.len() - j + i + 1);
    loop_bases.push(bases[j]);
    loop_bases.extend_from_slice(&bases[j + 1..]);
    loop_bases.extend_from_slice(&bases[..i]);
    loop_bases.push(bases[i]);
    let end = loop_bases.len() - 1;
    if boltzmann {
        model.hairpin_boltzmann_energy(&loop_bases, 0, end)
    } else {
        model.hairpin_energy(&loop_bases, 0, end)
    }
}

#[allow(clippy::too_many_arguments)]
fn circular_internal_energy(
    bases: &[u8],
    i: usize,
    j: usize,
    p: usize,
    q: usize,
    model: &EnergyModel,
    boltzmann: bool,
) -> f64 {
    let mut loop_bases = Vec::with_capacity((p - j) + (bases.len() - q + i) + 3);
    loop_bases.push(bases[j]);
    loop_bases.extend_from_slice(&bases[j + 1..p]);
    let inner_i = loop_bases.len();
    loop_bases.push(bases[p]);
    let inner_j = loop_bases.len();
    loop_bases.push(bases[q]);
    loop_bases.extend_from_slice(&bases[q + 1..]);
    loop_bases.extend_from_slice(&bases[..i]);
    loop_bases.push(bases[i]);
    let end = loop_bases.len() - 1;
    if boltzmann {
        model.internal_boltzmann_energy(&loop_bases, 0, end, inner_i, inner_j)
    } else {
        model.internal_energy(&loop_bases, 0, end, inner_i, inner_j)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::{PairConstraint, PairEnergy, SoftConstraintConfig};

    #[test]
    fn unfolded_circle_includes_the_standard_entropy_penalty() {
        let result = circular_fold("AAAAAA", 37.0, 3, 1.0, 0, 1.021).unwrap();
        assert_eq!(result.mfe_structure, "......");
        assert!((result.mfe_energy_kcal_mol - circular_exterior_penalty(6, 37.0)).abs() < 1e-12);
        assert!((result.ensemble_free_energy_kcal_mol - result.mfe_energy_kcal_mol).abs() < 1e-12);
    }

    #[test]
    fn reference_hairpin_case_matches_vienna_public_values() {
        let result = circular_fold("GGGAAACCCU", 37.0, 3, 1.0, 0, 1.021).unwrap();
        assert_eq!(result.mfe_structure, "..........");
        assert!((result.mfe_energy_kcal_mol - 4.83).abs() < 1e-9);
        assert!((result.ensemble_free_energy_kcal_mol - 4.81247326464).abs() < 2e-8);
    }

    #[test]
    fn hard_and_soft_constraints_flow_through_the_circular_root() {
        let sequence = "GGGAAACCCU";
        let hard = ConstraintConfig {
            force_pairs: vec![PairConstraint { i: 1, j: 7 }],
            force_unpaired: vec![2, 3, 4, 5, 6, 8, 9, 10],
            ..ConstraintConfig::default()
        };
        let constrained =
            circular_fold_with_constraints(sequence, 37.0, 3, 1.0, 0, 1.021, &hard).unwrap();
        assert_eq!(constrained.mfe_structure.as_bytes()[0], b'(');
        assert_eq!(constrained.mfe_structure.as_bytes()[6], b')');
        let forced_probability = constrained
            .pair_probabilities
            .iter()
            .find(|pair| pair.i == 1 && pair.j == 7)
            .unwrap()
            .probability;
        assert!((forced_probability - 1.0).abs() < 1.0e-10);
        for (index, &unpaired) in constrained.unpaired_probabilities.iter().enumerate() {
            let paired: f64 = constrained
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == index + 1 || pair.j == index + 1)
                .map(|pair| pair.probability)
                .sum();
            assert!((paired + unpaired - 1.0).abs() < 1.0e-8);
        }

        let soft = ConstraintConfig {
            soft: SoftConstraintConfig {
                pairs: vec![PairEnergy {
                    i: 1,
                    j: 7,
                    energy_kcal_mol: -2.0,
                }],
                ..SoftConstraintConfig::default()
            },
            ..ConstraintConfig::default()
        };
        let baseline = circular_fold(sequence, 37.0, 3, 1.0, 0, 1.021).unwrap();
        let biased =
            circular_fold_with_constraints(sequence, 37.0, 3, 1.0, 0, 1.021, &soft).unwrap();
        let baseline_pair = baseline
            .pair_probabilities
            .iter()
            .find(|pair| pair.i == 1 && pair.j == 7)
            .map_or(0.0, |pair| pair.probability);
        let biased_pair = biased
            .pair_probabilities
            .iter()
            .find(|pair| pair.i == 1 && pair.j == 7)
            .map_or(0.0, |pair| pair.probability);
        assert!(biased_pair > baseline_pair);
    }

    #[test]
    fn integrated_double_dangles_are_finite_and_normalized() {
        let result = circular_fold("GGGAAACCCUGGGAAACCCU", 37.0, 3, 1.0, 2, 1.021).unwrap();
        assert!(result.mfe_energy_kcal_mol.is_finite());
        assert!(result.log_partition_function.is_finite());
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
    fn every_requested_dangle_model_runs_through_the_circular_grammar() {
        for dangles in 0..=3 {
            let result =
                circular_fold("GGGAAACCCUGGGAAACCCU", 37.0, 3, 1.0, dangles, 1.021).unwrap();
            assert_eq!(result.dangles, dangles);
            assert!(result.mfe_energy_kcal_mol.is_finite());
            assert!(result.log_partition_function.is_finite());
        }
    }

    #[test]
    fn odd_models_evaluate_shared_dangles_and_coaxial_cycles_explicitly() {
        let sequence = "GAAACAGAAACAGAAAC";
        let constraints = ConstraintConfig {
            force_pairs: vec![
                PairConstraint { i: 1, j: 5 },
                PairConstraint { i: 7, j: 11 },
                PairConstraint { i: 13, j: 17 },
            ],
            ..ConstraintConfig::default()
        };
        let single =
            circular_fold_with_constraints(sequence, 37.0, 3, 1.0, 1, 1.021, &constraints).unwrap();
        let double =
            circular_fold_with_constraints(sequence, 37.0, 3, 1.0, 2, 1.021, &constraints).unwrap();
        let coaxial =
            circular_fold_with_constraints(sequence, 37.0, 3, 1.0, 3, 1.021, &constraints).unwrap();
        assert_eq!(single.mfe_structure, double.mfe_structure);
        assert_eq!(coaxial.mfe_structure, double.mfe_structure);
        assert_ne!(single.mfe_energy_kcal_mol, double.mfe_energy_kcal_mol);
        assert_ne!(coaxial.mfe_energy_kcal_mol, double.mfe_energy_kcal_mol);
        assert!(single.model.contains("single-dangle/coaxial"));
    }
}
