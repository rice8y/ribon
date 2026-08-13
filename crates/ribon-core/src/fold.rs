use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::energy::{CoaxialStack, EnergyModel, LoopEnergy};
use crate::structure::{normalize_sequence, pairs_to_dot_bracket, RnaError};
use serde::Serialize;

const INF: f64 = 1.0e100;

#[derive(Clone, Debug, Serialize)]
pub struct MfeResult {
    pub sequence: String,
    pub structure: String,
    pub energy_kcal_mol: f64,
    pub evaluated_energy_kcal_mol: f64,
    pub energy_difference_kcal_mol: f64,
    pub energy_difference_reason: Option<&'static str>,
    pub evaluated_coaxial_stacks: Vec<CoaxialStack>,
    pub evaluated_loop_energies: Vec<LoopEnergy>,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub salt_correction: bool,
    pub model: &'static str,
    pub constraints: ConstraintSummary,
}

#[derive(Clone, Copy, Debug)]
enum VChoice {
    Hairpin,
    Stack,
    Internal(usize, usize),
    Multiloop,
    Invalid,
}

#[derive(Clone, Copy, Debug)]
enum FChoice {
    Empty,
    Unpaired,
    Pair(usize),
}

#[derive(Clone, Copy, Debug)]
enum MChoice {
    Invalid,
    Unpaired,
    FirstPair(usize),
    AddPair(usize),
}

pub fn fold_mfe(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
) -> Result<MfeResult, RnaError> {
    let length = normalize_sequence(sequence)?.len();
    let constraints = ConstraintModel::unconstrained(length);
    fold_mfe_with_constraints(sequence, min_loop, model, &constraints)
}

pub fn fold_mfe_with_constraints(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<MfeResult, RnaError> {
    if model.dangles() % 2 == 1 {
        return fold_mfe_odd(sequence, min_loop, model, constraints);
    }
    fold_mfe_even(sequence, min_loop, model, constraints)
}

fn fold_mfe_even(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<MfeResult, RnaError> {
    let sequence = normalize_sequence(sequence)?;
    let bases = sequence.as_bytes();
    let n = bases.len();
    if n == 0 {
        return Err(RnaError::EmptySequence);
    }

    let mut f = vec![vec![0.0f64; n]; n];
    let mut v = vec![vec![INF; n]; n];
    let mut v_secure = vec![vec![INF; n]; n];
    let mut m1 = vec![vec![INF; n]; n];
    let mut m2 = vec![vec![INF; n]; n];
    let mut f_choice = vec![vec![FChoice::Empty; n]; n];
    let mut v_choice = vec![vec![VChoice::Invalid; n]; n];
    let mut m1_choice = vec![vec![MChoice::Invalid; n]; n];
    let mut m2_choice = vec![vec![MChoice::Invalid; n]; n];
    for (i, row) in f.iter_mut().enumerate() {
        row[i] = if constraints.allows_unpaired(i) {
            constraints.unpaired_energy(i)
        } else {
            INF
        };
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;

            if span > min_loop && constraints.allows_pair(bases, i, j, model) {
                let pair = model.pair_energy(bases[i], bases[j]) + constraints.pair_energy(i, j);
                let hairpin_soft = constraints.unpaired_range_energy(i + 1, j - 1);
                let mut best = hairpin_soft
                    .map(|soft| pair + model.hairpin_energy(bases, i, j) + soft)
                    .unwrap_or(INF);
                let mut choice = if best < INF / 2.0 {
                    VChoice::Hairpin
                } else {
                    VChoice::Invalid
                };

                if span >= 2 && v[i + 1][j - 1] < INF / 2.0 {
                    let candidate = pair
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1)
                        + v[i + 1][j - 1];
                    if candidate < best {
                        best = candidate;
                        choice = VChoice::Stack;
                    }
                }

                // Enumerate every geometrically possible bulge and internal
                // loop unless the caller explicitly restricts the model.
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    let remaining = internal_limit - left;
                    for right in 0..=remaining {
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
                        let candidate = pair
                            + model.internal_energy(bases, i, j, k, l)
                            + left_soft
                            + right_soft
                            + child;
                        if candidate < best {
                            best = candidate;
                            choice = VChoice::Internal(k, l);
                        }
                    }
                }

                // An affine multiloop needs at least two enclosed branches.
                if span >= 2 && m2[i + 1][j - 1] < INF / 2.0 {
                    let candidate = pair
                        + model.multiloop_closing()
                        + model.multiloop_closing_stem_energy(bases, i, j)
                        + constraints.terminal_pair_energy(i, j)
                        + m2[i + 1][j - 1];
                    if candidate < best {
                        best = candidate;
                        choice = VChoice::Multiloop;
                    }
                }
                v[i][j] = best;
                v_choice[i][j] = choice;
                if i + 1 < j && v[i + 1][j - 1] < INF / 2.0 {
                    v_secure[i][j] = pair
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1)
                        + v[i + 1][j - 1];
                }
            }

            let mut best = if constraints.allows_unpaired(j) {
                interval(&f, i, j.saturating_sub(1)) + constraints.unpaired_energy(j)
            } else {
                INF
            };
            let mut choice = if best < INF / 2.0 {
                FChoice::Unpaired
            } else {
                FChoice::Empty
            };
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branch_v = if constraints.no_lonely_pairs() {
                        v_secure[k][j]
                    } else {
                        v[k][j]
                    };
                    if branch_v >= INF / 2.0 {
                        continue;
                    }
                    let terminal = constraints.terminal_pair_energy(k, j);
                    let candidate = if k == i {
                        branch_v + model.exterior_stem_energy(bases, k, j) + terminal
                    } else {
                        f[i][k - 1] + branch_v + model.exterior_stem_energy(bases, k, j) + terminal
                    };
                    if candidate < best {
                        best = candidate;
                        choice = FChoice::Pair(k);
                    }
                }
            }
            f[i][j] = best;
            f_choice[i][j] = choice;

            // Affine multiloop segment with one-or-more branches.
            let mut best_m1 = if m1[i][j - 1] < INF / 2.0 && constraints.allows_unpaired(j) {
                m1[i][j - 1] + model.multiloop_unpaired() + constraints.unpaired_energy(j)
            } else {
                INF
            };
            let mut best_m1_choice = if best_m1 < INF / 2.0 {
                MChoice::Unpaired
            } else {
                MChoice::Invalid
            };
            let mut best_m2 = if m2[i][j - 1] < INF / 2.0 && constraints.allows_unpaired(j) {
                m2[i][j - 1] + model.multiloop_unpaired() + constraints.unpaired_energy(j)
            } else {
                INF
            };
            let mut best_m2_choice = if best_m2 < INF / 2.0 {
                MChoice::Unpaired
            } else {
                MChoice::Invalid
            };
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branch_v = if constraints.no_lonely_pairs() {
                        v_secure[k][j]
                    } else {
                        v[k][j]
                    };
                    if branch_v >= INF / 2.0 {
                        continue;
                    }
                    let branch = model.multiloop_stem_energy(bases, k, j)
                        + constraints.terminal_pair_energy(k, j)
                        + branch_v;
                    let leading_soft = if k == i {
                        Some(0.0)
                    } else {
                        constraints.unpaired_range_energy(i, k - 1)
                    };
                    if let Some(soft) = leading_soft {
                        let first = (k - i) as f64 * model.multiloop_unpaired() + soft + branch;
                        if first < best_m1 {
                            best_m1 = first;
                            best_m1_choice = MChoice::FirstPair(k);
                        }
                    }
                    if k > i && m1[i][k - 1] < INF / 2.0 {
                        let additional = m1[i][k - 1] + branch;
                        if additional < best_m1 {
                            best_m1 = additional;
                            best_m1_choice = MChoice::AddPair(k);
                        }
                        if additional < best_m2 {
                            best_m2 = additional;
                            best_m2_choice = MChoice::AddPair(k);
                        }
                    }
                }
            }
            m1[i][j] = best_m1;
            m2[i][j] = best_m2;
            m1_choice[i][j] = best_m1_choice;
            m2_choice[i][j] = best_m2_choice;
        }
    }

    let mut pairs = Vec::new();
    traceback_f(
        0,
        n - 1,
        &f_choice,
        &v_choice,
        &m1_choice,
        &m2_choice,
        constraints.no_lonely_pairs(),
        &mut pairs,
    );
    pairs.sort_unstable();
    let structure = pairs_to_dot_bracket(n, &pairs);
    if f[0][n - 1] >= INF / 2.0 {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    }
    let evaluated =
        model.evaluate_with_constraints(&sequence, &structure, min_loop, constraints)?;

    Ok(MfeResult {
        sequence,
        structure,
        energy_kcal_mol: f[0][n - 1],
        evaluated_energy_kcal_mol: evaluated.total_kcal_mol,
        energy_difference_kcal_mol: evaluated.total_kcal_mol - f[0][n - 1],
        energy_difference_reason: None,
        evaluated_coaxial_stacks: evaluated.coaxial_stacks,
        evaluated_loop_energies: evaluated.loop_energies,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        salt_correction: model.has_salt_correction(),
        model: model.model_name(),
        constraints: constraints.summary(),
    })
}

#[derive(Clone, Copy, Debug)]
enum OddVChoice {
    Hairpin,
    Stack,
    Internal(usize, usize),
    Multiloop { left: usize, right: usize },
    CoaxLeft(usize),
    CoaxRight(usize),
    Invalid,
}

#[derive(Clone, Copy, Debug)]
enum SegmentChoice {
    Invalid,
    UnpairedLeft,
    UnpairedRight,
    Stem(usize, usize),
    Split(usize),
    Coax(usize),
}

#[derive(Clone, Copy, Debug)]
enum SplitChoice {
    Invalid,
    Split(usize),
}

fn fold_mfe_odd(
    sequence: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<MfeResult, RnaError> {
    let sequence = normalize_sequence(sequence)?;
    let bases = sequence.as_bytes();
    let n = bases.len();
    if n == 0 {
        return Err(RnaError::EmptySequence);
    }

    let mut v = vec![vec![INF; n]; n];
    let mut v_secure = vec![vec![INF; n]; n];
    let mut fm = vec![vec![INF; n]; n];
    let mut fm2 = vec![vec![INF; n]; n];
    let mut fe = vec![vec![0.0; n]; n];
    let mut v_choice = vec![vec![OddVChoice::Invalid; n]; n];
    let mut fm_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    let mut fm2_choice = vec![vec![SplitChoice::Invalid; n]; n];
    let mut fe_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    for (i, row) in fe.iter_mut().enumerate() {
        row[i] = if constraints.allows_unpaired(i) {
            constraints.unpaired_energy(i)
        } else {
            INF
        };
    }

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;

            if span > min_loop && constraints.allows_pair(bases, i, j, model) {
                let pair = model.pair_energy(bases[i], bases[j]) + constraints.pair_energy(i, j);
                let mut best = constraints
                    .unpaired_range_energy(i + 1, j - 1)
                    .map(|soft| pair + model.hairpin_energy(bases, i, j) + soft)
                    .unwrap_or(INF);
                let mut choice = if best < INF / 2.0 {
                    OddVChoice::Hairpin
                } else {
                    OddVChoice::Invalid
                };

                if span >= 2 && v[i + 1][j - 1] < INF / 2.0 {
                    let candidate = pair
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1)
                        + v[i + 1][j - 1];
                    if candidate < best {
                        best = candidate;
                        choice = OddVChoice::Stack;
                    }
                }
                let internal_limit = model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + 1 + left;
                    if k >= j {
                        break;
                    }
                    for right in 0..=(internal_limit - left) {
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
                        let candidate = pair
                            + model.internal_energy(bases, i, j, k, l)
                            + left_soft
                            + right_soft
                            + child;
                        if candidate < best {
                            best = candidate;
                            choice = OddVChoice::Internal(k, l);
                        }
                    }
                }

                for left in 0..=1usize {
                    for right in 0..=1usize {
                        if i + 1 + left > j.saturating_sub(1 + right) {
                            continue;
                        }
                        let a = i + 1 + left;
                        let b = j - 1 - right;
                        if fm2[a][b] >= INF / 2.0 {
                            continue;
                        }
                        let Some(left_soft) = constraints
                            .unpaired_range_energy(i + 1, i + left)
                            .filter(|_| left == 1)
                            .or_else(|| (left == 0).then_some(0.0))
                        else {
                            continue;
                        };
                        let Some(right_soft) = constraints
                            .unpaired_range_energy(j - right, j - 1)
                            .filter(|_| right == 1)
                            .or_else(|| (right == 0).then_some(0.0))
                        else {
                            continue;
                        };
                        let candidate = pair
                            + model.multiloop_closing()
                            + (left + right) as f64 * model.multiloop_unpaired()
                            + left_soft
                            + right_soft
                            + model.multiloop_closing_stem_selected(
                                bases,
                                i,
                                j,
                                (right == 1).then_some(bases[j - 1]),
                                (left == 1).then_some(bases[i + 1]),
                            )
                            + constraints.terminal_pair_energy(i, j)
                            + fm2[a][b];
                        if candidate < best {
                            best = candidate;
                            choice = OddVChoice::Multiloop { left, right };
                        }
                    }
                }

                if model.dangles() == 3 && span >= 5 {
                    for k in (i + 1)..=(j - 2) {
                        let left_branch = if constraints.no_lonely_pairs() {
                            v_secure[i + 1][k]
                        } else {
                            v[i + 1][k]
                        };
                        if left_branch < INF / 2.0 && fm[k + 1][j - 1] < INF / 2.0 {
                            let candidate = pair
                                + left_branch
                                + fm[k + 1][j - 1]
                                + model.multiloop_closing()
                                + model.coaxial_energy(bases, i, j, true, i + 1, k, false);
                            if candidate < best {
                                best = candidate;
                                choice = OddVChoice::CoaxLeft(k);
                            }
                        }
                        let right_branch = if constraints.no_lonely_pairs() {
                            v_secure[k + 1][j - 1]
                        } else {
                            v[k + 1][j - 1]
                        };
                        if fm[i + 1][k] < INF / 2.0 && right_branch < INF / 2.0 {
                            let candidate = pair
                                + fm[i + 1][k]
                                + right_branch
                                + model.multiloop_closing()
                                + model.coaxial_energy(bases, i, j, true, k + 1, j - 1, false);
                            if candidate < best {
                                best = candidate;
                                choice = OddVChoice::CoaxRight(k);
                            }
                        }
                    }
                }
                v[i][j] = best;
                v_choice[i][j] = choice;
                if i + 1 < j && v[i + 1][j - 1] < INF / 2.0 {
                    v_secure[i][j] = pair
                        + model.internal_energy(bases, i, j, i + 1, j - 1)
                        + constraints.stack_energy(i, j, i + 1, j - 1)
                        + v[i + 1][j - 1];
                }
            }

            let mut best_fm = INF;
            let mut best_fm_choice = SegmentChoice::Invalid;
            if fm[i + 1][j] < INF / 2.0 && constraints.allows_unpaired(i) {
                best_fm =
                    fm[i + 1][j] + model.multiloop_unpaired() + constraints.unpaired_energy(i);
                best_fm_choice = SegmentChoice::UnpairedLeft;
            }
            if fm[i][j - 1] < INF / 2.0 && constraints.allows_unpaired(j) {
                let candidate =
                    fm[i][j - 1] + model.multiloop_unpaired() + constraints.unpaired_energy(j);
                if candidate < best_fm {
                    best_fm = candidate;
                    best_fm_choice = SegmentChoice::UnpairedRight;
                }
            }
            for left in 0..=1usize {
                for right in 0..=1usize {
                    let pi = i + left;
                    let Some(pj) = j.checked_sub(right) else {
                        continue;
                    };
                    let branch_v = if constraints.no_lonely_pairs() {
                        v_secure[pi][pj]
                    } else {
                        v[pi][pj]
                    };
                    if pi >= pj || branch_v >= INF / 2.0 {
                        continue;
                    }
                    if (left == 1 && !constraints.allows_unpaired(i))
                        || (right == 1 && !constraints.allows_unpaired(j))
                    {
                        continue;
                    }
                    let candidate = branch_v
                        + (left + right) as f64 * model.multiloop_unpaired()
                        + if left == 1 {
                            constraints.unpaired_energy(i)
                        } else {
                            0.0
                        }
                        + if right == 1 {
                            constraints.unpaired_energy(j)
                        } else {
                            0.0
                        }
                        + model.multiloop_stem_selected(
                            bases,
                            pi,
                            pj,
                            (left == 1).then_some(bases[i]),
                            (right == 1).then_some(bases[j]),
                        )
                        + constraints.terminal_pair_energy(pi, pj);
                    if candidate < best_fm {
                        best_fm = candidate;
                        best_fm_choice = SegmentChoice::Stem(pi, pj);
                    }
                }
            }
            for k in i..j {
                if fm[i][k] < INF / 2.0 && fm[k + 1][j] < INF / 2.0 {
                    let candidate = fm[i][k] + fm[k + 1][j];
                    if candidate < best_fm {
                        best_fm = candidate;
                        best_fm_choice = SegmentChoice::Split(k);
                    }
                }
                let left_branch = if constraints.no_lonely_pairs() {
                    v_secure[i][k]
                } else {
                    v[i][k]
                };
                let right_branch = if constraints.no_lonely_pairs() {
                    v_secure[k + 1][j]
                } else {
                    v[k + 1][j]
                };
                if model.dangles() == 3 && left_branch < INF / 2.0 && right_branch < INF / 2.0 {
                    let candidate = left_branch
                        + right_branch
                        + model.coaxial_energy(bases, i, k, false, k + 1, j, false);
                    if candidate < best_fm {
                        best_fm = candidate;
                        best_fm_choice = SegmentChoice::Coax(k);
                    }
                }
            }
            fm[i][j] = best_fm;
            fm_choice[i][j] = best_fm_choice;

            let mut best_two = INF;
            let mut best_two_choice = SplitChoice::Invalid;
            for k in i..j {
                if fm[i][k] < INF / 2.0 && fm[k + 1][j] < INF / 2.0 {
                    let candidate = fm[i][k] + fm[k + 1][j];
                    if candidate < best_two {
                        best_two = candidate;
                        best_two_choice = SplitChoice::Split(k);
                    }
                }
            }
            fm2[i][j] = best_two;
            fm2_choice[i][j] = best_two_choice;

            let mut best_fe = if constraints.allows_unpaired(i) {
                interval(&fe, i + 1, j) + constraints.unpaired_energy(i)
            } else {
                INF
            };
            let mut best_fe_choice = if best_fe < INF / 2.0 {
                SegmentChoice::UnpairedLeft
            } else {
                SegmentChoice::Invalid
            };
            let right_unpaired = if constraints.allows_unpaired(j) {
                interval(&fe, i, j - 1) + constraints.unpaired_energy(j)
            } else {
                INF
            };
            if right_unpaired < best_fe {
                best_fe = right_unpaired;
                best_fe_choice = SegmentChoice::UnpairedRight;
            }
            for left in 0..=1usize {
                for right in 0..=1usize {
                    let pi = i + left;
                    let Some(pj) = j.checked_sub(right) else {
                        continue;
                    };
                    let branch_v = if constraints.no_lonely_pairs() {
                        v_secure[pi][pj]
                    } else {
                        v[pi][pj]
                    };
                    if pi >= pj || branch_v >= INF / 2.0 {
                        continue;
                    }
                    if (left == 1 && !constraints.allows_unpaired(i))
                        || (right == 1 && !constraints.allows_unpaired(j))
                    {
                        continue;
                    }
                    let candidate = branch_v
                        + if left == 1 {
                            constraints.unpaired_energy(i)
                        } else {
                            0.0
                        }
                        + if right == 1 {
                            constraints.unpaired_energy(j)
                        } else {
                            0.0
                        }
                        + model.exterior_stem_selected(
                            bases,
                            pi,
                            pj,
                            (left == 1).then_some(bases[i]),
                            (right == 1).then_some(bases[j]),
                        )
                        + constraints.terminal_pair_energy(pi, pj);
                    if candidate < best_fe {
                        best_fe = candidate;
                        best_fe_choice = SegmentChoice::Stem(pi, pj);
                    }
                }
            }
            for k in i..j {
                let candidate = fe[i][k] + fe[k + 1][j];
                if candidate < best_fe {
                    best_fe = candidate;
                    best_fe_choice = SegmentChoice::Split(k);
                }
            }
            fe[i][j] = best_fe;
            fe_choice[i][j] = best_fe_choice;
        }
    }

    let mut pairs = Vec::new();
    traceback_odd_fe(
        0,
        n - 1,
        &fe_choice,
        &v_choice,
        &fm_choice,
        &fm2_choice,
        constraints.no_lonely_pairs(),
        &mut pairs,
    );
    pairs.sort_unstable();
    pairs.dedup();
    let structure = pairs_to_dot_bracket(n, &pairs);
    if fe[0][n - 1] >= INF / 2.0 {
        return Err(RnaError::InvalidOption(
            "constraints admit no valid secondary structure".into(),
        ));
    }
    let evaluated =
        model.evaluate_with_constraints(&sequence, &structure, min_loop, constraints)?;
    let energy_difference = evaluated.total_kcal_mol - fe[0][n - 1];
    let reason = if model.dangles() == 3 && energy_difference.abs() > 1.0e-9 {
        Some("dangles=3 uses distinct MFE and fixed-structure multiloop evaluation grammars")
    } else {
        None
    };
    Ok(MfeResult {
        sequence,
        structure,
        energy_kcal_mol: fe[0][n - 1],
        evaluated_energy_kcal_mol: evaluated.total_kcal_mol,
        energy_difference_kcal_mol: energy_difference,
        energy_difference_reason: reason,
        evaluated_coaxial_stacks: evaluated.coaxial_stacks,
        evaluated_loop_energies: evaluated.loop_energies,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        salt_correction: model.has_salt_correction(),
        model: model.model_name(),
        constraints: constraints.summary(),
    })
}

#[allow(clippy::too_many_arguments)]
fn traceback_odd_fe(
    i: usize,
    j: usize,
    fe_choice: &[Vec<SegmentChoice>],
    v_choice: &[Vec<OddVChoice>],
    fm_choice: &[Vec<SegmentChoice>],
    fm2_choice: &[Vec<SplitChoice>],
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    OddTrace {
        fe_choice,
        v_choice,
        fm_choice,
        fm2_choice,
        no_lonely_pairs,
    }
    .fe(i, j, pairs);
}

struct OddTrace<'a> {
    fe_choice: &'a [Vec<SegmentChoice>],
    v_choice: &'a [Vec<OddVChoice>],
    fm_choice: &'a [Vec<SegmentChoice>],
    fm2_choice: &'a [Vec<SplitChoice>],
    no_lonely_pairs: bool,
}

impl OddTrace<'_> {
    fn fe(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
        if i >= j || j >= self.fe_choice.len() {
            return;
        }
        match self.fe_choice[i][j] {
            SegmentChoice::Invalid => {}
            SegmentChoice::UnpairedLeft => self.fe(i + 1, j, pairs),
            SegmentChoice::UnpairedRight => self.fe(i, j - 1, pairs),
            SegmentChoice::Stem(k, l) => self.v_branch(k, l, pairs),
            SegmentChoice::Split(k) => {
                self.fe(i, k, pairs);
                self.fe(k + 1, j, pairs);
            }
            SegmentChoice::Coax(_) => unreachable!("exterior coaxial choice is not generated"),
        }
    }

    fn v_branch(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
        if self.no_lonely_pairs {
            pairs.push((i, j));
            self.v(i + 1, j - 1, pairs);
        } else {
            self.v(i, j, pairs);
        }
    }

    fn v(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
        if i >= j {
            return;
        }
        pairs.push((i, j));
        match self.v_choice[i][j] {
            OddVChoice::Hairpin | OddVChoice::Invalid => {}
            OddVChoice::Stack => self.v(i + 1, j - 1, pairs),
            OddVChoice::Internal(k, l) => self.v_branch(k, l, pairs),
            OddVChoice::Multiloop { left, right } => self.fm2(i + 1 + left, j - 1 - right, pairs),
            OddVChoice::CoaxLeft(k) => {
                self.v_branch(i + 1, k, pairs);
                self.fm(k + 1, j - 1, pairs);
            }
            OddVChoice::CoaxRight(k) => {
                self.fm(i + 1, k, pairs);
                self.v_branch(k + 1, j - 1, pairs);
            }
        }
    }

    fn fm(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
        if i >= j || j >= self.fm_choice.len() {
            return;
        }
        match self.fm_choice[i][j] {
            SegmentChoice::Invalid => {}
            SegmentChoice::UnpairedLeft => self.fm(i + 1, j, pairs),
            SegmentChoice::UnpairedRight => self.fm(i, j - 1, pairs),
            SegmentChoice::Stem(k, l) => self.v_branch(k, l, pairs),
            SegmentChoice::Split(k) => {
                self.fm(i, k, pairs);
                self.fm(k + 1, j, pairs);
            }
            SegmentChoice::Coax(k) => {
                self.v_branch(i, k, pairs);
                self.v_branch(k + 1, j, pairs);
            }
        }
    }

    fn fm2(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
        if i >= j || j >= self.fm2_choice.len() {
            return;
        }
        if let SplitChoice::Split(k) = self.fm2_choice[i][j] {
            self.fm(i, k, pairs);
            self.fm(k + 1, j, pairs);
        }
    }
}

fn interval(table: &[Vec<f64>], i: usize, j: usize) -> f64 {
    if i > j {
        0.0
    } else {
        table[i][j]
    }
}

#[allow(clippy::too_many_arguments)]
fn traceback_f(
    i: usize,
    j: usize,
    f_choice: &[Vec<FChoice>],
    v_choice: &[Vec<VChoice>],
    m1_choice: &[Vec<MChoice>],
    m2_choice: &[Vec<MChoice>],
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= f_choice.len() || j >= f_choice.len() || i >= j {
        return;
    }
    match f_choice[i][j] {
        FChoice::Empty => {}
        FChoice::Unpaired => traceback_f(
            i,
            j - 1,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
        FChoice::Pair(k) => {
            if k > i {
                traceback_f(
                    i,
                    k - 1,
                    f_choice,
                    v_choice,
                    m1_choice,
                    m2_choice,
                    no_lonely_pairs,
                    pairs,
                );
            }
            traceback_v_branch(
                k,
                j,
                f_choice,
                v_choice,
                m1_choice,
                m2_choice,
                no_lonely_pairs,
                pairs,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn traceback_v_branch(
    i: usize,
    j: usize,
    f_choice: &[Vec<FChoice>],
    v_choice: &[Vec<VChoice>],
    m1_choice: &[Vec<MChoice>],
    m2_choice: &[Vec<MChoice>],
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    if no_lonely_pairs {
        pairs.push((i, j));
        traceback_v(
            i + 1,
            j - 1,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        );
    } else {
        traceback_v(
            i,
            j,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn traceback_v(
    i: usize,
    j: usize,
    f_choice: &[Vec<FChoice>],
    v_choice: &[Vec<VChoice>],
    m1_choice: &[Vec<MChoice>],
    m2_choice: &[Vec<MChoice>],
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= j {
        return;
    }
    pairs.push((i, j));
    match v_choice[i][j] {
        VChoice::Hairpin | VChoice::Invalid => {}
        VChoice::Stack => traceback_v(
            i + 1,
            j - 1,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
        VChoice::Internal(k, l) => traceback_v_branch(
            k,
            l,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
        VChoice::Multiloop => traceback_m(
            i + 1,
            j - 1,
            true,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn traceback_m(
    i: usize,
    j: usize,
    at_least_two: bool,
    f_choice: &[Vec<FChoice>],
    v_choice: &[Vec<VChoice>],
    m1_choice: &[Vec<MChoice>],
    m2_choice: &[Vec<MChoice>],
    no_lonely_pairs: bool,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= j || j >= m1_choice.len() {
        return;
    }
    let choice = if at_least_two {
        m2_choice[i][j]
    } else {
        m1_choice[i][j]
    };
    match choice {
        MChoice::Invalid => {}
        MChoice::Unpaired => traceback_m(
            i,
            j - 1,
            at_least_two,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
        MChoice::FirstPair(k) => traceback_v_branch(
            k,
            j,
            f_choice,
            v_choice,
            m1_choice,
            m2_choice,
            no_lonely_pairs,
            pairs,
        ),
        MChoice::AddPair(k) => {
            traceback_m(
                i,
                k - 1,
                false,
                f_choice,
                v_choice,
                m1_choice,
                m2_choice,
                no_lonely_pairs,
                pairs,
            );
            traceback_v_branch(
                k,
                j,
                f_choice,
                v_choice,
                m1_choice,
                m2_choice,
                no_lonely_pairs,
                pairs,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::{ConstraintConfig, PairConstraint};

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
    fn predicts_a_stable_gc_hairpin() {
        let result = fold_mfe("GGGAAACCC", 3, &EnergyModel::default()).unwrap();
        assert_eq!(result.structure, "(((...)))");
        assert!(result.energy_kcal_mol < 0.0);
    }

    #[test]
    fn ambiguity_codes_remain_unpaired() {
        let result = fold_mfe("NNNNNN", 3, &EnergyModel::default()).unwrap();
        assert_eq!(result.structure, "......");
    }

    #[test]
    fn dynamic_program_matches_exhaustive_public_energy_for_short_sequences() {
        let model = EnergyModel::default();
        for sequence in ["GGGAAACCC", "GCGCUUCGCC"] {
            let predicted = fold_mfe(sequence, 3, &model).unwrap();
            let exhaustive =
                enumerate_matchings(sequence.as_bytes(), &model, 3, 0, sequence.len() - 1)
                    .into_iter()
                    .map(|pairs| {
                        let structure = pairs_to_dot_bracket(sequence.len(), &pairs);
                        model.evaluate(sequence, &structure).unwrap().total_kcal_mol
                    })
                    .fold(f64::INFINITY, f64::min);
            assert!((predicted.energy_kcal_mol - exhaustive).abs() < 1.0e-9);
        }
    }

    #[test]
    fn default_model_enumerates_internal_loops_beyond_thirty_nucleotides() {
        let mut bases = vec![b'A'; 44];
        bases[0] = b'G';
        bases[2] = b'G';
        bases[7] = b'C';
        bases[43] = b'C';
        let sequence = String::from_utf8(bases).unwrap();
        let config = ConstraintConfig {
            force_pairs: vec![
                PairConstraint { i: 1, j: 44 },
                PairConstraint { i: 3, j: 8 },
            ],
            force_unpaired: (1..=44)
                .filter(|position| !matches!(position, 1 | 3 | 8 | 44))
                .collect(),
            ..ConstraintConfig::default()
        };
        let constraints = ConstraintModel::compile(44, &config).unwrap();
        let complete =
            fold_mfe_with_constraints(&sequence, 3, &EnergyModel::default(), &constraints).unwrap();
        assert_eq!(
            complete.structure,
            "(.(....)...................................)"
        );

        let mut caller_limited = EnergyModel::default();
        caller_limited.max_internal_loop = Some(30);
        assert!(fold_mfe_with_constraints(&sequence, 3, &caller_limited, &constraints).is_err());
    }

    #[test]
    fn odd_dangle_models_fold_reference_sequences() {
        let common = [
            (
                "ACGAAUGGCCAUCUAGUGAUACUAUAAUCG",
                ".(((((((..(((....))).))))..)))",
                -2.7,
            ),
            (
                "GGGGGUGUUUCACUCGUGGCAGGAAGCAGC",
                "...(.((((((.((......)))))))).)",
                -6.2,
            ),
            (
                "CUCGAUUUGCGCACUCAUACCAGUCUGUCU",
                "........(((.(((......))).)))..",
                -1.0,
            ),
            (
                "CGAGUGAGACGCAUAGGGUAAUUGGCGAGG",
                ".........(((..(......)..)))...",
                -2.5,
            ),
            (
                "AUAAUUACGGCUAGCGACGGUCACCGCUGG",
                "..........((((((........))))))",
                -7.2,
            ),
        ];
        for (sequence, structure, reference) in common {
            for dangles in [1, 3] {
                let model = EnergyModel::with_dangles(37.0, dangles).unwrap();
                let predicted = fold_mfe(sequence, 3, &model).unwrap();
                assert_eq!(predicted.structure, structure);
                assert!((predicted.energy_kcal_mol - reference).abs() < 1.0e-9);
            }
        }
        let sequence = "GUCCCGGCCUCGAGACCUAUCCGGUUCGUCACGGAGCGCAGCCCGUGACGCGGGGUGACU";
        let single = fold_mfe(sequence, 3, &EnergyModel::with_dangles(37.0, 1).unwrap()).unwrap();
        assert_eq!(
            single.structure,
            "(((.((....)).))).(((((.(..((((((((.((...))))))))))).)))))..."
        );
        assert!((single.energy_kcal_mol + 24.0).abs() < 1.0e-9);

        // RNAstructure 6.6 `Fold --MFE --isolated` reports this structure and
        // -26.6 kcal/mol from the same official parameter bundle.
        let coaxial = fold_mfe(sequence, 3, &EnergyModel::with_dangles(37.0, 3).unwrap()).unwrap();
        assert_eq!(
            coaxial.structure,
            "(((...((((((.((((.....))))((((((((.((...)))))))))))))))))))."
        );
        assert!((coaxial.energy_kcal_mol + 26.6).abs() < 1.0e-9);
        assert!((coaxial.evaluated_energy_kcal_mol + 27.1).abs() < 1.0e-9);
        assert!((coaxial.energy_difference_kcal_mol + 0.5).abs() < 1.0e-9);
        assert!(coaxial.energy_difference_reason.is_some());
        assert_eq!(coaxial.evaluated_coaxial_stacks.len(), 1);
        let stack = &coaxial.evaluated_coaxial_stacks[0];
        assert_eq!((stack.first_i, stack.first_j), (27, 50));
        assert_eq!((stack.second_i, stack.second_j), (12, 51));
    }
}
