use crate::constraints::ConstraintModel;
use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{pairs_to_dot_bracket, RnaError};

#[derive(Clone, Copy, Debug)]
enum Choice {
    Empty,
    Unpaired,
    Pair(usize),
}

fn probability_matrix(length: usize, pairs: &[PairProbability]) -> Vec<Vec<f64>> {
    let mut matrix = vec![vec![0.0; length]; length];
    for pair in pairs {
        if pair.i > 0 && pair.j > 0 && pair.i <= length && pair.j <= length {
            matrix[pair.i - 1][pair.j - 1] = pair.probability;
        }
    }
    matrix
}

fn decode_with_scores(
    length: usize,
    min_loop: usize,
    pair_scores: &[Vec<f64>],
    unpaired_scores: &[f64],
) -> (String, f64) {
    if length == 0 {
        return (String::new(), 0.0);
    }
    let mut dp = vec![vec![0.0f64; length]; length];
    let mut choices = vec![vec![Choice::Empty; length]; length];
    for i in 0..length {
        dp[i][i] = unpaired_scores[i];
        choices[i][i] = Choice::Unpaired;
    }
    for span in 1..length {
        for i in 0..(length - span) {
            let j = i + span;
            let mut best = dp[i][j - 1] + unpaired_scores[j];
            let mut choice = Choice::Unpaired;
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let left = if k == i { 0.0 } else { dp[i][k - 1] };
                    let inner = if k + 1 >= j { 0.0 } else { dp[k + 1][j - 1] };
                    let candidate = left + inner + pair_scores[k][j];
                    if candidate > best + 1.0e-12 {
                        best = candidate;
                        choice = Choice::Pair(k);
                    }
                }
            }
            dp[i][j] = best;
            choices[i][j] = choice;
        }
    }

    fn trace(i: usize, j: usize, choices: &[Vec<Choice>], pairs: &mut Vec<(usize, usize)>) {
        if i > j || j >= choices.len() {
            return;
        }
        if i == j {
            return;
        }
        match choices[i][j] {
            Choice::Empty => {}
            Choice::Unpaired => trace(i, j - 1, choices, pairs),
            Choice::Pair(k) => {
                if k > i {
                    trace(i, k - 1, choices, pairs);
                }
                if k + 1 < j {
                    trace(k + 1, j - 1, choices, pairs);
                }
                pairs.push((k, j));
            }
        }
    }

    let mut pairs = Vec::new();
    trace(0, length - 1, &choices, &mut pairs);
    (pairs_to_dot_bracket(length, &pairs), dp[0][length - 1])
}

pub fn decode_centroid(length: usize, min_loop: usize, pairs: &[PairProbability]) -> (String, f64) {
    let probabilities = probability_matrix(length, pairs);
    let mut scores = vec![vec![0.0; length]; length];
    for i in 0..length {
        for j in (i + 1)..length {
            scores[i][j] = 2.0 * probabilities[i][j] - 1.0;
        }
    }
    decode_with_scores(length, min_loop, &scores, &vec![0.0; length])
}

pub fn decode_mea(
    length: usize,
    min_loop: usize,
    gamma: f64,
    pairs: &[PairProbability],
    unpaired: &[f64],
) -> (String, f64) {
    let probabilities = probability_matrix(length, pairs);
    let mut scores = vec![vec![0.0; length]; length];
    for i in 0..length {
        for j in (i + 1)..length {
            scores[i][j] = 2.0 * gamma * probabilities[i][j];
        }
    }
    decode_with_scores(length, min_loop, &scores, unpaired)
}

/// Constraint-aware centroid decoding. Unlike post-filtering, forced-paired
/// positions and the no-lonely-pair rule participate in the decoding DP.
pub fn decode_centroid_with_constraints(
    bases: &[u8],
    min_loop: usize,
    pairs: &[PairProbability],
    constraints: &ConstraintModel,
    model: &EnergyModel,
) -> Result<(String, f64), RnaError> {
    let probabilities = probability_matrix(bases.len(), pairs);
    let mut scores = vec![vec![0.0; bases.len()]; bases.len()];
    for i in 0..bases.len() {
        for j in (i + 1)..bases.len() {
            scores[i][j] = 2.0 * probabilities[i][j] - 1.0;
        }
    }
    decode_constrained(
        bases,
        min_loop,
        &scores,
        &vec![0.0; bases.len()],
        constraints,
        model,
    )
}

/// Constraint-aware maximum expected accuracy decoding.
pub fn decode_mea_with_constraints(
    bases: &[u8],
    min_loop: usize,
    gamma: f64,
    pairs: &[PairProbability],
    unpaired: &[f64],
    constraints: &ConstraintModel,
    model: &EnergyModel,
) -> Result<(String, f64), RnaError> {
    let probabilities = probability_matrix(bases.len(), pairs);
    let mut scores = vec![vec![0.0; bases.len()]; bases.len()];
    for i in 0..bases.len() {
        for j in (i + 1)..bases.len() {
            scores[i][j] = 2.0 * gamma * probabilities[i][j];
        }
    }
    decode_constrained(bases, min_loop, &scores, unpaired, constraints, model)
}

#[derive(Clone, Copy)]
enum ConstrainedVChoice {
    Invalid,
    Inside,
    Stack,
}

fn decode_constrained(
    bases: &[u8],
    min_loop: usize,
    pair_scores: &[Vec<f64>],
    unpaired_scores: &[f64],
    constraints: &ConstraintModel,
    model: &EnergyModel,
) -> Result<(String, f64), RnaError> {
    const SCORE_NEG_INF: f64 = f64::NEG_INFINITY;
    let length = bases.len();
    let mut dp = vec![vec![SCORE_NEG_INF; length]; length];
    let mut paired = vec![vec![SCORE_NEG_INF; length]; length];
    let mut secure = vec![vec![SCORE_NEG_INF; length]; length];
    let mut choices = vec![vec![Choice::Empty; length]; length];
    let mut paired_choices = vec![vec![ConstrainedVChoice::Invalid; length]; length];

    for i in 0..length {
        if constraints.allows_unpaired(i) {
            dp[i][i] = unpaired_scores[i];
            choices[i][i] = Choice::Unpaired;
        }
    }
    let empty_or = |table: &[Vec<f64>], i: usize, j: usize| {
        if i > j {
            0.0
        } else {
            table[i][j]
        }
    };

    for span in 1..length {
        for i in 0..(length - span) {
            let j = i + span;
            if span > min_loop && constraints.allows_pair(bases, i, j, model) {
                let inside = empty_or(&dp, i + 1, j - 1);
                if inside != SCORE_NEG_INF {
                    paired[i][j] = pair_scores[i][j] + inside;
                    paired_choices[i][j] = ConstrainedVChoice::Inside;
                }
                if i + 1 < j && paired[i + 1][j - 1] != SCORE_NEG_INF {
                    let stacked = pair_scores[i][j] + paired[i + 1][j - 1];
                    secure[i][j] = stacked;
                    if stacked > paired[i][j] {
                        paired[i][j] = stacked;
                        paired_choices[i][j] = ConstrainedVChoice::Stack;
                    }
                }
            }

            if constraints.allows_unpaired(j) && dp[i][j - 1] != SCORE_NEG_INF {
                dp[i][j] = dp[i][j - 1] + unpaired_scores[j];
                choices[i][j] = Choice::Unpaired;
            }
            if j > min_loop {
                for k in i..=(j - min_loop - 1) {
                    let branch = if constraints.no_lonely_pairs() {
                        secure[k][j]
                    } else {
                        paired[k][j]
                    };
                    if branch == SCORE_NEG_INF {
                        continue;
                    }
                    let left = if k == i { 0.0 } else { dp[i][k - 1] };
                    if left == SCORE_NEG_INF {
                        continue;
                    }
                    let candidate = left + branch;
                    if candidate > dp[i][j] {
                        dp[i][j] = candidate;
                        choices[i][j] = Choice::Pair(k);
                    }
                }
            }
        }
    }
    if dp[0][length - 1] == SCORE_NEG_INF {
        return Err(RnaError::InvalidOption(
            "constraints admit no centroid/MEA structure".into(),
        ));
    }

    struct Trace<'a> {
        choices: &'a [Vec<Choice>],
        paired_choices: &'a [Vec<ConstrainedVChoice>],
        no_lonely_pairs: bool,
    }
    impl Trace<'_> {
        fn interval(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
            if i > j || j >= self.choices.len() || i == j {
                return;
            }
            match self.choices[i][j] {
                Choice::Empty => {}
                Choice::Unpaired => self.interval(i, j - 1, pairs),
                Choice::Pair(k) => {
                    if k > i {
                        self.interval(i, k - 1, pairs);
                    }
                    self.branch(k, j, pairs);
                }
            }
        }
        fn branch(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
            if self.no_lonely_pairs {
                pairs.push((i, j));
                self.pair(i + 1, j - 1, pairs);
            } else {
                self.pair(i, j, pairs);
            }
        }
        fn pair(&self, i: usize, j: usize, pairs: &mut Vec<(usize, usize)>) {
            pairs.push((i, j));
            match self.paired_choices[i][j] {
                ConstrainedVChoice::Invalid => {}
                ConstrainedVChoice::Inside => {
                    if i + 1 < j {
                        self.interval(i + 1, j - 1, pairs);
                    }
                }
                ConstrainedVChoice::Stack => self.pair(i + 1, j - 1, pairs),
            }
        }
    }
    let mut result_pairs = Vec::new();
    Trace {
        choices: &choices,
        paired_choices: &paired_choices,
        no_lonely_pairs: constraints.no_lonely_pairs(),
    }
    .interval(0, length - 1, &mut result_pairs);
    Ok((
        pairs_to_dot_bracket(length, &result_pairs),
        dp[0][length - 1],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centroid_selects_pairs_above_one_half_without_crossing() {
        let pairs = vec![
            PairProbability {
                i: 1,
                j: 8,
                probability: 0.8,
            },
            PairProbability {
                i: 2,
                j: 7,
                probability: 0.7,
            },
            PairProbability {
                i: 3,
                j: 6,
                probability: 0.4,
            },
        ];
        let (structure, _) = decode_centroid(8, 2, &pairs);
        assert_eq!(structure, "((....))");
    }
}
