//! Exhaustive, unambiguous enumeration of pseudoknot-free structures.
//!
//! This is the correctness oracle for ensemble grammars whose local state
//! factorization is more involved than the ordinary McCaskill recurrence. It
//! deliberately has no sequence-length or structure-count guard: every finite
//! input is enumerated exactly, with resource exhaustion reported by the host
//! rather than silently changing the conformation space.

use crate::constraints::ConstraintModel;
use crate::energy::EnergyModel;
use crate::structure::RnaError;

pub(crate) fn for_each_noncrossing_structure(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    mut visit: impl FnMut(&[(usize, usize)]) -> Result<(), RnaError>,
) -> Result<(), RnaError> {
    let mut intervals = if bases.is_empty() {
        Vec::new()
    } else {
        vec![(0, bases.len() - 1)]
    };
    let mut pairs = Vec::new();
    enumerate(
        bases,
        min_loop,
        model,
        constraints,
        &mut intervals,
        &mut pairs,
        &mut visit,
    )
}

#[allow(clippy::too_many_arguments)]
fn enumerate(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    intervals: &mut Vec<(usize, usize)>,
    pairs: &mut Vec<(usize, usize)>,
    visit: &mut impl FnMut(&[(usize, usize)]) -> Result<(), RnaError>,
) -> Result<(), RnaError> {
    let Some((i, j)) = intervals.pop() else {
        return visit(pairs);
    };

    if i > j {
        let result = enumerate(bases, min_loop, model, constraints, intervals, pairs, visit);
        intervals.push((i, j));
        return result;
    }

    if constraints.allows_unpaired(i) {
        if i < j {
            intervals.push((i + 1, j));
        }
        enumerate(bases, min_loop, model, constraints, intervals, pairs, visit)?;
        if i < j {
            intervals.pop();
        }
    }

    if i + min_loop < j {
        for k in i + min_loop + 1..=j {
            if !constraints.allows_pair(bases, i, k, model) {
                continue;
            }
            pairs.push((i, k));
            let mut pushed = 0usize;
            if k < j {
                intervals.push((k + 1, j));
                pushed += 1;
            }
            if i + 1 < k {
                intervals.push((i + 1, k - 1));
                pushed += 1;
            }
            enumerate(bases, min_loop, model, constraints, intervals, pairs, visit)?;
            for _ in 0..pushed {
                intervals.pop();
            }
            pairs.pop();
        }
    }

    intervals.push((i, j));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::ConstraintConfig;
    use std::collections::HashSet;

    #[test]
    fn enumeration_is_unique_and_complete_for_a_short_sequence() {
        let bases = b"GCAAAUGC";
        let model = EnergyModel::with_dangles(37.0, 1).unwrap();
        let constraints =
            ConstraintModel::compile(bases.len(), &ConstraintConfig::default()).unwrap();
        let mut structures = HashSet::new();
        for_each_noncrossing_structure(bases, 3, &model, &constraints, |pairs| {
            assert!(structures.insert(pairs.to_vec()));
            Ok(())
        })
        .unwrap();
        assert!(structures.contains(&Vec::new()));
        assert!(structures.iter().any(|pairs| !pairs.is_empty()));
    }
}
