//! Exact secondary-structure energy-landscape paths.
//!
//! The complete pseudoknot-free state space admitted by the sequence, model,
//! minimum-loop rule, and constraints is enumerated. The minimax path is then
//! solved on the graph of single-base-pair insertions and deletions. No beam,
//! indirect-path width, or state-count limit is used.

use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::energy::EnergyModel;
use crate::exact_enumeration::for_each_noncrossing_structure;
use crate::structure::{
    is_pseudoknotted, normalize_sequence, pairs_to_dot_bracket, parse_structure, RnaError,
};
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

#[derive(Clone, Debug, Serialize)]
pub struct LandscapeState {
    pub structure: String,
    pub energy_kcal_mol: f64,
    pub degree: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct LandscapePathStep {
    pub step: usize,
    pub structure: String,
    pub energy_kcal_mol: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LandscapeResult {
    pub sequence: String,
    pub start_structure: String,
    pub target_structure: String,
    pub state_count: usize,
    pub edge_count: usize,
    pub local_minima: Vec<LandscapeState>,
    pub path: Vec<LandscapePathStep>,
    pub saddle_energy_kcal_mol: f64,
    pub start_barrier_kcal_mol: f64,
    pub target_barrier_kcal_mol: f64,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub constraints: ConstraintSummary,
    pub state_space_complete: bool,
    pub move_set: &'static str,
    pub algorithm: &'static str,
    pub time_complexity: &'static str,
    pub space_complexity: &'static str,
}

#[derive(Clone)]
struct State {
    structure: String,
    energy: f64,
    pairs: Vec<(usize, usize)>,
    partner: Vec<Option<usize>>,
}

#[derive(Clone, Copy)]
struct QueueEntry {
    saddle: f64,
    steps: usize,
    index: usize,
}

impl PartialEq for QueueEntry {
    fn eq(&self, other: &Self) -> bool {
        self.saddle.to_bits() == other.saddle.to_bits()
            && self.steps == other.steps
            && self.index == other.index
    }
}

impl Eq for QueueEntry {}

impl PartialOrd for QueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for QueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap is a max heap, so reverse every deterministic key.
        other
            .saddle
            .total_cmp(&self.saddle)
            .then_with(|| other.steps.cmp(&self.steps))
            .then_with(|| other.index.cmp(&self.index))
    }
}

pub fn exact_landscape_path(
    sequence: &str,
    start_structure: &str,
    target_structure: &str,
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
) -> Result<LandscapeResult, RnaError> {
    let sequence = normalize_sequence(sequence)?;
    let start = canonical_planar_structure(&sequence, start_structure)?;
    let target = canonical_planar_structure(&sequence, target_structure)?;

    // Validate canonical pairing, minimum-loop, noLP, hard constraints, and
    // all structure-dependent soft terms before constructing the graph.
    model.evaluate_with_constraints(&sequence, &start, min_loop, constraints)?;
    model.evaluate_with_constraints(&sequence, &target, min_loop, constraints)?;

    let bases = sequence.as_bytes();
    let mut states = Vec::new();
    for_each_noncrossing_structure(bases, min_loop, model, constraints, |pairs| {
        let mut sorted_pairs = pairs.to_vec();
        sorted_pairs.sort_unstable();
        let structure = pairs_to_dot_bracket(bases.len(), &sorted_pairs);
        let energy = model
            .evaluate_with_constraints(&sequence, &structure, min_loop, constraints)?
            .total_kcal_mol;
        let mut partner = vec![None; bases.len()];
        for &(i, j) in &sorted_pairs {
            partner[i] = Some(j);
            partner[j] = Some(i);
        }
        states.push(State {
            structure,
            energy,
            pairs: sorted_pairs,
            partner,
        });
        Ok(())
    })?;
    states.sort_by(|left, right| left.structure.cmp(&right.structure));

    let index: HashMap<String, usize> = states
        .iter()
        .enumerate()
        .map(|(position, state)| (state.structure.clone(), position))
        .collect();
    let start_index = *index.get(&start).ok_or_else(|| {
        RnaError::InvalidOption("start structure is outside the enumerated state space".into())
    })?;
    let target_index = *index.get(&target).ok_or_else(|| {
        RnaError::InvalidOption("target structure is outside the enumerated state space".into())
    })?;

    let adjacency = build_adjacency(bases, min_loop, model, constraints, &states, &index);
    let edge_count = adjacency.iter().map(Vec::len).sum::<usize>() / 2;
    let (saddle, parents) = minimax_path(&states, &adjacency, start_index, target_index)?;
    let path = recover_path(&parents, start_index, target_index)?
        .into_iter()
        .enumerate()
        .map(|(step, state_index)| LandscapePathStep {
            step,
            structure: states[state_index].structure.clone(),
            energy_kcal_mol: states[state_index].energy,
        })
        .collect();

    let mut local_minima = states
        .iter()
        .enumerate()
        .filter(|(state_index, state)| {
            adjacency[*state_index]
                .iter()
                .all(|&neighbor| states[neighbor].energy + 1.0e-12 >= state.energy)
        })
        .map(|(state_index, state)| LandscapeState {
            structure: state.structure.clone(),
            energy_kcal_mol: state.energy,
            degree: adjacency[state_index].len(),
        })
        .collect::<Vec<_>>();
    local_minima.sort_by(|left, right| {
        left.energy_kcal_mol
            .total_cmp(&right.energy_kcal_mol)
            .then_with(|| left.structure.cmp(&right.structure))
    });

    Ok(LandscapeResult {
        sequence,
        start_structure: start,
        target_structure: target,
        state_count: states.len(),
        edge_count,
        local_minima,
        path,
        saddle_energy_kcal_mol: saddle,
        start_barrier_kcal_mol: saddle - states[start_index].energy,
        target_barrier_kcal_mol: saddle - states[target_index].energy,
        temperature_celsius: model.temperature_celsius(),
        dangles: model.dangles(),
        salt_molar: model.salt_molar(),
        constraints: constraints.summary(),
        state_space_complete: true,
        move_set: "single base-pair insertion or deletion",
        algorithm: "complete planar state enumeration plus exact minimax Dijkstra",
        time_complexity: "exponential state space; polynomial in the explicit landscape graph",
        space_complexity: "exponential state space",
    })
}

fn canonical_planar_structure(sequence: &str, structure: &str) -> Result<String, RnaError> {
    let parsed = parse_structure(sequence, structure)?;
    if is_pseudoknotted(&parsed.pairs) {
        return Err(RnaError::PseudoknotUnsupported("exact landscape path"));
    }
    if parsed.pairs.iter().any(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(
            "landscape endpoints require canonical base pairs".into(),
        ));
    }
    let pairs = parsed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    Ok(pairs_to_dot_bracket(parsed.length, &pairs))
}

fn build_adjacency(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    states: &[State],
    index: &HashMap<String, usize>,
) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); states.len()];
    for (state_index, state) in states.iter().enumerate() {
        let mut neighbors = HashSet::new();
        for &(i, j) in &state.pairs {
            let mut candidate = state.structure.as_bytes().to_vec();
            candidate[i] = b'.';
            candidate[j] = b'.';
            if let Ok(candidate) = std::str::from_utf8(&candidate) {
                if let Some(&neighbor) = index.get(candidate) {
                    neighbors.insert(neighbor);
                }
            }
        }
        for i in 0..bases.len() {
            if state.partner[i].is_some() {
                continue;
            }
            for j in i + min_loop + 1..bases.len() {
                if state.partner[j].is_some()
                    || !constraints.allows_pair(bases, i, j, model)
                    || state
                        .pairs
                        .iter()
                        .any(|&(k, l)| (k < i && i < l && l < j) || (i < k && k < j && j < l))
                {
                    continue;
                }
                let mut candidate = state.structure.as_bytes().to_vec();
                candidate[i] = b'(';
                candidate[j] = b')';
                if let Ok(candidate) = std::str::from_utf8(&candidate) {
                    if let Some(&neighbor) = index.get(candidate) {
                        neighbors.insert(neighbor);
                    }
                }
            }
        }
        let mut neighbors = neighbors.into_iter().collect::<Vec<_>>();
        neighbors.sort_unstable();
        adjacency[state_index] = neighbors;
    }
    adjacency
}

fn minimax_path(
    states: &[State],
    adjacency: &[Vec<usize>],
    start: usize,
    target: usize,
) -> Result<(f64, Vec<Option<usize>>), RnaError> {
    let mut best = vec![f64::INFINITY; states.len()];
    let mut steps = vec![usize::MAX; states.len()];
    let mut parents = vec![None; states.len()];
    let mut queue = BinaryHeap::new();
    best[start] = states[start].energy;
    steps[start] = 0;
    queue.push(QueueEntry {
        saddle: best[start],
        steps: 0,
        index: start,
    });
    while let Some(current) = queue.pop() {
        if current.saddle > best[current.index] + 1.0e-12
            || ((current.saddle - best[current.index]).abs() <= 1.0e-12
                && current.steps > steps[current.index])
        {
            continue;
        }
        if current.index == target {
            return Ok((current.saddle, parents));
        }
        for &neighbor in &adjacency[current.index] {
            let candidate_saddle = current.saddle.max(states[neighbor].energy);
            let candidate_steps = current.steps + 1;
            let better = candidate_saddle + 1.0e-12 < best[neighbor]
                || ((candidate_saddle - best[neighbor]).abs() <= 1.0e-12
                    && (candidate_steps < steps[neighbor]
                        || (candidate_steps == steps[neighbor]
                            && parents[neighbor].is_none_or(|parent| current.index < parent))));
            if better {
                best[neighbor] = candidate_saddle;
                steps[neighbor] = candidate_steps;
                parents[neighbor] = Some(current.index);
                queue.push(QueueEntry {
                    saddle: candidate_saddle,
                    steps: candidate_steps,
                    index: neighbor,
                });
            }
        }
    }
    Err(RnaError::InvalidOption(
        "landscape endpoints are disconnected under single-pair moves".into(),
    ))
}

fn recover_path(
    parents: &[Option<usize>],
    start: usize,
    target: usize,
) -> Result<Vec<usize>, RnaError> {
    let mut path = vec![target];
    let mut current = target;
    while current != start {
        current = parents[current].ok_or_else(|| {
            RnaError::InvalidOption("failed to recover the exact landscape path".into())
        })?;
        path.push(current);
    }
    path.reverse();
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraints::ConstraintConfig;

    #[test]
    fn complete_landscape_returns_a_valid_global_minimax_path() {
        let sequence = "GGGAAACCC";
        let model = EnergyModel::with_dangles(37.0, 2).unwrap();
        let constraints =
            ConstraintModel::compile(sequence.len(), &ConstraintConfig::default()).unwrap();
        let result =
            exact_landscape_path(sequence, ".........", "(((...)))", 3, &model, &constraints)
                .unwrap();
        assert!(result.state_space_complete);
        assert!(result.state_count > result.path.len());
        assert_eq!(result.path.first().unwrap().structure, ".........");
        assert_eq!(result.path.last().unwrap().structure, "(((...)))");
        assert!(result.edge_count > 0);
        for adjacent in result.path.windows(2) {
            let changes = adjacent[0]
                .structure
                .bytes()
                .zip(adjacent[1].structure.bytes())
                .filter(|(left, right)| left != right)
                .count();
            assert_eq!(changes, 2);
        }
        let maximum = result
            .path
            .iter()
            .map(|step| step.energy_kcal_mol)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((maximum - result.saddle_energy_kcal_mol).abs() < 1.0e-12);
    }

    #[test]
    fn minimax_search_prefers_a_longer_path_with_a_lower_saddle() {
        let state = |energy| State {
            structure: String::new(),
            energy,
            pairs: Vec::new(),
            partner: Vec::new(),
        };
        // 0--1--3 is shorter but crosses energy 9. 0--2--4--3 has saddle 4.
        let states = vec![state(0.0), state(9.0), state(3.0), state(1.0), state(4.0)];
        let adjacency = vec![vec![1, 2], vec![0, 3], vec![0, 4], vec![1, 4], vec![2, 3]];
        let (saddle, parents) = minimax_path(&states, &adjacency, 0, 3).unwrap();
        assert_eq!(saddle, 4.0);
        assert_eq!(recover_path(&parents, 0, 3).unwrap(), vec![0, 2, 4, 3]);
    }
}
