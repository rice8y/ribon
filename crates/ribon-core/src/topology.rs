//! Fatgraph invariants for RNA chord diagrams.
//!
//! A secondary structure with at least one pair is represented by the
//! orientable one-vertex ribbon graph obtained by contracting the backbone
//! circle.  If `alpha` exchanges the two ends of every chord and `sigma` is
//! the cyclic order of chord ends along the backbone, the boundary components
//! are the cycles of `sigma alpha`.  Euler's formula then gives
//! `g = (1 + E - F) / 2`.

use crate::structure::{parse_structure, RnaError};
use serde::Serialize;
use std::collections::VecDeque;

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct FatgraphTopology {
    pub pair_count: usize,
    pub crossing_pair_count: usize,
    pub crossing_component_count: usize,
    pub maximum_crossing_degree: usize,
    pub boundary_components: usize,
    pub euler_characteristic: isize,
    pub genus: usize,
    pub planar: bool,
    pub signature: String,
    /// One-based pair indices grouped by connected component of the crossing
    /// graph. Isolated pairs are included as singleton components.
    pub crossing_components: Vec<Vec<usize>>,
    pub pairs: Vec<FatgraphPair>,
    pub method: &'static str,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct FatgraphPair {
    pub index: usize,
    pub i: usize,
    pub j: usize,
    pub crossing_degree: usize,
    pub component: usize,
}

/// Compute exact fatgraph invariants for an extended dot-bracket structure.
pub fn fatgraph_topology(sequence: &str, structure: &str) -> Result<FatgraphTopology, RnaError> {
    let parsed = parse_structure(sequence, structure)?;
    let mut pairs = parsed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    pairs.sort_unstable();
    Ok(fatgraph_from_pairs(&pairs))
}

pub(crate) fn fatgraph_from_pairs(pairs: &[(usize, usize)]) -> FatgraphTopology {
    let edge_count = pairs.len();
    if edge_count == 0 {
        return FatgraphTopology {
            pair_count: 0,
            crossing_pair_count: 0,
            crossing_component_count: 0,
            maximum_crossing_degree: 0,
            boundary_components: 1,
            euler_characteristic: 2,
            genus: 0,
            planar: true,
            signature: "g0/b1/e0".into(),
            crossing_components: Vec::new(),
            pairs: Vec::new(),
            method: "orientable one-vertex chord-diagram fatgraph",
        };
    }

    let mut endpoints = Vec::with_capacity(2 * edge_count);
    for (pair_index, &(i, j)) in pairs.iter().enumerate() {
        endpoints.push((i, pair_index));
        endpoints.push((j, pair_index));
    }
    endpoints.sort_unstable();

    let mut pair_darts = vec![[usize::MAX; 2]; edge_count];
    for (dart, &(_, pair_index)) in endpoints.iter().enumerate() {
        let slot = usize::from(pair_darts[pair_index][0] != usize::MAX);
        pair_darts[pair_index][slot] = dart;
    }
    let mut alpha = vec![0usize; 2 * edge_count];
    for darts in &pair_darts {
        alpha[darts[0]] = darts[1];
        alpha[darts[1]] = darts[0];
    }
    let mut seen = vec![false; 2 * edge_count];
    let mut boundary_components = 0usize;
    for start in 0..seen.len() {
        if seen[start] {
            continue;
        }
        boundary_components += 1;
        let mut dart = start;
        while !seen[dart] {
            seen[dart] = true;
            dart = (alpha[dart] + 1) % seen.len();
        }
    }

    let mut adjacency = vec![Vec::new(); edge_count];
    let mut crossing_pair_count = 0usize;
    for a in 0..edge_count {
        for b in (a + 1)..edge_count {
            if crosses(pairs[a], pairs[b]) {
                adjacency[a].push(b);
                adjacency[b].push(a);
                crossing_pair_count += 1;
            }
        }
    }
    let maximum_crossing_degree = adjacency.iter().map(Vec::len).max().unwrap_or(0);
    let mut visited = vec![false; edge_count];
    let mut crossing_components = Vec::new();
    for start in 0..edge_count {
        if visited[start] {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        visited[start] = true;
        let mut component = Vec::new();
        while let Some(current) = queue.pop_front() {
            component.push(current + 1);
            for &next in &adjacency[current] {
                if !visited[next] {
                    visited[next] = true;
                    queue.push_back(next);
                }
            }
        }
        component.sort_unstable();
        crossing_components.push(component);
    }
    crossing_components.sort();
    let crossing_component_count = crossing_components
        .iter()
        .filter(|component| {
            component
                .iter()
                .any(|&index| !adjacency[index - 1].is_empty())
        })
        .count();
    let mut component_of = vec![0usize; edge_count];
    let mut crossing_component = 0usize;
    for members in &crossing_components {
        if !members
            .iter()
            .any(|&member| !adjacency[member - 1].is_empty())
        {
            continue;
        }
        crossing_component += 1;
        for &member in members {
            component_of[member - 1] = crossing_component;
        }
    }
    let pair_summaries = pairs
        .iter()
        .enumerate()
        .map(|(index, &(i, j))| FatgraphPair {
            index: index + 1,
            i: i + 1,
            j: j + 1,
            crossing_degree: adjacency[index].len(),
            component: component_of[index],
        })
        .collect();

    let genus_numerator = 1isize + edge_count as isize - boundary_components as isize;
    debug_assert!(genus_numerator >= 0 && genus_numerator % 2 == 0);
    let genus = (genus_numerator / 2) as usize;
    let euler_characteristic = 2 - 2 * genus as isize;
    FatgraphTopology {
        pair_count: edge_count,
        crossing_pair_count,
        crossing_component_count,
        maximum_crossing_degree,
        boundary_components,
        euler_characteristic,
        genus,
        planar: genus == 0,
        signature: format!("g{genus}/b{boundary_components}/e{edge_count}"),
        crossing_components,
        pairs: pair_summaries,
        method: "orientable one-vertex chord-diagram fatgraph",
    }
}

fn crosses((i, j): (usize, usize), (k, l): (usize, usize)) -> bool {
    (i < k && k < j && j < l) || (k < i && i < l && l < j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planar_and_crossing_chord_diagrams_have_expected_genus() {
        let nested = fatgraph_from_pairs(&[(0, 5), (1, 4), (2, 3)]);
        assert_eq!((nested.genus, nested.boundary_components), (0, 4));
        assert!(nested.planar);

        let h_type = fatgraph_from_pairs(&[(0, 4), (2, 6)]);
        assert_eq!((h_type.genus, h_type.boundary_components), (1, 1));
        assert_eq!(h_type.crossing_pair_count, 1);
        assert_eq!(h_type.crossing_component_count, 1);
    }

    #[test]
    fn extended_dot_bracket_is_accepted() {
        let topology = fatgraph_topology("GCGCGCG", "(.[.).]").unwrap();
        assert_eq!(topology.genus, 1);
    }
}
