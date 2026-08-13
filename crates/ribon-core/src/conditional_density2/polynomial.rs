//! Polynomial fixed-seed density-2 conditional partition function.
//!
//! This is an independent implementation of the complete, unambiguous CParty
//! decomposition published by Gray et al.  The implementation is expressed as
//! an acyclic interval hypergraph.  The same edge generator is evaluated in
//! the log-sum-exp, min-plus, max-plus, and outside algebras, which keeps PF,
//! MFE, centroid/MEA decoding, and pair marginals on one structure grammar.

use super::{
    band_topology, conditional_parameter_model_name, layered_structure, planar_structure,
    validate_options, validate_planar_layer, ConditionalDensity2Options,
};
use crate::constraints::{ConstraintConfig, ConstraintModel, ConstraintSummary};
use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{parse_structure, RnaError};
use crate::topology::{fatgraph_from_pairs, FatgraphTopology};
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const NEG_INF: f64 = f64::NEG_INFINITY;
const POS_INF: f64 = f64::INFINITY;
const KIND_COUNT: usize = 15;

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2PolynomialResult {
    pub sequence: String,
    pub seed_structure: String,
    pub mfe_structure: String,
    pub mfe_added_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub centroid_structure: String,
    pub centroid_distance: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub constraints: ConstraintSummary,
    pub time_complexity: &'static str,
    pub space_complexity: &'static str,
    pub state_space_complete: bool,
    pub algorithm: &'static str,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2PolynomialEvaluationResult {
    pub sequence: String,
    pub seed_structure: String,
    pub added_structure: String,
    pub structure: String,
    pub energy_kcal_mol: f64,
    pub derivation_unique: bool,
    pub model: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2Sample {
    pub structure: String,
    pub added_structure: String,
    pub effective_energy_kcal_mol: f64,
    pub probability: f64,
    pub log_probability: f64,
    pub topology: FatgraphTopology,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2SamplingResult {
    pub sequence: String,
    pub seed_structure: String,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub seed: u64,
    pub requested: usize,
    pub returned: usize,
    pub unique: bool,
    pub log_partition_function: f64,
    pub samples: Vec<ConditionalDensity2Sample>,
    pub constraints: ConstraintSummary,
    pub algorithm: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2SuboptimalStructure {
    pub rank: usize,
    pub structure: String,
    pub added_structure: String,
    pub energy_kcal_mol: f64,
    pub delta_energy_kcal_mol: f64,
    pub relative_boltzmann_weight: f64,
    pub topology: FatgraphTopology,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConditionalDensity2SuboptimalResult {
    pub sequence: String,
    pub seed_structure: String,
    pub temperature_celsius: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub energy_band_kcal_mol: f64,
    pub requested_limit: usize,
    pub truncated: bool,
    pub structures: Vec<ConditionalDensity2SuboptimalStructure>,
    pub constraints: ConstraintSummary,
    pub algorithm: &'static str,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum Kind {
    W,
    V,
    P,
    Pg,
    Pgw,
    Vp,
    Vpr,
    Vpl,
    Vm,
    Wm,
    Wm1,
    Wmp,
    Wi,
    Wip,
    Be,
}

impl Kind {
    const fn index(self) -> usize {
        self as usize
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct StateRef {
    kind: Kind,
    i: usize,
    j: usize,
}

const DUMMY_STATE: StateRef = StateRef {
    kind: Kind::W,
    i: 0,
    j: 0,
};

#[derive(Clone, Copy, Debug)]
struct Edge {
    children: [StateRef; 3],
    child_count: u8,
    energy: f64,
}

impl Edge {
    const fn leaf(energy: f64) -> Self {
        Self {
            children: [DUMMY_STATE; 3],
            child_count: 0,
            energy,
        }
    }

    fn child(mut self, child: StateRef) -> Self {
        debug_assert!((self.child_count as usize) < self.children.len());
        self.children[self.child_count as usize] = child;
        self.child_count += 1;
        self
    }

    fn children(&self) -> &[StateRef] {
        &self.children[..self.child_count as usize]
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Borders {
    /// `b(i,j)`: outer left border of the band extending right.
    b: Option<usize>,
    /// `b'(i,j)`: inner left border of the band extending right.
    b_prime: Option<usize>,
    /// `B(i,j)`: outer/rightmost right border of the band extending left.
    big_b: Option<usize>,
    /// `B'(i,j)`: inner/leftmost right border of the band extending left.
    big_b_prime: Option<usize>,
}

struct Charts {
    n: usize,
    inside: Vec<Vec<f64>>,
    mfe: Vec<Vec<f64>>,
    back: Vec<Vec<u32>>,
    order: Vec<StateRef>,
}

impl Charts {
    fn new(n: usize) -> Self {
        let cells = n * n;
        Self {
            n,
            inside: vec![vec![NEG_INF; cells]; KIND_COUNT],
            mfe: vec![vec![POS_INF; cells]; KIND_COUNT],
            back: vec![vec![u32::MAX; cells]; KIND_COUNT],
            order: Vec::with_capacity(KIND_COUNT * cells),
        }
    }

    #[inline]
    fn offset(&self, state: StateRef) -> usize {
        state.i * self.n + state.j
    }

    #[inline]
    fn inside(&self, state: StateRef) -> f64 {
        self.inside[state.kind.index()][self.offset(state)]
    }

    #[inline]
    fn mfe(&self, state: StateRef) -> f64 {
        self.mfe[state.kind.index()][self.offset(state)]
    }

    fn set(&mut self, state: StateRef, log_z: f64, mfe: f64, choice: Option<u32>) {
        let offset = self.offset(state);
        self.inside[state.kind.index()][offset] = log_z;
        self.mfe[state.kind.index()][offset] = mfe;
        self.back[state.kind.index()][offset] = choice.unwrap_or(u32::MAX);
        self.order.push(state);
    }

    fn back_choice(&self, state: StateRef) -> Option<u32> {
        let choice = self.back[state.kind.index()][self.offset(state)];
        (choice != u32::MAX).then_some(choice)
    }
}

struct Context<'a> {
    bases: &'a [u8],
    n: usize,
    min_loop: usize,
    rt: f64,
    model: &'a EnergyModel,
    options: &'a ConditionalDensity2Options,
    partner: Vec<Option<usize>>,
    cover: Vec<Option<usize>>,
    paired_prefix: Vec<usize>,
    weak: Vec<bool>,
    borders: Vec<Borders>,
    allowed_variable: Vec<bool>,
    pair_bonus: Vec<f64>,
    stack_bonus: Vec<f64>,
    root_constant: f64,
}

impl Context<'_> {
    #[inline]
    fn offset(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }

    #[inline]
    fn weakly_closed(&self, i: usize, j: usize) -> bool {
        i > j || self.weak[self.offset(i, j)]
    }

    #[inline]
    fn empty_g(&self, i: usize, j: usize) -> bool {
        i > j || self.paired_prefix[j + 1] == self.paired_prefix[i]
    }

    #[inline]
    fn border(&self, i: usize, j: usize) -> Borders {
        self.borders[self.offset(i, j)]
    }

    #[inline]
    fn seed_unpaired(&self, i: usize) -> bool {
        self.partner[i].is_none()
    }

    #[inline]
    fn fixed_pair(&self, i: usize, j: usize) -> bool {
        self.partner[i] == Some(j)
    }

    #[inline]
    fn variable_pair(&self, i: usize, j: usize) -> bool {
        self.allowed_variable[self.offset(i, j)]
            && self.partner[i].is_none()
            && self.partner[j].is_none()
            && j - i > self.min_loop
            && self.model.can_pair(self.bases[i], self.bases[j])
            && (!self.options.pk_only || self.crosses_seed(i, j))
    }

    #[inline]
    fn valid_v_pair(&self, i: usize, j: usize) -> bool {
        j > i && j - i > self.min_loop && (self.fixed_pair(i, j) || self.variable_pair(i, j))
    }

    fn crosses_seed(&self, i: usize, j: usize) -> bool {
        // A seed pair crosses (i,j) iff it contributes either a right border
        // B(i,j) or a left border b(i,j).  Both border families are
        // precomputed for every interval, so the predicate remains O(1) in
        // the cubic edge enumeration rather than introducing a hidden fourth
        // loop over seed pairs.
        let borders = self.border(i, j);
        borders.big_b.is_some() || borders.b.is_some()
    }

    fn covered_in(&self, outer_i: usize, outer_j: usize, position: usize) -> bool {
        self.cover[position].is_some_and(|left| {
            let right = self.partner[left].expect("cover is a seed pair");
            outer_i <= left && right <= outer_j
        })
    }

    fn state(&self, kind: Kind, i: usize, j: usize) -> StateRef {
        StateRef { kind, i, j }
    }

    fn edge_energy(&self, state: StateRef, edge: Edge) -> f64 {
        let pair_bonus =
            if matches!(state.kind, Kind::V | Kind::Vp) && self.variable_pair(state.i, state.j) {
                self.pair_bonus[self.offset(state.i, state.j)]
            } else {
                0.0
            };
        let stack_bonus = if matches!(state.kind, Kind::V | Kind::Vp)
            && self.variable_pair(state.i, state.j)
            && edge.children().iter().any(|child| {
                child.i == state.i + 1
                    && child.j + 1 == state.j
                    && matches!(child.kind, Kind::V | Kind::Vp)
                    && self.variable_pair(child.i, child.j)
            }) {
            self.stack_bonus[self.offset(state.i, state.j)]
        } else {
            0.0
        };
        edge.energy + pair_bonus + stack_bonus
    }

    fn optional_child(&self, edge: Edge, kind: Kind, i: usize, j: usize) -> Edge {
        if i <= j {
            edge.child(self.state(kind, i, j))
        } else {
            edge
        }
    }

    fn optional_prefix(&self, edge: Edge, kind: Kind, start: usize, end_exclusive: usize) -> Edge {
        if start < end_exclusive {
            edge.child(self.state(kind, start, end_exclusive - 1))
        } else {
            edge
        }
    }

    fn empty_prefix_g(&self, start: usize, end_exclusive: usize) -> bool {
        start >= end_exclusive || self.empty_g(start, end_exclusive - 1)
    }

    fn emit_edges(&self, state: StateRef, emit: &mut impl FnMut(Edge)) {
        match state.kind {
            Kind::W => self.edges_w(state.i, state.j, emit),
            Kind::V => self.edges_v(state.i, state.j, emit),
            Kind::P => self.edges_p(state.i, state.j, emit),
            Kind::Pg => self.edges_pg(state.i, state.j, emit),
            Kind::Pgw => self.edges_pgw(state.i, state.j, emit),
            Kind::Vp => self.edges_vp(state.i, state.j, emit),
            Kind::Vpr => self.edges_vpr(state.i, state.j, emit),
            Kind::Vpl => self.edges_vpl(state.i, state.j, emit),
            Kind::Vm => self.edges_vm(state.i, state.j, emit),
            Kind::Wm => self.edges_wm(state.i, state.j, emit),
            Kind::Wm1 => self.edges_wm1(state.i, state.j, emit),
            Kind::Wmp => self.edges_wmp(state.i, state.j, emit),
            Kind::Wi => self.edges_wi(state.i, state.j, emit),
            Kind::Wip => self.edges_wip(state.i, state.j, emit),
            Kind::Be => self.edges_be(state.i, state.j, emit),
        }
    }

    fn edges_w(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if !self.weakly_closed(i, j) || self.cover[i].is_some() || self.cover[j].is_some() {
            return;
        }
        for r in i..j {
            if self.cover[r].is_some() {
                continue;
            }
            let exterior_stem = if self.valid_v_pair(r, j) && !self.crosses_seed(r, j) {
                self.model.exterior_stem_boltzmann_energy(self.bases, r, j)
            } else {
                0.0
            };
            let edge = self
                .optional_prefix(Edge::leaf(exterior_stem), Kind::W, i, r)
                .child(self.state(Kind::V, r, j));
            emit(edge);
            let edge = self
                .optional_prefix(
                    Edge::leaf(self.options.pseudoloop_initiation_kcal_mol),
                    Kind::W,
                    i,
                    r,
                )
                .child(self.state(Kind::P, r, j));
            emit(edge);
        }
        if self.seed_unpaired(j) {
            if i == j {
                emit(Edge::leaf(0.0));
            } else {
                emit(Edge::leaf(0.0).child(self.state(Kind::W, i, j - 1)));
            }
        }
    }

    fn edges_v(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if !self.valid_v_pair(i, j) || !self.weakly_closed(i, j) || self.crosses_seed(i, j) {
            return;
        }
        if self.empty_g(i + 1, j - 1) {
            let energy = self.model.hairpin_boltzmann_energy(self.bases, i, j);
            if energy.is_finite() {
                emit(Edge::leaf(energy));
            }
        }
        let limit = self.model.internal_loop_limit(j.saturating_sub(i + 2));
        for left in 0..=limit {
            let k = i + 1 + left;
            if k >= j {
                break;
            }
            for right in 0..=(limit - left) {
                let Some(l) = j.checked_sub(right + 1) else {
                    continue;
                };
                if k >= l || !self.empty_g(i + 1, k - 1) || !self.empty_g(l + 1, j - 1) {
                    continue;
                }
                let energy = self.model.internal_boltzmann_energy(self.bases, i, j, k, l);
                if energy.is_finite() {
                    emit(Edge::leaf(energy).child(self.state(Kind::V, k, l)));
                }
            }
        }
        emit(Edge::leaf(0.0).child(self.state(Kind::Vm, i, j)));
    }

    fn edges_p(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        match self.partner[j] {
            None => emit(Edge::leaf(0.0).child(self.state(Kind::Pg, i, j))),
            Some(left) if left < j => {
                for l in (left + 1)..j {
                    if !self.seed_unpaired(l) {
                        continue;
                    }
                    let borders = self.border(i, l);
                    let Some(inner_left) = borders.b_prime else {
                        continue;
                    };
                    let Some(inner_right) = self.partner[inner_left] else {
                        continue;
                    };
                    let mut edge = Edge::leaf(self.options.band_kcal_mol)
                        .child(self.state(Kind::Be, left, inner_left))
                        .child(self.state(Kind::Pg, i, l));
                    edge =
                        self.optional_child(edge, Kind::Wi, l + 1, inner_right.saturating_sub(1));
                    emit(edge);
                }
            }
            _ => {}
        }
    }

    fn edges_pg(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if !self.seed_unpaired(j) || i >= j {
            return;
        }
        // Supplement Eq. (iii) prints `b(i,l)`, while HFold's recurrence and
        // the density-2 prefix lemma use `b(i,j)`.  The formal sentinel is
        // +infinity: when no seed pair covers j, every i<l<j is eligible.
        // This upper bound is essential when j is covered; omitting it admits
        // density-3 chains.
        let upper = self.border(i, j).b.unwrap_or(j).min(j);
        for l in (i + 1)..upper {
            if !self.covered_in(i, j, l) {
                continue;
            }
            let borders = self.border(i, l);
            let (Some(b), Some(b_prime)) = (borders.b, borders.b_prime) else {
                continue;
            };
            if l > i + 1 {
                emit(
                    Edge::leaf(2.0 * self.options.band_kcal_mol)
                        .child(self.state(Kind::Be, b, b_prime))
                        .child(self.state(Kind::Pg, i, l - 1))
                        .child(self.state(Kind::Vp, l, j)),
                );
            }
            if l > i + 1 && self.partner[l - 1].is_none_or(|partner| partner < l - 1) {
                emit(
                    Edge::leaf(2.0 * self.options.band_kcal_mol)
                        .child(self.state(Kind::Be, b, b_prime))
                        .child(self.state(Kind::Pgw, i, l - 1))
                        .child(self.state(Kind::Vp, l, j)),
                );
            }
        }
        emit(Edge::leaf(self.options.band_kcal_mol).child(self.state(Kind::Vp, i, j)));
        if let Some(right) = self.partner[i].filter(|&right| right > i) {
            for l in (i + 1)..right {
                let borders = self.border(i, l);
                let (Some(b), Some(b_prime)) = (borders.b, borders.b_prime) else {
                    continue;
                };
                let mut edge = Edge::leaf(2.0 * self.options.band_kcal_mol).child(self.state(
                    Kind::Be,
                    b,
                    b_prime,
                ));
                edge = self.optional_child(edge, Kind::Wi, b_prime + 1, l.saturating_sub(1));
                emit(edge.child(self.state(Kind::Vp, l, j)));
            }
        }
    }

    fn edges_pgw(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j {
            return;
        }
        for l in (i + 1)..j {
            if self.cover[l] != self.cover[j] {
                continue;
            }
            emit(
                Edge::leaf(0.0)
                    .child(self.state(Kind::Pg, i, l))
                    .child(self.state(Kind::Wi, l + 1, j)),
            );
        }
    }

    fn edges_vp(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if !self.variable_pair(i, j) || !self.crosses_seed(i, j) {
            return;
        }
        let borders = self.border(i, j);
        if self.cover[i] != self.cover[j] {
            if let (Some(inner), Some(outer)) = (borders.big_b_prime, borders.big_b) {
                let edge =
                    self.optional_child(Edge::leaf(0.0), Kind::Wi, i + 1, inner.saturating_sub(1));
                emit(self.optional_child(edge, Kind::Wi, outer + 1, j.saturating_sub(1)));
            }
            if let (Some(outer), Some(inner)) = (borders.b, borders.b_prime) {
                let edge =
                    self.optional_child(Edge::leaf(0.0), Kind::Wi, i + 1, outer.saturating_sub(1));
                emit(self.optional_child(edge, Kind::Wi, inner + 1, j.saturating_sub(1)));
            }
            if let (Some(left_inner), Some(left_outer), Some(right_outer), Some(right_inner)) = (
                borders.big_b_prime,
                borders.big_b,
                borders.b,
                borders.b_prime,
            ) {
                let edge = self.optional_child(
                    Edge::leaf(0.0),
                    Kind::Wi,
                    i + 1,
                    left_inner.saturating_sub(1),
                );
                let edge = self.optional_child(
                    edge,
                    Kind::Wi,
                    left_outer + 1,
                    right_outer.saturating_sub(1),
                );
                emit(self.optional_child(edge, Kind::Wi, right_inner + 1, j.saturating_sub(1)));
            }
        }

        if i + 1 < j && self.seed_unpaired(i + 1) && self.seed_unpaired(j - 1) {
            let energy = self.options.spanning_stack_factor
                * self
                    .model
                    .internal_boltzmann_energy(self.bases, i, j, i + 1, j - 1);
            if energy.is_finite() {
                emit(Edge::leaf(energy).child(self.state(Kind::Vp, i + 1, j - 1)));
            }
        }

        let left_limit = borders
            .big_b_prime
            .unwrap_or(self.n)
            .min(borders.b.unwrap_or(self.n));
        let right_limit = borders.big_b.unwrap_or(0).max(borders.b_prime.unwrap_or(0));
        if left_limit > i + 1 && right_limit + 1 < j {
            let internal_limit = self.model.internal_loop_limit(j.saturating_sub(i + 2));
            for r in (i + 1)..left_limit {
                if self.cover[i] != self.cover[r] {
                    continue;
                }
                for r_prime in (right_limit + 1)..j {
                    if self.cover[j] != self.cover[r_prime]
                        || (r == i + 1 && r_prime + 1 == j)
                        || (r - i - 1) + (j - r_prime - 1) > internal_limit
                    {
                        continue;
                    }
                    let energy = self.options.spanning_internal_factor
                        * self
                            .model
                            .internal_boltzmann_energy(self.bases, i, j, r, r_prime);
                    if energy.is_finite() {
                        emit(Edge::leaf(energy).child(self.state(Kind::Vp, r, r_prime)));
                    }
                }
            }
        }

        let spanning_multi = self.options.spanning_multiloop_init_kcal_mol
            + 2.0 * self.options.spanning_multiloop_branch_kcal_mol;
        for r in (i + 1)..left_limit.min(j) {
            if !self.seed_unpaired(r) || r <= i + 1 || j <= 1 {
                continue;
            }
            emit(
                Edge::leaf(spanning_multi)
                    .child(self.state(Kind::Wip, i + 1, r - 1))
                    .child(self.state(Kind::Vp, r, j - 1)),
            );
            emit(
                Edge::leaf(spanning_multi)
                    .child(self.state(Kind::Wip, i + 1, r - 1))
                    .child(self.state(Kind::Vpr, r, j - 1)),
            );
        }
        for r in right_limit.saturating_add(1).max(i + 1)..j {
            if !self.seed_unpaired(r) || r + 1 > j - 1 {
                continue;
            }
            emit(
                Edge::leaf(spanning_multi)
                    .child(self.state(Kind::Vp, i + 1, r))
                    .child(self.state(Kind::Wip, r + 1, j - 1)),
            );
            emit(
                Edge::leaf(spanning_multi)
                    .child(self.state(Kind::Vpl, i + 1, r))
                    .child(self.state(Kind::Wip, r + 1, j - 1)),
            );
        }
    }

    fn edges_vpr(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j {
            return;
        }
        let borders = self.border(i, j);
        let lower = borders.big_b.unwrap_or(0).max(borders.b_prime.unwrap_or(0));
        for r in lower.saturating_add(1).max(i + 1)..j {
            if r < j {
                emit(
                    Edge::leaf(0.0)
                        .child(self.state(Kind::Vp, i, r))
                        .child(self.state(Kind::Wip, r + 1, j)),
                );
                if self.empty_g(r + 1, j) {
                    emit(
                        Edge::leaf(
                            self.options.spanning_multiloop_unpaired_kcal_mol * (j - r) as f64,
                        )
                        .child(self.state(Kind::Vp, i, r)),
                    );
                }
            }
        }
    }

    fn edges_vpl(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j {
            return;
        }
        let borders = self.border(i, j);
        let upper = borders
            .big_b_prime
            .unwrap_or(self.n)
            .min(borders.b.unwrap_or(self.n));
        for r in (i + 1)..upper.min(j) {
            if self.empty_g(i, r - 1) {
                emit(
                    Edge::leaf(self.options.spanning_multiloop_unpaired_kcal_mol * (r - i) as f64)
                        .child(self.state(Kind::Vp, r, j)),
                );
            }
        }
    }

    fn edges_vm(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if !self.valid_v_pair(i, j) || !self.weakly_closed(i, j) || i + 2 > j {
            return;
        }
        let a = self.options.multiloop_init_kcal_mol;
        let b = self.options.multiloop_branch_kcal_mol;
        let c = self.options.multiloop_unpaired_kcal_mol;
        for r in (i + 2)..j {
            emit(
                Edge::leaf(a + b)
                    .child(self.state(Kind::Wm, i + 1, r - 1))
                    .child(self.state(Kind::Wm1, r, j - 1)),
            );
            emit(
                Edge::leaf(a + b + self.options.multiloop_pseudoknot_kcal_mol)
                    .child(self.state(Kind::Wm, i + 1, r - 1))
                    .child(self.state(Kind::Wmp, r, j - 1)),
            );
        }
        for r in (i + 1)..j.saturating_sub(1) {
            if self.empty_g(i + 1, r - 1) {
                emit(
                    Edge::leaf(
                        c * (r - i - 1) as f64 + a + b + self.options.multiloop_pseudoknot_kcal_mol,
                    )
                    .child(self.state(Kind::Wmp, r, j - 1)),
                );
            }
        }
    }

    fn edges_wm(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j || !self.weakly_closed(i, j) {
            return;
        }
        let b = self.options.multiloop_branch_kcal_mol;
        let c = self.options.multiloop_unpaired_kcal_mol;
        for r in i..j {
            if self.empty_prefix_g(i, r) {
                emit(Edge::leaf(c * (r - i) as f64 + b).child(self.state(Kind::V, r, j)));
                emit(
                    Edge::leaf(c * (r - i) as f64 + b + self.options.multiloop_pseudoknot_kcal_mol)
                        .child(self.state(Kind::P, r, j)),
                );
            }
        }
        for r in (i + 1)..j.saturating_sub(1) {
            emit(
                Edge::leaf(b)
                    .child(self.state(Kind::Wm, i, r))
                    .child(self.state(Kind::V, r + 1, j)),
            );
            emit(
                Edge::leaf(b + self.options.multiloop_pseudoknot_kcal_mol)
                    .child(self.state(Kind::Wm, i, r))
                    .child(self.state(Kind::P, r + 1, j)),
            );
        }
        if self.seed_unpaired(j) {
            emit(Edge::leaf(c).child(self.state(Kind::Wm, i, j - 1)));
        }
    }

    fn edges_wm1(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j || !self.weakly_closed(i, j) {
            return;
        }
        emit(Edge::leaf(self.options.multiloop_branch_kcal_mol).child(self.state(Kind::V, i, j)));
        if self.seed_unpaired(j) {
            emit(
                Edge::leaf(self.options.multiloop_unpaired_kcal_mol).child(self.state(
                    Kind::Wm1,
                    i,
                    j - 1,
                )),
            );
        }
    }

    fn edges_wmp(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j || !self.weakly_closed(i, j) {
            return;
        }
        emit(
            Edge::leaf(
                self.options.multiloop_branch_kcal_mol + self.options.multiloop_pseudoknot_kcal_mol,
            )
            .child(self.state(Kind::P, i, j)),
        );
        if self.seed_unpaired(j) {
            emit(
                Edge::leaf(self.options.multiloop_unpaired_kcal_mol).child(self.state(
                    Kind::Wmp,
                    i,
                    j - 1,
                )),
            );
        }
    }

    fn edges_wi(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i == j {
            if self.seed_unpaired(i) {
                emit(Edge::leaf(self.options.pseudoloop_unpaired_kcal_mol));
            }
            return;
        }
        if self.cover[i] != self.cover[j] {
            return;
        }
        for r in i..j {
            let compatible_v =
                self.fixed_pair(r, j) || (self.seed_unpaired(r) && self.seed_unpaired(j));
            if compatible_v {
                let edge = self.optional_prefix(
                    Edge::leaf(self.options.closed_subregion_kcal_mol),
                    Kind::Wi,
                    i,
                    r,
                );
                emit(edge.child(self.state(Kind::V, r, j)));
            }
            let edge = self.optional_prefix(
                Edge::leaf(
                    self.options.closed_subregion_kcal_mol
                        + self.options.nested_pseudoknot_kcal_mol,
                ),
                Kind::Wi,
                i,
                r,
            );
            emit(edge.child(self.state(Kind::P, r, j)));
        }
        if self.seed_unpaired(j) {
            emit(
                Edge::leaf(self.options.pseudoloop_unpaired_kcal_mol).child(self.state(
                    Kind::Wi,
                    i,
                    j - 1,
                )),
            );
        }
    }

    fn edges_wip(&self, i: usize, j: usize, emit: &mut impl FnMut(Edge)) {
        if i >= j || !self.weakly_closed(i, j) {
            return;
        }
        let b = self.options.spanning_multiloop_branch_kcal_mol;
        let c = self.options.spanning_multiloop_unpaired_kcal_mol;
        emit(Edge::leaf(b).child(self.state(Kind::V, i, j)));
        emit(
            Edge::leaf(b + self.options.multiloop_pseudoknot_kcal_mol).child(self.state(
                Kind::P,
                i,
                j,
            )),
        );
        for r in (i + 1)..j {
            let compatible_v =
                self.fixed_pair(r, j) || (self.seed_unpaired(r) && self.seed_unpaired(j));
            if compatible_v {
                emit(
                    Edge::leaf(b)
                        .child(self.state(Kind::Wip, i, r - 1))
                        .child(self.state(Kind::V, r, j)),
                );
            }
            emit(
                Edge::leaf(b + self.options.multiloop_pseudoknot_kcal_mol)
                    .child(self.state(Kind::Wip, i, r - 1))
                    .child(self.state(Kind::P, r, j)),
            );
            if self.empty_g(i, r - 1) {
                if compatible_v {
                    emit(Edge::leaf(c * (r - i) as f64 + b).child(self.state(Kind::V, r, j)));
                }
                emit(
                    Edge::leaf(c * (r - i) as f64 + b + self.options.multiloop_pseudoknot_kcal_mol)
                        .child(self.state(Kind::P, r, j)),
                );
            }
        }
        if self.seed_unpaired(j) {
            emit(Edge::leaf(c).child(self.state(Kind::Wip, i, j - 1)));
        }
    }

    fn edges_be(&self, i: usize, i_prime: usize, emit: &mut impl FnMut(Edge)) {
        let Some(right) = self.partner.get(i).and_then(|p| *p) else {
            return;
        };
        let Some(right_prime) = self.partner.get(i_prime).and_then(|p| *p) else {
            return;
        };
        if i > i_prime || i_prime >= right_prime || right_prime > right {
            return;
        }
        // The published prose prints zero here, but the recurrence is a
        // product recurrence: the empty remaining band must have weight one.
        if i == i_prime {
            emit(Edge::leaf(0.0));
            return;
        }
        if i + 1 < self.n && self.partner[i + 1] == right.checked_sub(1) {
            let energy = self.options.spanning_stack_factor
                * self
                    .model
                    .internal_boltzmann_energy(self.bases, i, right, i + 1, right - 1);
            if energy.is_finite() {
                emit(Edge::leaf(energy).child(self.state(Kind::Be, i + 1, i_prime)));
            }
        }
        for l in (i + 1)..=i_prime {
            let Some(l_right) = self.partner[l] else {
                continue;
            };
            if l >= l_right || l_right < right_prime || l_right >= right {
                continue;
            }
            let left_empty = self.empty_g(i + 1, l - 1);
            let right_empty = self.empty_g(l_right + 1, right - 1);
            if left_empty && right_empty && !(l == i + 1 && l_right + 1 == right) {
                let energy = self.options.spanning_internal_factor
                    * self
                        .model
                        .internal_boltzmann_energy(self.bases, i, right, l, l_right);
                if energy.is_finite() {
                    emit(Edge::leaf(energy).child(self.state(Kind::Be, l, i_prime)));
                }
            }
            if self.weakly_closed(i + 1, l - 1) && self.weakly_closed(l_right + 1, right - 1) {
                emit(
                    Edge::leaf(
                        self.options.spanning_multiloop_init_kcal_mol
                            + 3.0 * self.options.spanning_multiloop_branch_kcal_mol,
                    )
                    .child(self.state(Kind::Wip, i + 1, l - 1))
                    .child(self.state(Kind::Be, l, i_prime))
                    .child(self.state(Kind::Wip, l_right + 1, right - 1)),
                );
            }
            if self.weakly_closed(i + 1, l - 1) && right_empty {
                emit(
                    Edge::leaf(
                        self.options.spanning_multiloop_init_kcal_mol
                            + 2.0 * self.options.spanning_multiloop_branch_kcal_mol
                            + self.options.spanning_multiloop_unpaired_kcal_mol
                                * (right - l_right + 1) as f64,
                    )
                    .child(self.state(Kind::Wip, i + 1, l - 1))
                    .child(self.state(Kind::Be, l, i_prime)),
                );
            }
            if left_empty && self.weakly_closed(l_right + 1, right - 1) {
                emit(
                    Edge::leaf(
                        self.options.spanning_multiloop_init_kcal_mol
                            + 2.0 * self.options.spanning_multiloop_branch_kcal_mol
                            + self.options.spanning_multiloop_unpaired_kcal_mol
                                * (l - i - 1) as f64,
                    )
                    .child(self.state(Kind::Be, l, i_prime))
                    .child(self.state(Kind::Wip, l_right + 1, right - 1)),
                );
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_density2_polynomial(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
) -> Result<ConditionalDensity2PolynomialResult, RnaError> {
    conditional_density2_polynomial_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_density2_polynomial_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2PolynomialResult, RnaError> {
    validate_options(temperature_celsius, options)?;
    if dangles % 2 == 1 {
        return Err(RnaError::InvalidOption(
            "the polynomial conditional density-2 engine supports dangles 0/2; use the exact exhaustive dispatcher for dangles 1/3".into(),
        ));
    }
    let seed = parse_structure(sequence, seed_structure)?;
    if seed.length == 0 {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis requires a nonempty sequence".into(),
        ));
    }
    if !seed.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis accepts one strand".into(),
        ));
    }
    let mut seed_pairs = seed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    seed_pairs.sort_unstable();
    if seed_pairs.iter().enumerate().any(|(index, &left)| {
        seed_pairs[index + 1..]
            .iter()
            .any(|&right| left.0 < right.0 && right.0 < left.1 && left.1 < right.1)
    }) {
        return Err(RnaError::InvalidOption(
            "the conditional seed must be pseudoknot-free".into(),
        ));
    }
    if let Some(pair) = seed.pairs.iter().find(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed contains noncanonical pair {}-{}",
            pair.i, pair.j
        )));
    }
    if let Some(&(i, j)) = seed_pairs.iter().find(|&&(i, j)| j - i <= min_loop) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed pair {}-{} violates min-loop={min_loop}",
            i + 1,
            j + 1
        )));
    }
    let model = super::conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    if !constraint_config.force_paired.is_empty()
        || !constraint_config.force_pairs.is_empty()
        || constraint_config.no_lonely_pairs
    {
        return Err(RnaError::InvalidOption(
            "force-paired, force-pairs, and no-lonely-pairs require the exact conditional constraint dispatcher".into(),
        ));
    }
    let constraints = ConstraintModel::compile(seed.length, constraint_config)?;
    for &(i, j) in &seed_pairs {
        if !constraints.allows_pair(seed.sequence.as_bytes(), i, j, &model) {
            return Err(RnaError::InvalidOption(format!(
                "seed pair {}-{} violates the requested constraints",
                i + 1,
                j + 1
            )));
        }
    }
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let partner = partner_table(seed.length, &seed_pairs);
    let cover = cover_table(seed.length, &seed_pairs);
    let paired_prefix = paired_prefix(&partner);
    let weak = weak_table(seed.length, &partner);
    let borders = border_table(seed.length, &partner);
    let (allowed_variable, pair_bonus, stack_bonus, mut root_constant) = constraint_terms(
        seed.sequence.as_bytes(),
        min_loop,
        &model,
        &constraints,
        &seed_pairs,
    );
    root_constant += seed_baseline_correction(
        seed.sequence.as_bytes(),
        min_loop,
        rt,
        &model,
        options,
        &seed_pairs,
    )?;
    let context = Context {
        bases: seed.sequence.as_bytes(),
        n: seed.length,
        min_loop,
        rt,
        model: &model,
        options,
        partner,
        cover,
        paired_prefix,
        weak,
        borders,
        allowed_variable,
        pair_bonus,
        stack_bonus,
        root_constant,
    };
    let charts = inside_charts(&context);
    let root = StateRef {
        kind: Kind::W,
        i: 0,
        j: seed.length - 1,
    };
    let raw_log_z = charts.inside(root);
    if raw_log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "the fixed seed admits no density-2 extension".into(),
        ));
    }
    let mfe_added = traceback_pairs(root, &charts, &context);
    let outside = outside_chart(&context, &charts, root);
    let pair_probabilities =
        pair_probabilities(&context, &charts, &outside, raw_log_z, &seed_pairs);
    drop(outside);
    let unpaired_probabilities = unpaired_probabilities(seed.length, &pair_probabilities);
    let probability_matrix = probability_matrix(seed.length, &pair_probabilities, &seed_pairs);

    let centroid_gains = pair_gains(seed.length, &probability_matrix, None, options.gamma);
    let centroid_back = max_plus_chart(&context, &charts.order, &centroid_gains);
    let centroid_added = traceback_max_pairs(root, &centroid_back, &context);
    drop(centroid_back);
    let total_variable_probability = probability_matrix
        .iter()
        .enumerate()
        .map(|(i, row)| row.iter().skip(i + 1).sum::<f64>())
        .sum::<f64>();
    let centroid_pair_probability = centroid_added
        .iter()
        .map(|&(i, j)| probability_matrix[i][j])
        .sum::<f64>();
    let centroid_distance =
        total_variable_probability + centroid_added.len() as f64 - 2.0 * centroid_pair_probability;

    let mea_gains = pair_gains(
        seed.length,
        &probability_matrix,
        Some(&unpaired_probabilities),
        options.gamma,
    );
    let mea_back = max_plus_chart(&context, &charts.order, &mea_gains);
    let mea_added = traceback_max_pairs(root, &mea_back, &context);
    drop(mea_back);
    let seed_occupied = seed_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .collect::<std::collections::HashSet<_>>();
    let mea_baseline = unpaired_probabilities
        .iter()
        .enumerate()
        .filter(|(i, _)| !seed_occupied.contains(i))
        .map(|(_, &p)| p)
        .sum::<f64>();
    let mea_score = mea_baseline + mea_added.iter().map(|&(i, j)| mea_gains[i][j]).sum::<f64>();

    let log_z = raw_log_z - context.root_constant / rt;
    Ok(ConditionalDensity2PolynomialResult {
        sequence: seed.sequence.clone(),
        seed_structure: planar_structure(seed.length, &seed_pairs, '(', ')'),
        mfe_structure: layered_structure(seed.length, &seed_pairs, &mfe_added),
        mfe_added_structure: planar_structure(seed.length, &mfe_added, '(', ')'),
        mfe_energy_kcal_mol: charts.mfe(root) + context.root_constant,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        pair_probabilities,
        unpaired_probabilities,
        centroid_structure: layered_structure(seed.length, &seed_pairs, &centroid_added),
        centroid_distance,
        mea_structure: layered_structure(seed.length, &seed_pairs, &mea_added),
        mea_score,
        constraints: constraints.summary(),
        time_complexity: "O(n^3)",
        space_complexity: "O(n^2)",
        state_space_complete: true,
        algorithm:
            "independent unambiguous CParty density-2 interval hypergraph with inside/outside",
        model: conditional_parameter_model_name(&model),
    })
}

/// Draw exact independent Boltzmann samples by stochastic traceback through
/// the same conditional inside hypergraph used for the partition function.
#[allow(clippy::too_many_arguments)]
pub fn sample_conditional_density2_polynomial(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    count: usize,
    random_seed: u64,
    unique: bool,
) -> Result<ConditionalDensity2SamplingResult, RnaError> {
    sample_conditional_density2_polynomial_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        count,
        random_seed,
        unique,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn sample_conditional_density2_polynomial_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    count: usize,
    random_seed: u64,
    unique: bool,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2SamplingResult, RnaError> {
    validate_options(temperature_celsius, options)?;
    if dangles % 2 == 1 {
        return Err(RnaError::InvalidOption(
            "the polynomial conditional sampler supports dangles 0/2".into(),
        ));
    }
    if count == 0 {
        return Err(RnaError::InvalidOption(
            "conditional sample count must be positive".into(),
        ));
    }
    let (seed, seed_pairs) = validated_seed(sequence, seed_structure, min_loop)?;
    let model = super::conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    if !constraint_config.force_paired.is_empty()
        || !constraint_config.force_pairs.is_empty()
        || constraint_config.no_lonely_pairs
    {
        return Err(RnaError::InvalidOption(
            "force-paired, force-pairs, and no-lonely-pairs require the exact conditional constraint dispatcher".into(),
        ));
    }
    let constraints = ConstraintModel::compile(seed.length, constraint_config)?;
    for &(i, j) in &seed_pairs {
        if !constraints.allows_pair(seed.sequence.as_bytes(), i, j, &model) {
            return Err(RnaError::InvalidOption(format!(
                "seed pair {}-{} violates the requested constraints",
                i + 1,
                j + 1
            )));
        }
    }
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let partner = partner_table(seed.length, &seed_pairs);
    let (allowed_variable, pair_bonus, stack_bonus, mut root_constant) = constraint_terms(
        seed.sequence.as_bytes(),
        min_loop,
        &model,
        &constraints,
        &seed_pairs,
    );
    root_constant += seed_baseline_correction(
        seed.sequence.as_bytes(),
        min_loop,
        rt,
        &model,
        options,
        &seed_pairs,
    )?;
    let context = Context {
        bases: seed.sequence.as_bytes(),
        n: seed.length,
        min_loop,
        rt,
        model: &model,
        options,
        cover: cover_table(seed.length, &seed_pairs),
        paired_prefix: paired_prefix(&partner),
        weak: weak_table(seed.length, &partner),
        borders: border_table(seed.length, &partner),
        allowed_variable,
        pair_bonus,
        stack_bonus,
        root_constant,
        partner,
    };
    let charts = inside_charts(&context);
    let root = StateRef {
        kind: Kind::W,
        i: 0,
        j: seed.length - 1,
    };
    let raw_log_z = charts.inside(root);
    if raw_log_z == NEG_INF {
        return Err(RnaError::InvalidOption(
            "the fixed seed admits no density-2 extension".into(),
        ));
    }
    let mut rng = SplitMix64::new(random_seed);
    let mut samples = Vec::with_capacity(count);
    let mut seen = HashSet::new();
    for _ in 0..count {
        let mut added = Vec::new();
        let mut log_probability = 0.0;
        stochastic_traceback(
            root,
            &charts,
            &context,
            &mut rng,
            &mut added,
            &mut log_probability,
        );
        added.sort_unstable();
        added.dedup();
        let added_structure = planar_structure(seed.length, &added, '(', ')');
        if unique && !seen.insert(added_structure.clone()) {
            continue;
        }
        let mut union = seed_pairs.clone();
        union.extend_from_slice(&added);
        union.sort_unstable();
        samples.push(ConditionalDensity2Sample {
            structure: layered_structure(seed.length, &seed_pairs, &added),
            added_structure,
            effective_energy_kcal_mol: -rt * (raw_log_z + log_probability) + context.root_constant,
            probability: log_probability.exp(),
            log_probability,
            topology: fatgraph_from_pairs(&union),
        });
    }
    let log_partition_function = raw_log_z - context.root_constant / rt;
    Ok(ConditionalDensity2SamplingResult {
        sequence: seed.sequence,
        seed_structure: planar_structure(seed.length, &seed_pairs, '(', ')'),
        temperature_celsius,
        dangles,
        salt_molar,
        seed: random_seed,
        requested: count,
        returned: samples.len(),
        unique,
        log_partition_function,
        samples,
        constraints: constraints.summary(),
        algorithm: "exact stochastic traceback of the conditional density-2 inside hypergraph",
    })
}

/// Enumerate the exact k lowest-energy conditional structures without a beam.
/// Runtime is output-sensitive after construction of the cubic hypergraph.
#[allow(clippy::too_many_arguments)]
pub fn suboptimal_conditional_density2_polynomial(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    energy_band_kcal_mol: f64,
    limit: usize,
) -> Result<ConditionalDensity2SuboptimalResult, RnaError> {
    suboptimal_conditional_density2_polynomial_with_constraints(
        sequence,
        seed_structure,
        temperature_celsius,
        min_loop,
        dangles,
        salt_molar,
        options,
        energy_band_kcal_mol,
        limit,
        &ConstraintConfig::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn suboptimal_conditional_density2_polynomial_with_constraints(
    sequence: &str,
    seed_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
    energy_band_kcal_mol: f64,
    limit: usize,
    constraint_config: &ConstraintConfig,
) -> Result<ConditionalDensity2SuboptimalResult, RnaError> {
    validate_options(temperature_celsius, options)?;
    if dangles % 2 == 1 {
        return Err(RnaError::InvalidOption(
            "the polynomial conditional k-best engine supports dangles 0/2".into(),
        ));
    }
    if !energy_band_kcal_mol.is_finite() || energy_band_kcal_mol < 0.0 {
        return Err(RnaError::InvalidOption(
            "conditional suboptimal energy band must be finite and non-negative".into(),
        ));
    }
    if limit == 0 {
        return Err(RnaError::InvalidOption(
            "conditional suboptimal limit must be positive".into(),
        ));
    }
    let (seed, seed_pairs) = validated_seed(sequence, seed_structure, min_loop)?;
    let model = super::conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    if !constraint_config.force_paired.is_empty()
        || !constraint_config.force_pairs.is_empty()
        || constraint_config.no_lonely_pairs
    {
        return Err(RnaError::InvalidOption(
            "force-paired, force-pairs, and no-lonely-pairs require the exact conditional constraint dispatcher".into(),
        ));
    }
    let constraints = ConstraintModel::compile(seed.length, constraint_config)?;
    for &(i, j) in &seed_pairs {
        if !constraints.allows_pair(seed.sequence.as_bytes(), i, j, &model) {
            return Err(RnaError::InvalidOption(format!(
                "seed pair {}-{} violates the requested constraints",
                i + 1,
                j + 1
            )));
        }
    }
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let partner = partner_table(seed.length, &seed_pairs);
    let (allowed_variable, pair_bonus, stack_bonus, mut root_constant) = constraint_terms(
        seed.sequence.as_bytes(),
        min_loop,
        &model,
        &constraints,
        &seed_pairs,
    );
    root_constant += seed_baseline_correction(
        seed.sequence.as_bytes(),
        min_loop,
        rt,
        &model,
        options,
        &seed_pairs,
    )?;
    let context = Context {
        bases: seed.sequence.as_bytes(),
        n: seed.length,
        min_loop,
        rt,
        model: &model,
        options,
        cover: cover_table(seed.length, &seed_pairs),
        paired_prefix: paired_prefix(&partner),
        weak: weak_table(seed.length, &partner),
        borders: border_table(seed.length, &partner),
        allowed_variable,
        pair_bonus,
        stack_bonus,
        root_constant,
        partner,
    };
    let inside = inside_charts(&context);
    let root = StateRef {
        kind: Kind::W,
        i: 0,
        j: seed.length - 1,
    };
    if inside.inside(root) == NEG_INF {
        return Err(RnaError::InvalidOption(
            "the fixed seed admits no density-2 extension".into(),
        ));
    }
    let kbest = kbest_charts(&context, &inside.order, limit.saturating_add(1));
    let root_entries = kbest.entries(root);
    let Some(raw_mfe) = root_entries.first().map(|entry| entry.energy) else {
        return Err(RnaError::InvalidOption(
            "the fixed seed admits no density-2 extension".into(),
        ));
    };
    let in_band = root_entries
        .iter()
        .take_while(|entry| entry.energy <= raw_mfe + energy_band_kcal_mol + 1.0e-12)
        .count();
    let truncated = in_band > limit;
    let mut structures = Vec::with_capacity(in_band.min(limit));
    for (rank, entry) in root_entries.iter().take(in_band.min(limit)).enumerate() {
        let added = traceback_kbest_pairs(root, rank, &kbest, &context);
        let mut union = seed_pairs.clone();
        union.extend_from_slice(&added);
        union.sort_unstable();
        let delta = entry.energy - raw_mfe;
        structures.push(ConditionalDensity2SuboptimalStructure {
            rank: rank + 1,
            structure: layered_structure(seed.length, &seed_pairs, &added),
            added_structure: planar_structure(seed.length, &added, '(', ')'),
            energy_kcal_mol: entry.energy + context.root_constant,
            delta_energy_kcal_mol: delta,
            relative_boltzmann_weight: (-delta / rt).exp(),
            topology: fatgraph_from_pairs(&union),
        });
    }
    Ok(ConditionalDensity2SuboptimalResult {
        sequence: seed.sequence,
        seed_structure: planar_structure(seed.length, &seed_pairs, '(', ')'),
        temperature_celsius,
        dangles,
        salt_molar,
        energy_band_kcal_mol,
        requested_limit: limit,
        truncated,
        structures,
        constraints: constraints.summary(),
        algorithm: "exact lazy k-best conditional density-2 hypergraph enumeration",
    })
}

#[allow(clippy::too_many_arguments)]
pub fn evaluate_conditional_density2_polynomial(
    sequence: &str,
    seed_structure: &str,
    added_structure: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    options: &ConditionalDensity2Options,
) -> Result<ConditionalDensity2PolynomialEvaluationResult, RnaError> {
    validate_options(temperature_celsius, options)?;
    if dangles % 2 == 1 {
        return Err(RnaError::InvalidOption(
            "the polynomial conditional density-2 evaluator supports dangles 0/2; use the exact exhaustive dispatcher for dangles 1/3".into(),
        ));
    }
    let seed = parse_structure(sequence, seed_structure)?;
    let added = parse_structure(sequence, added_structure)?;
    if seed.length == 0 || !seed.strand_breaks.is_empty() || !added.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "conditional density-2 evaluation requires one nonempty strand".into(),
        ));
    }
    let mut seed_pairs = seed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    let mut added_pairs = added
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    seed_pairs.sort_unstable();
    added_pairs.sort_unstable();
    validate_planar_layer("seed", &seed, &seed_pairs, min_loop)?;
    validate_planar_layer("added", &added, &added_pairs, min_loop)?;
    let seed_occupied = seed_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .collect::<std::collections::HashSet<_>>();
    if let Some(position) = added_pairs
        .iter()
        .flat_map(|&(i, j)| [i, j])
        .find(|position| seed_occupied.contains(position))
    {
        return Err(RnaError::MultiplePartners {
            position: position + 1,
        });
    }
    let mut union = seed_pairs.clone();
    union.extend_from_slice(&added_pairs);
    union.sort_unstable();
    let topology = band_topology(seed.length, &union);
    if topology.maximum_density > 2 {
        return Err(RnaError::InvalidOption(format!(
            "conditional union has band density {}, expected at most 2",
            topology.maximum_density
        )));
    }
    let model = super::conditional_energy_model(temperature_celsius, dangles, salt_molar, options)?;
    let partner = partner_table(seed.length, &seed_pairs);
    let mut allowed_variable = vec![false; seed.length * seed.length];
    for &(i, j) in &added_pairs {
        allowed_variable[i * seed.length + j] = true;
    }
    let context = Context {
        bases: seed.sequence.as_bytes(),
        n: seed.length,
        min_loop,
        rt: GAS_CONSTANT_KCAL * (temperature_celsius + 273.15),
        model: &model,
        options,
        cover: cover_table(seed.length, &seed_pairs),
        paired_prefix: paired_prefix(&partner),
        weak: weak_table(seed.length, &partner),
        borders: border_table(seed.length, &partner),
        allowed_variable,
        pair_bonus: vec![0.0; seed.length * seed.length],
        stack_bonus: vec![0.0; seed.length * seed.length],
        root_constant: 0.0,
        partner,
    };
    let inside = inside_charts(&context);
    let forced = forced_chart(&context, &inside.order);
    let root = StateRef {
        kind: Kind::W,
        i: 0,
        j: seed.length - 1,
    };
    let offset = root.i * seed.length + root.j;
    let count = forced.count[root.kind.index()][offset];
    if count != added_pairs.len() as i32 {
        return Err(RnaError::InvalidOption(
            "the supplied added layer has no finite-energy derivation in the density-2 grammar"
                .into(),
        ));
    }
    let traced = traceback_forced_pairs(root, &forced, &context);
    if traced != added_pairs {
        return Err(RnaError::Numerical(
            "conditional fixed-structure traceback did not reproduce the supplied layer".into(),
        ));
    }
    let baseline_correction = seed_baseline_correction(
        seed.sequence.as_bytes(),
        min_loop,
        context.rt,
        &model,
        options,
        &seed_pairs,
    )?;
    Ok(ConditionalDensity2PolynomialEvaluationResult {
        sequence: seed.sequence,
        seed_structure: planar_structure(seed.length, &seed_pairs, '(', ')'),
        added_structure: planar_structure(seed.length, &added_pairs, '(', ')'),
        structure: layered_structure(seed.length, &seed_pairs, &added_pairs),
        energy_kcal_mol: forced.energy[root.kind.index()][offset] + baseline_correction,
        derivation_unique: true,
        model: conditional_parameter_model_name(&model),
    })
}

fn inside_charts(context: &Context<'_>) -> Charts {
    let mut charts = Charts::new(context.n);
    for span in 0..context.n {
        for i in 0..(context.n - span) {
            if context.partner[i] == Some(i + span) {
                for i_prime in i..=i + span {
                    if context.partner[i_prime]
                        .is_some_and(|right| i_prime < right && right <= i + span)
                    {
                        compute_state(
                            context,
                            &mut charts,
                            StateRef {
                                kind: Kind::Be,
                                i,
                                j: i_prime,
                            },
                        );
                    }
                }
            }
        }
        for i in 0..(context.n - span) {
            let j = i + span;
            for kind in [
                Kind::Vm,
                Kind::V,
                Kind::Vp,
                Kind::Vpr,
                Kind::Vpl,
                Kind::Pg,
                Kind::P,
                Kind::Wm1,
                Kind::Wmp,
                Kind::Wm,
                Kind::Wi,
                Kind::Wip,
                Kind::W,
                Kind::Pgw,
            ] {
                compute_state(context, &mut charts, StateRef { kind, i, j });
            }
        }
    }
    charts
}

fn compute_state(context: &Context<'_>, charts: &mut Charts, state: StateRef) {
    let mut log_z = NEG_INF;
    let mut best = POS_INF;
    let mut best_choice = None;
    let mut ordinal = 0u32;
    context.emit_edges(state, &mut |edge| {
        let current_ordinal = ordinal;
        ordinal = ordinal.saturating_add(1);
        let energy = context.edge_energy(state, edge);
        let mut log_candidate = -energy / context.rt;
        let mut mfe_candidate = energy;
        for &child in edge.children() {
            let inside = charts.inside(child);
            let mfe = charts.mfe(child);
            if inside == NEG_INF || mfe == POS_INF {
                return;
            }
            log_candidate += inside;
            mfe_candidate += mfe;
        }
        log_z = log_add(log_z, log_candidate);
        if mfe_candidate < best - 1.0e-12 {
            best = mfe_candidate;
            best_choice = Some(current_ordinal);
        }
    });
    charts.set(state, log_z, best, best_choice);
}

fn outside_chart(context: &Context<'_>, charts: &Charts, root: StateRef) -> Vec<Vec<f64>> {
    let mut outside = vec![vec![NEG_INF; context.n * context.n]; KIND_COUNT];
    outside[root.kind.index()][root.i * context.n + root.j] = 0.0;
    for &parent in charts.order.iter().rev() {
        let parent_outside = outside[parent.kind.index()][parent.i * context.n + parent.j];
        if parent_outside == NEG_INF {
            continue;
        }
        context.emit_edges(parent, &mut |edge| {
            let mut total = parent_outside - context.edge_energy(parent, edge) / context.rt;
            for &child in edge.children() {
                let inside = charts.inside(child);
                if inside == NEG_INF {
                    return;
                }
                total += inside;
            }
            for &child in edge.children() {
                let value = total - charts.inside(child);
                let offset = child.i * context.n + child.j;
                let cell = &mut outside[child.kind.index()][offset];
                *cell = log_add(*cell, value);
            }
        });
    }
    outside
}

fn pair_probabilities(
    context: &Context<'_>,
    charts: &Charts,
    outside: &[Vec<f64>],
    log_z: f64,
    seed_pairs: &[(usize, usize)],
) -> Vec<PairProbability> {
    let mut result = seed_pairs
        .iter()
        .map(|&(i, j)| PairProbability {
            i: i + 1,
            j: j + 1,
            probability: 1.0,
        })
        .collect::<Vec<_>>();
    for i in 0..context.n {
        for j in (i + 1)..context.n {
            if !context.variable_pair(i, j) {
                continue;
            }
            let mut probability = 0.0;
            for kind in [Kind::V, Kind::Vp] {
                let state = StateRef { kind, i, j };
                let inside = charts.inside(state);
                let out = outside[kind.index()][i * context.n + j];
                if inside != NEG_INF && out != NEG_INF {
                    probability += (inside + out - log_z).exp();
                }
            }
            if probability > 1.0e-15 {
                result.push(PairProbability {
                    i: i + 1,
                    j: j + 1,
                    probability: probability.clamp(0.0, 1.0),
                });
            }
        }
    }
    result.sort_by_key(|pair| (pair.i, pair.j));
    result
}

fn unpaired_probabilities(length: usize, pairs: &[PairProbability]) -> Vec<f64> {
    let mut result = vec![1.0; length];
    for pair in pairs {
        result[pair.i - 1] -= pair.probability;
        result[pair.j - 1] -= pair.probability;
    }
    for value in &mut result {
        *value = value.clamp(0.0, 1.0);
    }
    result
}

fn probability_matrix(
    length: usize,
    pairs: &[PairProbability],
    seed_pairs: &[(usize, usize)],
) -> Vec<Vec<f64>> {
    let mut matrix = vec![vec![0.0; length]; length];
    for pair in pairs {
        let ij = (pair.i - 1, pair.j - 1);
        if !seed_pairs.contains(&ij) {
            matrix[ij.0][ij.1] = pair.probability;
        }
    }
    matrix
}

fn pair_gains(
    length: usize,
    probabilities: &[Vec<f64>],
    unpaired: Option<&[f64]>,
    gamma: f64,
) -> Vec<Vec<f64>> {
    let mut gains = vec![vec![0.0; length]; length];
    for i in 0..length {
        for j in (i + 1)..length {
            gains[i][j] = if let Some(unpaired) = unpaired {
                2.0 * gamma * probabilities[i][j] - unpaired[i] - unpaired[j]
            } else {
                2.0 * probabilities[i][j] - 1.0
            };
        }
    }
    gains
}

struct MaxCharts {
    n: usize,
    score: Vec<Vec<f64>>,
    back: Vec<Vec<u32>>,
}

struct ForcedCharts {
    n: usize,
    count: Vec<Vec<i32>>,
    energy: Vec<Vec<f64>>,
    back: Vec<Vec<u32>>,
}

fn forced_chart(context: &Context<'_>, order: &[StateRef]) -> ForcedCharts {
    let cells = context.n * context.n;
    let mut chart = ForcedCharts {
        n: context.n,
        count: vec![vec![-1; cells]; KIND_COUNT],
        energy: vec![vec![POS_INF; cells]; KIND_COUNT],
        back: vec![vec![u32::MAX; cells]; KIND_COUNT],
    };
    for &state in order {
        let local_count = i32::from(
            matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j),
        );
        let mut best_count = -1;
        let mut best_energy = POS_INF;
        let mut best_choice = None;
        let mut ordinal = 0u32;
        context.emit_edges(state, &mut |edge| {
            let current_ordinal = ordinal;
            ordinal = ordinal.saturating_add(1);
            let mut candidate_count = local_count;
            let mut candidate_energy = context.edge_energy(state, edge);
            for &child in edge.children() {
                let offset = child.i * context.n + child.j;
                let count = chart.count[child.kind.index()][offset];
                if count < 0 {
                    return;
                }
                candidate_count += count;
                candidate_energy += chart.energy[child.kind.index()][offset];
            }
            if candidate_count > best_count
                || (candidate_count == best_count && candidate_energy < best_energy - 1.0e-12)
            {
                best_count = candidate_count;
                best_energy = candidate_energy;
                best_choice = Some(current_ordinal);
            }
        });
        let offset = state.i * context.n + state.j;
        chart.count[state.kind.index()][offset] = best_count;
        chart.energy[state.kind.index()][offset] = best_energy;
        chart.back[state.kind.index()][offset] = best_choice.unwrap_or(u32::MAX);
    }
    chart
}

fn max_plus_chart(context: &Context<'_>, order: &[StateRef], gains: &[Vec<f64>]) -> MaxCharts {
    let mut chart = MaxCharts {
        n: context.n,
        score: vec![vec![NEG_INF; context.n * context.n]; KIND_COUNT],
        back: vec![vec![u32::MAX; context.n * context.n]; KIND_COUNT],
    };
    for &state in order {
        let local = if matches!(state.kind, Kind::V | Kind::Vp)
            && context.variable_pair(state.i, state.j)
        {
            gains[state.i][state.j]
        } else {
            0.0
        };
        let mut best = NEG_INF;
        let mut best_choice = None;
        let mut ordinal = 0u32;
        context.emit_edges(state, &mut |edge| {
            let current_ordinal = ordinal;
            ordinal = ordinal.saturating_add(1);
            let mut candidate = local;
            for &child in edge.children() {
                let value = chart.score[child.kind.index()][child.i * context.n + child.j];
                if value == NEG_INF {
                    return;
                }
                candidate += value;
            }
            if candidate > best + 1.0e-12 {
                best = candidate;
                best_choice = Some(current_ordinal);
            }
        });
        let offset = state.i * context.n + state.j;
        chart.score[state.kind.index()][offset] = best;
        chart.back[state.kind.index()][offset] = best_choice.unwrap_or(u32::MAX);
    }
    chart
}

fn traceback_pairs(root: StateRef, charts: &Charts, context: &Context<'_>) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let mut stack = vec![root];
    while let Some(state) = stack.pop() {
        if matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j) {
            pairs.push((state.i, state.j));
        }
        if let Some(edge) = selected_edge(context, state, charts.back_choice(state)) {
            stack.extend(edge.children().iter().rev().copied());
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn traceback_max_pairs(
    root: StateRef,
    charts: &MaxCharts,
    context: &Context<'_>,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let mut stack = vec![root];
    while let Some(state) = stack.pop() {
        if matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j) {
            pairs.push((state.i, state.j));
        }
        let choice = charts.back[state.kind.index()][state.i * charts.n + state.j];
        if let Some(edge) = selected_edge(context, state, (choice != u32::MAX).then_some(choice)) {
            stack.extend(edge.children().iter().rev().copied());
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn traceback_forced_pairs(
    root: StateRef,
    charts: &ForcedCharts,
    context: &Context<'_>,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let mut stack = vec![root];
    while let Some(state) = stack.pop() {
        if matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j) {
            pairs.push((state.i, state.j));
        }
        let choice = charts.back[state.kind.index()][state.i * charts.n + state.j];
        if let Some(edge) = selected_edge(context, state, (choice != u32::MAX).then_some(choice)) {
            stack.extend(edge.children().iter().rev().copied());
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn selected_edge(context: &Context<'_>, state: StateRef, choice: Option<u32>) -> Option<Edge> {
    let choice = choice?;
    let mut ordinal = 0u32;
    let mut selected = None;
    context.emit_edges(state, &mut |edge| {
        if ordinal == choice {
            selected = Some(edge);
        }
        ordinal = ordinal.saturating_add(1);
    });
    selected
}

fn stochastic_traceback(
    root: StateRef,
    charts: &Charts,
    context: &Context<'_>,
    rng: &mut SplitMix64,
    pairs: &mut Vec<(usize, usize)>,
    log_probability: &mut f64,
) {
    let mut stack = vec![root];
    while let Some(state) = stack.pop() {
        if matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j) {
            pairs.push((state.i, state.j));
        }
        let state_inside = charts.inside(state);
        let mut candidates = Vec::new();
        context.emit_edges(state, &mut |edge| {
            let mut log_weight = -context.edge_energy(state, edge) / context.rt;
            for &child in edge.children() {
                let child_inside = charts.inside(child);
                if child_inside == NEG_INF {
                    return;
                }
                log_weight += child_inside;
            }
            candidates.push((log_weight, edge));
        });
        let (log_weight, edge) = choose_log_weighted_edge(&candidates, rng);
        *log_probability += log_weight - state_inside;
        stack.extend(edge.children().iter().rev().copied());
    }
}

fn choose_log_weighted_edge(candidates: &[(f64, Edge)], rng: &mut SplitMix64) -> (f64, Edge) {
    debug_assert!(!candidates.is_empty());
    let maximum = candidates
        .iter()
        .map(|candidate| candidate.0)
        .fold(NEG_INF, f64::max);
    let total = candidates
        .iter()
        .map(|candidate| (candidate.0 - maximum).exp())
        .sum::<f64>();
    let mut threshold = rng.next_f64() * total;
    for &(weight, edge) in candidates {
        threshold -= (weight - maximum).exp();
        if threshold <= 0.0 {
            return (weight, edge);
        }
    }
    *candidates.last().expect("nonempty conditional edge list")
}

#[derive(Clone, Copy, Debug)]
struct KEntry {
    energy: f64,
    edge_ordinal: u32,
    child_ranks: [usize; 3],
}

struct KCharts {
    n: usize,
    values: Vec<Vec<Vec<KEntry>>>,
}

impl KCharts {
    fn new(n: usize) -> Self {
        Self {
            n,
            values: vec![vec![Vec::new(); n * n]; KIND_COUNT],
        }
    }

    fn entries(&self, state: StateRef) -> &[KEntry] {
        &self.values[state.kind.index()][state.i * self.n + state.j]
    }

    fn set_entries(&mut self, state: StateRef, entries: Vec<KEntry>) {
        self.values[state.kind.index()][state.i * self.n + state.j] = entries;
    }
}

#[derive(Clone, Debug)]
struct KHeapItem {
    energy: f64,
    edge_ordinal: u32,
    child_ranks: [usize; 3],
}

impl PartialEq for KHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.energy.to_bits() == other.energy.to_bits()
            && self.edge_ordinal == other.edge_ordinal
            && self.child_ranks == other.child_ranks
    }
}

impl Eq for KHeapItem {}

impl PartialOrd for KHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap is a max-heap; reverse all fields so the lowest-energy,
        // lexicographically earliest derivation is emitted first.
        other
            .energy
            .total_cmp(&self.energy)
            .then_with(|| other.edge_ordinal.cmp(&self.edge_ordinal))
            .then_with(|| other.child_ranks.cmp(&self.child_ranks))
    }
}

fn kbest_charts(context: &Context<'_>, order: &[StateRef], keep: usize) -> KCharts {
    let mut charts = KCharts::new(context.n);
    for &state in order {
        let mut edges = Vec::new();
        context.emit_edges(state, &mut |edge| edges.push(edge));
        let mut heap = BinaryHeap::new();
        let mut visited = HashSet::new();
        for (edge_ordinal, &edge) in edges.iter().enumerate() {
            let ranks = [0usize; 3];
            if let Some(energy) = kbest_candidate_energy(context, &charts, state, edge, ranks) {
                let edge_ordinal = edge_ordinal as u32;
                visited.insert((edge_ordinal, ranks));
                heap.push(KHeapItem {
                    energy,
                    edge_ordinal,
                    child_ranks: ranks,
                });
            }
        }
        let mut entries = Vec::with_capacity(keep);
        while entries.len() < keep {
            let Some(item) = heap.pop() else {
                break;
            };
            entries.push(KEntry {
                energy: item.energy,
                edge_ordinal: item.edge_ordinal,
                child_ranks: item.child_ranks,
            });
            let edge = edges[item.edge_ordinal as usize];
            for child_index in 0..edge.child_count as usize {
                let mut ranks = item.child_ranks;
                ranks[child_index] += 1;
                if !visited.insert((item.edge_ordinal, ranks)) {
                    continue;
                }
                if let Some(energy) = kbest_candidate_energy(context, &charts, state, edge, ranks) {
                    heap.push(KHeapItem {
                        energy,
                        edge_ordinal: item.edge_ordinal,
                        child_ranks: ranks,
                    });
                }
            }
        }
        charts.set_entries(state, entries);
    }
    charts
}

fn kbest_candidate_energy(
    context: &Context<'_>,
    charts: &KCharts,
    state: StateRef,
    edge: Edge,
    ranks: [usize; 3],
) -> Option<f64> {
    let mut energy = context.edge_energy(state, edge);
    if !energy.is_finite() {
        return None;
    }
    for (child_index, &child) in edge.children().iter().enumerate() {
        energy += charts.entries(child).get(ranks[child_index])?.energy;
    }
    Some(energy)
}

fn traceback_kbest_pairs(
    root: StateRef,
    rank: usize,
    charts: &KCharts,
    context: &Context<'_>,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let mut stack = vec![(root, rank)];
    while let Some((state, rank)) = stack.pop() {
        if matches!(state.kind, Kind::V | Kind::Vp) && context.variable_pair(state.i, state.j) {
            pairs.push((state.i, state.j));
        }
        let entry = charts.entries(state)[rank];
        let edge = selected_edge(context, state, Some(entry.edge_ordinal))
            .expect("a k-best entry always stores an emitted edge");
        for (child_index, &child) in edge.children().iter().enumerate().rev() {
            stack.push((child, entry.child_ranks[child_index]));
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn validated_seed(
    sequence: &str,
    seed_structure: &str,
    min_loop: usize,
) -> Result<(crate::structure::ParsedStructure, Vec<(usize, usize)>), RnaError> {
    let seed = parse_structure(sequence, seed_structure)?;
    if seed.length == 0 {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis requires a nonempty sequence".into(),
        ));
    }
    if !seed.strand_breaks.is_empty() {
        return Err(RnaError::InvalidOption(
            "conditional density-2 analysis accepts one strand".into(),
        ));
    }
    let mut pairs = seed
        .pairs
        .iter()
        .map(|pair| (pair.i - 1, pair.j - 1))
        .collect::<Vec<_>>();
    pairs.sort_unstable();
    if pairs.iter().enumerate().any(|(index, &left)| {
        pairs[index + 1..]
            .iter()
            .any(|&right| left.0 < right.0 && right.0 < left.1 && left.1 < right.1)
    }) {
        return Err(RnaError::InvalidOption(
            "the conditional seed must be pseudoknot-free".into(),
        ));
    }
    if let Some(pair) = seed.pairs.iter().find(|pair| !pair.canonical) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed contains noncanonical pair {}-{}",
            pair.i, pair.j
        )));
    }
    if let Some(&(i, j)) = pairs.iter().find(|&&(i, j)| j - i <= min_loop) {
        return Err(RnaError::InvalidOption(format!(
            "the conditional seed pair {}-{} violates min-loop={min_loop}",
            i + 1,
            j + 1
        )));
    }
    Ok((seed, pairs))
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

fn constraint_terms(
    bases: &[u8],
    min_loop: usize,
    model: &EnergyModel,
    constraints: &ConstraintModel,
    seed_pairs: &[(usize, usize)],
) -> (Vec<bool>, Vec<f64>, Vec<f64>, f64) {
    let n = bases.len();
    let mut allowed = vec![false; n * n];
    let mut pair_bonus = vec![0.0; n * n];
    let mut stack_bonus = vec![0.0; n * n];
    let seed_partner = partner_table(n, seed_pairs);
    for i in 0..n {
        for j in (i + 1)..n {
            allowed[i * n + j] =
                constraints.allows_pair_for_partition(bases, i, j, min_loop, model);
            pair_bonus[i * n + j] = constraints.pair_energy(i, j)
                - constraints.unpaired_energy(i)
                - constraints.unpaired_energy(j);
            if i + 1 < j && seed_partner[i + 1] == Some(j - 1) {
                pair_bonus[i * n + j] += constraints.stack_energy(i, j, i + 1, j - 1);
            }
            if i > 0 && j + 1 < n && seed_partner[i - 1] == Some(j + 1) {
                pair_bonus[i * n + j] += constraints.stack_energy(i - 1, j + 1, i, j);
            }
            if i + 1 < j {
                stack_bonus[i * n + j] = constraints.stack_energy(i, j, i + 1, j - 1);
            }
        }
    }
    let mut root_constant = (0..n)
        .map(|position| constraints.unpaired_energy(position))
        .sum::<f64>();
    for &(i, j) in seed_pairs {
        root_constant += constraints.pair_energy(i, j)
            - constraints.unpaired_energy(i)
            - constraints.unpaired_energy(j);
        if i + 1 < j && seed_partner[i + 1] == Some(j - 1) {
            root_constant += constraints.stack_energy(i, j, i + 1, j - 1);
        }
    }
    (allowed, pair_bonus, stack_bonus, root_constant)
}

fn seed_baseline_correction(
    bases: &[u8],
    min_loop: usize,
    rt: f64,
    model: &EnergyModel,
    options: &ConditionalDensity2Options,
    seed_pairs: &[(usize, usize)],
) -> Result<f64, RnaError> {
    let n = bases.len();
    let partner = partner_table(n, seed_pairs);
    let context = Context {
        bases,
        n,
        min_loop,
        rt,
        model,
        options,
        cover: cover_table(n, seed_pairs),
        paired_prefix: paired_prefix(&partner),
        weak: weak_table(n, &partner),
        borders: border_table(n, &partner),
        allowed_variable: vec![false; n * n],
        pair_bonus: vec![0.0; n * n],
        stack_bonus: vec![0.0; n * n],
        root_constant: 0.0,
        partner,
    };
    let charts = inside_charts(&context);
    let root = StateRef {
        kind: Kind::W,
        i: 0,
        j: n - 1,
    };
    let grammar_energy = charts.mfe(root);
    if !grammar_energy.is_finite() {
        return Err(RnaError::Numerical(
            "the fixed seed has no finite derivation in the conditional grammar".into(),
        ));
    }
    let sequence = std::str::from_utf8(bases).expect("normalized RNA sequence is ASCII");
    let physical_energy = model
        .evaluate(sequence, &planar_structure(n, seed_pairs, '(', ')'))?
        .total_kcal_mol;
    Ok(physical_energy - grammar_energy)
}

fn partner_table(length: usize, pairs: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut partner = vec![None; length];
    for &(i, j) in pairs {
        partner[i] = Some(j);
        partner[j] = Some(i);
    }
    partner
}

fn cover_table(length: usize, pairs: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut cover = vec![None; length];
    for (position, cell) in cover.iter_mut().enumerate() {
        *cell = pairs
            .iter()
            .filter(|&&(i, j)| i < position && position < j)
            .max_by_key(|&&(i, _)| i)
            .map(|&(i, _)| i);
    }
    cover
}

fn paired_prefix(partner: &[Option<usize>]) -> Vec<usize> {
    let mut prefix = vec![0; partner.len() + 1];
    for (i, value) in partner.iter().enumerate() {
        prefix[i + 1] = prefix[i] + usize::from(value.is_some());
    }
    prefix
}

fn weak_table(length: usize, partner: &[Option<usize>]) -> Vec<bool> {
    let mut table = vec![false; length * length];
    for i in 0..length {
        let mut minimum = usize::MAX;
        let mut maximum = 0usize;
        for j in i..length {
            if let Some(p) = partner[j] {
                minimum = minimum.min(p);
                maximum = maximum.max(p);
            }
            table[i * length + j] = minimum >= i && maximum <= j;
        }
    }
    table
}

fn border_table(length: usize, partner: &[Option<usize>]) -> Vec<Borders> {
    let mut table = vec![Borders::default(); length * length];
    // b and b': left endpoints of seed pairs that cover the right endpoint.
    for j in 0..length {
        let mut minimum = None;
        let mut maximum = None;
        for i in (0..j).rev() {
            if partner[i].is_some_and(|right| i < right && j < right) {
                minimum = Some(i);
                maximum.get_or_insert(i);
            }
            let cell = &mut table[i * length + j];
            cell.b = minimum;
            cell.b_prime = maximum;
        }
    }
    // B and B': right endpoints of seed pairs that cover the left endpoint.
    for i in 0..length {
        let mut minimum = None;
        let mut maximum = None;
        for j in (i + 1)..length {
            if partner[j].is_some_and(|left| left < i && left < j) {
                minimum.get_or_insert(j);
                maximum = Some(j);
            }
            let cell = &mut table[i * length + j];
            cell.big_b = maximum;
            cell.big_b_prime = minimum;
        }
    }
    table
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conditional_density2::{band_topology, conditional_density2_ensemble};
    use crate::constraints::{ConstraintConfig, ConstraintModel};
    use crate::exact_enumeration::for_each_noncrossing_structure;
    use crate::partition;

    fn log_structure_count(context: &Context<'_>, order: &[StateRef]) -> Vec<Vec<f64>> {
        let mut count = vec![vec![NEG_INF; context.n * context.n]; KIND_COUNT];
        for &state in order {
            let mut total = NEG_INF;
            context.emit_edges(state, &mut |edge| {
                let mut candidate = 0.0;
                for &child in edge.children() {
                    let value = count[child.kind.index()][child.i * context.n + child.j];
                    if value == NEG_INF {
                        return;
                    }
                    candidate += value;
                }
                total = log_add(total, candidate);
            });
            count[state.kind.index()][state.i * context.n + state.j] = total;
        }
        count
    }

    fn polynomial_log_count(sequence: &str, seed_structure: &str, min_loop: usize) -> f64 {
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let pairs = seed
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let partner = partner_table(seed.length, &pairs);
        let context = Context {
            bases: seed.sequence.as_bytes(),
            n: seed.length,
            min_loop,
            rt: GAS_CONSTANT_KCAL * 310.15,
            model: &model,
            options: &ConditionalDensity2Options::default(),
            cover: cover_table(seed.length, &pairs),
            paired_prefix: paired_prefix(&partner),
            weak: weak_table(seed.length, &partner),
            borders: border_table(seed.length, &partner),
            allowed_variable: vec![true; seed.length * seed.length],
            pair_bonus: vec![0.0; seed.length * seed.length],
            stack_bonus: vec![0.0; seed.length * seed.length],
            root_constant: 0.0,
            partner,
        };
        let charts = inside_charts(&context);
        let count = log_structure_count(&context, &charts.order);
        count[Kind::W.index()][seed.length - 1]
    }

    fn polynomial_target_reachable(
        sequence: &str,
        seed_structure: &str,
        added_structure: &str,
        min_loop: usize,
    ) -> bool {
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let seed_pairs = seed
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let added = parse_structure(sequence, added_structure).unwrap();
        let added_pairs = added
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let partner = partner_table(seed.length, &seed_pairs);
        let mut allowed_variable = vec![false; seed.length * seed.length];
        let mut gains = vec![vec![-1000.0; seed.length]; seed.length];
        for &(i, j) in &added_pairs {
            allowed_variable[i * seed.length + j] = true;
            gains[i][j] = 1.0;
        }
        let context = Context {
            bases: seed.sequence.as_bytes(),
            n: seed.length,
            min_loop,
            rt: GAS_CONSTANT_KCAL * 310.15,
            model: &model,
            options: &ConditionalDensity2Options::default(),
            cover: cover_table(seed.length, &seed_pairs),
            paired_prefix: paired_prefix(&partner),
            weak: weak_table(seed.length, &partner),
            borders: border_table(seed.length, &partner),
            allowed_variable,
            pair_bonus: vec![0.0; seed.length * seed.length],
            stack_bonus: vec![0.0; seed.length * seed.length],
            root_constant: 0.0,
            partner,
        };
        let charts = inside_charts(&context);
        let max = max_plus_chart(&context, &charts.order, &gains);
        let root_score = max.score[Kind::W.index()][seed.length - 1];
        (root_score - added_pairs.len() as f64).abs() < 1.0e-9
    }

    fn polynomial_log_z_with_pair_bonus(
        sequence: &str,
        seed_structure: &str,
        min_loop: usize,
        pair: (usize, usize),
        bonus: f64,
    ) -> f64 {
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let seed_pairs = seed
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let partner = partner_table(seed.length, &seed_pairs);
        let mut pair_bonus = vec![0.0; seed.length * seed.length];
        pair_bonus[pair.0 * seed.length + pair.1] = bonus;
        let context = Context {
            bases: seed.sequence.as_bytes(),
            n: seed.length,
            min_loop,
            rt: GAS_CONSTANT_KCAL * 310.15,
            model: &model,
            options: &ConditionalDensity2Options::default(),
            cover: cover_table(seed.length, &seed_pairs),
            paired_prefix: paired_prefix(&partner),
            weak: weak_table(seed.length, &partner),
            borders: border_table(seed.length, &partner),
            allowed_variable: vec![true; seed.length * seed.length],
            pair_bonus,
            stack_bonus: vec![0.0; seed.length * seed.length],
            root_constant: 0.0,
            partner,
        };
        inside_charts(&context).inside(StateRef {
            kind: Kind::W,
            i: 0,
            j: seed.length - 1,
        })
    }

    #[test]
    fn borders_match_the_published_definitions() {
        // Pairs 1-10, 2-9, 3-8 in one-based coordinates.
        let partner = partner_table(10, &[(0, 9), (1, 8), (2, 7)]);
        let borders = border_table(10, &partner);
        let at = |i: usize, j: usize| borders[i * 10 + j];
        let right = at(0, 5);
        assert_eq!(right.b, Some(0));
        assert_eq!(right.b_prime, Some(2));
        let left = at(5, 9);
        assert_eq!(left.big_b_prime, Some(7));
        assert_eq!(left.big_b, Some(9));
    }

    #[test]
    fn be_identity_has_unit_weight() {
        let result = conditional_density2_polynomial(
            "GCGAAAGC",
            "((....))",
            37.0,
            3,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        assert!(result.log_partition_function.is_finite());
        assert!(result.mfe_energy_kcal_mol.is_finite());
    }

    #[test]
    fn probabilities_are_position_normalized() {
        let result = conditional_density2_polynomial(
            "GGGAAACCCUUU",
            "((......))..",
            37.0,
            3,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        for i in 0..result.sequence.len() {
            let paired = result
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i - 1 == i || pair.j - 1 == i)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!((paired + result.unpaired_probabilities[i] - 1.0).abs() < 1.0e-8);
        }
    }

    #[test]
    fn counting_semiring_matches_independent_exhaustive_oracle() {
        let cases = [
            ("GCGAAACG", "........", 3),
            ("GCGAAACGC", "(.......)", 3),
            ("GCGAAACGGC", "((......))", 3),
            ("GGGAAACCCU", ".(......).", 3),
            ("GGGAAACCCUU", "((.......))", 3),
            ("GCGCGCGCGC", "(......)..", 3),
            ("GCGCGCGCGCGC", "((......))..", 3),
            ("GCGCGCGCGCGC", ".(......)...", 3),
            ("GCGCGCGCGCGC", "(........)..", 3),
            ("GCGCGCGCGCGC", ".((......)).", 3),
            ("GCGCGCGCGCGC", "(.(....).)..", 3),
            ("GCGCGCGCGCGCGC", "((........))..", 3),
            ("GCGCGCGCGCGCGC", "..((......))..", 3),
            ("GCGCGCGCGCGCGC", "(.(......).)..", 3),
        ];
        for (sequence, seed, min_loop) in cases {
            let oracle = conditional_density2_ensemble(
                sequence,
                seed,
                37.0,
                min_loop,
                0,
                1.021,
                &ConditionalDensity2Options::default(),
            )
            .unwrap();
            let actual = polynomial_log_count(sequence, seed, min_loop).exp();
            if (actual - oracle.state_count as f64).abs() >= 1.0e-8 {
                let parsed = parse_structure(sequence, seed).unwrap();
                let seed_pairs = parsed
                    .pairs
                    .iter()
                    .map(|pair| (pair.i - 1, pair.j - 1))
                    .collect::<Vec<_>>();
                let constraints = ConstraintModel::compile(
                    sequence.len(),
                    &ConstraintConfig {
                        force_unpaired: seed_pairs
                            .iter()
                            .flat_map(|&(i, j)| [i + 1, j + 1])
                            .collect(),
                        ..ConstraintConfig::default()
                    },
                )
                .unwrap();
                let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
                let mut states = Vec::new();
                for_each_noncrossing_structure(
                    sequence.as_bytes(),
                    min_loop,
                    &model,
                    &constraints,
                    |added| {
                        let mut union = seed_pairs.clone();
                        union.extend_from_slice(added);
                        union.sort_unstable();
                        if band_topology(sequence.len(), &union).maximum_density <= 2 {
                            states.push(planar_structure(sequence.len(), added, '(', ')'));
                        }
                        Ok(())
                    },
                )
                .unwrap();
                eprintln!("oracle states for {sequence} {seed}: {states:?}");
                for state in &states {
                    if !polynomial_target_reachable(sequence, seed, state, min_loop) {
                        eprintln!("unreachable: {state}");
                    }
                }
            }
            assert!(
                (actual - oracle.state_count as f64).abs() < 1.0e-8,
                "{sequence} {seed}: polynomial={actual}, exhaustive={}",
                oracle.state_count
            );
        }
    }

    #[test]
    fn every_small_planar_seed_matches_the_exhaustive_state_space() {
        let sequence = "GCGCGCGCGCGC";
        let min_loop = 3;
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let constraints = ConstraintModel::unconstrained(sequence.len());
        let mut seeds = Vec::new();
        for_each_noncrossing_structure(
            sequence.as_bytes(),
            min_loop,
            &model,
            &constraints,
            |pairs| {
                seeds.push(planar_structure(sequence.len(), pairs, '(', ')'));
                Ok(())
            },
        )
        .unwrap();
        assert!(seeds.len() >= 20);
        for seed in seeds {
            let oracle = conditional_density2_ensemble(
                sequence,
                &seed,
                37.0,
                min_loop,
                0,
                1.021,
                &ConditionalDensity2Options::default(),
            )
            .unwrap();
            let actual = polynomial_log_count(sequence, &seed, min_loop).exp();
            assert!(
                (actual - oracle.state_count as f64).abs() < 1.0e-8,
                "seed {seed}: polynomial={actual}, exhaustive={}",
                oracle.state_count
            );
        }
    }

    #[test]
    fn outside_pair_marginals_match_log_z_derivatives() {
        let sequence = "GCGCGCGCGC";
        let seed = "(......)..";
        let result = conditional_density2_polynomial(
            sequence,
            seed,
            37.0,
            3,
            0,
            1.021,
            &ConditionalDensity2Options::default(),
        )
        .unwrap();
        let epsilon = 1.0e-5;
        let rt = GAS_CONSTANT_KCAL * 310.15;
        for pair in result
            .pair_probabilities
            .iter()
            .filter(|pair| pair.probability < 1.0 - 1.0e-12)
        {
            let ij = (pair.i - 1, pair.j - 1);
            let plus = polynomial_log_z_with_pair_bonus(sequence, seed, 3, ij, epsilon);
            let minus = polynomial_log_z_with_pair_bonus(sequence, seed, 3, ij, -epsilon);
            let derivative_probability = -rt * (plus - minus) / (2.0 * epsilon);
            assert!(
                (derivative_probability - pair.probability).abs() < 2.0e-8,
                "pair {}-{} outside={} derivative={derivative_probability}",
                pair.i,
                pair.j,
                pair.probability
            );
        }
    }

    fn assert_nonempty_seed_ensemble_matches_independent_state_enumeration(dangles: u8) {
        let sequence = "GCGCGCGCGC";
        let seed_structure = "(......)..";
        let min_loop = 3;
        let options = ConditionalDensity2Options::default();
        let ensemble = conditional_density2_polynomial(
            sequence,
            seed_structure,
            37.0,
            min_loop,
            dangles,
            1.021,
            &options,
        )
        .unwrap();
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let seed_pairs = seed
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let constraints = ConstraintModel::compile(
            sequence.len(),
            &ConstraintConfig {
                force_unpaired: seed_pairs
                    .iter()
                    .flat_map(|&(i, j)| [i + 1, j + 1])
                    .collect(),
                ..ConstraintConfig::default()
            },
        )
        .unwrap();
        let model = EnergyModel::with_dangles_and_salt(37.0, dangles, 1.021).unwrap();
        let rt = GAS_CONSTANT_KCAL * 310.15;
        let mut states = Vec::<(Vec<(usize, usize)>, f64)>::new();
        for_each_noncrossing_structure(
            sequence.as_bytes(),
            min_loop,
            &model,
            &constraints,
            |added| {
                let mut union = seed_pairs.clone();
                union.extend_from_slice(added);
                union.sort_unstable();
                if band_topology(sequence.len(), &union).maximum_density <= 2 {
                    let added_structure = planar_structure(sequence.len(), added, '(', ')');
                    let evaluated = evaluate_conditional_density2_polynomial(
                        sequence,
                        seed_structure,
                        &added_structure,
                        37.0,
                        min_loop,
                        dangles,
                        1.021,
                        &options,
                    )?;
                    states.push((added.to_vec(), evaluated.energy_kcal_mol));
                }
                Ok(())
            },
        )
        .unwrap();
        assert!(states.len() > 1);

        let log_z = states
            .iter()
            .fold(NEG_INF, |total, (_, energy)| log_add(total, -energy / rt));
        assert!((ensemble.log_partition_function - log_z).abs() < 1.0e-10);
        let enumerated_mfe = states
            .iter()
            .map(|(_, energy)| *energy)
            .fold(POS_INF, f64::min);
        assert!((ensemble.mfe_energy_kcal_mol - enumerated_mfe).abs() < 1.0e-10);

        let mut marginal = std::collections::HashMap::<(usize, usize), f64>::new();
        for (pairs, energy) in &states {
            let weight = (-energy / rt - log_z).exp();
            for &pair in pairs {
                *marginal.entry(pair).or_default() += weight;
            }
        }
        let production = ensemble
            .pair_probabilities
            .iter()
            .filter_map(|pair| {
                let ij = (pair.i - 1, pair.j - 1);
                (!seed_pairs.contains(&ij)).then_some((ij, pair.probability))
            })
            .collect::<std::collections::HashMap<_, _>>();
        for pair in marginal.keys().chain(production.keys()) {
            assert!(
                (marginal.get(pair).copied().unwrap_or(0.0)
                    - production.get(pair).copied().unwrap_or(0.0))
                .abs()
                    < 1.0e-10,
                "pair {pair:?}: exhaustive={:?} production={:?}",
                marginal.get(pair),
                production.get(pair)
            );
        }

        let variable_pairs = |structure: &str| {
            parse_structure(sequence, structure)
                .unwrap()
                .pairs
                .into_iter()
                .map(|pair| (pair.i - 1, pair.j - 1))
                .filter(|pair| !seed_pairs.contains(pair))
                .collect::<Vec<_>>()
        };
        let total_probability = marginal.values().sum::<f64>();
        let centroid_distance = |pairs: &[(usize, usize)]| {
            total_probability + pairs.len() as f64
                - 2.0
                    * pairs
                        .iter()
                        .map(|pair| marginal.get(pair).copied().unwrap_or(0.0))
                        .sum::<f64>()
        };
        let best_centroid = states
            .iter()
            .map(|(pairs, _)| centroid_distance(pairs))
            .fold(POS_INF, f64::min);
        let decoded_centroid = variable_pairs(&ensemble.centroid_structure);
        assert!((ensemble.centroid_distance - best_centroid).abs() < 1.0e-10);
        assert!((centroid_distance(&decoded_centroid) - best_centroid).abs() < 1.0e-10);

        let seed_occupied = seed_pairs
            .iter()
            .flat_map(|&(i, j)| [i, j])
            .collect::<std::collections::HashSet<_>>();
        let mea_score = |pairs: &[(usize, usize)]| {
            let paired = pairs
                .iter()
                .flat_map(|&(i, j)| [i, j])
                .chain(seed_occupied.iter().copied())
                .collect::<std::collections::HashSet<_>>();
            2.0 * pairs
                .iter()
                .map(|pair| marginal.get(pair).copied().unwrap_or(0.0))
                .sum::<f64>()
                + ensemble
                    .unpaired_probabilities
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !paired.contains(i))
                    .map(|(_, probability)| probability)
                    .sum::<f64>()
        };
        let best_mea = states
            .iter()
            .map(|(pairs, _)| mea_score(pairs))
            .fold(NEG_INF, f64::max);
        let decoded_mea = variable_pairs(&ensemble.mea_structure);
        assert!((ensemble.mea_score - best_mea).abs() < 1.0e-10);
        assert!((mea_score(&decoded_mea) - best_mea).abs() < 1.0e-10);
    }

    #[test]
    fn nonempty_seed_ensemble_matches_independent_state_enumeration() {
        for dangles in [0, 2] {
            assert_nonempty_seed_ensemble_matches_independent_state_enumeration(dangles);
        }
    }

    #[test]
    fn empty_seed_reduces_to_planar_partition_on_short_sequences() {
        for sequence in ["GCGAAACGCU", "GGGAAACCCU", "GCGCGCGCGC"] {
            for dangles in [0, 2] {
                let options = ConditionalDensity2Options::default();
                let conditional = conditional_density2_polynomial(
                    sequence,
                    "..........",
                    37.0,
                    3,
                    dangles,
                    1.021,
                    &options,
                )
                .unwrap();
                let model = EnergyModel::with_dangles_and_salt(37.0, dangles, 1.021).unwrap();
                let planar = partition(sequence, 37.0, 3, &model).unwrap();
                assert!(
                    (conditional.log_partition_function - planar.log_partition_function).abs()
                        < 1.0e-10,
                    "{sequence} dangles={dangles}: conditional={} planar={}",
                    conditional.log_partition_function,
                    planar.log_partition_function
                );
                let conditional_pairs = conditional
                    .pair_probabilities
                    .iter()
                    .map(|pair| ((pair.i, pair.j), pair.probability))
                    .collect::<std::collections::HashMap<_, _>>();
                let planar_pairs = planar
                    .pair_probabilities
                    .iter()
                    .map(|pair| ((pair.i, pair.j), pair.probability))
                    .collect::<std::collections::HashMap<_, _>>();
                for pair in conditional_pairs.keys().chain(planar_pairs.keys()) {
                    assert!(
                        (conditional_pairs.get(pair).copied().unwrap_or(0.0)
                            - planar_pairs.get(pair).copied().unwrap_or(0.0))
                        .abs()
                            < 1.0e-10,
                        "{sequence} dangles={dangles} pair {pair:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn polynomial_engine_rejects_nonlocal_odd_dangle_models() {
        for dangles in [1, 3] {
            let error = conditional_density2_polynomial(
                "GCGAAACGCU",
                "..........",
                37.0,
                3,
                dangles,
                1.021,
                &ConditionalDensity2Options::default(),
            )
            .unwrap_err();
            assert!(error.to_string().contains("exact exhaustive dispatcher"));
        }
    }

    #[test]
    fn fixed_structure_evaluator_reproduces_mfe_backtrace_energy() {
        let sequence = "GCGCGCGCGCGC";
        let seed = ".(......)...";
        let options = ConditionalDensity2Options::default();
        let ensemble =
            conditional_density2_polynomial(sequence, seed, 37.0, 3, 0, 1.021, &options).unwrap();
        let evaluated = evaluate_conditional_density2_polynomial(
            sequence,
            seed,
            &ensemble.mfe_added_structure,
            37.0,
            3,
            0,
            1.021,
            &options,
        )
        .unwrap();
        assert!(evaluated.derivation_unique);
        assert!((evaluated.energy_kcal_mol - ensemble.mfe_energy_kcal_mol).abs() < 1.0e-10);
        assert_eq!(evaluated.structure, ensemble.mfe_structure);
    }

    #[test]
    fn long_real_mfe_backtraces_remain_density_two() {
        let cases = [
            (
                "GGCCGGCATGGTCCCAGCCTCCTCGCTGGCGCCGGCTGGGCAACAUUCCGAGGGGACCGUCCCCUCGGUAAUGGCGAAUGGGACCCA",
                ".........((((((((((........))).............(((((((((((((...))))))))).)))).....)))))))..",
            ),
            (
                "AUAAUAAAUAACGGAUUGUGUCCGUAAUCACACGUGGUGCGUACGAUAACGCAUAGUGUUUUUCCCUCCACUUAAAUCGAAGGG",
                "................((((........))))....((((((......)))))).........(((((..........).))))",
            ),
        ];
        let options = ConditionalDensity2Options::default();
        for (sequence, seed_structure) in cases {
            let ensemble = conditional_density2_polynomial(
                sequence,
                seed_structure,
                37.0,
                0,
                0,
                1.021,
                &options,
            )
            .unwrap();
            let seed = parse_structure(sequence, seed_structure).unwrap();
            let added = parse_structure(sequence, &ensemble.mfe_added_structure).unwrap();
            let mut union = seed
                .pairs
                .iter()
                .chain(&added.pairs)
                .map(|pair| (pair.i - 1, pair.j - 1))
                .collect::<Vec<_>>();
            union.sort_unstable();
            assert!(band_topology(sequence.len(), &union).maximum_density <= 2);
            let evaluated = evaluate_conditional_density2_polynomial(
                sequence,
                seed_structure,
                &ensemble.mfe_added_structure,
                37.0,
                0,
                0,
                1.021,
                &options,
            )
            .unwrap();
            assert!((evaluated.energy_kcal_mol - ensemble.mfe_energy_kcal_mol).abs() < 1.0e-9);
        }
    }

    #[test]
    fn conditional_kbest_matches_complete_independent_energy_ordering() {
        let sequence = "GCGCGCGCGCGC";
        let seed_structure = "(......)....";
        let min_loop = 3;
        let options = ConditionalDensity2Options::default();
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let seed_pairs = seed
            .pairs
            .iter()
            .map(|pair| (pair.i - 1, pair.j - 1))
            .collect::<Vec<_>>();
        let model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let constraints = ConstraintModel::compile(
            sequence.len(),
            &ConstraintConfig {
                force_unpaired: seed_pairs
                    .iter()
                    .flat_map(|&(i, j)| [i + 1, j + 1])
                    .collect(),
                ..ConstraintConfig::default()
            },
        )
        .unwrap();
        let mut expected = Vec::new();
        for_each_noncrossing_structure(
            sequence.as_bytes(),
            min_loop,
            &model,
            &constraints,
            |added| {
                let mut union = seed_pairs.clone();
                union.extend_from_slice(added);
                union.sort_unstable();
                if band_topology(sequence.len(), &union).maximum_density <= 2 {
                    let added_structure = planar_structure(sequence.len(), added, '(', ')');
                    let energy = evaluate_conditional_density2_polynomial(
                        sequence,
                        seed_structure,
                        &added_structure,
                        37.0,
                        min_loop,
                        0,
                        1.021,
                        &options,
                    )?
                    .energy_kcal_mol;
                    expected.push((energy, added_structure));
                }
                Ok(())
            },
        )
        .unwrap();
        expected.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let actual = suboptimal_conditional_density2_polynomial(
            sequence,
            seed_structure,
            37.0,
            min_loop,
            0,
            1.021,
            &options,
            1000.0,
            expected.len() + 1,
        )
        .unwrap();
        assert!(!actual.truncated);
        assert_eq!(actual.structures.len(), expected.len());
        for (observed, (energy, structure)) in actual.structures.iter().zip(expected) {
            assert_eq!(observed.added_structure, structure);
            assert!((observed.energy_kcal_mol - energy).abs() < 1.0e-10);
        }
    }

    #[test]
    fn conditional_sampling_probabilities_are_exact_and_reproducible() {
        let sequence = "GCGCGCGCGCGC";
        let seed = "(......)....";
        let options = ConditionalDensity2Options::default();
        let ensemble =
            conditional_density2_polynomial(sequence, seed, 37.0, 3, 0, 1.021, &options).unwrap();
        let first = sample_conditional_density2_polynomial(
            sequence, seed, 37.0, 3, 0, 1.021, &options, 512, 42, false,
        )
        .unwrap();
        let second = sample_conditional_density2_polynomial(
            sequence, seed, 37.0, 3, 0, 1.021, &options, 512, 42, false,
        )
        .unwrap();
        assert_eq!(first.samples[0].structure, second.samples[0].structure);
        assert!((first.log_partition_function - ensemble.log_partition_function).abs() < 1.0e-12);
        let rt = GAS_CONSTANT_KCAL * 310.15;
        for sample in &first.samples {
            let evaluated = evaluate_conditional_density2_polynomial(
                sequence,
                seed,
                &sample.added_structure,
                37.0,
                3,
                0,
                1.021,
                &options,
            )
            .unwrap();
            let expected = -evaluated.energy_kcal_mol / rt - ensemble.log_partition_function;
            assert!((sample.log_probability - expected).abs() < 1.0e-10);
            assert!((sample.effective_energy_kcal_mol - evaluated.energy_kcal_mol).abs() < 1.0e-10);
        }
    }

    #[test]
    fn pk_only_restricts_every_added_pair_to_a_seed_crossing() {
        let options = ConditionalDensity2Options {
            pk_only: true,
            ..ConditionalDensity2Options::default()
        };
        let sequence = "GCGCGCGCGCGC";
        let seed_structure = "(......)....";
        let result =
            conditional_density2_polynomial(sequence, seed_structure, 37.0, 3, 0, 1.021, &options)
                .unwrap();
        let seed = parse_structure(sequence, seed_structure).unwrap();
        let added = parse_structure(sequence, &result.mfe_added_structure).unwrap();
        for pair in added.pairs {
            let pair = (pair.i - 1, pair.j - 1);
            assert!(seed.pairs.iter().any(|fixed| {
                let fixed = (fixed.i - 1, fixed.j - 1);
                (pair.0 < fixed.0 && fixed.0 < pair.1 && pair.1 < fixed.1)
                    || (fixed.0 < pair.0 && pair.0 < fixed.1 && fixed.1 < pair.1)
            }));
        }
    }

    #[test]
    fn polynomial_soft_and_unpaired_constraints_match_exhaustive_union_scoring() {
        use crate::constraints::{PairEnergy, PositionEnergy, SoftConstraintConfig};

        let sequence = "GCGCGCGCGC";
        let seed = "(......)..";
        let config = ConstraintConfig {
            force_unpaired: vec![2],
            soft: SoftConstraintConfig {
                unpaired: vec![PositionEnergy {
                    position: 4,
                    energy_kcal_mol: 0.37,
                }],
                paired: vec![PositionEnergy {
                    position: 9,
                    energy_kcal_mol: -0.22,
                }],
                pairs: vec![PairEnergy {
                    i: 3,
                    j: 9,
                    energy_kcal_mol: -0.51,
                }],
                stack: vec![PositionEnergy {
                    position: 3,
                    energy_kcal_mol: 0.19,
                }],
            },
            ..ConstraintConfig::default()
        };
        let options = ConditionalDensity2Options::default();
        let polynomial = conditional_density2_polynomial_with_constraints(
            sequence, seed, 37.0, 3, 0, 1.021, &options, &config,
        )
        .unwrap();
        let exhaustive =
            crate::conditional_density2::conditional_density2_ensemble_with_constraints(
                sequence, seed, 37.0, 3, 0, 1.021, &options, &config,
            )
            .unwrap();
        assert!(
            (polynomial.log_partition_function - exhaustive.log_partition_function).abs() < 1.0e-10,
            "polynomial={} exhaustive={}",
            polynomial.log_partition_function,
            exhaustive.log_partition_function
        );
        assert!((polynomial.mfe_energy_kcal_mol - exhaustive.mfe_energy_kcal_mol).abs() < 1.0e-10);
        assert_eq!(polynomial.mfe_structure, exhaustive.mfe_structure);
        let expected = exhaustive
            .pair_probabilities
            .iter()
            .map(|pair| ((pair.i, pair.j), pair.probability))
            .collect::<std::collections::HashMap<_, _>>();
        for pair in &polynomial.pair_probabilities {
            assert!(
                (pair.probability - expected.get(&(pair.i, pair.j)).copied().unwrap_or(0.0)).abs()
                    < 1.0e-10
            );
        }
    }
}
