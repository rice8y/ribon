//! Alignment-aware covariation-guided consensus folding.

use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::energy::EnergyModel;
use crate::ensemble::summarize;
use crate::partition::{PairProbability, PartitionResult};
use crate::structure::{pairs_to_dot_bracket, RnaError};
use crate::{AnalysisResult, ModelDescription};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const INF: f64 = 1.0e100;
const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct ComparativeOptions {
    pub covariance_weight_kcal_mol: f64,
    pub incompatible_penalty: f64,
    pub minimum_pair_occupancy: f64,
}

impl Default for ComparativeOptions {
    fn default() -> Self {
        Self {
            covariance_weight_kcal_mol: 1.0,
            incompatible_penalty: 1.0,
            minimum_pair_occupancy: 0.5,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct CovariationEntry {
    pub i: usize,
    pub j: usize,
    pub score: f64,
    pub mutual_information_bits: f64,
    pub canonical_fraction: f64,
    pub compensatory_changes: usize,
    pub consistent_changes: usize,
    pub incompatible_sequences: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct ComparativeResult {
    pub alignment: Vec<String>,
    pub sequence_count: usize,
    pub alignment_length: usize,
    pub consensus_sequence: String,
    pub consensus_structure: String,
    pub consensus_energy_kcal_mol: f64,
    pub analysis: AnalysisResult,
    pub covariation: Vec<CovariationEntry>,
    pub covariance_weight_kcal_mol: f64,
    pub model: &'static str,
}

#[allow(clippy::too_many_arguments)]
pub fn comparative_fold(
    alignment: &[String],
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
    options: &ComparativeOptions,
) -> Result<ComparativeResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    comparative_fold_with_model(alignment, min_loop, gamma, &model, options)
}

pub fn comparative_fold_with_model(
    alignment: &[String],
    min_loop: usize,
    gamma: f64,
    model: &EnergyModel,
    options: &ComparativeOptions,
) -> Result<ComparativeResult, RnaError> {
    validate_options(options)?;
    let alignment = normalize_alignment(alignment)?;
    let length = alignment[0].len();
    let consensus = consensus_sequence(&alignment)?;
    let mut covariation = Vec::new();
    let mut covariance_energies = vec![vec![0.0; length]; length];
    let mut allowed_pairs = vec![vec![false; length]; length];
    for i in 0..length {
        for j in i + min_loop + 1..length {
            let entry = covariation_entry(&alignment, i, j, options.incompatible_penalty);
            if entry.canonical_fraction >= options.minimum_pair_occupancy {
                allowed_pairs[i][j] = true;
                covariance_energies[i][j] = -options.covariance_weight_kcal_mol * entry.score;
            }
            if entry.canonical_fraction > 0.0 || entry.score != 0.0 {
                covariation.push(entry);
            }
        }
    }
    let analysis = analyze_alignment(
        &alignment,
        &consensus,
        &allowed_pairs,
        &covariance_energies,
        min_loop,
        gamma,
        model,
    )?;
    let display_sequence = |sequence: &str| {
        if model.nucleic_acid() == crate::energy::NucleicAcid::Dna {
            sequence.replace('U', "T")
        } else {
            sequence.to_owned()
        }
    };
    Ok(ComparativeResult {
        sequence_count: alignment.len(),
        alignment_length: length,
        consensus_sequence: display_sequence(&consensus),
        consensus_structure: analysis.mfe_structure.clone(),
        consensus_energy_kcal_mol: analysis.mfe_energy_kcal_mol,
        analysis,
        alignment: alignment.iter().map(|row| display_sequence(row)).collect(),
        covariation,
        covariance_weight_kcal_mol: options.covariance_weight_kcal_mol,
        model: if model.parameter_profile_name().is_some() {
            "alignment-column custom loop-energy average with covariation pseudo-energy"
        } else {
            match model.nucleic_acid() {
                crate::energy::NucleicAcid::Rna => {
                    "alignment-column RNAstructure 6.6 RNA loop-energy average with covariation pseudo-energy"
                }
                crate::energy::NucleicAcid::Dna => {
                    "alignment-column RNAstructure 6.6 DNA loop-energy average with covariation pseudo-energy"
                }
            }
        },
    })
}

#[derive(Clone)]
struct AlignmentRow {
    bases: Vec<u8>,
    column_to_position: Vec<Option<usize>>,
}

struct AlignmentScorer<'a> {
    rows: Vec<AlignmentRow>,
    model: &'a EnergyModel,
    allowed_pairs: &'a [Vec<bool>],
    covariance_energies: &'a [Vec<f64>],
}

impl AlignmentScorer<'_> {
    fn allows_pair(&self, i: usize, j: usize) -> bool {
        self.allowed_pairs[i][j]
    }

    fn covariance(&self, i: usize, j: usize) -> f64 {
        self.covariance_energies[i][j]
    }

    fn average_pair_energy(
        &self,
        i: usize,
        j: usize,
        mut evaluate: impl FnMut(&EnergyModel, &[u8], usize, usize) -> f64,
    ) -> f64 {
        let mut total = 0.0;
        let mut supported = 0usize;
        for row in &self.rows {
            let (Some(a), Some(b)) = (row.column_to_position[i], row.column_to_position[j]) else {
                continue;
            };
            if a >= b || !self.model.can_pair(row.bases[a], row.bases[b]) {
                continue;
            }
            let energy = evaluate(self.model, &row.bases, a, b);
            if energy.is_finite() {
                total += energy;
                supported += 1;
            }
        }
        if supported == 0 {
            INF
        } else {
            total / self.rows.len() as f64
        }
    }

    fn average_internal_energy(
        &self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        boltzmann: bool,
    ) -> f64 {
        let mut total = 0.0;
        let mut supported = 0usize;
        for row in &self.rows {
            let (Some(a), Some(b), Some(c), Some(d)) = (
                row.column_to_position[i],
                row.column_to_position[j],
                row.column_to_position[k],
                row.column_to_position[l],
            ) else {
                continue;
            };
            if !(a < c
                && c < d
                && d < b
                && self.model.can_pair(row.bases[a], row.bases[b])
                && self.model.can_pair(row.bases[c], row.bases[d]))
            {
                continue;
            }
            let energy = if boltzmann {
                self.model.internal_boltzmann_energy(&row.bases, a, b, c, d)
            } else {
                self.model.internal_energy(&row.bases, a, b, c, d)
            };
            if energy.is_finite() {
                total += energy;
                supported += 1;
            }
        }
        if supported == 0 {
            INF
        } else {
            total / self.rows.len() as f64
        }
    }

    fn hairpin(&self, i: usize, j: usize, boltzmann: bool) -> f64 {
        self.average_pair_energy(i, j, |model, bases, a, b| {
            if boltzmann {
                model.hairpin_boltzmann_energy(bases, a, b)
            } else {
                model.hairpin_energy(bases, a, b)
            }
        })
    }

    fn internal(&self, i: usize, j: usize, k: usize, l: usize, boltzmann: bool) -> f64 {
        self.average_internal_energy(i, j, k, l, boltzmann)
    }

    fn stem(&self, i: usize, j: usize, context: StemContext, boltzmann: bool) -> f64 {
        self.average_pair_energy(i, j, |model, bases, a, b| match (context, boltzmann) {
            (StemContext::Exterior, false) => model.exterior_stem_energy(bases, a, b),
            (StemContext::Exterior, true) => model.exterior_stem_boltzmann_energy(bases, a, b),
            (StemContext::Multiloop, false) => model.multiloop_stem_energy(bases, a, b),
            (StemContext::Multiloop, true) => model.multiloop_stem_boltzmann_energy(bases, a, b),
            (StemContext::Closing, false) => model.multiloop_closing_stem_energy(bases, a, b),
            (StemContext::Closing, true) => {
                model.multiloop_closing_stem_boltzmann_energy(bases, a, b)
            }
        })
    }

    fn unpaired_multiloop(&self, column: usize, boltzmann: bool) -> f64 {
        let occupancy = self
            .rows
            .iter()
            .filter(|row| row.column_to_position[column].is_some())
            .count() as f64
            / self.rows.len() as f64;
        occupancy
            * if boltzmann {
                self.model.multiloop_unpaired_boltzmann()
            } else {
                self.model.multiloop_unpaired()
            }
    }
}

#[derive(Clone, Copy)]
enum StemContext {
    Exterior,
    Multiloop,
    Closing,
}

#[derive(Clone, Copy)]
enum PairChoice {
    Invalid,
    Hairpin,
    Internal(usize, usize),
    Multiloop,
}

#[derive(Clone, Copy)]
enum SegmentChoice {
    Invalid,
    Unpaired,
    FirstPair(usize),
    AddPair(usize),
}

#[derive(Clone, Copy)]
enum ExteriorChoice {
    Invalid,
    Unpaired,
    Pair(usize),
}

struct AlignmentMfe {
    energy: f64,
    structure: String,
}

#[allow(clippy::too_many_arguments)]
fn analyze_alignment(
    alignment: &[String],
    consensus: &str,
    allowed_pairs: &[Vec<bool>],
    covariance_energies: &[Vec<f64>],
    min_loop: usize,
    gamma: f64,
    model: &EnergyModel,
) -> Result<AnalysisResult, RnaError> {
    let temperature_celsius = model.temperature_celsius();
    let dangles = model.dangles();
    let salt_molar = model.salt_molar();
    let rows = alignment_rows(alignment);
    let scorer = AlignmentScorer {
        rows,
        model,
        allowed_pairs,
        covariance_energies,
    };
    let mfe = alignment_mfe(consensus.len(), min_loop, &scorer)?;
    let partition = alignment_partition(consensus.len(), temperature_celsius, min_loop, &scorer)?;
    let (centroid_structure, centroid_score) =
        crate::decode::decode_centroid(consensus.len(), min_loop, &partition.pair_probabilities);
    let (mea_structure, mea_score) = crate::decode::decode_mea(
        consensus.len(),
        min_loop,
        gamma,
        &partition.pair_probabilities,
        &partition.unpaired_probabilities,
    );
    let ensemble = summarize(&partition);
    Ok(AnalysisResult {
        sequence: if model.nucleic_acid() == crate::energy::NucleicAcid::Dna {
            consensus.replace('U', "T")
        } else {
            consensus.to_string()
        },
        length: consensus.len(),
        temperature_celsius,
        model: ModelDescription {
            parameter_set: if model.parameter_profile_name().is_some() {
                "Ribon custom thermodynamic overlay / alignment-column average"
            } else {
                match model.nucleic_acid() {
                    crate::energy::NucleicAcid::Rna => {
                        "Ribon RNAstructure 6.6 RNA / alignment-column average"
                    }
                    crate::energy::NucleicAcid::Dna => {
                        "Ribon RNAstructure 6.6 DNA / alignment-column average"
                    }
                }
            },
            parameter_profile_name: model.parameter_profile_name().map(str::to_owned),
            parameter_fingerprint_sha256: model
                .parameter_profile_fingerprint()
                .map(str::to_owned),
            mfe: if model.parameter_profile_name().is_some() {
                "alignment-averaged custom loop grammar with covariance"
            } else {
                match model.nucleic_acid() {
                    crate::energy::NucleicAcid::Rna => {
                        "alignment-averaged RNAstructure 6.6 RNA loop grammar with covariance"
                    }
                    crate::energy::NucleicAcid::Dna => {
                        "alignment-averaged RNAstructure 6.6 DNA loop grammar with covariance"
                    }
                }
            },
            ensemble: if model.parameter_profile_name().is_some() {
                "alignment-averaged custom McCaskill grammar with covariance"
            } else {
                match model.nucleic_acid() {
                    crate::energy::NucleicAcid::Rna => {
                        "alignment-averaged RNA McCaskill grammar with covariance"
                    }
                    crate::energy::NucleicAcid::Dna => {
                        "alignment-averaged DNA McCaskill grammar with covariance"
                    }
                }
            },
            dangles,
            salt_correction: model.has_salt_correction(),
            salt_molar,
            note: "Each loop contribution is evaluated on every ungapped alignment row, averaged across all rows, and augmented by the column-pair covariance pseudo-energy.",
        },
        mfe_structure: mfe.structure,
        mfe_energy_kcal_mol: mfe.energy,
        ensemble_free_energy_kcal_mol: partition.ensemble_free_energy_kcal_mol,
        partition_function: partition.partition_function,
        log_partition_function: partition.log_partition_function,
        centroid_structure,
        centroid_score,
        mea_structure,
        mea_score,
        mea_gamma: gamma,
        pair_probabilities: partition.pair_probabilities,
        unpaired_probabilities: partition.unpaired_probabilities,
        ensemble,
        constraints: ConstraintModel::unconstrained(consensus.len()).summary(),
    })
}

fn alignment_rows(alignment: &[String]) -> Vec<AlignmentRow> {
    alignment
        .iter()
        .map(|row| {
            let mut bases = Vec::new();
            let mut column_to_position = Vec::with_capacity(row.len());
            for &base in row.as_bytes() {
                if base == b'-' {
                    column_to_position.push(None);
                } else {
                    column_to_position.push(Some(bases.len()));
                    bases.push(base);
                }
            }
            AlignmentRow {
                bases,
                column_to_position,
            }
        })
        .collect()
}

fn alignment_mfe(
    n: usize,
    min_loop: usize,
    scorer: &AlignmentScorer<'_>,
) -> Result<AlignmentMfe, RnaError> {
    let mut q = vec![vec![0.0; n]; n];
    let mut qb = vec![vec![INF; n]; n];
    let mut m1 = vec![vec![INF; n]; n];
    let mut m2 = vec![vec![INF; n]; n];
    let mut q_choice = vec![vec![ExteriorChoice::Invalid; n]; n];
    let mut qb_choice = vec![vec![PairChoice::Invalid; n]; n];
    let mut m1_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    let mut m2_choice = vec![vec![SegmentChoice::Invalid; n]; n];
    for (i, row) in q_choice.iter_mut().enumerate() {
        row[i] = ExteriorChoice::Unpaired;
    }

    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if span > min_loop && scorer.allows_pair(i, j) {
                let covariance = scorer.covariance(i, j);
                let mut best = scorer.hairpin(i, j, false) + covariance;
                let mut choice = if best.is_finite() {
                    PairChoice::Hairpin
                } else {
                    PairChoice::Invalid
                };
                let internal_limit = scorer.model.internal_loop_limit(j.saturating_sub(i + 2));
                for left in 0..=internal_limit {
                    let k = i + left + 1;
                    if k >= j {
                        break;
                    }
                    for right in 0..=internal_limit - left {
                        let Some(l) = j.checked_sub(right + 1) else {
                            continue;
                        };
                        if k >= l || l - k <= min_loop || qb[k][l] >= INF / 2.0 {
                            continue;
                        }
                        let energy = scorer.internal(i, j, k, l, false);
                        if energy.is_finite() && covariance + energy + qb[k][l] < best {
                            best = covariance + energy + qb[k][l];
                            choice = PairChoice::Internal(k, l);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] < INF / 2.0 {
                    let energy = covariance
                        + scorer.model.multiloop_closing()
                        + scorer.stem(i, j, StemContext::Closing, false)
                        + m2[i + 1][j - 1];
                    if energy < best {
                        best = energy;
                        choice = PairChoice::Multiloop;
                    }
                }
                qb[i][j] = best;
                qb_choice[i][j] = choice;
            }

            let unpaired = scorer.unpaired_multiloop(j, false);
            let mut best1 = if m1[i][j - 1] < INF / 2.0 {
                m1[i][j - 1] + unpaired
            } else {
                INF
            };
            let mut choice1 = if best1 < INF / 2.0 {
                SegmentChoice::Unpaired
            } else {
                SegmentChoice::Invalid
            };
            let mut best2 = if m2[i][j - 1] < INF / 2.0 {
                m2[i][j - 1] + unpaired
            } else {
                INF
            };
            let mut choice2 = if best2 < INF / 2.0 {
                SegmentChoice::Unpaired
            } else {
                SegmentChoice::Invalid
            };
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    if qb_row[j] >= INF / 2.0 {
                        continue;
                    }
                    let branch = qb_row[j] + scorer.stem(k, j, StemContext::Multiloop, false);
                    let leading = (i..k)
                        .map(|column| scorer.unpaired_multiloop(column, false))
                        .sum::<f64>();
                    if leading + branch < best1 {
                        best1 = leading + branch;
                        choice1 = SegmentChoice::FirstPair(k);
                    }
                    if k > i && m1[i][k - 1] < INF / 2.0 {
                        let additional = m1[i][k - 1] + branch;
                        if additional < best1 {
                            best1 = additional;
                            choice1 = SegmentChoice::AddPair(k);
                        }
                        if additional < best2 {
                            best2 = additional;
                            choice2 = SegmentChoice::AddPair(k);
                        }
                    }
                }
            }
            m1[i][j] = best1;
            m2[i][j] = best2;
            m1_choice[i][j] = choice1;
            m2_choice[i][j] = choice2;

            let mut best_q = q[i][j - 1];
            let mut choice_q = ExteriorChoice::Unpaired;
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    if qb_row[j] >= INF / 2.0 {
                        continue;
                    }
                    let left = if k == i { 0.0 } else { q[i][k - 1] };
                    let candidate =
                        left + qb_row[j] + scorer.stem(k, j, StemContext::Exterior, false);
                    if candidate < best_q {
                        best_q = candidate;
                        choice_q = ExteriorChoice::Pair(k);
                    }
                }
            }
            q[i][j] = best_q;
            q_choice[i][j] = choice_q;
        }
    }

    let mut pairs = Vec::new();
    trace_alignment_q(
        0,
        n - 1,
        &q_choice,
        &qb_choice,
        &m1_choice,
        &m2_choice,
        &mut pairs,
    );
    pairs.sort_unstable();
    Ok(AlignmentMfe {
        energy: q[0][n - 1],
        structure: pairs_to_dot_bracket(n, &pairs),
    })
}

fn trace_alignment_q(
    i: usize,
    j: usize,
    q_choice: &[Vec<ExteriorChoice>],
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= j {
        return;
    }
    match q_choice[i][j] {
        ExteriorChoice::Invalid => {}
        ExteriorChoice::Unpaired => {
            trace_alignment_q(i, j - 1, q_choice, qb_choice, m1_choice, m2_choice, pairs)
        }
        ExteriorChoice::Pair(k) => {
            if k > i {
                trace_alignment_q(i, k - 1, q_choice, qb_choice, m1_choice, m2_choice, pairs);
            }
            trace_alignment_pair(k, j, qb_choice, m1_choice, m2_choice, pairs);
        }
    }
}

fn trace_alignment_pair(
    i: usize,
    j: usize,
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    pairs: &mut Vec<(usize, usize)>,
) {
    pairs.push((i, j));
    match qb_choice[i][j] {
        PairChoice::Invalid | PairChoice::Hairpin => {}
        PairChoice::Internal(k, l) => {
            trace_alignment_pair(k, l, qb_choice, m1_choice, m2_choice, pairs)
        }
        PairChoice::Multiloop => {
            trace_alignment_segment(i + 1, j - 1, false, qb_choice, m1_choice, m2_choice, pairs)
        }
    }
}

fn trace_alignment_segment(
    i: usize,
    j: usize,
    one: bool,
    qb_choice: &[Vec<PairChoice>],
    m1_choice: &[Vec<SegmentChoice>],
    m2_choice: &[Vec<SegmentChoice>],
    pairs: &mut Vec<(usize, usize)>,
) {
    if i > j {
        return;
    }
    let choice = if one {
        m1_choice[i][j]
    } else {
        m2_choice[i][j]
    };
    match choice {
        SegmentChoice::Invalid => {}
        SegmentChoice::Unpaired => {
            trace_alignment_segment(i, j - 1, one, qb_choice, m1_choice, m2_choice, pairs)
        }
        SegmentChoice::FirstPair(k) => {
            trace_alignment_pair(k, j, qb_choice, m1_choice, m2_choice, pairs)
        }
        SegmentChoice::AddPair(k) => {
            trace_alignment_segment(i, k - 1, true, qb_choice, m1_choice, m2_choice, pairs);
            trace_alignment_pair(k, j, qb_choice, m1_choice, m2_choice, pairs);
        }
    }
}

fn alignment_partition(
    n: usize,
    temperature_celsius: f64,
    min_loop: usize,
    scorer: &AlignmentScorer<'_>,
) -> Result<PartitionResult, RnaError> {
    let rt = GAS_CONSTANT_KCAL * (temperature_celsius + 273.15);
    let mut q = vec![vec![NEG_INF; n]; n];
    let mut qb = vec![vec![NEG_INF; n]; n];
    let mut m1 = vec![vec![NEG_INF; n]; n];
    let mut m2 = vec![vec![NEG_INF; n]; n];
    for (i, row) in q.iter_mut().enumerate() {
        row[i] = 0.0;
    }
    for span in 1..n {
        for i in 0..n - span {
            let j = i + span;
            if span > min_loop && scorer.allows_pair(i, j) {
                let covariance = scorer.covariance(i, j);
                let hairpin = scorer.hairpin(i, j, true);
                if hairpin.is_finite() {
                    qb[i][j] = -(hairpin + covariance) / rt;
                }
                let internal_limit = scorer.model.internal_loop_limit(j.saturating_sub(i + 2));
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
                        let energy = scorer.internal(i, j, k, l, true);
                        if energy.is_finite() {
                            let child = qb[k][l];
                            log_update(&mut qb[i][j], child - (energy + covariance) / rt);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] != NEG_INF {
                    let energy = covariance
                        + scorer.model.multiloop_closing_boltzmann()
                        + scorer.stem(i, j, StemContext::Closing, true);
                    log_update(&mut qb[i][j], m2[i + 1][j - 1] - energy / rt);
                }
            }

            let log_unpaired = -scorer.unpaired_multiloop(j, true) / rt;
            if m1[i][j - 1] != NEG_INF {
                m1[i][j] = m1[i][j - 1] + log_unpaired;
            }
            if m2[i][j - 1] != NEG_INF {
                m2[i][j] = m2[i][j - 1] + log_unpaired;
            }
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    if qb_row[j] == NEG_INF {
                        continue;
                    }
                    let branch = qb_row[j] - scorer.stem(k, j, StemContext::Multiloop, true) / rt;
                    let leading = (i..k)
                        .map(|column| scorer.unpaired_multiloop(column, true))
                        .sum::<f64>();
                    log_update(&mut m1[i][j], branch - leading / rt);
                    if k > i && m1[i][k - 1] != NEG_INF {
                        let additional = m1[i][k - 1] + branch;
                        log_update(&mut m1[i][j], additional);
                        log_update(&mut m2[i][j], additional);
                    }
                }
            }

            q[i][j] = q[i][j - 1];
            if j > min_loop {
                for (k, qb_row) in qb.iter().enumerate().take(j - min_loop).skip(i) {
                    if qb_row[j] == NEG_INF {
                        continue;
                    }
                    let left = if k == i { 0.0 } else { q[i][k - 1] };
                    let stem = -scorer.stem(k, j, StemContext::Exterior, true) / rt;
                    log_update(&mut q[i][j], left + stem + qb_row[j]);
                }
            }
        }
    }

    let log_z = q[0][n - 1];
    if log_z == NEG_INF {
        return Err(RnaError::Numerical(
            "comparative partition function has zero weight".into(),
        ));
    }
    let mut oq = vec![vec![NEG_INF; n]; n];
    let mut oqb = vec![vec![NEG_INF; n]; n];
    let mut om1 = vec![vec![NEG_INF; n]; n];
    let mut om2 = vec![vec![NEG_INF; n]; n];
    oq[0][n - 1] = 0.0;

    for span in (1..n).rev() {
        for i in 0..n - span {
            let j = i + span;
            let q_parent = oq[i][j];
            if q_parent != NEG_INF {
                log_update(&mut oq[i][j - 1], q_parent);
                if j > min_loop {
                    for k in i..=j - min_loop - 1 {
                        if qb[k][j] == NEG_INF {
                            continue;
                        }
                        let left = if k == i { 0.0 } else { q[i][k - 1] };
                        let stem = -scorer.stem(k, j, StemContext::Exterior, true) / rt;
                        log_update(&mut oqb[k][j], q_parent + left + stem);
                        if k > i {
                            log_update(&mut oq[i][k - 1], q_parent + stem + qb[k][j]);
                        }
                    }
                }
            }

            let parent2 = om2[i][j];
            if parent2 != NEG_INF {
                let unpaired = -scorer.unpaired_multiloop(j, true) / rt;
                if m2[i][j - 1] != NEG_INF {
                    log_update(&mut om2[i][j - 1], parent2 + unpaired);
                }
                if j > min_loop {
                    for k in i + 1..=j - min_loop - 1 {
                        if m1[i][k - 1] == NEG_INF || qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -scorer.stem(k, j, StemContext::Multiloop, true) / rt;
                        log_update(&mut om1[i][k - 1], parent2 + stem + qb[k][j]);
                        log_update(&mut oqb[k][j], parent2 + m1[i][k - 1] + stem);
                    }
                }
            }

            let parent1 = om1[i][j];
            if parent1 != NEG_INF {
                let unpaired = -scorer.unpaired_multiloop(j, true) / rt;
                if m1[i][j - 1] != NEG_INF {
                    log_update(&mut om1[i][j - 1], parent1 + unpaired);
                }
                if j > min_loop {
                    for k in i..=j - min_loop - 1 {
                        if qb[k][j] == NEG_INF {
                            continue;
                        }
                        let stem = -scorer.stem(k, j, StemContext::Multiloop, true) / rt;
                        let leading = (i..k)
                            .map(|column| scorer.unpaired_multiloop(column, true))
                            .sum::<f64>();
                        log_update(&mut oqb[k][j], parent1 + stem - leading / rt);
                        if k > i && m1[i][k - 1] != NEG_INF {
                            log_update(&mut om1[i][k - 1], parent1 + stem + qb[k][j]);
                            log_update(&mut oqb[k][j], parent1 + m1[i][k - 1] + stem);
                        }
                    }
                }
            }

            let pair_parent = oqb[i][j];
            if pair_parent != NEG_INF && qb[i][j] != NEG_INF {
                let covariance = scorer.covariance(i, j);
                let internal_limit = scorer.model.internal_loop_limit(j.saturating_sub(i + 2));
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
                        let energy = scorer.internal(i, j, k, l, true);
                        if energy.is_finite() {
                            log_update(&mut oqb[k][l], pair_parent - (energy + covariance) / rt);
                        }
                    }
                }
                if i + 1 < j && m2[i + 1][j - 1] != NEG_INF {
                    let energy = covariance
                        + scorer.model.multiloop_closing_boltzmann()
                        + scorer.stem(i, j, StemContext::Closing, true);
                    log_update(&mut om2[i + 1][j - 1], pair_parent - energy / rt);
                }
            }
        }
    }

    let mut pair_probabilities = Vec::new();
    let mut paired_mass = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
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
                paired_mass[i] += probability;
                paired_mass[j] += probability;
            }
        }
    }
    let unpaired_probabilities = paired_mass
        .into_iter()
        .map(|mass| (1.0 - mass).clamp(0.0, 1.0))
        .collect();
    Ok(PartitionResult {
        partition_function: if log_z < f64::MAX.ln() {
            log_z.exp()
        } else {
            f64::MAX
        },
        log_partition_function: log_z,
        ensemble_free_energy_kcal_mol: -rt * log_z,
        pair_probabilities,
        unpaired_probabilities,
        model: "alignment-averaged McCaskill grammar with covariance",
        salt_molar: scorer.model.salt_molar(),
        constraints: ConstraintSummary {
            enabled: false,
            hard_constraints: false,
            soft_constraints: true,
            no_gu: false,
            no_lonely_pairs: false,
            max_span: None,
            probing_kind: None,
            probing_method: None,
            probing_reactivities: Vec::new(),
        },
    })
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

fn validate_options(options: &ComparativeOptions) -> Result<(), RnaError> {
    if !options.covariance_weight_kcal_mol.is_finite()
        || options.covariance_weight_kcal_mol < 0.0
        || !options.incompatible_penalty.is_finite()
        || options.incompatible_penalty < 0.0
        || !options.minimum_pair_occupancy.is_finite()
        || !(0.0..=1.0).contains(&options.minimum_pair_occupancy)
    {
        return Err(RnaError::InvalidOption(
            "comparative-fold options are outside their valid ranges".into(),
        ));
    }
    Ok(())
}

fn normalize_alignment(alignment: &[String]) -> Result<Vec<String>, RnaError> {
    if alignment.is_empty() {
        return Err(RnaError::InvalidOption(
            "comparative folding requires at least one aligned sequence".into(),
        ));
    }
    let mut result = Vec::with_capacity(alignment.len());
    let mut expected = None;
    for (row_index, row) in alignment.iter().enumerate() {
        let mut normalized = String::new();
        for raw in row.chars().filter(|symbol| !symbol.is_whitespace()) {
            let base = raw.to_ascii_uppercase();
            let base = if base == 'T' { 'U' } else { base };
            if !matches!(base, 'A' | 'C' | 'G' | 'U' | '-' | '.') {
                return Err(RnaError::InvalidOption(format!(
                    "invalid alignment symbol {raw:?} in row {}",
                    row_index + 1
                )));
            }
            normalized.push(if base == '.' { '-' } else { base });
        }
        if normalized.is_empty() {
            return Err(RnaError::InvalidOption(format!(
                "alignment row {} is empty",
                row_index + 1
            )));
        }
        if let Some(length) = expected {
            if normalized.len() != length {
                return Err(RnaError::LengthMismatch {
                    sequence: length,
                    structure: normalized.len(),
                });
            }
        } else {
            expected = Some(normalized.len());
        }
        result.push(normalized);
    }
    let length = result[0].len();
    let keep = (0..length)
        .filter(|&column| result.iter().any(|row| row.as_bytes()[column] != b'-'))
        .collect::<Vec<_>>();
    if keep.is_empty() {
        return Err(RnaError::InvalidOption(
            "alignment contains no nucleotide columns".into(),
        ));
    }
    if keep.len() != length {
        result = result
            .into_iter()
            .map(|row| {
                let bytes = row.as_bytes();
                keep.iter().map(|&column| bytes[column] as char).collect()
            })
            .collect();
    }
    Ok(result)
}

fn consensus_sequence(alignment: &[String]) -> Result<String, RnaError> {
    let length = alignment[0].len();
    let rows = alignment
        .iter()
        .map(|row| row.as_bytes())
        .collect::<Vec<_>>();
    let mut consensus = String::with_capacity(length);
    for column in 0..length {
        let mut counts = [0usize; 4];
        for row in &rows {
            if let Some(index) = base_index(row[column]) {
                counts[index] += 1;
            }
        }
        let Some((index, count)) = counts
            .iter()
            .enumerate()
            .max_by_key(|(index, count)| (**count, std::cmp::Reverse(*index)))
        else {
            unreachable!();
        };
        if *count == 0 {
            return Err(RnaError::InvalidOption(format!(
                "alignment column {} contains only gaps",
                column + 1
            )));
        }
        consensus.push([b'A', b'C', b'G', b'U'][index] as char);
    }
    Ok(consensus)
}

fn covariation_entry(
    alignment: &[String],
    i: usize,
    j: usize,
    incompatible_penalty: f64,
) -> CovariationEntry {
    let observations = alignment
        .iter()
        .map(|row| {
            let bytes = row.as_bytes();
            (bytes[i], bytes[j])
        })
        .collect::<Vec<_>>();
    let canonical = observations
        .iter()
        .filter(|&&(a, b)| pair_type(a, b).is_some())
        .count();
    let incompatible = observations
        .iter()
        .filter(|&&(a, b)| pair_type(a, b).is_none() && (a, b) != (b'-', b'-'))
        .count();
    let mut compensatory = 0usize;
    let mut consistent = 0usize;
    let canonical_observations = observations
        .iter()
        .copied()
        .filter(|&(a, b)| pair_type(a, b).is_some())
        .collect::<Vec<_>>();
    for (index, &(a, b)) in canonical_observations.iter().enumerate() {
        for &(c, d) in &canonical_observations[index + 1..] {
            if a != c && b != d {
                compensatory += 1;
            } else if a != c || b != d {
                consistent += 1;
            }
        }
    }
    // RNAalifold 2002 conservation term: sum the two Hamming distances for
    // every pair of sequences that can both realize the base pair, normalized
    // by the number of alignment rows (not by the number of row pairs).
    let conservation = canonical_observations
        .iter()
        .enumerate()
        .map(|(index, &(a, b))| {
            canonical_observations[index + 1..]
                .iter()
                .map(|&(c, d)| usize::from(a != c) + usize::from(b != d))
                .sum::<usize>()
        })
        .sum::<usize>() as f64
        / alignment.len() as f64;
    // The published inconsistency term assigns 0 to canonical pairs, 0.25 to
    // a double gap, and 1 to every other non-pairing row.
    let inconsistency = observations
        .iter()
        .map(|&(a, b)| {
            if pair_type(a, b).is_some() {
                0.0
            } else if (a, b) == (b'-', b'-') {
                0.25
            } else {
                1.0
            }
        })
        .sum::<f64>()
        / alignment.len() as f64;
    let nongap_observations = observations
        .iter()
        .copied()
        .filter(|&(a, b)| a != b'-' && b != b'-')
        .collect::<Vec<_>>();
    CovariationEntry {
        i: i + 1,
        j: j + 1,
        score: conservation - incompatible_penalty * inconsistency,
        mutual_information_bits: mutual_information(&nongap_observations),
        canonical_fraction: canonical as f64 / observations.len() as f64,
        compensatory_changes: compensatory,
        consistent_changes: consistent,
        incompatible_sequences: incompatible,
    }
}

fn mutual_information(observations: &[(u8, u8)]) -> f64 {
    if observations.is_empty() {
        return 0.0;
    }
    let mut joint = HashMap::<(u8, u8), usize>::new();
    let mut left = HashMap::<u8, usize>::new();
    let mut right = HashMap::<u8, usize>::new();
    for &(a, b) in observations {
        *joint.entry((a, b)).or_default() += 1;
        *left.entry(a).or_default() += 1;
        *right.entry(b).or_default() += 1;
    }
    let n = observations.len() as f64;
    joint
        .into_iter()
        .map(|((a, b), count)| {
            let p = count as f64 / n;
            let pa = left[&a] as f64 / n;
            let pb = right[&b] as f64 / n;
            p * (p / (pa * pb)).log2()
        })
        .sum()
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'U' => Some(3),
        _ => None,
    }
}

fn pair_type(a: u8, b: u8) -> Option<u8> {
    match (a, b) {
        (b'C', b'G') => Some(0),
        (b'G', b'C') => Some(1),
        (b'G', b'U') => Some(2),
        (b'U', b'G') => Some(3),
        (b'A', b'U') => Some(4),
        (b'U', b'A') => Some(5),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_sequence_reduces_exactly_to_ordinary_fold() {
        let alignment = vec!["GGGAAACCC".to_string()];
        let comparative = comparative_fold(
            &alignment,
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ComparativeOptions::default(),
        )
        .unwrap();
        let ordinary = crate::analyze("GGGAAACCC", 37.0, 3, 1.0).unwrap();
        assert_eq!(comparative.consensus_structure, ordinary.mfe_structure);
        assert!(
            (comparative.consensus_energy_kcal_mol - ordinary.mfe_energy_kcal_mol).abs() < 1.0e-12
        );
        assert!(
            (comparative.analysis.log_partition_function - ordinary.log_partition_function).abs()
                < 1.0e-12
        );
        assert_eq!(
            comparative.analysis.pair_probabilities.len(),
            ordinary.pair_probabilities.len()
        );
        for (actual, expected) in comparative
            .analysis
            .pair_probabilities
            .iter()
            .zip(&ordinary.pair_probabilities)
        {
            assert_eq!((actual.i, actual.j), (expected.i, expected.j));
            assert!((actual.probability - expected.probability).abs() < 1.0e-12);
        }
    }

    #[test]
    fn dna_consensus_and_alignment_preserve_thymine() {
        let alignment = vec!["GGGTTTCCC".to_string()];
        let model =
            EnergyModel::with_parameter_family(37.0, 0, 1.021, crate::energy::NucleicAcid::Dna)
                .unwrap();
        let result =
            comparative_fold_with_model(&alignment, 3, 1.0, &model, &ComparativeOptions::default())
                .unwrap();
        assert_eq!(result.consensus_sequence, "GGGTTTCCC");
        assert_eq!(result.analysis.sequence, "GGGTTTCCC");
        assert_eq!(result.alignment, alignment);
    }

    #[test]
    fn compensatory_change_scores_more_than_consistent_change() {
        let compensatory = vec!["GAAAC".to_string(), "CAAAG".to_string()];
        let consistent = vec!["GAAAC".to_string(), "GAAAU".to_string()];
        let a = covariation_entry(&normalize_alignment(&compensatory).unwrap(), 0, 4, 1.0);
        let b = covariation_entry(&normalize_alignment(&consistent).unwrap(), 0, 4, 1.0);
        assert!(a.score > b.score);
        assert_eq!(a.compensatory_changes, 1);
        assert_eq!(b.consistent_changes, 1);
    }

    #[test]
    fn covariance_uses_the_published_rnaalifold_2002_normalization() {
        let alignment = vec!["AU".into(), "GU".into(), "GC".into(), "--".into()];
        let entry = covariation_entry(&alignment, 0, 1, 1.0);
        // Canonical row pairs contribute Hamming sums 1, 2, and 1; the
        // double-gap row contributes an inconsistency of 0.25. Both are
        // normalized by N=4.
        assert!((entry.score - (4.0 / 4.0 - 0.25 / 4.0)).abs() < 1e-12);
        assert_eq!(entry.compensatory_changes, 1);
        assert_eq!(entry.consistent_changes, 2);
    }

    #[test]
    fn gaps_are_removed_per_row_before_loop_energy_evaluation() {
        let alignment = vec!["GGG-AAACCC".to_string(), "GGG-AAACCC".to_string()];
        let comparative = comparative_fold(
            &alignment,
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ComparativeOptions::default(),
        )
        .unwrap();
        let ordinary = crate::analyze("GGGAAACCC", 37.0, 3, 1.0).unwrap();
        assert_eq!(comparative.analysis.mfe_structure, ordinary.mfe_structure);
        assert!(
            (comparative.analysis.mfe_energy_kcal_mol - ordinary.mfe_energy_kcal_mol).abs()
                < 1.0e-12
        );
    }
}
