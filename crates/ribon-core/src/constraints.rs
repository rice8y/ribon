use crate::energy::EnergyModel;
use crate::structure::RnaError;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct ConstraintConfig {
    pub force_unpaired: Vec<usize>,
    pub force_paired: Vec<usize>,
    pub force_pairs: Vec<PairConstraint>,
    pub forbid_pairs: Vec<PairConstraint>,
    pub max_span: Option<usize>,
    pub no_gu: bool,
    pub no_lonely_pairs: bool,
    pub soft: SoftConstraintConfig,
    pub probing: Option<ProbingConfig>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct PairConstraint {
    pub i: usize,
    pub j: usize,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct SoftConstraintConfig {
    /// Pseudo-energy added whenever the nucleotide is unpaired.
    pub unpaired: Vec<PositionEnergy>,
    /// Pseudo-energy added whenever the nucleotide is paired.
    pub paired: Vec<PositionEnergy>,
    /// Pseudo-energy added for one particular base pair.
    pub pairs: Vec<PairEnergy>,
    /// Per-nucleotide pseudo-energy added whenever that nucleotide
    /// participates in a nearest-neighbor stack.
    pub stack: Vec<PositionEnergy>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct PositionEnergy {
    pub position: usize,
    pub energy_kcal_mol: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
pub struct PairEnergy {
    pub i: usize,
    pub j: usize,
    pub energy_kcal_mol: f64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
pub struct ProbingConfig {
    /// `shape` and `dms` are metadata labels; both use the selected
    /// pseudo-energy conversion.
    pub kind: String,
    /// `deigan` or `zarringhalam`.
    pub method: String,
    /// One entry per nucleotide. JSON `null` and negative values mean missing.
    pub reactivities: Vec<Option<f64>>,
    pub slope: f64,
    pub intercept: f64,
    pub beta: f64,
    /// Zarringhalam conversion: `O`, `L`, `M`, `C`, or `S`.
    pub conversion: String,
    pub default_probability: f64,
}

impl Default for ProbingConfig {
    fn default() -> Self {
        Self {
            kind: "shape".into(),
            method: "deigan".into(),
            reactivities: Vec::new(),
            slope: 1.8,
            intercept: -0.6,
            beta: 0.89,
            conversion: "O".into(),
            default_probability: 0.5,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ConstraintSummary {
    pub enabled: bool,
    pub hard_constraints: bool,
    pub soft_constraints: bool,
    pub no_gu: bool,
    pub no_lonely_pairs: bool,
    pub max_span: Option<usize>,
    pub probing_kind: Option<String>,
    pub probing_method: Option<String>,
    pub probing_reactivities: Vec<Option<f64>>,
}

#[derive(Clone, Debug)]
pub struct ConstraintModel {
    length: usize,
    force_unpaired: Vec<bool>,
    force_paired: Vec<bool>,
    forced_partner: Vec<Option<usize>>,
    forbidden_pairs: Vec<Vec<bool>>,
    max_span: Option<usize>,
    no_gu: bool,
    no_lonely_pairs: bool,
    unpaired_energy: Vec<f64>,
    unpaired_prefix: Vec<f64>,
    paired_energy: Vec<f64>,
    pair_energy: HashMap<(usize, usize), f64>,
    stack_energy: Vec<f64>,
    context_stack_energy: HashMap<(usize, usize, usize, usize), f64>,
    terminal_pair_energy: HashMap<(usize, usize), f64>,
    summary: ConstraintSummary,
}

impl ConstraintModel {
    pub fn unconstrained(length: usize) -> Self {
        Self::compile(length, &ConstraintConfig::default())
            .expect("the empty constraint set is valid")
    }

    pub fn compile(length: usize, config: &ConstraintConfig) -> Result<Self, RnaError> {
        if length == 0 {
            return Err(RnaError::EmptySequence);
        }
        let mut force_unpaired = vec![false; length];
        let mut force_paired = vec![false; length];
        let mut forced_partner = vec![None; length];
        let mut forbidden_pairs = vec![vec![false; length]; length];
        let mut unpaired_energy = vec![0.0; length];
        let mut paired_energy = vec![0.0; length];
        let mut pair_energy = HashMap::new();
        let mut stack_energy = vec![0.0; length];

        for &position in &config.force_unpaired {
            force_unpaired[index(position, length, "force_unpaired")?] = true;
        }
        for &position in &config.force_paired {
            force_paired[index(position, length, "force_paired")?] = true;
        }
        for pair in &config.force_pairs {
            let (i, j) = pair_indices(*pair, length, "force_pairs")?;
            if let Some(existing) = forced_partner[i] {
                if existing != j {
                    return Err(invalid(format!(
                        "nucleotide {} has conflicting forced partners {} and {}",
                        i + 1,
                        existing + 1,
                        j + 1
                    )));
                }
            }
            if let Some(existing) = forced_partner[j] {
                if existing != i {
                    return Err(invalid(format!(
                        "nucleotide {} has conflicting forced partners {} and {}",
                        j + 1,
                        existing + 1,
                        i + 1
                    )));
                }
            }
            forced_partner[i] = Some(j);
            forced_partner[j] = Some(i);
            force_paired[i] = true;
            force_paired[j] = true;
        }
        for pair in &config.forbid_pairs {
            let (i, j) = pair_indices(*pair, length, "forbid_pairs")?;
            forbidden_pairs[i][j] = true;
            forbidden_pairs[j][i] = true;
        }

        for position in 0..length {
            if force_unpaired[position] && force_paired[position] {
                return Err(invalid(format!(
                    "nucleotide {} is forced both paired and unpaired",
                    position + 1
                )));
            }
        }
        for i in 0..length {
            if let Some(j) = forced_partner[i] {
                if forbidden_pairs[i][j] {
                    return Err(invalid(format!(
                        "base pair {}-{} is both forced and forbidden",
                        i + 1,
                        j + 1
                    )));
                }
            }
        }
        let forced_pairs: Vec<(usize, usize)> = forced_partner
            .iter()
            .enumerate()
            .filter_map(|(i, partner)| partner.filter(|&j| i < j).map(|j| (i, j)))
            .collect();
        for (position, &(i, j)) in forced_pairs.iter().enumerate() {
            for &(k, l) in &forced_pairs[position + 1..] {
                if (i < k && k < j && j < l) || (k < i && i < l && l < j) {
                    return Err(RnaError::PseudoknotUnsupported("forced constraints"));
                }
            }
        }

        for entry in &config.soft.unpaired {
            let i = index(entry.position, length, "soft.unpaired")?;
            unpaired_energy[i] +=
                quantize_soft(finite(entry.energy_kcal_mol, "soft.unpaired energy")?);
        }
        for entry in &config.soft.paired {
            let i = index(entry.position, length, "soft.paired")?;
            paired_energy[i] += quantize_soft(finite(entry.energy_kcal_mol, "soft.paired energy")?);
        }
        for entry in &config.soft.pairs {
            let (i, j) = pair_indices(
                PairConstraint {
                    i: entry.i,
                    j: entry.j,
                },
                length,
                "soft.pairs",
            )?;
            *pair_energy.entry((i, j)).or_insert(0.0) +=
                quantize_soft(finite(entry.energy_kcal_mol, "soft pair energy")?);
        }
        for entry in &config.soft.stack {
            let i = index(entry.position, length, "soft.stack")?;
            stack_energy[i] += quantize_soft(finite(entry.energy_kcal_mol, "soft.stack energy")?);
        }

        let mut probing_kind = None;
        let mut probing_method = None;
        let mut probing_reactivities = Vec::new();
        if let Some(probing) = &config.probing {
            if probing.reactivities.len() != length {
                return Err(invalid(format!(
                    "probing reactivity length {} does not match sequence length {length}",
                    probing.reactivities.len()
                )));
            }
            let kind = probing.kind.to_ascii_lowercase();
            if !matches!(kind.as_str(), "shape" | "dms") {
                return Err(invalid("probing kind must be `shape` or `dms`"));
            }
            let method = probing.method.to_ascii_lowercase();
            match method.as_str() {
                "deigan" => apply_deigan(probing, &mut stack_energy)?,
                "zarringhalam" => {
                    apply_zarringhalam(probing, &mut unpaired_energy, &mut paired_energy)?
                }
                _ => return Err(invalid("probing method must be `deigan` or `zarringhalam`")),
            }
            probing_kind = Some(kind);
            probing_method = Some(method);
            probing_reactivities = probing.reactivities.clone();
        }

        if let Some(max_span) = config.max_span {
            if max_span == 0 {
                return Err(invalid("max_span must be positive"));
            }
            for &(i, j) in &forced_pairs {
                if j - i >= max_span {
                    return Err(invalid(format!(
                        "forced pair {}-{} exceeds max_span={max_span}",
                        i + 1,
                        j + 1
                    )));
                }
            }
        }

        let mut unpaired_prefix = Vec::with_capacity(length + 1);
        unpaired_prefix.push(0.0);
        for &energy in &unpaired_energy {
            unpaired_prefix.push(unpaired_prefix.last().copied().unwrap_or(0.0) + energy);
        }
        let hard_constraints = force_unpaired.iter().any(|&value| value)
            || force_paired.iter().any(|&value| value)
            || forbidden_pairs
                .iter()
                .any(|row| row.iter().any(|&value| value))
            || config.max_span.is_some()
            || config.no_gu
            || config.no_lonely_pairs;
        let soft_constraints = unpaired_energy.iter().any(|energy| energy.abs() > 0.0)
            || paired_energy.iter().any(|energy| energy.abs() > 0.0)
            || pair_energy.values().any(|energy| energy.abs() > 0.0)
            || stack_energy.iter().any(|energy| energy.abs() > 0.0)
            || config.probing.is_some();
        let summary = ConstraintSummary {
            enabled: hard_constraints || soft_constraints,
            hard_constraints,
            soft_constraints,
            no_gu: config.no_gu,
            no_lonely_pairs: config.no_lonely_pairs,
            max_span: config.max_span,
            probing_kind,
            probing_method,
            probing_reactivities,
        };

        Ok(Self {
            length,
            force_unpaired,
            force_paired,
            forced_partner,
            forbidden_pairs,
            max_span: config.max_span,
            no_gu: config.no_gu,
            no_lonely_pairs: config.no_lonely_pairs,
            unpaired_energy,
            unpaired_prefix,
            paired_energy,
            pair_energy,
            stack_energy,
            context_stack_energy: HashMap::new(),
            terminal_pair_energy: HashMap::new(),
            summary,
        })
    }

    pub fn summary(&self) -> ConstraintSummary {
        self.summary.clone()
    }

    pub fn no_lonely_pairs(&self) -> bool {
        self.no_lonely_pairs
    }

    pub fn allows_unpaired(&self, position: usize) -> bool {
        !self.force_paired[position]
    }

    pub fn unpaired_energy(&self, position: usize) -> f64 {
        self.unpaired_energy[position]
    }

    pub fn unpaired_range_energy(&self, start: usize, end: usize) -> Option<f64> {
        if start > end {
            return Some(0.0);
        }
        if self.force_paired[start..=end].iter().any(|&forced| forced) {
            None
        } else {
            Some(self.unpaired_prefix[end + 1] - self.unpaired_prefix[start])
        }
    }

    pub fn allows_pair(&self, bases: &[u8], i: usize, j: usize, model: &EnergyModel) -> bool {
        if i >= j || j >= self.length || !model.can_pair(bases[i], bases[j]) {
            return false;
        }
        if self.force_unpaired[i]
            || self.force_unpaired[j]
            || self.forbidden_pairs[i][j]
            || self.max_span.is_some_and(|span| j - i >= span)
            || (self.no_gu && matches!((bases[i], bases[j]), (b'G', b'U') | (b'U', b'G')))
        {
            return false;
        }
        if self.forced_partner[i].is_some_and(|partner| partner != j)
            || self.forced_partner[j].is_some_and(|partner| partner != i)
        {
            return false;
        }
        true
    }

    /// ViennaRNA's partition-function noLP path removes pair types that can
    /// never participate in a stack, while its MFE grammar enforces noLP per
    /// structure. Keep that stage-specific convention explicit.
    pub fn allows_pair_for_partition(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        min_loop: usize,
        model: &EnergyModel,
    ) -> bool {
        if !self.allows_pair(bases, i, j, model) {
            return false;
        }
        if !self.no_lonely_pairs {
            return true;
        }
        // ViennaRNA computes this default noLP mask from the raw pair matrix.
        // In particular, the adjacent supporting pair is allowed to be GU even
        // with noGU enabled; noGU is applied to the candidate pair itself.
        // User-installed hard constraints likewise do not retroactively change
        // this precomputed stack-capability mask.
        let adjacent_pair_allowed = |a: usize, b: usize| {
            model.can_pair(bases[a], bases[b]) && self.max_span.is_none_or(|span| b - a < span)
        };
        let outer = i > 0 && j + 1 < self.length && adjacent_pair_allowed(i - 1, j + 1);
        let inner = i + 1 < j && j - i - 2 > min_loop && adjacent_pair_allowed(i + 1, j - 1);
        outer || inner
    }

    pub fn pair_energy(&self, i: usize, j: usize) -> f64 {
        self.paired_energy[i]
            + self.paired_energy[j]
            + self
                .pair_energy
                .get(&(i.min(j), i.max(j)))
                .copied()
                .unwrap_or(0.0)
    }

    pub fn stack_energy(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        if k == i + 1 && l + 1 == j {
            self.stack_energy[i]
                + self.stack_energy[j]
                + self.stack_energy[k]
                + self.stack_energy[l]
                + self
                    .context_stack_energy
                    .get(&(i, j, k, l))
                    .copied()
                    .unwrap_or(0.0)
        } else {
            0.0
        }
    }

    pub(crate) fn terminal_pair_energy(&self, i: usize, j: usize) -> f64 {
        self.terminal_pair_energy
            .get(&(i.min(j), i.max(j)))
            .copied()
            .unwrap_or(0.0)
    }

    pub(crate) fn add_context_stack_energy(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        energy: f64,
    ) {
        if energy != 0.0 {
            *self.context_stack_energy.entry((i, j, k, l)).or_insert(0.0) += energy;
            self.summary.enabled = true;
            self.summary.soft_constraints = true;
        }
    }

    pub(crate) fn add_terminal_pair_energy(&mut self, i: usize, j: usize, energy: f64) {
        if energy != 0.0 {
            *self
                .terminal_pair_energy
                .entry((i.min(j), i.max(j)))
                .or_insert(0.0) += energy;
            self.summary.enabled = true;
            self.summary.soft_constraints = true;
        }
    }

    pub fn validate_structure(
        &self,
        bases: &[u8],
        partner: &[Option<usize>],
        min_loop: usize,
        model: &EnergyModel,
    ) -> Result<(), RnaError> {
        for (position, &actual_partner) in partner.iter().enumerate().take(self.length) {
            if self.force_unpaired[position] && actual_partner.is_some() {
                return Err(invalid(format!(
                    "structure pairs nucleotide {}, which is forced unpaired",
                    position + 1
                )));
            }
            if self.force_paired[position] && actual_partner.is_none() {
                return Err(invalid(format!(
                    "structure leaves nucleotide {} unpaired, but it is forced paired",
                    position + 1
                )));
            }
            if let Some(required) = self.forced_partner[position] {
                if actual_partner != Some(required) {
                    return Err(invalid(format!(
                        "structure does not contain forced pair {}-{}",
                        position + 1,
                        required + 1
                    )));
                }
            }
        }
        for (i, &actual_partner) in partner.iter().enumerate().take(self.length) {
            if let Some(j) = actual_partner.filter(|&j| i < j) {
                if j - i <= min_loop || !self.allows_pair(bases, i, j, model) {
                    return Err(invalid(format!(
                        "structure contains disallowed pair {}-{}",
                        i + 1,
                        j + 1
                    )));
                }
                if self.no_lonely_pairs {
                    let inner = i + 1 < j && partner[i + 1] == Some(j - 1);
                    let outer = i > 0 && j + 1 < self.length && partner[i - 1] == Some(j + 1);
                    if !inner && !outer {
                        return Err(invalid(format!(
                            "structure contains lonely pair {}-{} while no_lonely_pairs is enabled",
                            i + 1,
                            j + 1
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    pub fn structure_energy(&self, partner: &[Option<usize>]) -> f64 {
        let mut energy = 0.0;
        for (i, paired) in partner.iter().enumerate() {
            if paired.is_none() {
                energy += self.unpaired_energy[i];
            }
        }
        for i in 0..self.length {
            if let Some(j) = partner[i].filter(|&j| i < j) {
                energy += self.pair_energy(i, j);
                if i + 1 < j && partner[i + 1] == Some(j - 1) {
                    energy += self.stack_energy(i, j, i + 1, j - 1);
                }
            }
        }
        energy
    }
}

fn apply_deigan(config: &ProbingConfig, stack: &mut [f64]) -> Result<(), RnaError> {
    let slope = finite(config.slope, "Deigan slope")?;
    let intercept = finite(config.intercept, "Deigan intercept")?;
    for (index, value) in config.reactivities.iter().enumerate() {
        if let Some(value) = value.filter(|value| *value >= 0.0) {
            stack[index] += quantize_soft(slope * value.ln_1p() + intercept);
        }
    }
    Ok(())
}

fn apply_zarringhalam(
    config: &ProbingConfig,
    unpaired: &mut [f64],
    paired: &mut [f64],
) -> Result<(), RnaError> {
    let beta = finite(config.beta, "Zarringhalam beta")?;
    let default_probability = finite(
        config.default_probability,
        "Zarringhalam default_probability",
    )?;
    if !(0.0..=1.0).contains(&default_probability) {
        return Err(invalid(
            "Zarringhalam default_probability must be between 0 and 1",
        ));
    }
    let observed_max = config
        .reactivities
        .iter()
        .flatten()
        .copied()
        .filter(|value| value.is_finite() && *value >= 0.0)
        .fold(0.0f64, f64::max);
    for (index, value) in config.reactivities.iter().enumerate() {
        let probability = match value {
            Some(value) if value.is_finite() && *value >= 0.0 => {
                zarringhalam_probability(*value, observed_max, &config.conversion)?
            }
            _ => default_probability,
        };
        unpaired[index] += quantize_soft(beta * (probability - 1.0).abs());
        paired[index] += quantize_soft(beta * probability.abs());
    }
    Ok(())
}

fn zarringhalam_probability(value: f64, maximum: f64, conversion: &str) -> Result<f64, RnaError> {
    let code = conversion
        .chars()
        .next()
        .unwrap_or('O')
        .to_ascii_uppercase();
    let probability = match code {
        'S' => value,
        'C' => (value >= 0.25) as u8 as f64,
        'L' => (value - 0.2) / 0.68,
        'O' => {
            if value <= 0.0 {
                0.0
            } else {
                (value.ln() + 2.29) / 1.6
            }
        }
        'M' => {
            let points = [(0.0, 0.0), (0.25, 0.35), (0.3, 0.55), (0.7, 0.85)];
            if value <= 0.0 {
                0.0
            } else if value >= maximum || maximum <= 0.7 {
                1.0
            } else if value > 0.7 {
                0.85 + (value - 0.7) / (maximum - 0.7) * 0.15
            } else {
                let mut mapped = 0.0;
                for window in points.windows(2) {
                    let (x0, y0) = window[0];
                    let (x1, y1) = window[1];
                    if value <= x1 {
                        mapped = y0 + (value - x0) / (x1 - x0) * (y1 - y0);
                        break;
                    }
                }
                mapped
            }
        }
        _ => {
            return Err(invalid(
                "Zarringhalam conversion must start with O, L, M, C, or S",
            ))
        }
    };
    Ok(probability.clamp(0.0, 1.0))
}

fn index(position: usize, length: usize, field: &str) -> Result<usize, RnaError> {
    if position == 0 || position > length {
        Err(invalid(format!(
            "{field} position {position} is outside 1..={length}"
        )))
    } else {
        Ok(position - 1)
    }
}

fn pair_indices(
    pair: PairConstraint,
    length: usize,
    field: &str,
) -> Result<(usize, usize), RnaError> {
    let i = index(pair.i, length, field)?;
    let j = index(pair.j, length, field)?;
    if i == j {
        return Err(invalid(format!(
            "{field} cannot pair a position with itself"
        )));
    }
    Ok((i.min(j), i.max(j)))
}

fn finite(value: f64, name: &str) -> Result<f64, RnaError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("{name} must be finite")))
    }
}

fn quantize_soft(value: f64) -> f64 {
    // ViennaRNA's static soft-constraint arrays store deca-calorie integers
    // and use roundf(value * 100) at the public API boundary.
    ((value as f32 * 100.0).round() / 100.0) as f64
}

fn invalid(message: impl Into<String>) -> RnaError {
    RnaError::InvalidOption(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_contradictory_and_crossing_hard_constraints() {
        let contradictory = ConstraintConfig {
            force_unpaired: vec![2],
            force_paired: vec![2],
            ..ConstraintConfig::default()
        };
        assert!(ConstraintModel::compile(8, &contradictory).is_err());

        let crossing = ConstraintConfig {
            force_pairs: vec![PairConstraint { i: 1, j: 6 }, PairConstraint { i: 3, j: 8 }],
            ..ConstraintConfig::default()
        };
        assert!(matches!(
            ConstraintModel::compile(8, &crossing),
            Err(RnaError::PseudoknotUnsupported(_))
        ));
    }

    #[test]
    fn deigan_and_zarringhalam_compile_to_shared_pseudo_energies() {
        let deigan = ConstraintConfig {
            probing: Some(ProbingConfig {
                reactivities: vec![Some(0.0), Some(1.0)],
                ..ProbingConfig::default()
            }),
            ..ConstraintConfig::default()
        };
        let compiled = ConstraintModel::compile(2, &deigan).unwrap();
        assert!((compiled.stack_energy[0] + 0.6).abs() < 1.0e-6);
        assert!((compiled.stack_energy[1] - 0.65).abs() < 1.0e-6);

        let zarringhalam = ConstraintConfig {
            probing: Some(ProbingConfig {
                method: "zarringhalam".into(),
                conversion: "S".into(),
                reactivities: vec![Some(0.0), Some(1.0)],
                ..ProbingConfig::default()
            }),
            ..ConstraintConfig::default()
        };
        let compiled = ConstraintModel::compile(2, &zarringhalam).unwrap();
        assert!((compiled.unpaired_energy[0] - 0.89).abs() < 1.0e-6);
        assert!((compiled.paired_energy[1] - 0.89).abs() < 1.0e-6);
    }
}
