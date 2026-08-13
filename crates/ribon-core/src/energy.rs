use crate::constraints::{ConstraintModel, ConstraintSummary};
use crate::salt::{SaltCorrections, STANDARD_MOLAR};
use crate::structure::{is_pseudoknotted, parse_structure, ParsedStructure, RnaError};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

mod turner2004 {
    include!(concat!(env!("OUT_DIR"), "/turner2004_generated.rs"));
}
mod mathews2004_dna {
    include!(concat!(env!("OUT_DIR"), "/mathews2004_dna_generated.rs"));
}

/// Nucleic-acid parameter family used by the thermodynamic engine.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum NucleicAcid {
    #[default]
    Rna,
    Dna,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SpecialLoopParameter {
    pub sequence: String,
    pub free_energy_37_centi_kcal_mol: i32,
    pub enthalpy_centi_kcal_mol: i32,
}

/// Optional normalized table replacements for a custom thermodynamic model.
/// Integer values use 0.01 kcal/mol and retain `10_000_000` as infinity.
/// Omitted fields inherit the selected RNA or DNA base family exactly.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct ThermodynamicParameterOverrides {
    pub schema_version: u32,
    pub name: String,
    pub fingerprint_sha256: String,
    pub stack_37: Option<Vec<i32>>,
    pub stack_dh: Option<Vec<i32>>,
    pub mismatch_h_37: Option<Vec<i32>>,
    pub mismatch_h_dh: Option<Vec<i32>>,
    pub mismatch_i_37: Option<Vec<i32>>,
    pub mismatch_i_dh: Option<Vec<i32>>,
    pub mismatch_1n_37: Option<Vec<i32>>,
    pub mismatch_1n_dh: Option<Vec<i32>>,
    pub mismatch_23_37: Option<Vec<i32>>,
    pub mismatch_23_dh: Option<Vec<i32>>,
    pub mismatch_m_37: Option<Vec<i32>>,
    pub mismatch_m_dh: Option<Vec<i32>>,
    pub mismatch_ext_37: Option<Vec<i32>>,
    pub mismatch_ext_dh: Option<Vec<i32>>,
    pub dangle5_37: Option<Vec<i32>>,
    pub dangle5_dh: Option<Vec<i32>>,
    pub dangle3_37: Option<Vec<i32>>,
    pub dangle3_dh: Option<Vec<i32>>,
    pub int11_37: Option<Vec<i32>>,
    pub int11_dh: Option<Vec<i32>>,
    pub int21_37: Option<Vec<i32>>,
    pub int21_dh: Option<Vec<i32>>,
    pub int22_37: Option<Vec<i32>>,
    pub int22_dh: Option<Vec<i32>>,
    pub hairpin_37: Option<Vec<i32>>,
    pub hairpin_dh: Option<Vec<i32>>,
    pub bulge_37: Option<Vec<i32>>,
    pub bulge_dh: Option<Vec<i32>>,
    pub internal_37: Option<Vec<i32>>,
    pub internal_dh: Option<Vec<i32>>,
    pub ml_params: Option<Vec<i32>>,
    pub ninio: Option<Vec<i32>>,
    pub misc: Option<Vec<i32>>,
    pub lxc_37: Option<f64>,
    pub triloops: Option<Vec<SpecialLoopParameter>>,
    pub tetraloops: Option<Vec<SpecialLoopParameter>>,
    pub hexaloops: Option<Vec<SpecialLoopParameter>>,
}

pub fn validate_parameter_overrides(
    profile: &ThermodynamicParameterOverrides,
) -> Result<(), RnaError> {
    if profile.schema_version != 1 {
        return Err(RnaError::InvalidOption(format!(
            "unsupported normalized parameter schema {}; expected 1",
            profile.schema_version
        )));
    }
    if profile.name.trim().is_empty() {
        return Err(RnaError::InvalidOption(
            "custom parameter profile name must not be empty".into(),
        ));
    }
    if profile.fingerprint_sha256.len() != 64
        || !profile
            .fingerprint_sha256
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit())
    {
        return Err(RnaError::InvalidOption(
            "custom parameter fingerprint_sha256 must contain 64 hexadecimal characters".into(),
        ));
    }
    fn table(values: &Option<Vec<i32>>, name: &str, expected: usize) -> Result<(), RnaError> {
        if let Some(values) = values {
            if values.len() != expected {
                return Err(RnaError::InvalidOption(format!(
                    "custom parameter table {name} has length {}; expected {expected}",
                    values.len()
                )));
            }
            if values
                .iter()
                .any(|value| i64::from(*value).abs() > i64::from(INF_PARAMETER))
            {
                return Err(RnaError::InvalidOption(format!(
                    "custom parameter table {name} contains a value outside the normalized range"
                )));
            }
        }
        Ok(())
    }
    macro_rules! tables {
        ($(($field:ident, $length:expr)),+ $(,)?) => {
            $(table(&profile.$field, stringify!($field), $length)?;)+
        };
    }
    tables!(
        (stack_37, 49),
        (stack_dh, 49),
        (mismatch_h_37, 175),
        (mismatch_h_dh, 175),
        (mismatch_i_37, 175),
        (mismatch_i_dh, 175),
        (mismatch_1n_37, 175),
        (mismatch_1n_dh, 175),
        (mismatch_23_37, 175),
        (mismatch_23_dh, 175),
        (mismatch_m_37, 175),
        (mismatch_m_dh, 175),
        (mismatch_ext_37, 175),
        (mismatch_ext_dh, 175),
        (dangle5_37, 35),
        (dangle5_dh, 35),
        (dangle3_37, 35),
        (dangle3_dh, 35),
        (int11_37, 1225),
        (int11_dh, 1225),
        (int21_37, 6125),
        (int21_dh, 6125),
        (int22_37, 9216),
        (int22_dh, 9216),
        (hairpin_37, 31),
        (hairpin_dh, 31),
        (bulge_37, 31),
        (bulge_dh, 31),
        (internal_37, 31),
        (internal_dh, 31),
        (ml_params, 6),
        (ninio, 3),
        (misc, 4),
    );
    if profile.lxc_37.is_some_and(|value| !value.is_finite()) {
        return Err(RnaError::InvalidOption(
            "custom lxc_37 must be finite".into(),
        ));
    }
    for (name, expected_length, values) in [
        ("triloops", 5, &profile.triloops),
        ("tetraloops", 6, &profile.tetraloops),
        ("hexaloops", 8, &profile.hexaloops),
    ] {
        if let Some(values) = values {
            let mut seen = std::collections::HashSet::new();
            for entry in values {
                let sequence = entry.sequence.to_ascii_uppercase().replace('T', "U");
                if sequence.len() != expected_length
                    || !sequence
                        .bytes()
                        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'U'))
                    || !seen.insert(sequence)
                    || i64::from(entry.free_energy_37_centi_kcal_mol).abs()
                        > i64::from(INF_PARAMETER)
                    || i64::from(entry.enthalpy_centi_kcal_mol).abs() > i64::from(INF_PARAMETER)
                {
                    return Err(RnaError::InvalidOption(format!(
                        "custom {name} contains an invalid or duplicate motif"
                    )));
                }
            }
        }
    }
    Ok(())
}

macro_rules! params {
    ($model:expr, $name:ident, $field:ident) => {
        if let Some(custom) = &$model.parameter_overrides {
            if let Some(values) = &custom.$field {
                values.as_slice()
            } else {
                match $model.nucleic_acid {
                    NucleicAcid::Rna => &turner2004::$name[..],
                    NucleicAcid::Dna => &mathews2004_dna::$name[..],
                }
            }
        } else {
            match $model.nucleic_acid {
                NucleicAcid::Rna => &turner2004::$name[..],
                NucleicAcid::Dna => &mathews2004_dna::$name[..],
            }
        }
    };
}

macro_rules! lxc_parameter {
    ($model:expr) => {
        if let Some(value) = $model
            .parameter_overrides
            .as_ref()
            .and_then(|custom| custom.lxc_37)
        {
            value
        } else {
            match $model.nucleic_acid {
                NucleicAcid::Rna => turner2004::LXC_37,
                NucleicAcid::Dna => mathews2004_dna::LXC_37,
            }
        }
    };
}

const INF_PARAMETER: i32 = 10_000_000;
const T_MEASURE_KELVIN: f64 = 310.15;
pub const CUSTOM_MODEL_ID: &str = "ribon-custom-thermodynamic-v1";

/// Turner 2004 nearest-neighbor parameters at a selected temperature.
///
/// The embedded source tables are generated from the RNAstructure 6.6
/// `data_tables/rna.*` bundle. MFE folding and fixed-structure evaluation
/// support dangle models 0 through 3. The partition ensemble uses the
/// standard polynomial grammars for models 0 and 2. Models 1 and 3 use exact
/// fixed-structure enumeration wherever their shared-dangle/coaxial state does
/// not admit the ordinary local factorization.
#[derive(Clone, Debug)]
pub struct EnergyModel {
    temperature_celsius: f64,
    dangles: u8,
    salt: SaltCorrections,
    nucleic_acid: NucleicAcid,
    parameter_overrides: Option<Arc<ThermodynamicParameterOverrides>>,
    /// Optional caller-selected internal-loop size restriction.  `None`
    /// enumerates every geometrically possible loop; the Turner logarithmic
    /// extension is defined beyond the 30-entry tabulated range.
    pub max_internal_loop: Option<usize>,
}

#[derive(Clone, Copy)]
enum SpecialLoopKind {
    Tri,
    Tetra,
    Hexa,
}

impl Default for EnergyModel {
    fn default() -> Self {
        Self {
            temperature_celsius: 37.0,
            dangles: 2,
            salt: SaltCorrections::new(STANDARD_MOLAR, 37.0),
            nucleic_acid: NucleicAcid::Rna,
            parameter_overrides: None,
            max_internal_loop: None,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct EnergyBreakdown {
    pub total_kcal_mol: f64,
    pub pair_kcal_mol: f64,
    pub stack_kcal_mol: f64,
    pub hairpin_kcal_mol: f64,
    pub internal_loop_kcal_mol: f64,
    pub multiloop_kcal_mol: f64,
    pub exterior_kcal_mol: f64,
    pub constraint_kcal_mol: f64,
    pub model: &'static str,
    pub temperature_celsius: f64,
    pub dangle_model: u8,
    pub salt_molar: f64,
    pub salt_correction: bool,
    /// Coaxial helix pairs selected by the fixed-structure evaluator.
    pub coaxial_stacks: Vec<CoaxialStack>,
    /// Per-loop terms suitable for differential diagnostics and annotation.
    pub loop_energies: Vec<LoopEnergy>,
    pub constraints: ConstraintSummary,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct CoaxialStack {
    /// One-based closing pair of the multiloop containing this interaction.
    pub loop_i: usize,
    pub loop_j: usize,
    /// One-based base-pair coordinates of the first helix around the loop.
    pub first_i: usize,
    pub first_j: usize,
    /// One-based base-pair coordinates of the adjacent helix.
    pub second_i: usize,
    pub second_j: usize,
    /// Full replacement contribution: stack + 2 * MLintern[GC].
    pub energy_kcal_mol: f64,
    /// Change relative to scoring the two stems independently.
    pub stabilization_kcal_mol: f64,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct LoopEnergy {
    pub kind: &'static str,
    /// One-based closing pair; zero denotes the exterior loop.
    pub i: usize,
    pub j: usize,
    pub energy_kcal_mol: f64,
}

impl EnergyBreakdown {
    fn new(
        temperature_celsius: f64,
        dangles: u8,
        salt_molar: f64,
        model: &'static str,
        constraints: ConstraintSummary,
    ) -> Self {
        Self {
            total_kcal_mol: 0.0,
            pair_kcal_mol: 0.0,
            stack_kcal_mol: 0.0,
            hairpin_kcal_mol: 0.0,
            internal_loop_kcal_mol: 0.0,
            multiloop_kcal_mol: 0.0,
            exterior_kcal_mol: 0.0,
            constraint_kcal_mol: 0.0,
            model,
            temperature_celsius,
            dangle_model: dangles,
            salt_molar,
            salt_correction: (salt_molar - STANDARD_MOLAR).abs() >= 1.0e-12,
            coaxial_stacks: Vec::new(),
            loop_energies: Vec::new(),
            constraints,
        }
    }

    fn finish(mut self) -> Self {
        self.total_kcal_mol = self.pair_kcal_mol
            + self.stack_kcal_mol
            + self.hairpin_kcal_mol
            + self.internal_loop_kcal_mol
            + self.multiloop_kcal_mol
            + self.exterior_kcal_mol
            + self.constraint_kcal_mol;
        self
    }
}

impl EnergyModel {
    /// Construct the standard default model (`dangles=2`).
    pub fn new(temperature_celsius: f64) -> Result<Self, RnaError> {
        Self::with_dangles(temperature_celsius, 2)
    }

    /// Construct a model with an explicit dangle convention.
    pub fn with_dangles(temperature_celsius: f64, dangles: u8) -> Result<Self, RnaError> {
        Self::with_dangles_and_salt(temperature_celsius, dangles, STANDARD_MOLAR)
    }

    /// Construct a model with explicit dangles and monovalent salt molarity.
    pub fn with_dangles_and_salt(
        temperature_celsius: f64,
        dangles: u8,
        salt_molar: f64,
    ) -> Result<Self, RnaError> {
        Self::with_parameter_family(temperature_celsius, dangles, salt_molar, NucleicAcid::Rna)
    }

    /// Construct an RNA or DNA model from the independently generated,
    /// embedded RNAstructure 6.6 parameter bundle.
    pub fn with_parameter_family(
        temperature_celsius: f64,
        dangles: u8,
        salt_molar: f64,
        nucleic_acid: NucleicAcid,
    ) -> Result<Self, RnaError> {
        if !temperature_celsius.is_finite() || temperature_celsius <= -273.15 {
            return Err(RnaError::InvalidOption(
                "temperature must be finite and above absolute zero".into(),
            ));
        }
        if dangles > 3 {
            return Err(RnaError::InvalidOption(format!(
                "dangles must be one of 0, 1, 2, or 3, got {dangles}"
            )));
        }
        if !salt_molar.is_finite() || salt_molar <= 0.0 {
            return Err(RnaError::InvalidOption(
                "salt molarity must be finite and positive".into(),
            ));
        }
        Ok(Self {
            temperature_celsius,
            dangles,
            salt: SaltCorrections::new(salt_molar, temperature_celsius),
            nucleic_acid,
            parameter_overrides: None,
            max_internal_loop: None,
        })
    }

    /// Construct a validated custom profile by overlaying normalized table
    /// replacements on one complete built-in family.
    pub fn with_parameter_overrides(
        temperature_celsius: f64,
        dangles: u8,
        salt_molar: f64,
        nucleic_acid: NucleicAcid,
        overrides: ThermodynamicParameterOverrides,
    ) -> Result<Self, RnaError> {
        validate_parameter_overrides(&overrides)?;
        let mut model =
            Self::with_parameter_family(temperature_celsius, dangles, salt_molar, nucleic_acid)?;
        model.parameter_overrides = Some(Arc::new(overrides));
        Ok(model)
    }

    pub fn temperature_celsius(&self) -> f64 {
        self.temperature_celsius
    }

    pub fn dangles(&self) -> u8 {
        self.dangles
    }

    pub fn with_dangle_model(&self, dangles: u8) -> Result<Self, RnaError> {
        if dangles > 3 {
            return Err(RnaError::InvalidOption(format!(
                "dangles must be one of 0, 1, 2, or 3, got {dangles}"
            )));
        }
        let mut model = self.clone();
        model.dangles = dangles;
        Ok(model)
    }

    pub fn nucleic_acid(&self) -> NucleicAcid {
        self.nucleic_acid
    }

    pub fn parameter_model_id(&self) -> &'static str {
        if self.parameter_overrides.is_some() {
            return CUSTOM_MODEL_ID;
        }
        match self.nucleic_acid {
            NucleicAcid::Rna => "ribon-turner-2004",
            NucleicAcid::Dna => "ribon-mathews-dna-2004",
        }
    }

    pub fn parameter_profile_name(&self) -> Option<&str> {
        self.parameter_overrides
            .as_ref()
            .map(|profile| profile.name.as_str())
    }

    pub fn parameter_profile_fingerprint(&self) -> Option<&str> {
        self.parameter_overrides
            .as_ref()
            .map(|profile| profile.fingerprint_sha256.as_str())
    }

    pub fn salt_molar(&self) -> f64 {
        self.salt.molar
    }

    pub fn with_salt_molar(&self, salt_molar: f64) -> Result<Self, RnaError> {
        if !salt_molar.is_finite() || salt_molar <= 0.0 {
            return Err(RnaError::InvalidOption(
                "salt molarity must be finite and positive".into(),
            ));
        }
        let mut model = self.clone();
        model.salt = SaltCorrections::new(salt_molar, self.temperature_celsius);
        Ok(model)
    }

    pub fn has_salt_correction(&self) -> bool {
        (self.salt.molar - STANDARD_MOLAR).abs() >= 1.0e-12
    }

    pub fn supports_partition(&self) -> bool {
        true
    }

    /// Return the effective internal-loop bound for one finite interval.
    ///
    /// The default is the complete geometric range.  A finite value is only
    /// used when a Rust caller explicitly opts into a restricted model.
    pub fn internal_loop_limit(&self, available: usize) -> usize {
        self.max_internal_loop.unwrap_or(available).min(available)
    }

    pub fn model_name(&self) -> &'static str {
        if self.parameter_overrides.is_some() {
            return "custom normalized nearest-neighbor model";
        }
        energy_model_name(self.nucleic_acid, self.dangles)
    }

    pub fn ensemble_model_name(&self) -> &'static str {
        if self.parameter_overrides.is_some() {
            return if self.dangles % 2 == 0 {
                "custom normalized nearest-neighbor McCaskill ensemble"
            } else {
                "custom normalized exact fixed-structure ensemble"
            };
        }
        ensemble_model_name(self.nucleic_acid, self.dangles)
    }

    pub fn can_pair(&self, a: u8, b: u8) -> bool {
        pair_type(a, b).is_some()
    }

    /// Pair formation has no context-free contribution in the Turner model.
    pub fn pair_energy(&self, a: u8, b: u8) -> f64 {
        if self.can_pair(a, b) {
            0.0
        } else {
            f64::INFINITY
        }
    }

    fn special_loop_energy(&self, kind: SpecialLoopKind, sequence: &str) -> Option<(i32, i32)> {
        if let Some(custom) = &self.parameter_overrides {
            let replacement = match kind {
                SpecialLoopKind::Tri => &custom.triloops,
                SpecialLoopKind::Tetra => &custom.tetraloops,
                SpecialLoopKind::Hexa => &custom.hexaloops,
            };
            if let Some(replacement) = replacement {
                return replacement.iter().find_map(|entry| {
                    let normalized = entry.sequence.to_ascii_uppercase().replace('T', "U");
                    normalized.eq(sequence).then_some((
                        entry.free_energy_37_centi_kcal_mol,
                        entry.enthalpy_centi_kcal_mol,
                    ))
                });
            }
        }
        let table = match (self.nucleic_acid, kind) {
            (NucleicAcid::Rna, SpecialLoopKind::Tri) => turner2004::TRILOOPS,
            (NucleicAcid::Rna, SpecialLoopKind::Tetra) => turner2004::TETRALOOPS,
            (NucleicAcid::Rna, SpecialLoopKind::Hexa) => turner2004::HEXALOOPS,
            (NucleicAcid::Dna, SpecialLoopKind::Tri) => mathews2004_dna::TRILOOPS,
            (NucleicAcid::Dna, SpecialLoopKind::Tetra) => mathews2004_dna::TETRALOOPS,
            (NucleicAcid::Dna, SpecialLoopKind::Hexa) => mathews2004_dna::HEXALOOPS,
        };
        find_special(table, sequence)
    }

    fn temperature_ratio(&self) -> f64 {
        (self.temperature_celsius + 273.15) / T_MEASURE_KELVIN
    }

    fn scaled(&self, g37: i32, enthalpy: i32) -> f64 {
        if g37 >= INF_PARAMETER || enthalpy >= INF_PARAMETER {
            return f64::INFINITY;
        }
        // ViennaRNA's default non-smooth parameter path truncates scaled
        // centi-kcal values toward zero before MFE and Boltzmann conversion.
        let value = enthalpy as f64 - (enthalpy as f64 - g37 as f64) * self.temperature_ratio();
        value.trunc() / 100.0
    }

    fn scaled_centi_raw(&self, g37: i32, enthalpy: i32) -> f64 {
        enthalpy as f64 - (enthalpy as f64 - g37 as f64) * self.temperature_ratio()
    }

    fn scaled_boltzmann(&self, g37: i32, enthalpy: i32) -> f64 {
        if g37 >= INF_PARAMETER || enthalpy >= INF_PARAMETER {
            return f64::INFINITY;
        }
        self.scaled_centi_raw(g37, enthalpy) / 100.0
    }

    fn table(&self, g37: &[i32], enthalpy: &[i32], index: usize) -> f64 {
        self.scaled(g37[index], enthalpy[index])
    }

    fn table_boltzmann(&self, g37: &[i32], enthalpy: &[i32], index: usize) -> f64 {
        self.scaled_boltzmann(g37[index], enthalpy[index])
    }

    pub fn stack_energy(&self, a: u8, b: u8, c: u8, d: u8) -> f64 {
        let Some(outer) = pair_type(a, b) else {
            return f64::INFINITY;
        };
        let Some(inner) = pair_type(d, c) else {
            return f64::INFINITY;
        };
        self.table(
            params!(self, STACK_37, stack_37),
            params!(self, STACK_DH, stack_dh),
            outer * 7 + inner,
        ) + self.salt.stack_centi as f64 / 100.0
    }

    fn loop_initiation(&self, g37: &[i32], enthalpy: &[i32], size: usize) -> f64 {
        if size <= 30 {
            self.table(g37, enthalpy, size)
        } else {
            let extension_centi = (lxc_parameter!(self)
                * 100.0
                * self.temperature_ratio()
                * (size as f64 / 30.0).ln())
            .trunc();
            self.table(g37, enthalpy, 30) + extension_centi / 100.0
        }
    }

    fn loop_initiation_boltzmann(&self, g37: &[i32], enthalpy: &[i32], size: usize) -> f64 {
        if size <= 30 {
            self.table_boltzmann(g37, enthalpy, size)
        } else {
            let extension_centi =
                lxc_parameter!(self) * 100.0 * self.temperature_ratio() * (size as f64 / 30.0).ln();
            self.table_boltzmann(g37, enthalpy, 30) + extension_centi / 100.0
        }
    }

    pub fn hairpin_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        let size = j - i - 1;
        let salt = self.salt.loop_centi(size + 1, self.temperature_celsius) as f64 / 100.0;
        let initiation = self.loop_initiation(
            params!(self, HAIRPIN_37, hairpin_37),
            params!(self, HAIRPIN_DH, hairpin_dh),
            size,
        );
        if !initiation.is_finite() || size < 3 {
            return initiation + salt;
        }

        let motif = std::str::from_utf8(&bases[i..=j]).ok();
        let special = match size {
            3 => {
                motif.and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Tri, sequence))
            }
            4 => motif
                .and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Tetra, sequence)),
            6 => {
                motif.and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Hexa, sequence))
            }
            _ => None,
        };
        if let Some((g37, enthalpy)) = special {
            return self.scaled(g37, enthalpy) + salt;
        }
        if size == 3 {
            return initiation + self.terminal_au_for_type(pair) + salt;
        }

        let first = base_index(bases[i + 1]);
        let last = base_index(bases[j - 1]);
        initiation
            + self.table(
                params!(self, MISMATCH_H_37, mismatch_h_37),
                params!(self, MISMATCH_H_DH, mismatch_h_dh),
                index_mismatch(pair, first, last),
            )
            + salt
    }

    pub(crate) fn hairpin_boltzmann_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        let size = j - i - 1;
        let salt = self.salt.loop_centi(size + 1, self.temperature_celsius) as f64 / 100.0;
        let initiation = self.loop_initiation_boltzmann(
            params!(self, HAIRPIN_37, hairpin_37),
            params!(self, HAIRPIN_DH, hairpin_dh),
            size,
        );
        if !initiation.is_finite() || size < 3 {
            return initiation + salt;
        }

        let motif = std::str::from_utf8(&bases[i..=j]).ok();
        let special = match size {
            3 => {
                motif.and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Tri, sequence))
            }
            4 => motif
                .and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Tetra, sequence)),
            6 => {
                motif.and_then(|sequence| self.special_loop_energy(SpecialLoopKind::Hexa, sequence))
            }
            _ => None,
        };
        if let Some((g37, enthalpy)) = special {
            return self.scaled_boltzmann(g37, enthalpy) + salt;
        }
        if size == 3 {
            return initiation + self.terminal_au_for_type_boltzmann(pair) + salt;
        }

        let first = base_index(bases[i + 1]);
        let last = base_index(bases[j - 1]);
        initiation
            + self.table_boltzmann(
                params!(self, MISMATCH_H_37, mismatch_h_37),
                params!(self, MISMATCH_H_DH, mismatch_h_dh),
                index_mismatch(pair, first, last),
            )
            + salt
    }

    pub fn internal_energy(&self, bases: &[u8], i: usize, j: usize, k: usize, l: usize) -> f64 {
        if !(i < k && k < l && l < j && j < bases.len()) {
            return f64::INFINITY;
        }
        let Some(outer) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        let Some(inner) = pair_type(bases[l], bases[k]) else {
            return f64::INFINITY;
        };
        let n1 = k - i - 1;
        let n2 = j - l - 1;
        let small = n1.min(n2);
        let large = n1.max(n2);
        if large == 0 {
            return self.stack_energy(bases[i], bases[j], bases[k], bases[l]);
        }
        let salt = self.salt.loop_centi(n1 + n2 + 2, self.temperature_celsius) as f64 / 100.0;

        if small == 0 {
            let mut energy = self.loop_initiation(
                params!(self, BULGE_37, bulge_37),
                params!(self, BULGE_DH, bulge_dh),
                large,
            );
            if large == 1 {
                energy += self.table(
                    params!(self, STACK_37, stack_37),
                    params!(self, STACK_DH, stack_dh),
                    outer * 7 + inner,
                );
            } else {
                energy += self.terminal_au_for_type(outer);
                energy += self.terminal_au_for_type(inner);
            }
            return energy + salt;
        }

        let si1 = base_index(bases[i + 1]);
        let sj1 = base_index(bases[j - 1]);
        let sp1 = base_index(bases[k - 1]);
        let sq1 = base_index(bases[l + 1]);
        let energy = match (small, large) {
            (1, 1) => self.table(
                params!(self, INT11_37, int11_37),
                params!(self, INT11_DH, int11_dh),
                index_int11(outer, inner, si1, sj1),
            ),
            (1, 2) if n1 == 1 => self.table(
                params!(self, INT21_37, int21_37),
                params!(self, INT21_DH, int21_dh),
                index_int21(outer, inner, si1, sq1, sj1),
            ),
            (1, 2) => self.table(
                params!(self, INT21_37, int21_37),
                params!(self, INT21_DH, int21_dh),
                index_int21(inner, outer, sq1, si1, sp1),
            ),
            (2, 2) => self.table(
                params!(self, INT22_37, int22_37),
                params!(self, INT22_DH, int22_dh),
                index_int22(outer, inner, si1, sp1, sq1, sj1),
            ),
            (2, 3) => {
                self.loop_initiation(
                    params!(self, INTERNAL_37, internal_37),
                    params!(self, INTERNAL_DH, internal_dh),
                    5,
                ) + self.ninio()
                    + self.mismatch_23(outer, si1, sj1)
                    + self.mismatch_23(inner, sq1, sp1)
            }
            (1, _) => {
                self.generic_internal_initiation(n1, n2)
                    + self.mismatch_1n(outer, si1, sj1)
                    + self.mismatch_1n(inner, sq1, sp1)
            }
            _ => {
                self.generic_internal_initiation(n1, n2)
                    + self.mismatch_internal(outer, si1, sj1)
                    + self.mismatch_internal(inner, sq1, sp1)
            }
        };
        energy + salt
    }

    pub(crate) fn internal_boltzmann_energy(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        k: usize,
        l: usize,
    ) -> f64 {
        if !(i < k && k < l && l < j && j < bases.len()) {
            return f64::INFINITY;
        }
        let Some(outer) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        let Some(inner) = pair_type(bases[l], bases[k]) else {
            return f64::INFINITY;
        };
        let n1 = k - i - 1;
        let n2 = j - l - 1;
        let small = n1.min(n2);
        let large = n1.max(n2);
        if large == 0 {
            return self.stack_boltzmann_energy(bases[i], bases[j], bases[k], bases[l]);
        }
        let salt = self.salt.loop_centi(n1 + n2 + 2, self.temperature_celsius) as f64 / 100.0;

        if small == 0 {
            let mut energy = self.loop_initiation_boltzmann(
                params!(self, BULGE_37, bulge_37),
                params!(self, BULGE_DH, bulge_dh),
                large,
            );
            if large == 1 {
                energy += self.table_boltzmann(
                    params!(self, STACK_37, stack_37),
                    params!(self, STACK_DH, stack_dh),
                    outer * 7 + inner,
                );
            } else {
                energy += self.terminal_au_for_type_boltzmann(outer);
                energy += self.terminal_au_for_type_boltzmann(inner);
            }
            return energy + salt;
        }

        let si1 = base_index(bases[i + 1]);
        let sj1 = base_index(bases[j - 1]);
        let sp1 = base_index(bases[k - 1]);
        let sq1 = base_index(bases[l + 1]);
        let energy = match (small, large) {
            (1, 1) => self.table_boltzmann(
                params!(self, INT11_37, int11_37),
                params!(self, INT11_DH, int11_dh),
                index_int11(outer, inner, si1, sj1),
            ),
            (1, 2) if n1 == 1 => self.table_boltzmann(
                params!(self, INT21_37, int21_37),
                params!(self, INT21_DH, int21_dh),
                index_int21(outer, inner, si1, sq1, sj1),
            ),
            (1, 2) => self.table_boltzmann(
                params!(self, INT21_37, int21_37),
                params!(self, INT21_DH, int21_dh),
                index_int21(inner, outer, sq1, si1, sp1),
            ),
            (2, 2) => self.table_boltzmann(
                params!(self, INT22_37, int22_37),
                params!(self, INT22_DH, int22_dh),
                index_int22(outer, inner, si1, sp1, sq1, sj1),
            ),
            (2, 3) => {
                self.loop_initiation_boltzmann(
                    params!(self, INTERNAL_37, internal_37),
                    params!(self, INTERNAL_DH, internal_dh),
                    5,
                ) + self.ninio_boltzmann()
                    + self.mismatch_23_boltzmann(outer, si1, sj1)
                    + self.mismatch_23_boltzmann(inner, sq1, sp1)
            }
            (1, _) => {
                self.generic_internal_initiation_boltzmann(n1, n2)
                    + self.mismatch_1n_boltzmann(outer, si1, sj1)
                    + self.mismatch_1n_boltzmann(inner, sq1, sp1)
            }
            _ => {
                self.generic_internal_initiation_boltzmann(n1, n2)
                    + self.mismatch_internal_boltzmann(outer, si1, sj1)
                    + self.mismatch_internal_boltzmann(inner, sq1, sp1)
            }
        };
        energy + salt
    }

    fn stack_boltzmann_energy(&self, a: u8, b: u8, c: u8, d: u8) -> f64 {
        let Some(outer) = pair_type(a, b) else {
            return f64::INFINITY;
        };
        let Some(inner) = pair_type(d, c) else {
            return f64::INFINITY;
        };
        self.table_boltzmann(
            params!(self, STACK_37, stack_37),
            params!(self, STACK_DH, stack_dh),
            outer * 7 + inner,
        ) + self.salt.stack_centi as f64 / 100.0
    }

    fn generic_internal_initiation_boltzmann(&self, n1: usize, n2: usize) -> f64 {
        let total = n1 + n2;
        let asymmetry = (n1.abs_diff(n2) as f64 * self.ninio_boltzmann()).min(self.max_ninio());
        self.loop_initiation_boltzmann(
            params!(self, INTERNAL_37, internal_37),
            params!(self, INTERNAL_DH, internal_dh),
            total,
        ) + asymmetry
    }

    fn mismatch_internal_boltzmann(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table_boltzmann(
            params!(self, MISMATCH_I_37, mismatch_i_37),
            params!(self, MISMATCH_I_DH, mismatch_i_dh),
            index_mismatch(pair, left, right),
        )
    }

    fn mismatch_1n_boltzmann(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table_boltzmann(
            params!(self, MISMATCH_1N_37, mismatch_1n_37),
            params!(self, MISMATCH_1N_DH, mismatch_1n_dh),
            index_mismatch(pair, left, right),
        )
    }

    fn mismatch_23_boltzmann(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table_boltzmann(
            params!(self, MISMATCH_23_37, mismatch_23_37),
            params!(self, MISMATCH_23_DH, mismatch_23_dh),
            index_mismatch(pair, left, right),
        )
    }

    fn generic_internal_initiation(&self, n1: usize, n2: usize) -> f64 {
        let total = n1 + n2;
        let asymmetry = (n1.abs_diff(n2) as f64 * self.ninio()).min(self.max_ninio());
        self.loop_initiation(
            params!(self, INTERNAL_37, internal_37),
            params!(self, INTERNAL_DH, internal_dh),
            total,
        ) + asymmetry
    }

    fn mismatch_internal(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table(
            params!(self, MISMATCH_I_37, mismatch_i_37),
            params!(self, MISMATCH_I_DH, mismatch_i_dh),
            index_mismatch(pair, left, right),
        )
    }

    fn mismatch_1n(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table(
            params!(self, MISMATCH_1N_37, mismatch_1n_37),
            params!(self, MISMATCH_1N_DH, mismatch_1n_dh),
            index_mismatch(pair, left, right),
        )
    }

    fn mismatch_23(&self, pair: usize, left: usize, right: usize) -> f64 {
        self.table(
            params!(self, MISMATCH_23_37, mismatch_23_37),
            params!(self, MISMATCH_23_DH, mismatch_23_dh),
            index_mismatch(pair, left, right),
        )
    }

    pub fn multiloop_closing(&self) -> f64 {
        self.scaled(
            params!(self, ML_PARAMS, ml_params)[2],
            params!(self, ML_PARAMS, ml_params)[3],
        ) + (self.salt.ml_base_centi + self.salt.ml_closing_centi) as f64 / 100.0
    }

    pub(crate) fn multiloop_closing_boltzmann(&self) -> f64 {
        self.scaled_boltzmann(
            params!(self, ML_PARAMS, ml_params)[2],
            params!(self, ML_PARAMS, ml_params)[3],
        ) + (self.salt.ml_base_centi + self.salt.ml_closing_centi) as f64 / 100.0
    }

    pub fn multiloop_unpaired(&self) -> f64 {
        self.scaled(
            params!(self, ML_PARAMS, ml_params)[0],
            params!(self, ML_PARAMS, ml_params)[1],
        ) + self.salt.ml_base_centi as f64 / 100.0
    }

    pub(crate) fn multiloop_unpaired_boltzmann(&self) -> f64 {
        self.scaled_boltzmann(
            params!(self, ML_PARAMS, ml_params)[0],
            params!(self, ML_PARAMS, ml_params)[1],
        ) + self.salt.ml_base_centi as f64 / 100.0
    }

    /// Context-free affine branch term used when a non-canonical, terminal
    /// substructure (currently a G-quadruplex) is embedded in a multiloop.
    /// Canonical stems add terminal-AU and dangle terms on top of this value.
    pub(crate) fn multiloop_branch(&self) -> f64 {
        self.scaled(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
    }

    pub(crate) fn multiloop_branch_boltzmann(&self) -> f64 {
        self.scaled_boltzmann(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
    }

    fn multiloop_stem_for_type(
        &self,
        pair: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        self.scaled(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
            + self.terminal_au_for_type(pair)
            + self.dangle_context(pair, five_prime, three_prime, true)
    }

    pub(crate) fn multiloop_stem_selected(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = bases
            .get(i)
            .zip(bases.get(j))
            .and_then(|(&a, &b)| pair_type(a, b))
        else {
            return f64::INFINITY;
        };
        self.multiloop_stem_for_type(pair, five_prime, three_prime)
    }

    pub(crate) fn multiloop_closing_stem_selected(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = bases
            .get(j)
            .zip(bases.get(i))
            .and_then(|(&a, &b)| pair_type(a, b))
        else {
            return f64::INFINITY;
        };
        self.multiloop_stem_for_type(pair, five_prime, three_prime)
    }

    pub(crate) fn exterior_stem_selected(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = bases
            .get(i)
            .zip(bases.get(j))
            .and_then(|(&a, &b)| pair_type(a, b))
        else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type(pair) + self.dangle_context(pair, five_prime, three_prime, false)
    }

    pub(crate) fn exterior_stem_boltzmann_selected(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = bases
            .get(i)
            .zip(bases.get(j))
            .and_then(|(&a, &b)| pair_type(a, b))
        else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type_boltzmann(pair)
            + self.dangle_context_boltzmann(pair, five_prime, three_prime, false)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn coaxial_energy(
        &self,
        bases: &[u8],
        first_i: usize,
        first_j: usize,
        first_reversed: bool,
        second_i: usize,
        second_j: usize,
        second_reversed: bool,
    ) -> f64 {
        let first = if first_reversed {
            pair_type(bases[first_j], bases[first_i])
        } else {
            pair_type(bases[first_i], bases[first_j])
        };
        let second = if second_reversed {
            pair_type(bases[second_j], bases[second_i])
        } else {
            pair_type(bases[second_i], bases[second_j])
        };
        let (Some(first), Some(second)) = (first, second) else {
            return f64::INFINITY;
        };
        self.table(
            params!(self, STACK_37, stack_37),
            params!(self, STACK_DH, stack_dh),
            reverse_pair_type(first) * 7 + reverse_pair_type(second),
        ) + 2.0
            * (self.scaled(
                params!(self, ML_PARAMS, ml_params)[4],
                params!(self, ML_PARAMS, ml_params)[5],
            ) + self.salt.ml_base_centi as f64 / 100.0)
    }

    /// Energy of a stem branching from a multiloop. The adjacent nucleotides
    /// are outside the pair `(i, j)`, matching `vrna_E_multibranch_stem()`.
    pub fn multiloop_stem_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        self.multiloop_stem_for_type(
            pair,
            i.checked_sub(1).map(|index| bases[index]),
            bases.get(j + 1).copied(),
        )
    }

    pub(crate) fn multiloop_stem_boltzmann_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        self.scaled_boltzmann(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
            + self.terminal_au_for_type_boltzmann(pair)
            + self.dangle_context_boltzmann(
                pair,
                i.checked_sub(1).map(|index| bases[index]),
                bases.get(j + 1).copied(),
                true,
            )
    }

    /// Stem contribution for the pair that closes a multiloop. The closing
    /// pair has reversed orientation, and its dangles lie inside the loop.
    pub fn multiloop_closing_stem_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[j], bases[i]) else {
            return f64::INFINITY;
        };
        self.multiloop_stem_for_type(
            pair,
            (j > i + 1).then_some(bases[j - 1]),
            (i + 1 < j).then_some(bases[i + 1]),
        )
    }

    pub(crate) fn multiloop_closing_stem_boltzmann_energy(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
    ) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[j], bases[i]) else {
            return f64::INFINITY;
        };
        self.scaled_boltzmann(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
            + self.terminal_au_for_type_boltzmann(pair)
            + self.dangle_context_boltzmann(
                pair,
                (j > i + 1).then_some(bases[j - 1]),
                (i + 1 < j).then_some(bases[i + 1]),
                true,
            )
    }

    pub fn exterior_stem_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type(pair)
            + self.dangle_context(
                pair,
                i.checked_sub(1).map(|index| bases[index]),
                bases.get(j + 1).copied(),
                false,
            )
    }

    pub(crate) fn exterior_stem_boltzmann_energy(&self, bases: &[u8], i: usize, j: usize) -> f64 {
        if i >= j || j >= bases.len() {
            return f64::INFINITY;
        }
        let Some(pair) = pair_type(bases[i], bases[j]) else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type_boltzmann(pair)
            + self.dangle_context_boltzmann(
                pair,
                i.checked_sub(1).map(|index| bases[index]),
                bases.get(j + 1).copied(),
                false,
            )
    }

    /// Exterior stem contribution for an explicitly oriented pair. This is
    /// used at the inner end of an intermolecular helix, where the strand cut
    /// reverses the ordinary `(i,j)` orientation.
    pub(crate) fn oriented_exterior_stem_energy(
        &self,
        five_base: u8,
        three_base: u8,
        five_dangle: Option<u8>,
        three_dangle: Option<u8>,
    ) -> f64 {
        let Some(pair) = pair_type(five_base, three_base) else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type(pair)
            + self.dangle_context(pair, five_dangle, three_dangle, false)
    }

    pub(crate) fn oriented_exterior_stem_boltzmann_energy(
        &self,
        five_base: u8,
        three_base: u8,
        five_dangle: Option<u8>,
        three_dangle: Option<u8>,
    ) -> f64 {
        let Some(pair) = pair_type(five_base, three_base) else {
            return f64::INFINITY;
        };
        self.terminal_au_for_type_boltzmann(pair)
            + self.dangle_context_boltzmann(pair, five_dangle, three_dangle, false)
    }

    pub(crate) fn cofold_exterior_stem_energy(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        cut: usize,
    ) -> f64 {
        self.oriented_exterior_stem_energy(
            bases[i],
            bases[j],
            i.checked_sub(1)
                .filter(|_| i != cut)
                .map(|position| bases[position]),
            (j + 1 < bases.len() && j + 1 != cut).then(|| bases[j + 1]),
        )
    }

    pub(crate) fn cofold_exterior_stem_boltzmann_energy(
        &self,
        bases: &[u8],
        i: usize,
        j: usize,
        cut: usize,
    ) -> f64 {
        self.oriented_exterior_stem_boltzmann_energy(
            bases[i],
            bases[j],
            i.checked_sub(1)
                .filter(|_| i != cut)
                .map(|position| bases[position]),
            (j + 1 < bases.len() && j + 1 != cut).then(|| bases[j + 1]),
        )
    }

    fn dangle_context(
        &self,
        pair: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
        multiloop: bool,
    ) -> f64 {
        if self.dangles == 0 {
            return 0.0;
        }

        match (five_prime, three_prime) {
            (Some(five), Some(three)) => {
                let (g37, dh) = if multiloop {
                    (
                        params!(self, MISMATCH_M_37, mismatch_m_37),
                        params!(self, MISMATCH_M_DH, mismatch_m_dh),
                    )
                } else {
                    (
                        params!(self, MISMATCH_EXT_37, mismatch_ext_37),
                        params!(self, MISMATCH_EXT_DH, mismatch_ext_dh),
                    )
                };
                self.table(
                    g37,
                    dh,
                    index_mismatch(pair, base_index(five), base_index(three)),
                )
                .min(0.0)
            }
            (Some(five), None) => self
                .table(
                    params!(self, DANGLE5_37, dangle5_37),
                    params!(self, DANGLE5_DH, dangle5_dh),
                    pair * 5 + base_index(five),
                )
                .min(0.0),
            (None, Some(three)) => self
                .table(
                    params!(self, DANGLE3_37, dangle3_37),
                    params!(self, DANGLE3_DH, dangle3_dh),
                    pair * 5 + base_index(three),
                )
                .min(0.0),
            (None, None) => 0.0,
        }
    }

    // The historical RNAduplex recurrence always supplies the nucleotides
    // adjacent to both interaction ends to vrna_E_exterior_stem(), independent
    // of the global dangle model. Keep that API-specific convention isolated
    // from the ordinary single-strand exterior-loop grammar.
    pub(crate) fn duplex_end_energy(
        &self,
        a: u8,
        b: u8,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = pair_type(a, b) else {
            return f64::INFINITY;
        };
        let context = match (five_prime, three_prime) {
            // vrna_params() zeros mismatchExt for the dangles=0 parameter
            // object used by RNAduplex, while leaving single-dangle tables
            // available to this odd-dangle recurrence.
            (Some(_), Some(_)) => 0.0,
            (Some(five), None) => self.table(
                params!(self, DANGLE5_37, dangle5_37),
                params!(self, DANGLE5_DH, dangle5_dh),
                pair * 5 + base_index(five),
            ),
            (None, Some(three)) => self.table(
                params!(self, DANGLE3_37, dangle3_37),
                params!(self, DANGLE3_DH, dangle3_dh),
                pair * 5 + base_index(three),
            ),
            (None, None) => 0.0,
        };
        self.terminal_au_for_type(pair) + context
    }

    pub(crate) fn duplex_end_boltzmann_energy(
        &self,
        a: u8,
        b: u8,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
    ) -> f64 {
        let Some(pair) = pair_type(a, b) else {
            return f64::INFINITY;
        };
        let (g37, dh) = match (five_prime, three_prime) {
            (Some(_), Some(_)) => {
                return self.terminal_au_for_type_boltzmann(pair);
            }
            (Some(five), None) => {
                let index = pair * 5 + base_index(five);
                (
                    params!(self, DANGLE5_37, dangle5_37)[index],
                    params!(self, DANGLE5_DH, dangle5_dh)[index],
                )
            }
            (None, Some(three)) => {
                let index = pair * 5 + base_index(three);
                (
                    params!(self, DANGLE3_37, dangle3_37)[index],
                    params!(self, DANGLE3_DH, dangle3_dh)[index],
                )
            }
            (None, None) => {
                return self.terminal_au_for_type_boltzmann(pair);
            }
        };
        let favorable_centi = -self.scaled_centi_raw(g37, dh);
        self.terminal_au_for_type_boltzmann(pair) - smooth_favorable_centi(favorable_centi) / 100.0
    }

    /// Effective energy whose ordinary Boltzmann factor equals ViennaRNA's
    /// smoothed PF dangle/mismatch factor (`pf_smooth=1`).
    fn dangle_context_boltzmann(
        &self,
        pair: usize,
        five_prime: Option<u8>,
        three_prime: Option<u8>,
        multiloop: bool,
    ) -> f64 {
        if self.dangles == 0 {
            return 0.0;
        }

        let (g37, dh) = match (five_prime, three_prime) {
            (Some(five), Some(three)) => {
                let (g37, dh) = if multiloop {
                    (
                        params!(self, MISMATCH_M_37, mismatch_m_37),
                        params!(self, MISMATCH_M_DH, mismatch_m_dh),
                    )
                } else {
                    (
                        params!(self, MISMATCH_EXT_37, mismatch_ext_37),
                        params!(self, MISMATCH_EXT_DH, mismatch_ext_dh),
                    )
                };
                let index = index_mismatch(pair, base_index(five), base_index(three));
                (g37[index], dh[index])
            }
            (Some(five), None) => {
                let index = pair * 5 + base_index(five);
                (
                    params!(self, DANGLE5_37, dangle5_37)[index],
                    params!(self, DANGLE5_DH, dangle5_dh)[index],
                )
            }
            (None, Some(three)) => {
                let index = pair * 5 + base_index(three);
                (
                    params!(self, DANGLE3_37, dangle3_37)[index],
                    params!(self, DANGLE3_DH, dangle3_dh)[index],
                )
            }
            (None, None) => return 0.0,
        };

        let favorable_centi = -self.scaled_centi_raw(g37, dh);
        -smooth_favorable_centi(favorable_centi) / 100.0
    }

    pub fn terminal_au(&self) -> f64 {
        self.scaled(params!(self, MISC, misc)[2], params!(self, MISC, misc)[3])
    }

    pub fn terminal_pair_energy(&self, a: u8, b: u8) -> f64 {
        pair_type(a, b)
            .map(|pair| self.terminal_au_for_type(pair))
            .unwrap_or(f64::INFINITY)
    }

    pub fn duplex_initialization_energy(&self) -> f64 {
        self.scaled(params!(self, MISC, misc)[0], params!(self, MISC, misc)[1])
            + self.salt.duplex_init_centi as f64 / 100.0
    }

    pub(crate) fn duplex_initialization_boltzmann_energy(&self) -> f64 {
        self.scaled_boltzmann(params!(self, MISC, misc)[0], params!(self, MISC, misc)[1])
            + self.salt.duplex_init_centi as f64 / 100.0
    }

    fn terminal_au_for_type(&self, pair: usize) -> f64 {
        if pair > 1 {
            self.terminal_au()
        } else {
            0.0
        }
    }

    fn terminal_au_for_type_boltzmann(&self, pair: usize) -> f64 {
        if pair > 1 {
            self.scaled_boltzmann(params!(self, MISC, misc)[2], params!(self, MISC, misc)[3])
        } else {
            0.0
        }
    }

    fn ninio(&self) -> f64 {
        self.scaled(
            params!(self, NINIO, ninio)[0],
            params!(self, NINIO, ninio)[1],
        )
    }

    fn ninio_boltzmann(&self) -> f64 {
        self.scaled_boltzmann(
            params!(self, NINIO, ninio)[0],
            params!(self, NINIO, ninio)[1],
        )
    }

    fn max_ninio(&self) -> f64 {
        params!(self, NINIO, ninio)[2] as f64 / 100.0
    }

    pub fn evaluate(&self, sequence: &str, structure: &str) -> Result<EnergyBreakdown, RnaError> {
        let parsed = parse_structure(sequence, structure)?;
        let constraints = ConstraintModel::unconstrained(parsed.length);
        self.evaluate_parsed(&parsed, 0, &constraints)
    }

    /// Evaluate a pseudoknot-free two-strand structure with one explicit cut.
    /// The face incident to the strand ends is an exterior face rather than a
    /// hairpin; all other closed faces retain the ordinary Turner loop rules.
    pub(crate) fn evaluate_cofold(
        &self,
        sequence_a: &str,
        sequence_b: &str,
        structure: &str,
    ) -> Result<EnergyBreakdown, RnaError> {
        let sequence = format!("{sequence_a}&{sequence_b}");
        let parsed = parse_structure(&sequence, structure)?;
        if is_pseudoknotted(&parsed.pairs) {
            return Err(RnaError::PseudoknotUnsupported("cofold energy evaluation"));
        }
        if parsed.pairs.iter().any(|pair| !pair.canonical) {
            return Err(RnaError::InvalidOption(
                "cofold energy evaluation requires canonical pairs".into(),
            ));
        }
        let cut = sequence_a.len();
        let mut result = EnergyBreakdown::new(
            self.temperature_celsius,
            self.dangles,
            self.salt.molar,
            self.model_name(),
            ConstraintModel::unconstrained(parsed.length).summary(),
        );
        let top_level = parsed
            .pairs
            .iter()
            .filter_map(|pair| {
                let i = pair.i - 1;
                let j = pair.j - 1;
                (!parsed
                    .pairs
                    .iter()
                    .any(|outer| outer.i - 1 < i && j < outer.j - 1))
                .then_some((i, j))
            })
            .collect::<Vec<_>>();
        let mut exterior = if self.dangles % 2 == 1 {
            let stems = top_level
                .iter()
                .map(|&(i, j)| LoopStem {
                    i,
                    j,
                    reversed: false,
                    five: i
                        .checked_sub(1)
                        .filter(|_| i != cut)
                        .filter(|&p| parsed.partner[p].is_none()),
                    three: (j + 1 < parsed.length && j + 1 != cut)
                        .then_some(j + 1)
                        .filter(|&p| parsed.partner[p].is_none()),
                })
                .collect::<Vec<_>>();
            self.optimize_odd_stems(parsed.sequence.as_bytes(), &stems, false, false)
        } else {
            top_level
                .iter()
                .map(|&(i, j)| {
                    self.cofold_exterior_stem_energy(parsed.sequence.as_bytes(), i, j, cut)
                })
                .sum()
        };
        if top_level.iter().any(|&(i, j)| i < cut && cut <= j) {
            exterior += self.duplex_initialization_energy();
        }
        result.exterior_kcal_mol += exterior;
        result.loop_energies.push(LoopEnergy {
            kind: "intermolecular-exterior",
            i: 0,
            j: 0,
            energy_kcal_mol: exterior,
        });
        for (i, j) in top_level {
            if i < cut && cut <= j {
                self.evaluate_cofold_pair(&parsed, i, j, cut, &mut result);
            } else {
                self.evaluate_pair(&parsed, i, j, &mut result);
            }
        }
        Ok(result.finish())
    }

    fn evaluate_cofold_pair(
        &self,
        parsed: &ParsedStructure,
        i: usize,
        j: usize,
        cut: usize,
        result: &mut EnergyBreakdown,
    ) {
        let bases = parsed.sequence.as_bytes();
        let mut children = Vec::new();
        let mut cursor = i + 1;
        while cursor < j {
            if let Some(partner) = parsed.partner[cursor] {
                if partner > cursor && partner < j {
                    children.push((cursor, partner));
                    cursor = partner + 1;
                    continue;
                }
            }
            cursor += 1;
        }
        let spanning = children
            .iter()
            .position(|&(left, right)| left < cut && cut <= right);
        match (spanning, children.as_slice()) {
            (None, _) => {
                let mut stems = Vec::with_capacity(children.len() + 1);
                stems.push(LoopStem {
                    i,
                    j,
                    reversed: true,
                    five: (j > cut && parsed.partner[j - 1].is_none()).then_some(j - 1),
                    three: (i + 1 < cut && parsed.partner[i + 1].is_none()).then_some(i + 1),
                });
                stems.extend(children.iter().map(|&(k, l)| {
                    LoopStem {
                        i: k,
                        j: l,
                        reversed: false,
                        five: k
                            .checked_sub(1)
                            .filter(|_| k != cut)
                            .filter(|&p| parsed.partner[p].is_none()),
                        three: (l + 1 < parsed.length && l + 1 != cut)
                            .then_some(l + 1)
                            .filter(|&p| parsed.partner[p].is_none()),
                    }
                }));
                let energy = if self.dangles % 2 == 1 {
                    self.optimize_odd_stems(bases, &stems, false, false)
                } else {
                    self.oriented_exterior_stem_energy(
                        bases[j],
                        bases[i],
                        (j > cut).then_some(bases[j - 1]),
                        (i + 1 < cut).then_some(bases[i + 1]),
                    ) + children
                        .iter()
                        .map(|&(k, l)| self.cofold_exterior_stem_energy(bases, k, l, cut))
                        .sum::<f64>()
                };
                result.exterior_kcal_mol += energy;
                result.loop_energies.push(LoopEnergy {
                    kind: "intermolecular-inner-exterior",
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
                for &(k, l) in &children {
                    self.evaluate_pair(parsed, k, l, result);
                }
            }
            (Some(_), &[(k, l)]) => {
                let energy = self.internal_energy(bases, i, j, k, l);
                if k == i + 1 && l + 1 == j {
                    result.stack_kcal_mol += energy;
                } else {
                    result.internal_loop_kcal_mol += energy;
                }
                result.loop_energies.push(LoopEnergy {
                    kind: if k == i + 1 && l + 1 == j {
                        "stack"
                    } else {
                        "internal"
                    },
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
                self.evaluate_cofold_pair(parsed, k, l, cut, result);
            }
            (Some(spanning_index), _) => {
                let occupied: usize = children.iter().map(|&(k, l)| l - k + 1).sum();
                let unpaired = (j - i - 1).saturating_sub(occupied);
                let mut energy =
                    self.multiloop_closing() + self.multiloop_unpaired() * unpaired as f64;
                if self.dangles % 2 == 1 {
                    let evaluated = self.evaluate_odd_multiloop(parsed, i, j, &children);
                    energy += evaluated.energy;
                    result.coaxial_stacks.extend(evaluated.coaxial_stacks);
                } else {
                    energy += self.multiloop_closing_stem_energy(bases, i, j);
                    energy += children
                        .iter()
                        .map(|&(k, l)| self.multiloop_stem_energy(bases, k, l))
                        .sum::<f64>();
                }
                result.multiloop_kcal_mol += energy;
                result.loop_energies.push(LoopEnergy {
                    kind: "multiloop",
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
                for (index, &(k, l)) in children.iter().enumerate() {
                    if index == spanning_index {
                        self.evaluate_cofold_pair(parsed, k, l, cut, result);
                    } else {
                        self.evaluate_pair(parsed, k, l, result);
                    }
                }
            }
        }
    }

    /// Evaluate a supplied structure under the same hard and soft constraints
    /// used by constrained MFE/partition prediction.
    pub fn evaluate_with_constraints(
        &self,
        sequence: &str,
        structure: &str,
        min_loop: usize,
        constraints: &ConstraintModel,
    ) -> Result<EnergyBreakdown, RnaError> {
        let parsed = parse_structure(sequence, structure)?;
        self.evaluate_parsed(&parsed, min_loop, constraints)
    }

    fn evaluate_parsed(
        &self,
        parsed: &ParsedStructure,
        min_loop: usize,
        constraints: &ConstraintModel,
    ) -> Result<EnergyBreakdown, RnaError> {
        if is_pseudoknotted(&parsed.pairs) {
            return Err(RnaError::PseudoknotUnsupported("energy evaluation"));
        }
        if let Some(pair) = parsed.pairs.iter().find(|pair| !pair.canonical) {
            return Err(RnaError::InvalidOption(format!(
                "energy evaluation does not parameterize the noncanonical pair {}-{}",
                pair.i, pair.j
            )));
        }
        constraints.validate_structure(
            parsed.sequence.as_bytes(),
            &parsed.partner,
            min_loop,
            self,
        )?;
        let mut result = EnergyBreakdown::new(
            self.temperature_celsius,
            self.dangles,
            self.salt.molar,
            self.model_name(),
            constraints.summary(),
        );
        let mut top_level = Vec::new();
        for pair in &parsed.pairs {
            let i = pair.i - 1;
            let j = pair.j - 1;
            let enclosed = parsed
                .pairs
                .iter()
                .any(|outer| outer.i - 1 < i && j < outer.j - 1);
            if !enclosed {
                top_level.push((i, j));
            }
        }
        let exterior = if self.dangles % 2 == 1 {
            self.evaluate_odd_exterior(parsed, &top_level)
        } else {
            top_level
                .iter()
                .map(|&(i, j)| self.exterior_stem_energy(parsed.sequence.as_bytes(), i, j))
                .sum()
        };
        result.exterior_kcal_mol += exterior;
        result.loop_energies.push(LoopEnergy {
            kind: "exterior",
            i: 0,
            j: 0,
            energy_kcal_mol: exterior,
        });
        for (i, j) in top_level {
            self.evaluate_pair(parsed, i, j, &mut result);
        }
        result.constraint_kcal_mol = constraints.structure_energy(&parsed.partner);
        Ok(result.finish())
    }

    fn evaluate_pair(
        &self,
        parsed: &ParsedStructure,
        i: usize,
        j: usize,
        result: &mut EnergyBreakdown,
    ) {
        let bases = parsed.sequence.as_bytes();
        let mut children = Vec::new();
        let mut cursor = i + 1;
        while cursor < j {
            if let Some(partner) = parsed.partner[cursor] {
                if partner > cursor && partner < j {
                    children.push((cursor, partner));
                    cursor = partner + 1;
                    continue;
                }
            }
            cursor += 1;
        }

        match children.as_slice() {
            [] => {
                let energy = self.hairpin_energy(bases, i, j);
                result.hairpin_kcal_mol += energy;
                result.loop_energies.push(LoopEnergy {
                    kind: "hairpin",
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
            }
            &[(k, l)] => {
                let energy = self.internal_energy(bases, i, j, k, l);
                if k == i + 1 && l + 1 == j {
                    result.stack_kcal_mol += energy;
                } else {
                    result.internal_loop_kcal_mol += energy;
                }
                result.loop_energies.push(LoopEnergy {
                    kind: if k == i + 1 && l + 1 == j {
                        "stack"
                    } else {
                        "internal"
                    },
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
                self.evaluate_pair(parsed, k, l, result);
            }
            _ => {
                let occupied: usize = children.iter().map(|&(k, l)| l - k + 1).sum();
                let unpaired = (j - i - 1).saturating_sub(occupied);
                let mut energy =
                    self.multiloop_closing() + self.multiloop_unpaired() * unpaired as f64;
                if self.dangles % 2 == 1 {
                    let evaluated = self.evaluate_odd_multiloop(parsed, i, j, &children);
                    energy += evaluated.energy;
                    result.coaxial_stacks.extend(evaluated.coaxial_stacks);
                } else {
                    energy += self.multiloop_closing_stem_energy(bases, i, j);
                    for &(k, l) in &children {
                        energy += self.multiloop_stem_energy(bases, k, l);
                    }
                }
                result.multiloop_kcal_mol += energy;
                result.loop_energies.push(LoopEnergy {
                    kind: "multiloop",
                    i: i + 1,
                    j: j + 1,
                    energy_kcal_mol: energy,
                });
                for &(k, l) in &children {
                    self.evaluate_pair(parsed, k, l, result);
                }
            }
        }
    }

    fn evaluate_odd_exterior(&self, parsed: &ParsedStructure, stems: &[(usize, usize)]) -> f64 {
        if stems.is_empty() {
            return 0.0;
        }
        let bases = parsed.sequence.as_bytes();
        let loop_stems: Vec<LoopStem> = stems
            .iter()
            .map(|&(i, j)| LoopStem {
                i,
                j,
                reversed: false,
                five: i.checked_sub(1).filter(|&p| parsed.partner[p].is_none()),
                three: (j + 1 < bases.len() && parsed.partner[j + 1].is_none()).then_some(j + 1),
            })
            .collect();
        self.optimize_odd_stems(bases, &loop_stems, false, false)
    }

    /// Stem contribution for a circular multiloop without a distinguished
    /// closing pair. Neighbor lookup wraps across the sequence cut and odd
    /// models optimize shared single-dangles/coaxial stacks on the cycle.
    pub(crate) fn circular_multiloop_stems_energy(
        &self,
        bases: &[u8],
        partner: &[Option<usize>],
        stems: &[(usize, usize)],
    ) -> f64 {
        if stems.is_empty() {
            return 0.0;
        }
        let n = bases.len();
        let loop_stems = stems
            .iter()
            .map(|&(i, j)| {
                let before = (i + n - 1) % n;
                let after = (j + 1) % n;
                LoopStem {
                    i,
                    j,
                    reversed: false,
                    five: partner[before].is_none().then_some(before),
                    three: partner[after].is_none().then_some(after),
                }
            })
            .collect::<Vec<_>>();
        match self.dangles {
            1 => self.optimize_odd_stems(bases, &loop_stems, true, true),
            3 if stems.len() >= 3 => self.evaluate_d3_cyclic_multiloop(bases, &loop_stems).energy,
            _ => loop_stems
                .iter()
                .map(|stem| {
                    self.multiloop_stem_selected(
                        bases,
                        stem.i,
                        stem.j,
                        stem.five.map(|position| bases[position]),
                        stem.three.map(|position| bases[position]),
                    )
                })
                .sum(),
        }
    }

    fn evaluate_odd_multiloop(
        &self,
        parsed: &ParsedStructure,
        i: usize,
        j: usize,
        children: &[(usize, usize)],
    ) -> OddLoopEvaluation {
        let bases = parsed.sequence.as_bytes();
        let mut stems = Vec::with_capacity(children.len() + 1);
        stems.push(LoopStem {
            i,
            j,
            reversed: true,
            five: (j > i + 1 && parsed.partner[j - 1].is_none()).then_some(j - 1),
            three: (i + 1 < j && parsed.partner[i + 1].is_none()).then_some(i + 1),
        });
        stems.extend(children.iter().map(|&(k, l)| LoopStem {
            i: k,
            j: l,
            reversed: false,
            five: (k > i + 1 && parsed.partner[k - 1].is_none()).then_some(k - 1),
            three: (l + 1 < j && parsed.partner[l + 1].is_none()).then_some(l + 1),
        }));
        if self.dangles == 3 {
            self.evaluate_d3_cyclic_multiloop(bases, &stems)
        } else {
            OddLoopEvaluation {
                energy: self.optimize_odd_stems(bases, &stems, true, true),
                coaxial_stacks: Vec::new(),
            }
        }
    }

    /// ViennaRNA's dangles=3 fixed-structure evaluator walks a multiloop
    /// twice. Each pass chooses a different cut in the cycle so every legal
    /// coaxial matching is considered while the closing edge of a pass is
    /// deliberately excluded. This reproduces `energy_of_ml_pt()` including
    /// its `ld5` bookkeeping for a shared single-nucleotide dangle.
    fn evaluate_d3_cyclic_multiloop(&self, bases: &[u8], stems: &[LoopStem]) -> OddLoopEvaluation {
        debug_assert!(stems.len() >= 3);
        let first = self.evaluate_d3_cyclic_pass(bases, stems, 0);
        let second = self.evaluate_d3_cyclic_pass(bases, stems, 1);
        if second.energy < first.energy {
            second
        } else {
            first
        }
    }

    fn evaluate_d3_cyclic_pass(
        &self,
        bases: &[u8],
        stems: &[LoopStem],
        start: usize,
    ) -> OddLoopEvaluation {
        let count = stems.len();
        let before_start = (start + count - 1) % count;
        let mut ld5 = self.odd_dangle5(bases, &stems[start]);
        if loop_gap(&stems[before_start], &stems[start]) == 1
            && self.odd_dangle3(bases, &stems[before_start]) < ld5
        {
            ld5 = 0.0;
        }

        let mut energy = D3EvaluationState::default();
        let mut coaxial: Option<D3EvaluationState> = None;
        let mut previous = start;

        // Include every stem exactly once as `current`; the final transition
        // returns to `start`, and any new coax created there is not accepted.
        for step in 1..=count {
            let current = (start + step) % count;
            // `energy_of_ml_pt()` deliberately uses P->MLintern[type]
            // directly for dangles=3. Unlike vrna_E_multibranch_stem(), this
            // path does not add a terminal AU/GU penalty.
            let current_base = self.multiloop_intern_base();
            energy.energy += current_base;
            if let Some(state) = &mut coaxial {
                state.energy += current_base;
            }

            let gap = loop_gap(&stems[previous], &stems[current]);
            let dangle5 = self.odd_dangle5(bases, &stems[current]);
            let dangle3 = self.odd_dangle3(bases, &stems[previous]);

            let next_coaxial = match gap {
                0 => {
                    let replacement = self.coaxial_energy(
                        bases,
                        stems[previous].i,
                        stems[previous].j,
                        stems[previous].reversed,
                        stems[current].i,
                        stems[current].j,
                        stems[current].reversed,
                    );
                    let stabilization = replacement - self.multiloop_intern_base() - current_base;
                    let mut state = energy.clone();
                    state.energy += stabilization - ld5;
                    state.coaxial_stacks.push(CoaxialStack {
                        loop_i: stems[0].i + 1,
                        loop_j: stems[0].j + 1,
                        first_i: stems[previous].i + 1,
                        first_j: stems[previous].j + 1,
                        second_i: stems[current].i + 1,
                        second_j: stems[current].j + 1,
                        energy_kcal_mol: replacement,
                        stabilization_kcal_mol: stabilization,
                    });
                    ld5 = 0.0;
                    if let Some(previous_coaxial) = coaxial.take() {
                        if previous_coaxial.energy < energy.energy {
                            energy = previous_coaxial;
                        }
                    }
                    Some(state)
                }
                1 => {
                    let dangling = dangle3.min(dangle5);
                    energy.energy += dangling;
                    ld5 = dangling - dangle3;
                    if let Some(mut previous_coaxial) = coaxial.take() {
                        previous_coaxial.energy += dangle5;
                        if previous_coaxial.energy < energy.energy {
                            ld5 = dangle5;
                            energy = previous_coaxial;
                        }
                    }
                    None
                }
                _ => {
                    energy.energy += dangle5 + dangle3;
                    if let Some(mut previous_coaxial) = coaxial.take() {
                        previous_coaxial.energy += dangle5;
                        if previous_coaxial.energy < energy.energy {
                            energy = previous_coaxial;
                        }
                    }
                    ld5 = dangle5;
                    None
                }
            };

            coaxial = next_coaxial;
            previous = current;
        }

        OddLoopEvaluation {
            energy: energy.energy,
            coaxial_stacks: energy.coaxial_stacks,
        }
    }

    fn odd_dangle5(&self, bases: &[u8], stem: &LoopStem) -> f64 {
        let Some(index) = stem.five else {
            return 0.0;
        };
        let pair = self.odd_stem_pair_type(bases, stem);
        self.dangle_context(pair, Some(bases[index]), None, true)
    }

    fn multiloop_intern_base(&self) -> f64 {
        self.scaled(
            params!(self, ML_PARAMS, ml_params)[4],
            params!(self, ML_PARAMS, ml_params)[5],
        ) + self.salt.ml_base_centi as f64 / 100.0
    }

    fn odd_dangle3(&self, bases: &[u8], stem: &LoopStem) -> f64 {
        let Some(index) = stem.three else {
            return 0.0;
        };
        let pair = self.odd_stem_pair_type(bases, stem);
        self.dangle_context(pair, None, Some(bases[index]), true)
    }

    fn odd_stem_pair_type(&self, bases: &[u8], stem: &LoopStem) -> usize {
        if stem.reversed {
            pair_type(bases[stem.j], bases[stem.i])
        } else {
            pair_type(bases[stem.i], bases[stem.j])
        }
        .expect("canonical stems were validated before evaluation")
    }

    fn optimize_odd_stems(
        &self,
        bases: &[u8],
        stems: &[LoopStem],
        cyclic: bool,
        multiloop: bool,
    ) -> f64 {
        if stems.is_empty() {
            return 0.0;
        }
        let states: Vec<Vec<OddStemState>> = stems
            .iter()
            .map(|stem| {
                let mut options = Vec::new();
                for mask in 0..4u8 {
                    let five = mask & 1 != 0;
                    let three = mask & 2 != 0;
                    if (five && stem.five.is_none()) || (three && stem.three.is_none()) {
                        continue;
                    }
                    let energy = self.odd_stem_energy(bases, stem, five, three, multiloop);
                    options.push(OddStemState {
                        five,
                        three,
                        energy,
                    });
                }
                options
            })
            .collect();

        let mut best = f64::INFINITY;
        for &first in &states[0] {
            let mut previous = vec![(first, first.energy)];
            for index in 1..stems.len() {
                let mut next_values = Vec::new();
                for &current in &states[index] {
                    let mut value = f64::INFINITY;
                    for &(prior, prior_energy) in &previous {
                        if odd_states_compatible(prior, current, &stems[index - 1], &stems[index]) {
                            value = value.min(prior_energy + current.energy);
                        }
                    }
                    next_values.push((current, value));
                }
                previous = next_values;
            }
            for &(last, value) in &previous {
                if !cyclic || odd_states_compatible(last, first, stems.last().unwrap(), &stems[0]) {
                    best = best.min(value);
                }
            }
        }
        best
    }

    fn odd_stem_base(&self, bases: &[u8], stem: &LoopStem, multiloop: bool) -> f64 {
        let pair = if stem.reversed {
            pair_type(bases[stem.j], bases[stem.i])
        } else {
            pair_type(bases[stem.i], bases[stem.j])
        }
        .expect("canonical stems were validated before evaluation");
        let intern = if multiloop {
            self.scaled(
                params!(self, ML_PARAMS, ml_params)[4],
                params!(self, ML_PARAMS, ml_params)[5],
            ) + self.salt.ml_base_centi as f64 / 100.0
        } else {
            0.0
        };
        intern + self.terminal_au_for_type(pair)
    }

    fn odd_stem_energy(
        &self,
        bases: &[u8],
        stem: &LoopStem,
        use_five: bool,
        use_three: bool,
        multiloop: bool,
    ) -> f64 {
        let pair = if stem.reversed {
            pair_type(bases[stem.j], bases[stem.i])
        } else {
            pair_type(bases[stem.i], bases[stem.j])
        }
        .expect("canonical stems were validated before evaluation");
        let five = use_five.then(|| bases[stem.five.unwrap()]);
        let three = use_three.then(|| bases[stem.three.unwrap()]);
        let context = if self.dangles == 3 && multiloop {
            self.dangle_context(pair, five, None, multiloop)
                + self.dangle_context(pair, None, three, multiloop)
        } else {
            self.dangle_context(pair, five, three, multiloop)
        };
        self.odd_stem_base(bases, stem, multiloop) + context
    }
}

#[derive(Clone, Copy)]
struct LoopStem {
    i: usize,
    j: usize,
    reversed: bool,
    five: Option<usize>,
    three: Option<usize>,
}

#[derive(Clone, Copy)]
struct OddStemState {
    five: bool,
    three: bool,
    energy: f64,
}

struct OddLoopEvaluation {
    energy: f64,
    coaxial_stacks: Vec<CoaxialStack>,
}

#[derive(Clone, Default)]
struct D3EvaluationState {
    energy: f64,
    coaxial_stacks: Vec<CoaxialStack>,
}

fn odd_states_compatible(
    previous: OddStemState,
    current: OddStemState,
    previous_stem: &LoopStem,
    current_stem: &LoopStem,
) -> bool {
    !(previous.three
        && current.five
        && previous_stem.three.is_some()
        && previous_stem.three == current_stem.five)
}

fn loop_gap(previous: &LoopStem, current: &LoopStem) -> usize {
    if previous.reversed {
        current.i.saturating_sub(previous.i + 1)
    } else if current.reversed {
        current.j.saturating_sub(previous.j + 1)
    } else {
        current.i.saturating_sub(previous.j + 1)
    }
}

fn energy_model_name(nucleic_acid: NucleicAcid, dangles: u8) -> &'static str {
    match (nucleic_acid, dangles) {
        (NucleicAcid::Rna, 0) => "Turner 2004 RNA nearest-neighbor, dangles=0",
        (NucleicAcid::Rna, 1) => "Turner 2004 RNA nearest-neighbor, dangles=1",
        (NucleicAcid::Rna, 2) => "Turner 2004 RNA nearest-neighbor, dangles=2",
        (NucleicAcid::Rna, 3) => {
            "Turner 2004 RNA nearest-neighbor, dangles=3 with coaxial stacking"
        }
        (NucleicAcid::Dna, 0) => "Mathews 2004 DNA nearest-neighbor, dangles=0",
        (NucleicAcid::Dna, 1) => "Mathews 2004 DNA nearest-neighbor, dangles=1",
        (NucleicAcid::Dna, 2) => "Mathews 2004 DNA nearest-neighbor, dangles=2",
        (NucleicAcid::Dna, 3) => {
            "Mathews 2004 DNA nearest-neighbor, dangles=3 with coaxial stacking"
        }
        _ => unreachable!("validated dangle model"),
    }
}

fn ensemble_model_name(nucleic_acid: NucleicAcid, dangles: u8) -> &'static str {
    match (nucleic_acid, dangles) {
        (NucleicAcid::Rna, 0) => "Turner 2004 RNA McCaskill, dangles=0",
        (NucleicAcid::Rna, 2) => "Turner 2004 RNA McCaskill, dangles=2",
        (NucleicAcid::Rna, 1) => {
            "Turner 2004 RNA exact fixed-structure ensemble, exclusive single dangles"
        }
        (NucleicAcid::Rna, 3) => {
            "Turner 2004 RNA exact fixed-structure ensemble, single dangles and coaxial stacking"
        }
        (NucleicAcid::Dna, 0) => "Mathews 2004 DNA McCaskill, dangles=0",
        (NucleicAcid::Dna, 2) => "Mathews 2004 DNA McCaskill, dangles=2",
        (NucleicAcid::Dna, 1) => {
            "Mathews 2004 DNA exact fixed-structure ensemble, exclusive single dangles"
        }
        (NucleicAcid::Dna, 3) => {
            "Mathews 2004 DNA exact fixed-structure ensemble, single dangles and coaxial stacking"
        }
        _ => unreachable!("validated dangle model"),
    }
}

fn smooth_favorable_centi(value: f64) -> f64 {
    const SCALE: f64 = 10.0;
    let scaled = value / SCALE;
    if scaled < -1.228_369_7 {
        0.0
    } else if scaled > 0.866_025_4 {
        value
    } else {
        let wave = (scaled - 0.342_426_63).sin() + 1.0;
        SCALE * 0.384_900_18 * wave * wave
    }
}

fn find_special(table: &[(&str, i32, i32)], sequence: &str) -> Option<(i32, i32)> {
    table
        .iter()
        .find(|(motif, _, _)| *motif == sequence)
        .map(|(_, energy, enthalpy)| (*energy, *enthalpy))
}

fn pair_type(a: u8, b: u8) -> Option<usize> {
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

fn reverse_pair_type(pair: usize) -> usize {
    match pair {
        0 => 1,
        1 => 0,
        2 => 3,
        3 => 2,
        4 => 5,
        5 => 4,
        _ => pair,
    }
}

fn base_index(base: u8) -> usize {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 3,
        b'U' => 4,
        _ => 0,
    }
}

fn index_mismatch(pair: usize, left: usize, right: usize) -> usize {
    (pair * 5 + left) * 5 + right
}

fn index_int11(outer: usize, inner: usize, left: usize, right: usize) -> usize {
    ((outer * 7 + inner) * 5 + left) * 5 + right
}

fn index_int21(outer: usize, inner: usize, first: usize, second: usize, third: usize) -> usize {
    (((outer * 7 + inner) * 5 + first) * 5 + second) * 5 + third
}

fn index_int22(
    outer: usize,
    inner: usize,
    first: usize,
    second: usize,
    third: usize,
    fourth: usize,
) -> usize {
    let [first, second, third, fourth] = [first, second, third, fourth].map(|base| base - 1);
    (((((outer * 6 + inner) * 4 + first) * 4 + second) * 4 + third) * 4) + fourth
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_gc_helix_is_stabilized_by_stacking() {
        let model = EnergyModel::default();
        let single = model.evaluate("GAAAC", "(...)").unwrap();
        let helix = model.evaluate("GGAAACC", "((...))").unwrap();
        assert!(helix.total_kcal_mol < single.total_kcal_mol);
        assert!(helix.stack_kcal_mol < 0.0);
    }

    #[test]
    fn turner_special_tetraloop_uses_absolute_energy() {
        let model = EnergyModel::default();
        // CUUCGG is a Turner 2004 special tetraloop with 3.7 kcal/mol.
        assert!((model.hairpin_energy(b"CUUCGG", 0, 5) - 3.7).abs() < 1.0e-12);
    }

    #[test]
    fn temperature_rescaling_uses_enthalpy_tables() {
        let cold = EnergyModel::new(10.0).unwrap();
        let warm = EnergyModel::new(70.0).unwrap();
        assert_ne!(
            cold.stack_energy(b'G', b'C', b'G', b'C'),
            warm.stack_energy(b'G', b'C', b'G', b'C')
        );
    }

    #[test]
    fn mathews_2004_dna_bundle_is_independent_and_matches_reference_energy() {
        let rna =
            EnergyModel::with_parameter_family(37.0, 0, STANDARD_MOLAR, NucleicAcid::Rna).unwrap();
        let dna =
            EnergyModel::with_parameter_family(37.0, 0, STANDARD_MOLAR, NucleicAcid::Dna).unwrap();
        let dna_energy = dna.evaluate("GGGTTTCCC", "(((...)))").unwrap();
        // CParty 1.0 with dna_Matthews04.par and all positions constrained
        // reports -0.2 kcal/mol for this structure.
        assert!((dna_energy.total_kcal_mol + 0.2).abs() < 1.0e-12);
        assert_ne!(
            dna.stack_energy(b'G', b'C', b'G', b'C'),
            rna.stack_energy(b'G', b'C', b'G', b'C')
        );
        assert_eq!(dna.parameter_model_id(), "ribon-mathews-dna-2004");
    }

    #[test]
    fn normalized_custom_overlay_changes_only_the_requested_table() {
        let baseline = EnergyModel::with_dangles(37.0, 0).unwrap();
        let profile = ThermodynamicParameterOverrides {
            schema_version: 1,
            name: "terminal-au-test".into(),
            fingerprint_sha256: "ab".repeat(32),
            misc: Some(vec![410, 360, 500, 370]),
            ..ThermodynamicParameterOverrides::default()
        };
        let custom = EnergyModel::with_parameter_overrides(
            37.0,
            0,
            STANDARD_MOLAR,
            NucleicAcid::Rna,
            profile,
        )
        .unwrap();
        let ordinary = baseline.evaluate("AAAAU", "(...)").unwrap();
        let changed = custom.evaluate("AAAAU", "(...)").unwrap();
        assert!(
            (changed.total_kcal_mol - ordinary.total_kcal_mol - 9.0).abs() < 1.0e-12,
            "ordinary={} changed={} delta={}",
            ordinary.total_kcal_mol,
            changed.total_kcal_mol,
            changed.total_kcal_mol - ordinary.total_kcal_mol,
        );
        assert_eq!(custom.parameter_model_id(), CUSTOM_MODEL_ID);
        assert_eq!(custom.parameter_profile_name(), Some("terminal-au-test"));
    }

    #[test]
    fn malformed_custom_overlay_is_rejected_before_analysis() {
        let invalid = ThermodynamicParameterOverrides {
            schema_version: 1,
            name: "bad".into(),
            fingerprint_sha256: "00".repeat(32),
            stack_37: Some(vec![0; 48]),
            ..ThermodynamicParameterOverrides::default()
        };
        let error = EnergyModel::with_parameter_overrides(
            37.0,
            2,
            STANDARD_MOLAR,
            NucleicAcid::Rna,
            invalid,
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("stack_37 has length 48; expected 49"));
    }

    #[test]
    fn custom_overlay_validation_handles_minimum_i32_without_panicking() {
        let invalid = ThermodynamicParameterOverrides {
            schema_version: 1,
            name: "out-of-range".into(),
            fingerprint_sha256: "00".repeat(32),
            misc: Some(vec![i32::MIN, 0, 0, 0]),
            ..ThermodynamicParameterOverrides::default()
        };
        let error = EnergyModel::with_parameter_overrides(
            37.0,
            2,
            STANDARD_MOLAR,
            NucleicAcid::Rna,
            invalid,
        )
        .unwrap_err();
        assert!(error.to_string().contains("outside the normalized range"));
    }

    #[test]
    fn custom_dna_special_loops_normalize_t_to_u_for_lookup() {
        let profile = ThermodynamicParameterOverrides {
            schema_version: 1,
            name: "dna-tetraloop".into(),
            fingerprint_sha256: "cd".repeat(32),
            tetraloops: Some(vec![SpecialLoopParameter {
                sequence: "CTTCGG".into(),
                free_energy_37_centi_kcal_mol: 123,
                enthalpy_centi_kcal_mol: 123,
            }]),
            ..ThermodynamicParameterOverrides::default()
        };
        let model = EnergyModel::with_parameter_overrides(
            37.0,
            0,
            STANDARD_MOLAR,
            NucleicAcid::Dna,
            profile,
        )
        .unwrap();
        assert!((model.hairpin_energy(b"CUUCGG", 0, 5) - 1.23).abs() < 1.0e-12);
    }

    #[test]
    fn noncanonical_pairs_fail_before_serialization() {
        let model = EnergyModel::default();
        let error = model.evaluate("AAAAA", "(...)").unwrap_err();
        assert!(error.to_string().contains("noncanonical pair 1-5"));
    }

    #[test]
    fn dangle_model_is_explicit_and_validated() {
        assert_eq!(EnergyModel::default().dangles(), 2);
        for supported in 0..=3 {
            assert_eq!(
                EnergyModel::with_dangles(37.0, supported)
                    .unwrap()
                    .dangles(),
                supported
            );
        }
        assert!(EnergyModel::with_dangles(37.0, 4).is_err());
        assert!(EnergyModel::with_dangles(37.0, 1)
            .unwrap()
            .supports_partition());
    }

    #[test]
    fn exterior_context_changes_supplied_structure_energy() {
        let none = EnergyModel::with_dangles(37.0, 0)
            .unwrap()
            .evaluate("AGGGAAACCCU", ".(((...))).")
            .unwrap();
        let double = EnergyModel::with_dangles(37.0, 2)
            .unwrap()
            .evaluate("AGGGAAACCCU", ".(((...))).")
            .unwrap();
        assert!(double.exterior_kcal_mol < none.exterior_kcal_mol);
        assert_eq!(double.dangle_model, 2);
    }

    #[test]
    fn odd_dangle_evaluator_is_deterministic_under_ribon_parameters() {
        for (sequence, structure, d1, d3) in [
            ("GGAAACAGAAACC", "((...).(...))", 16.6, 13.3),
            ("GAGAAACGAAACC", "(.(...)(...))", 17.2, 14.8),
            ("GAGAAACAGAAACAC", "(.(...).(...).)", 13.8, 13.8),
        ] {
            let single = EnergyModel::with_dangles(37.0, 1)
                .unwrap()
                .evaluate(sequence, structure)
                .unwrap();
            let coaxial = EnergyModel::with_dangles(37.0, 3)
                .unwrap()
                .evaluate(sequence, structure)
                .unwrap();
            assert!((single.total_kcal_mol - d1).abs() < 1.0e-9);
            assert!((coaxial.total_kcal_mol - d3).abs() < 1.0e-9);
        }

        let sequence = "GUCCCGGCCUCGAGACCUAUCCGGUUCGUCACGGAGCGCAGCCCGUGACGCGGGGUGACU";
        let structure = "(((...((((((.((((.....))))((((((((.((...))))))))))))))))))).";
        let evaluated = EnergyModel::with_dangles(37.0, 3)
            .unwrap()
            .evaluate(sequence, structure)
            .unwrap();
        assert!((evaluated.total_kcal_mol + 27.1).abs() < 1.0e-9);
        assert!((evaluated.exterior_kcal_mol + 1.2).abs() < 1.0e-9);
        let multiloop = evaluated
            .loop_energies
            .iter()
            .find(|entry| entry.kind == "multiloop" && entry.i == 12)
            .unwrap();
        assert!((multiloop.energy_kcal_mol - 3.8).abs() < 1.0e-9);
        assert_eq!(evaluated.coaxial_stacks.len(), 1);
    }
}
