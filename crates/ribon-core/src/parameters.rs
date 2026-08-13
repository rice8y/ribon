//! Reproducible thermodynamic-model metadata and analysis profiles.

use crate::constraints::ConstraintConfig;
use crate::structure::RnaError;
use crate::{analyze_with_model, AnalysisResult, EnergyModel, NucleicAcid};
use serde::{Deserialize, Serialize};

pub const MODEL_ID: &str = "ribon-turner-2004";
pub const DNA_MODEL_ID: &str = "ribon-mathews-dna-2004";
pub const PARAMETER_BUNDLE_SHA256: &str =
    "0c00a31400f1dedbe9a3e161b2f9b1b74cde54941144ee988f48173d33bbcd7b";
pub const SOURCE_ARCHIVE_SHA256: &str =
    "8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8";
pub const DNA_PARAMETER_BUNDLE_SHA256: &str =
    "019ad1d5c3dac421df37e0a5aeded6d3da50da03deecc23ba0ae5a6d5d06b977";

/// Serializable configuration for a complete, reproducible analysis.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct ParameterProfile {
    pub schema_version: u32,
    pub name: String,
    pub model_id: String,
    pub parameter_bundle_sha256: String,
    pub temperature_celsius: f64,
    pub min_loop: usize,
    pub mea_gamma: f64,
    pub dangles: u8,
    pub salt_molar: f64,
    pub constraints: ConstraintConfig,
}

impl Default for ParameterProfile {
    fn default() -> Self {
        Self {
            schema_version: 1,
            name: "ribon-standard".into(),
            model_id: MODEL_ID.into(),
            parameter_bundle_sha256: PARAMETER_BUNDLE_SHA256.into(),
            temperature_celsius: 37.0,
            min_loop: 3,
            mea_gamma: 1.0,
            dangles: 2,
            salt_molar: 1.021,
            constraints: ConstraintConfig::default(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ParameterManifest {
    pub schema_version: u32,
    pub model_id: &'static str,
    pub model_name: &'static str,
    pub parameterization: &'static str,
    pub data_source: &'static str,
    pub source_version: &'static str,
    pub source_archive_sha256: &'static str,
    pub parameter_bundle_sha256: &'static str,
    pub parameter_file_count: usize,
    pub license: &'static str,
    pub units: &'static str,
    pub supported_dangles_mfe: [u8; 4],
    pub supported_dangles_partition: [u8; 4],
}

pub fn parameter_manifest() -> ParameterManifest {
    ParameterManifest {
        schema_version: 1,
        model_id: MODEL_ID,
        model_name: "Ribon Turner 2004 nearest-neighbor model",
        parameterization: "Turner 2004 RNA nearest-neighbor free energies and enthalpies",
        data_source: "RNAstructure data_tables/rna.*",
        source_version: "RNAstructure 6.6",
        source_archive_sha256: SOURCE_ARCHIVE_SHA256,
        parameter_bundle_sha256: PARAMETER_BUNDLE_SHA256,
        parameter_file_count: 34,
        license: "GPL-2.0-only",
        units: "source kcal/mol; generated centi-kcal/mol tables; API kcal/mol",
        supported_dangles_mfe: [0, 1, 2, 3],
        supported_dangles_partition: [0, 1, 2, 3],
    }
}

pub fn dna_parameter_manifest() -> ParameterManifest {
    ParameterManifest {
        schema_version: 1,
        model_id: DNA_MODEL_ID,
        model_name: "Ribon Mathews 2004 DNA nearest-neighbor model",
        parameterization: "Mathews 2004 DNA nearest-neighbor free energies and enthalpies",
        data_source: "RNAstructure data_tables/dna.*",
        source_version: "RNAstructure 6.6",
        source_archive_sha256: SOURCE_ARCHIVE_SHA256,
        parameter_bundle_sha256: DNA_PARAMETER_BUNDLE_SHA256,
        parameter_file_count: 33,
        license: "GPL-2.0-only",
        units: "source kcal/mol; generated centi-kcal/mol tables; API kcal/mol",
        supported_dangles_mfe: [0, 1, 2, 3],
        supported_dangles_partition: [0, 1, 2, 3],
    }
}

pub fn parameter_catalog() -> Vec<ParameterManifest> {
    vec![parameter_manifest(), dna_parameter_manifest()]
}

pub fn validate_parameter_profile(profile: &ParameterProfile) -> Result<(), RnaError> {
    if profile.schema_version != 1 {
        return Err(RnaError::InvalidOption(format!(
            "unsupported parameter-profile schema version {}; expected 1",
            profile.schema_version
        )));
    }
    if profile.name.trim().is_empty() {
        return Err(RnaError::InvalidOption(
            "parameter profile name must not be empty".into(),
        ));
    }
    if profile.model_id != MODEL_ID && profile.model_id != DNA_MODEL_ID {
        return Err(RnaError::InvalidOption(format!(
            "unknown model {:?}; expected {MODEL_ID:?} or {DNA_MODEL_ID:?}",
            profile.model_id
        )));
    }
    let expected_fingerprint = match profile.model_id.as_str() {
        DNA_MODEL_ID => DNA_PARAMETER_BUNDLE_SHA256,
        _ => PARAMETER_BUNDLE_SHA256,
    };
    if profile.parameter_bundle_sha256 != expected_fingerprint {
        return Err(RnaError::InvalidOption(
            "thermodynamic parameter bundle fingerprint mismatch".into(),
        ));
    }
    if !profile.temperature_celsius.is_finite() || profile.temperature_celsius <= -273.15 {
        return Err(RnaError::InvalidOption(
            "profile temperature must be finite and above absolute zero".into(),
        ));
    }
    if !profile.mea_gamma.is_finite() || profile.mea_gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "profile MEA gamma must be finite and positive".into(),
        ));
    }
    if profile.dangles > 3 {
        return Err(RnaError::InvalidOption(
            "profile dangles must be one of 0, 1, 2, or 3".into(),
        ));
    }
    if !profile.salt_molar.is_finite() || profile.salt_molar <= 0.0 {
        return Err(RnaError::InvalidOption(
            "profile salt molarity must be finite and positive".into(),
        ));
    }
    Ok(())
}

pub fn analyze_with_profile(
    sequence: &str,
    profile: &ParameterProfile,
) -> Result<AnalysisResult, RnaError> {
    validate_parameter_profile(profile)?;
    let family = if profile.model_id == DNA_MODEL_ID {
        NucleicAcid::Dna
    } else {
        NucleicAcid::Rna
    };
    let model = EnergyModel::with_parameter_family(
        profile.temperature_celsius,
        profile.dangles,
        profile.salt_molar,
        family,
    )?;
    analyze_with_model(
        sequence.into(),
        profile.min_loop,
        profile.mea_gamma,
        model,
        &profile.constraints,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_profile_is_the_default_integrated_analysis() {
        let profile = ParameterProfile::default();
        let profiled = analyze_with_profile("GGGAAACCC", &profile).unwrap();
        let direct = crate::analyze("GGGAAACCC", 37.0, 3, 1.0).unwrap();
        assert_eq!(profiled.mfe_structure, direct.mfe_structure);
        assert_eq!(
            profile.parameter_bundle_sha256,
            parameter_manifest().parameter_bundle_sha256
        );
    }

    #[test]
    fn fingerprint_mismatch_is_rejected() {
        let profile = ParameterProfile {
            parameter_bundle_sha256: "00".repeat(32),
            ..ParameterProfile::default()
        };
        assert!(validate_parameter_profile(&profile).is_err());
    }

    #[test]
    fn dna_profile_selects_the_dna_bundle() {
        let profile = ParameterProfile {
            name: "dna".into(),
            model_id: DNA_MODEL_ID.into(),
            parameter_bundle_sha256: DNA_PARAMETER_BUNDLE_SHA256.into(),
            dangles: 0,
            ..ParameterProfile::default()
        };
        let result = analyze_with_profile("GGGTTTCCC", &profile).unwrap();
        assert_eq!(result.mfe_structure, "(((...)))");
        assert!((result.mfe_energy_kcal_mol + 0.2).abs() < 1.0e-12);
        assert!(result.model.parameter_set.contains("DNA"));
    }
}
