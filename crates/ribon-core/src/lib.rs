//! Deterministic RNA secondary-structure analysis and layout core.
//!
//! The crate is deliberately independent from Typst. The `ribon-plugin` crate
//! supplies the byte-oriented WebAssembly boundary, while Typst performs all
//! final vector rendering.

pub mod accessibility;
mod circular_standard;
pub mod cofold;
mod cofold_standard;
pub mod comparative;
pub mod conditional_density2;
pub mod constraints;
pub mod decode;
pub mod design;
pub mod duplex;
pub mod energy;
pub mod ensemble;
mod exact_enumeration;
pub mod extended;
pub mod fold;
mod gquad_integrated;
pub mod landscape;
pub mod layout;
pub mod ligand;
pub mod local;
mod modified_parameters;
mod naview;
pub mod parameters;
pub mod partition;
pub mod pseudoknot;
mod puzzler;
mod salt;
pub mod sampling;
pub mod structure;
pub mod suboptimal;
pub mod topology;
mod turtle;
use serde::Serialize;

pub use accessibility::{
    accessibility, accessibility_with_model, AccessibilityEntry, AccessibilityResult,
    AccessibilityWindow,
};
pub use circular_standard::{
    circular_fold, circular_fold_with_constraints, circular_fold_with_model, CircularResult,
};
pub use cofold::{
    cofold, cofold_with_model, CofoldComplex, CofoldMonomer, CofoldResult,
    EquilibriumConcentrations,
};
pub use comparative::{
    comparative_fold, comparative_fold_with_model, ComparativeOptions, ComparativeResult,
    CovariationEntry,
};
pub use conditional_density2::{
    conditional_density2, conditional_density2_ensemble,
    conditional_density2_ensemble_with_constraints, conditional_density2_polynomial,
    conditional_density2_polynomial_with_constraints, evaluate_conditional_density2,
    evaluate_conditional_density2_polynomial, evaluate_conditional_density2_structure,
    evaluate_conditional_density2_structure_with_constraints, sample_conditional_density2,
    sample_conditional_density2_polynomial,
    sample_conditional_density2_polynomial_with_constraints,
    sample_conditional_density2_with_constraints, suboptimal_conditional_density2,
    suboptimal_conditional_density2_polynomial,
    suboptimal_conditional_density2_polynomial_with_constraints,
    suboptimal_conditional_density2_with_constraints, ConditionalDensity2Analysis,
    ConditionalDensity2Energy, ConditionalDensity2Evaluation, ConditionalDensity2EvaluationResult,
    ConditionalDensity2Options, ConditionalDensity2PolynomialEvaluationResult,
    ConditionalDensity2PolynomialResult, ConditionalDensity2Result, ConditionalDensity2Sample,
    ConditionalDensity2SamplingResult, ConditionalDensity2SuboptimalResult,
    ConditionalDensity2SuboptimalStructure,
};
pub use constraints::{
    ConstraintConfig, ConstraintModel, ConstraintSummary, PairConstraint, PairEnergy,
    PositionEnergy, ProbingConfig, SoftConstraintConfig,
};
pub use decode::{
    decode_centroid, decode_centroid_with_constraints, decode_mea, decode_mea_with_constraints,
};
pub use design::{inverse_fold_exact, DesignCandidate, InverseDesignOptions, InverseDesignResult};
pub use duplex::{duplex, duplex_with_model, DuplexResult};
pub use energy::{
    validate_parameter_overrides, CoaxialStack, EnergyBreakdown, EnergyModel, LoopEnergy,
    NucleicAcid, SpecialLoopParameter, ThermodynamicParameterOverrides, CUSTOM_MODEL_ID,
};
pub use ensemble::{
    ensemble_defect, summarize as summarize_ensemble, EnsembleDefectResult, EnsembleSummary,
};
pub use extended::{
    fold_gquad, fold_gquad_with_model, fold_modified, fold_modified_with_model, GQuadProbability,
    GQuadResult, GQuadruplex, ModifiedBase, ModifiedFoldResult, ModifiedParameterUse,
};
pub use fold::{fold_mfe, fold_mfe_with_constraints, MfeResult};
pub use landscape::{exact_landscape_path, LandscapePathStep, LandscapeResult, LandscapeState};
pub use layout::{layout_structure, LayoutKind, LayoutResult};
pub use ligand::{ligand_ensemble_exact, LigandEnsembleResult, LigandMotif, LigandOccupancy};
pub use local::{
    local_fold, local_fold_with_model, LocalAccessibility, LocalFoldResult, LocalPairProbability,
    LocalWindowSummary,
};
pub use modified_parameters::ModifiedBaseKind;
pub use parameters::{
    analyze_with_profile, dna_parameter_manifest, parameter_catalog, parameter_manifest,
    validate_parameter_profile, ParameterManifest, ParameterProfile, DNA_MODEL_ID,
    DNA_PARAMETER_BUNDLE_SHA256, MODEL_ID, PARAMETER_BUNDLE_SHA256, REFERENCE_ARCHIVE_SHA256,
};
pub use partition::{partition, partition_with_constraints, PairProbability, PartitionResult};
pub use pseudoknot::{
    decode_pseudoknot, decode_pseudoknot_with_model, evaluate_pseudoknot_structure,
    evaluate_pseudoknot_structure_with_model, predict_pseudoknot, predict_pseudoknot_with_model,
    ExactArbitraryEnsemble, PseudoknotEnergyBreakdown, PseudoknotEvaluationResult,
    PseudoknotOptions, PseudoknotPair, PseudoknotResult,
};
pub use sampling::{sample_boltzmann, BoltzmannSample, SamplingResult};
pub use structure::{normalize_sequence, parse_structure, Pair, ParsedStructure, RnaError};
pub use suboptimal::{suboptimal_structures, SuboptimalResult, SuboptimalStructure};
pub use topology::{fatgraph_topology, FatgraphPair, FatgraphTopology};

/// Complete single-sequence prediction output returned to the Typst wrapper.
#[derive(Clone, Debug, Serialize)]
pub struct AnalysisResult {
    pub sequence: String,
    pub length: usize,
    pub temperature_celsius: f64,
    pub model: ModelDescription,
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub partition_function: f64,
    pub log_partition_function: f64,
    pub centroid_structure: String,
    pub centroid_score: f64,
    pub mea_structure: String,
    pub mea_score: f64,
    pub mea_gamma: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
    pub ensemble: EnsembleSummary,
    pub constraints: ConstraintSummary,
}

#[derive(Clone, Debug, Serialize)]
pub struct ModelDescription {
    pub parameter_set: &'static str,
    pub parameter_profile_name: Option<String>,
    pub parameter_fingerprint_sha256: Option<String>,
    pub mfe: &'static str,
    pub ensemble: &'static str,
    pub dangles: u8,
    pub salt_correction: bool,
    pub salt_molar: f64,
    pub note: &'static str,
}

/// Run all prediction stages once so the result can directly feed plotting,
/// dot plots, and annotation layers.
pub fn analyze(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
) -> Result<AnalysisResult, RnaError> {
    analyze_with_dangles(sequence, temperature_celsius, min_loop, gamma, 2)
}

/// Run all prediction stages with an explicit nearest-neighbor dangle convention.
///
/// PF-derived quantities for requested models `1` and `3` enumerate the exact
/// fixed-structure single-dangle/coaxial energy over the complete planar state
/// space. Models `0` and `2` use the polynomial inside/outside grammar.
pub fn analyze_with_dangles(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
) -> Result<AnalysisResult, RnaError> {
    analyze_with_constraints(
        sequence,
        temperature_celsius,
        min_loop,
        gamma,
        dangles,
        &ConstraintConfig::default(),
    )
}

/// Run MFE, partition, centroid, and MEA under a shared constraint model.
pub fn analyze_with_constraints(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    constraint_config: &ConstraintConfig,
) -> Result<AnalysisResult, RnaError> {
    analyze_with_options(
        sequence,
        temperature_celsius,
        min_loop,
        gamma,
        dangles,
        1.021,
        constraint_config,
    )
}

/// Run integrated prediction with explicit dangles, salt molarity, and constraints.
#[allow(clippy::too_many_arguments)]
pub fn analyze_with_options(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
    constraint_config: &ConstraintConfig,
) -> Result<AnalysisResult, RnaError> {
    if sequence.chars().any(|symbol| symbol == '&') {
        return Err(RnaError::InvalidOption(
            "prediction accepts one strand; use supplied structures for multi-strand drawing"
                .into(),
        ));
    }
    if !temperature_celsius.is_finite() || temperature_celsius <= -273.15 {
        return Err(RnaError::InvalidOption(
            "temperature must be finite and above absolute zero".into(),
        ));
    }
    if !gamma.is_finite() || gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "MEA gamma must be finite and positive".into(),
        ));
    }

    let sequence = normalize_sequence(sequence)?;
    let energy_model =
        EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    analyze_with_model(sequence, min_loop, gamma, energy_model, constraint_config)
}

/// Run the integrated analysis with an already selected thermodynamic
/// parameter family. This is the common entry point for built-in RNA, DNA,
/// and future validated parameter profiles.
pub fn analyze_with_model(
    sequence: String,
    min_loop: usize,
    gamma: f64,
    energy_model: EnergyModel,
    constraint_config: &ConstraintConfig,
) -> Result<AnalysisResult, RnaError> {
    if !gamma.is_finite() || gamma <= 0.0 {
        return Err(RnaError::InvalidOption(
            "MEA gamma must be finite and positive".into(),
        ));
    }
    let sequence = normalize_sequence(&sequence)?;
    let temperature_celsius = energy_model.temperature_celsius();
    let constraints = ConstraintModel::compile(sequence.len(), constraint_config)?;
    analyze_with_compiled_constraints(
        sequence,
        temperature_celsius,
        min_loop,
        gamma,
        energy_model,
        constraints,
    )
}

pub fn analyze_with_compiled_constraints(
    sequence: String,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    energy_model: EnergyModel,
    constraints: ConstraintModel,
) -> Result<AnalysisResult, RnaError> {
    let mfe = fold_mfe_with_constraints(&sequence, min_loop, &energy_model, &constraints)?;
    let ensemble = partition_with_constraints(
        &sequence,
        temperature_celsius,
        min_loop,
        &energy_model,
        &constraints,
    )?;
    let (centroid_structure, centroid_score) = decode_centroid_with_constraints(
        sequence.as_bytes(),
        min_loop,
        &ensemble.pair_probabilities,
        &constraints,
        &energy_model,
    )?;
    let (mea_structure, mea_score) = decode_mea_with_constraints(
        sequence.as_bytes(),
        min_loop,
        gamma,
        &ensemble.pair_probabilities,
        &ensemble.unpaired_probabilities,
        &constraints,
        &energy_model,
    )?;

    let ensemble_summary = summarize_ensemble(&ensemble);
    Ok(AnalysisResult {
        length: sequence.len(),
        sequence,
        temperature_celsius,
        model: ModelDescription {
            parameter_set: if energy_model.parameter_profile_name().is_some() {
                "Ribon normalized custom thermodynamic parameter overlay"
            } else {
                match energy_model.nucleic_acid() {
                    NucleicAcid::Rna => "Ribon RNAstructure 6.6 standard RNA parameter family",
                    NucleicAcid::Dna => "Ribon RNAstructure 6.6 standard DNA parameter family",
                }
            },
            parameter_profile_name: energy_model.parameter_profile_name().map(str::to_owned),
            parameter_fingerprint_sha256: energy_model
                .parameter_profile_fingerprint()
                .map(str::to_owned),
            mfe: energy_model.model_name(),
            ensemble: energy_model.ensemble_model_name(),
            dangles: energy_model.dangles(),
            salt_correction: energy_model.has_salt_correction(),
            salt_molar: energy_model.salt_molar(),
            note: model_note(
                energy_model.nucleic_acid(),
                energy_model.dangles(),
                energy_model.has_salt_correction(),
            ),
        },
        mfe_structure: mfe.structure,
        mfe_energy_kcal_mol: mfe.energy_kcal_mol,
        ensemble_free_energy_kcal_mol: ensemble.ensemble_free_energy_kcal_mol,
        partition_function: ensemble.partition_function,
        log_partition_function: ensemble.log_partition_function,
        centroid_structure,
        centroid_score,
        mea_structure,
        mea_score,
        mea_gamma: gamma,
        pair_probabilities: ensemble.pair_probabilities,
        unpaired_probabilities: ensemble.unpaired_probabilities,
        ensemble: ensemble_summary,
        constraints: constraints.summary(),
    })
}

/// MFE-only prediction for dangle models 0, 1, 2, and 3.
///
/// Integrated analysis accepts odd requested dangle models too; their
/// PF-derived quantities use exact fixed-structure enumeration.
pub fn fold_sequence(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
) -> Result<MfeResult, RnaError> {
    fold_sequence_with_constraints(
        sequence,
        temperature_celsius,
        min_loop,
        dangles,
        &ConstraintConfig::default(),
    )
}

/// MFE prediction under hard, soft, and probing-derived constraints.
pub fn fold_sequence_with_constraints(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    constraint_config: &ConstraintConfig,
) -> Result<MfeResult, RnaError> {
    fold_sequence_with_options(
        sequence,
        temperature_celsius,
        min_loop,
        dangles,
        1.021,
        constraint_config,
    )
}

/// MFE prediction with explicit monovalent salt molarity.
pub fn fold_sequence_with_options(
    sequence: &str,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    constraint_config: &ConstraintConfig,
) -> Result<MfeResult, RnaError> {
    if sequence.chars().any(|symbol| symbol == '&') {
        return Err(RnaError::InvalidOption(
            "fold accepts one strand; use supplied structures for multi-strand drawing".into(),
        ));
    }
    let sequence = normalize_sequence(sequence)?;
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    let constraints = ConstraintModel::compile(sequence.len(), constraint_config)?;
    fold_mfe_with_constraints(&sequence, min_loop, &model, &constraints)
}

fn model_note(nucleic_acid: NucleicAcid, dangles: u8, salt_correction: bool) -> &'static str {
    if nucleic_acid == NucleicAcid::Dna {
        debug_assert!(!salt_correction);
        return "RNAstructure 6.6 DNA free-energy and enthalpy tables with explicit temperature interpolation and no additional salt correction.";
    }
    match (dangles, salt_correction) {
        (0, false) => "MFE, supplied-structure evaluation, and partition function share the RNAstructure 6.6 standard RNA tables. The grammar is single-strand, pseudoknot-free, and uses dangles=0 at the 1.021 M salt-correction reference concentration.",
        (2, false) => "MFE, supplied-structure evaluation, and partition function share the RNAstructure 6.6 standard RNA tables and the double-dangle/mismatch convention at the 1.021 M salt-correction reference concentration.",
        (0, true) => "RNAstructure 6.6 RNA dangles=0 with published monovalent-salt corrections on stacks and hairpin, internal, and multibranch loops.",
        (2, true) => "RNAstructure 6.6 RNA dangles=2 with published monovalent-salt corrections on stacks and hairpin, internal, and multibranch loops.",
        (1, false) => "MFE, supplied-structure evaluation, and the exact fixed-structure ensemble use exclusive single dangling ends at the 1.021 M reference salt concentration.",
        (3, false) => "MFE and supplied-structure evaluation include coaxial stacking; the exact fixed-structure ensemble sums the same evaluated structure energies at the 1.021 M reference salt concentration.",
        (1, true) => "Single-dangle MFE and exact fixed-structure ensemble with published monovalent-salt corrections.",
        (3, true) => "Coaxial-stacking MFE and exact fixed-structure ensemble with published monovalent-salt corrections.",
        _ => unreachable!("validated dangle model"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_analysis_accepts_odd_requested_dangle_models() {
        let folded = fold_sequence("GGGAAACCC", 37.0, 3, 3).unwrap();
        assert_eq!(folded.structure, "(((...)))");
        assert_eq!(folded.dangles, 3);
        let analysis = analyze_with_dangles("GGGAAACCC", 37.0, 3, 1.0, 3).unwrap();
        assert_eq!(analysis.model.dangles, 3);
        assert!(analysis.model.ensemble.contains("exact fixed-structure"));
    }

    #[test]
    fn hard_constraints_flow_through_mfe_partition_and_decoding() {
        let config = ConstraintConfig {
            force_pairs: vec![PairConstraint { i: 1, j: 9 }],
            force_unpaired: vec![4, 5, 6],
            ..ConstraintConfig::default()
        };
        let folded = fold_sequence_with_constraints("GGGAAACCC", 37.0, 3, 2, &config).unwrap();
        assert_eq!(&folded.structure[0..1], "(");
        let model = EnergyModel::default();
        let compiled = ConstraintModel::compile(9, &config).unwrap();
        partition_with_constraints("GGGAAACCC", 37.0, 3, &model, &compiled).unwrap();
        let result = analyze_with_constraints("GGGAAACCC", 37.0, 3, 1.0, 2, &config).unwrap();
        assert_eq!(&result.mfe_structure[0..1], "(");
        assert_eq!(&result.mfe_structure[8..9], ")");
        let forced = result
            .pair_probabilities
            .iter()
            .find(|pair| pair.i == 1 && pair.j == 9)
            .unwrap();
        assert!((forced.probability - 1.0).abs() < 1.0e-10);
        assert_eq!(&result.centroid_structure[0..1], "(");
        assert_eq!(&result.mea_structure[8..9], ")");
        assert!(result.constraints.hard_constraints);
    }

    #[test]
    fn soft_pair_bonus_changes_both_mfe_and_ensemble() {
        let config = ConstraintConfig {
            soft: SoftConstraintConfig {
                pairs: vec![PairEnergy {
                    i: 1,
                    j: 5,
                    energy_kcal_mol: -10.0,
                }],
                ..SoftConstraintConfig::default()
            },
            ..ConstraintConfig::default()
        };
        let result = analyze_with_constraints("GAAAC", 37.0, 3, 1.0, 0, &config).unwrap();
        assert_eq!(result.mfe_structure, "(...)");
        let probability = result
            .pair_probabilities
            .iter()
            .find(|pair| pair.i == 1 && pair.j == 5)
            .unwrap()
            .probability;
        assert!(probability > 0.99);
        assert!(result.constraints.soft_constraints);
    }

    #[test]
    fn no_lonely_pairs_is_enforced_by_mfe_partition_and_evaluation() {
        let config = ConstraintConfig {
            no_lonely_pairs: true,
            soft: SoftConstraintConfig {
                pairs: vec![PairEnergy {
                    i: 1,
                    j: 5,
                    energy_kcal_mol: -10.0,
                }],
                ..SoftConstraintConfig::default()
            },
            ..ConstraintConfig::default()
        };
        let result = analyze_with_constraints("GAAAC", 37.0, 3, 1.0, 0, &config).unwrap();
        assert_eq!(result.mfe_structure, ".....");
        assert!(result.pair_probabilities.is_empty());

        let odd = fold_sequence_with_constraints("GAAAC", 37.0, 3, 1, &config).unwrap();
        assert_eq!(odd.structure, ".....");

        let compiled = ConstraintModel::compile(5, &config).unwrap();
        let error = EnergyModel::with_dangles(37.0, 0)
            .unwrap()
            .evaluate_with_constraints("GAAAC", "(...)", 3, &compiled)
            .unwrap_err();
        assert!(error.to_string().contains("lonely pair"));
    }
}
