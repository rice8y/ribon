use ribon_core::{
    accessibility_with_model, analyze_with_model, circular_fold_with_model, cofold_with_model,
    comparative_fold_with_model, conditional_density2, conditional_density2_ensemble,
    duplex_with_model, evaluate_conditional_density2, evaluate_conditional_density2_structure,
    evaluate_pseudoknot_structure_with_model, exact_landscape_path, fatgraph_topology,
    fold_gquad_with_model, fold_modified_with_model, inverse_fold_exact, layout_structure,
    ligand_ensemble_exact, local_fold_with_model, parameter_catalog, parameter_manifest,
    parse_structure, predict_pseudoknot_with_model, sample_boltzmann,
    sample_conditional_density2_with_constraints, suboptimal_conditional_density2_with_constraints,
    suboptimal_structures, AccessibilityWindow, ComparativeOptions, ConditionalDensity2Options,
    ConstraintConfig, ConstraintModel, EnergyModel, InverseDesignOptions, LayoutKind, LigandMotif,
    ModifiedBase, NucleicAcid, PseudoknotOptions, ThermodynamicParameterOverrides,
};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::str::FromStr;

const SCHEMA_VERSION: u8 = 1;

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u8,
    #[serde(default)]
    id: Option<String>,
    operation: String,
    input: Value,
    #[serde(default)]
    model: Model,
    #[serde(default)]
    constraints: ConstraintConfig,
    #[serde(default)]
    options: Value,
    #[serde(default)]
    execution: ExecutionPolicy,
}

/// Resource policy for document-time execution. Expensive exact algorithms
/// remain available, but require an explicit opt-in beyond conservative
/// limits so an accidental request cannot stall an entire Typst build.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields, rename_all = "kebab-case")]
struct ExecutionPolicy {
    allow_expensive: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
struct Model {
    id: String,
    parameter_base: NucleicAcid,
    parameter_overrides: Option<ThermodynamicParameterOverrides>,
    temperature_celsius: f64,
    min_loop: usize,
    dangles: u8,
    salt_molar: f64,
    mea_gamma: f64,
}

impl Default for Model {
    fn default() -> Self {
        Self {
            id: "ribon-turner-2004".into(),
            parameter_base: NucleicAcid::Rna,
            parameter_overrides: None,
            temperature_celsius: 37.0,
            min_loop: 3,
            dangles: 2,
            salt_molar: 1.021,
            mea_gamma: 1.0,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
struct Engine {
    name: &'static str,
    backend: &'static str,
    api: &'static str,
    model_id: &'static str,
    parameter_bundle_sha256: &'static str,
}

#[derive(Clone, Debug, Serialize)]
struct Payload {
    kind: String,
    data: Value,
}

#[derive(Clone, Debug, Serialize)]
struct ErrorBody {
    code: String,
    message: String,
}

#[derive(Clone, Debug, Serialize)]
struct Response {
    schema_version: u8,
    id: Option<String>,
    ok: bool,
    engine: Engine,
    operation: Option<String>,
    model: Option<Model>,
    constraints: Option<ConstraintConfig>,
    execution: Option<ExecutionPolicy>,
    result: Option<Payload>,
    error: Option<ErrorBody>,
}

#[derive(Debug)]
struct ApiError {
    code: &'static str,
    message: String,
}

impl ApiError {
    fn request(message: impl Into<String>) -> Self {
        Self {
            code: "invalid_request",
            message: message.into(),
        }
    }

    fn operation(message: impl Into<String>) -> Self {
        Self {
            code: "unsupported_operation",
            message: message.into(),
        }
    }

    fn analysis(message: impl Into<String>) -> Self {
        Self {
            code: "analysis_failed",
            message: message.into(),
        }
    }

    fn resource(message: impl Into<String>) -> Self {
        Self {
            code: "resource_limit",
            message: message.into(),
        }
    }
}

impl Engine {
    fn current() -> Self {
        Self {
            name: "ribon-thermo",
            backend: "pure-rust",
            api: "ribon.analysis/1",
            model_id: ribon_core::MODEL_ID,
            parameter_bundle_sha256: ribon_core::PARAMETER_BUNDLE_SHA256,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SequenceInput {
    sequence: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StructureInput {
    sequence: String,
    structure: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ConditionalStructureInput {
    sequence: String,
    seed_structure: String,
    added_structure: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DuplexInput {
    sequence_a: String,
    sequence_b: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AlignmentInput {
    alignment: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LandscapeInput {
    sequence: String,
    start_structure: String,
    target_structure: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InverseDesignInput {
    target_structure: String,
    template: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LigandInput {
    sequence: String,
    motifs: Vec<LigandMotif>,
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct LayoutOptions {
    method: String,
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct SamplingOptions {
    count: usize,
    seed: u64,
    unique: bool,
}

impl Default for SamplingOptions {
    fn default() -> Self {
        Self {
            count: 100,
            seed: 0,
            unique: false,
        }
    }
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct AccessibilityOptions {
    windows: Vec<AccessibilityWindow>,
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct SuboptimalOptions {
    energy_band_kcal_mol: f64,
    limit: usize,
}

impl Default for SuboptimalOptions {
    fn default() -> Self {
        Self {
            energy_band_kcal_mol: 5.0,
            limit: 100,
        }
    }
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct LocalOptions {
    window_size: usize,
    max_pair_span: usize,
    max_unpaired: usize,
}

impl Default for LocalOptions {
    fn default() -> Self {
        Self {
            window_size: 150,
            max_pair_span: 100,
            max_unpaired: 30,
        }
    }
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct CofoldOptions {
    concentration_a_molar: Option<f64>,
    concentration_b_molar: Option<f64>,
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ModifiedOptions {
    modifications: Vec<ModifiedBase>,
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ComparativeApiOptions {
    #[serde(flatten)]
    folding: ComparativeOptions,
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct PseudoknotApiOptions {
    #[serde(flatten)]
    decoding: PseudoknotOptions,
}

#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ConditionalDensity2ApiOptions {
    #[serde(flatten)]
    ensemble: ConditionalDensity2Options,
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ConditionalSamplingOptions {
    count: usize,
    seed: u64,
    unique: bool,
    #[serde(flatten)]
    ensemble: ConditionalDensity2Options,
}

impl Default for ConditionalSamplingOptions {
    fn default() -> Self {
        Self {
            count: 100,
            seed: 0,
            unique: false,
            ensemble: ConditionalDensity2Options::default(),
        }
    }
}

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ConditionalSuboptimalOptions {
    energy_band_kcal_mol: f64,
    limit: usize,
    #[serde(flatten)]
    ensemble: ConditionalDensity2Options,
}

impl Default for ConditionalSuboptimalOptions {
    fn default() -> Self {
        Self {
            energy_band_kcal_mol: 5.0,
            limit: 100,
            ensemble: ConditionalDensity2Options::default(),
        }
    }
}

fn from_value<T: DeserializeOwned>(value: Value, field: &str) -> Result<T, ApiError> {
    serde_json::from_value(value)
        .map_err(|error| ApiError::request(format!("invalid {field}: {error}")))
}

fn options<T: DeserializeOwned + Default>(value: Value) -> Result<T, ApiError> {
    if value.is_null() {
        Ok(T::default())
    } else {
        from_value(value, "options")
    }
}

fn value<T: Serialize>(value: T) -> Result<Value, ApiError> {
    serde_json::to_value(value)
        .map_err(|error| ApiError::analysis(format!("result serialization failed: {error}")))
}

const MAX_STANDARD_SEQUENCE_NT: usize = 500;
const MAX_LAYOUT_NT: usize = 20_000;
const MAX_LOCAL_SEQUENCE_NT: usize = 5_000;
const MAX_LOCAL_WINDOW_NT: usize = 200;
const MAX_DUPLEX_TOTAL_NT: usize = 400;
const MAX_DUPLEX_LENGTH_PRODUCT: usize = 40_000;
const MAX_ALIGNMENT_ROWS: usize = 128;
const MAX_SAMPLING_DRAWS: usize = 2_000;
const MAX_SUBOPTIMAL_STRUCTURES: usize = 500;
const MAX_ACCESSIBILITY_WINDOWS: usize = 2_000;
const MAX_ODD_DANGLE_EXACT_NT: usize = 24;
const MAX_EXACT_ARBITRARY_NT: usize = 18;
const MAX_LANDSCAPE_NT: usize = 14;
const MAX_LIGAND_NT: usize = 18;
const MAX_INVERSE_CANDIDATES: u128 = 262_144;

fn nucleotide_count(sequence: &str) -> usize {
    sequence
        .chars()
        .filter(|symbol| !symbol.is_whitespace() && *symbol != '&')
        .count()
}

fn input_sequence(request: &Request) -> Option<&str> {
    request.input.get("sequence").and_then(Value::as_str)
}

fn inverse_candidate_upper_bound(template: &str) -> u128 {
    template.chars().fold(1_u128, |count, symbol| {
        let choices = match symbol.to_ascii_uppercase() {
            'N' => 4,
            'B' | 'D' | 'H' | 'V' => 3,
            'R' | 'Y' | 'S' | 'W' | 'K' | 'M' => 2,
            'A' | 'C' | 'G' | 'U' | 'T' => 1,
            _ => 4,
        };
        count.saturating_mul(choices)
    })
}

fn resource_error(operation: &str, detail: impl AsRef<str>) -> ApiError {
    ApiError::resource(format!(
        "{operation}: {}; pass execution.allow_expensive=true only when this exact document-time computation is intentional",
        detail.as_ref()
    ))
}

fn enforce_resource_policy(request: &Request) -> Result<(), ApiError> {
    if request.execution.allow_expensive {
        return Ok(());
    }
    let operation = request.operation.as_str();
    let n = input_sequence(request).map(nucleotide_count).unwrap_or(0);

    if matches!(operation, "layout" | "validate" | "fatgraph-topology") {
        if n > MAX_LAYOUT_NT {
            return Err(resource_error(
                operation,
                format!("{n} nt exceeds the document layout limit of {MAX_LAYOUT_NT} nt"),
            ));
        }
    } else if operation == "local" {
        let options: LocalOptions = options(request.options.clone())?;
        if n > MAX_LOCAL_SEQUENCE_NT || options.window_size > MAX_LOCAL_WINDOW_NT {
            return Err(resource_error(
                operation,
                format!(
                    "{n} nt with window {} exceeds the safe local limits ({MAX_LOCAL_SEQUENCE_NT} nt, window {MAX_LOCAL_WINDOW_NT})",
                    options.window_size
                ),
            ));
        }
    } else if operation == "duplex" || operation == "cofold" {
        let input: DuplexInput = from_value(request.input.clone(), "input")?;
        let a = nucleotide_count(&input.sequence_a);
        let b = nucleotide_count(&input.sequence_b);
        if a + b > MAX_DUPLEX_TOTAL_NT || a.saturating_mul(b) > MAX_DUPLEX_LENGTH_PRODUCT {
            return Err(resource_error(
                operation,
                format!(
                    "strand lengths {a}+{b} exceed the safe total/product limits ({MAX_DUPLEX_TOTAL_NT}, {MAX_DUPLEX_LENGTH_PRODUCT})"
                ),
            ));
        }
    } else if operation == "comparative" {
        let input: AlignmentInput = from_value(request.input.clone(), "input")?;
        let width = input
            .alignment
            .first()
            .map_or(0, |row| nucleotide_count(row));
        if input.alignment.len() > MAX_ALIGNMENT_ROWS || width > MAX_STANDARD_SEQUENCE_NT {
            return Err(resource_error(
                operation,
                format!(
                    "alignment {}x{width} exceeds the safe limits ({MAX_ALIGNMENT_ROWS} rows, {MAX_STANDARD_SEQUENCE_NT} columns)",
                    input.alignment.len()
                ),
            ));
        }
    } else if n > MAX_STANDARD_SEQUENCE_NT {
        return Err(resource_error(
            operation,
            format!("{n} nt exceeds the document analysis limit of {MAX_STANDARD_SEQUENCE_NT} nt"),
        ));
    }

    if operation == "sample" {
        let options: SamplingOptions = options(request.options.clone())?;
        if options.count > MAX_SAMPLING_DRAWS {
            return Err(resource_error(
                operation,
                format!(
                    "{} draws exceed the safe limit of {MAX_SAMPLING_DRAWS}",
                    options.count
                ),
            ));
        }
    }
    if operation == "suboptimal" || operation == "conditional-density2-suboptimal" {
        let limit = request
            .options
            .get("limit")
            .and_then(Value::as_u64)
            .unwrap_or(100) as usize;
        if limit > MAX_SUBOPTIMAL_STRUCTURES {
            return Err(resource_error(
                operation,
                format!("{limit} structures exceed the safe limit of {MAX_SUBOPTIMAL_STRUCTURES}"),
            ));
        }
    }
    if operation == "accessibility" {
        let options: AccessibilityOptions = options(request.options.clone())?;
        if options.windows.len() > MAX_ACCESSIBILITY_WINDOWS {
            return Err(resource_error(
                operation,
                format!(
                    "{} windows exceed the safe limit of {MAX_ACCESSIBILITY_WINDOWS}",
                    options.windows.len()
                ),
            ));
        }
    }
    if operation == "conditional-density2-sample" {
        let options: ConditionalSamplingOptions = options(request.options.clone())?;
        if options.count > MAX_SAMPLING_DRAWS {
            return Err(resource_error(
                operation,
                format!(
                    "{} draws exceed the safe limit of {MAX_SAMPLING_DRAWS}",
                    options.count
                ),
            ));
        }
    }

    let odd_ensemble = request.model.dangles % 2 == 1
        && matches!(
            operation,
            "analyze"
                | "sample"
                | "accessibility"
                | "suboptimal"
                | "cofold"
                | "local"
                | "circular"
                | "modified"
                | "gquad"
                | "pseudoknot"
                | "conditional-density2"
                | "conditional-density2-sample"
                | "conditional-density2-suboptimal"
        );
    if odd_ensemble && n > MAX_ODD_DANGLE_EXACT_NT {
        return Err(resource_error(
            operation,
            format!(
                "the exact odd-dangle ensemble is exponential and {n} nt exceeds its safe limit of {MAX_ODD_DANGLE_EXACT_NT} nt"
            ),
        ));
    }

    if matches!(
        operation,
        "conditional-density2-oracle" | "evaluate-conditional-density2-oracle"
    ) && n > MAX_EXACT_ARBITRARY_NT
    {
        return Err(resource_error(
            operation,
            format!("the exhaustive oracle is limited to {MAX_EXACT_ARBITRARY_NT} nt"),
        ));
    }
    if operation == "pseudoknot"
        && request
            .options
            .get("exact-arbitrary-ensemble")
            .and_then(Value::as_bool)
            .unwrap_or(false)
        && n > MAX_EXACT_ARBITRARY_NT
    {
        return Err(resource_error(
            operation,
            format!("the arbitrary matching ensemble is limited to {MAX_EXACT_ARBITRARY_NT} nt"),
        ));
    }
    if operation == "landscape" && n > MAX_LANDSCAPE_NT {
        return Err(resource_error(
            operation,
            format!("complete landscape enumeration is limited to {MAX_LANDSCAPE_NT} nt"),
        ));
    }
    if operation == "ligand" && n > MAX_LIGAND_NT {
        return Err(resource_error(
            operation,
            format!("the exact joint ligand ensemble is limited to {MAX_LIGAND_NT} nt"),
        ));
    }
    if operation == "inverse-design" {
        let input: InverseDesignInput = from_value(request.input.clone(), "input")?;
        let candidates = inverse_candidate_upper_bound(&input.template);
        if candidates > MAX_INVERSE_CANDIDATES {
            return Err(resource_error(
                operation,
                format!(
                    "the IUPAC template permits up to {candidates} sequences, above the safe limit of {MAX_INVERSE_CANDIDATES}"
                ),
            ));
        }
    }
    Ok(())
}

fn validate_model(model: &Model) -> Result<(), ApiError> {
    if model.id != ribon_core::MODEL_ID
        && model.id != ribon_core::DNA_MODEL_ID
        && model.id != ribon_core::CUSTOM_MODEL_ID
    {
        return Err(ApiError::request(format!(
            "unknown model {:?}; expected {:?}, {:?}, or {:?}",
            model.id,
            ribon_core::MODEL_ID,
            ribon_core::DNA_MODEL_ID,
            ribon_core::CUSTOM_MODEL_ID,
        )));
    }
    if (model.id == ribon_core::CUSTOM_MODEL_ID) != model.parameter_overrides.is_some() {
        return Err(ApiError::request(
            "custom model id and model.parameter_overrides must be supplied together",
        ));
    }
    if !model.temperature_celsius.is_finite() || model.temperature_celsius <= -273.15 {
        return Err(ApiError::request(
            "model.temperature_celsius must be finite and above absolute zero",
        ));
    }
    if !matches!(model.dangles, 0..=3) {
        return Err(ApiError::request(
            "model.dangles must be one of 0, 1, 2, or 3",
        ));
    }
    if !model.salt_molar.is_finite() || model.salt_molar <= 0.0 {
        return Err(ApiError::request(
            "model.salt_molar must be finite and positive",
        ));
    }
    if !model.mea_gamma.is_finite() || model.mea_gamma <= 0.0 {
        return Err(ApiError::request(
            "model.mea_gamma must be finite and positive",
        ));
    }
    Ok(())
}

fn energy_model(model: &Model) -> Result<EnergyModel, ApiError> {
    energy_model_with_dangles(model, model.dangles)
}

fn energy_model_with_dangles(model: &Model, dangles: u8) -> Result<EnergyModel, ApiError> {
    let family = selected_family(model);
    if let Some(overrides) = &model.parameter_overrides {
        return EnergyModel::with_parameter_overrides(
            model.temperature_celsius,
            dangles,
            model.salt_molar,
            family,
            overrides.clone(),
        )
        .map_err(|error| ApiError::analysis(error.to_string()));
    }
    EnergyModel::with_parameter_family(model.temperature_celsius, dangles, model.salt_molar, family)
        .map_err(|error| ApiError::analysis(error.to_string()))
}

fn selected_family(model: &Model) -> NucleicAcid {
    if model.id == ribon_core::CUSTOM_MODEL_ID {
        model.parameter_base
    } else if model.id == ribon_core::DNA_MODEL_ID {
        NucleicAcid::Dna
    } else {
        NucleicAcid::Rna
    }
}

fn require_rna(model: &Model, operation: &str) -> Result<(), ApiError> {
    if selected_family(model) == NucleicAcid::Dna {
        return Err(ApiError::request(format!(
            "operation {operation:?} is RNA-chemistry-specific; select {:?}",
            ribon_core::MODEL_ID
        )));
    }
    Ok(())
}

fn configure_conditional_model(options: &mut ConditionalDensity2Options, model: &Model) {
    options.nucleic_acid = selected_family(model);
    options.parameter_overrides = model.parameter_overrides.clone();
}

fn dispatch(request: &Request) -> Result<Payload, ApiError> {
    let model = &request.model;
    let map_error = |error: ribon_core::RnaError| ApiError::analysis(error.to_string());
    let (kind, data) = match request.operation.as_str() {
        "analyze" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let result = analyze_with_model(
                input.sequence,
                model.min_loop,
                model.mea_gamma,
                energy_model(model)?,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("analysis", value(result)?)
        }
        "fold" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let sequence = ribon_core::normalize_sequence(&input.sequence).map_err(map_error)?;
            let constraints = ConstraintModel::compile(sequence.len(), &request.constraints)
                .map_err(map_error)?;
            let energy = energy_model(model)?;
            let result = ribon_core::fold_mfe_with_constraints(
                &sequence,
                model.min_loop,
                &energy,
                &constraints,
            )
            .map_err(map_error)?;
            ("mfe", value(result)?)
        }
        "evaluate" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let parsed = parse_structure(&input.sequence, &input.structure).map_err(map_error)?;
            let constraints =
                ConstraintModel::compile(parsed.length, &request.constraints).map_err(map_error)?;
            let result = energy_model(model)?
                .evaluate_with_constraints(
                    &input.sequence,
                    &input.structure,
                    model.min_loop,
                    &constraints,
                )
                .map_err(map_error)?;
            let mut data = value(result)?;
            if let Some(object) = data.as_object_mut() {
                object.insert("sequence".into(), Value::String(input.sequence));
                object.insert("structure".into(), Value::String(input.structure));
            }
            ("energy", data)
        }
        "validate" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            (
                "structure",
                value(parse_structure(&input.sequence, &input.structure).map_err(map_error)?)?,
            )
        }
        "layout" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let options: LayoutOptions = options(request.options.clone())?;
            let method = if options.method.is_empty() {
                "naview"
            } else {
                &options.method
            };
            let kind = LayoutKind::from_str(method).map_err(map_error)?;
            (
                "layout",
                value(
                    layout_structure(&input.sequence, &input.structure, kind).map_err(map_error)?,
                )?,
            )
        }
        "sample" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: SamplingOptions = options(request.options.clone())?;
            let sequence = ribon_core::normalize_sequence(&input.sequence).map_err(map_error)?;
            let energy = energy_model(model)?;
            let constraints = ConstraintModel::compile(sequence.len(), &request.constraints)
                .map_err(map_error)?;
            let result = sample_boltzmann(
                &sequence,
                model.temperature_celsius,
                model.min_loop,
                &energy,
                &constraints,
                options.count,
                options.seed,
                options.unique,
            )
            .map_err(map_error)?;
            ("samples", value(result)?)
        }
        "accessibility" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: AccessibilityOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = accessibility_with_model(
                &input.sequence,
                model.min_loop,
                &energy,
                &request.constraints,
                &options.windows,
            )
            .map_err(map_error)?;
            ("accessibility", value(result)?)
        }
        "suboptimal" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: SuboptimalOptions = options(request.options.clone())?;
            let sequence = ribon_core::normalize_sequence(&input.sequence).map_err(map_error)?;
            let energy = energy_model(model)?;
            let constraints = ConstraintModel::compile(sequence.len(), &request.constraints)
                .map_err(map_error)?;
            let result = suboptimal_structures(
                &sequence,
                model.min_loop,
                &energy,
                &constraints,
                options.energy_band_kcal_mol,
                options.limit,
            )
            .map_err(map_error)?;
            ("suboptimal", value(result)?)
        }
        "duplex" => {
            let input: DuplexInput = from_value(request.input.clone(), "input")?;
            let energy = energy_model_with_dangles(model, 0)?;
            let result = duplex_with_model(&input.sequence_a, &input.sequence_b, &energy)
                .map_err(map_error)?;
            ("duplex", value(result)?)
        }
        "cofold" => {
            let input: DuplexInput = from_value(request.input.clone(), "input")?;
            let options: CofoldOptions = options(request.options.clone())?;
            let concentrations =
                match (options.concentration_a_molar, options.concentration_b_molar) {
                    (None, None) => None,
                    (Some(a), Some(b)) => Some((a, b)),
                    _ => {
                        return Err(ApiError::request(
                            "both cofold concentrations must be supplied together",
                        ))
                    }
                };
            let result = cofold_with_model(
                &input.sequence_a,
                &input.sequence_b,
                model.min_loop,
                model.mea_gamma,
                energy_model(model)?,
                concentrations,
            )
            .map_err(map_error)?;
            ("cofold", value(result)?)
        }
        "local" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: LocalOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = local_fold_with_model(
                &input.sequence,
                model.min_loop,
                &energy,
                options.window_size,
                options.max_pair_span,
                options.max_unpaired,
            )
            .map_err(map_error)?;
            ("local", value(result)?)
        }
        "circular" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let energy = energy_model(model)?;
            let result = circular_fold_with_model(
                &input.sequence,
                model.min_loop,
                model.mea_gamma,
                &energy,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("circular", value(result)?)
        }
        "modified" => {
            require_rna(model, "modified")?;
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: ModifiedOptions = options(request.options.clone())?;
            let result = fold_modified_with_model(
                &input.sequence,
                &options.modifications,
                model.min_loop,
                model.mea_gamma,
                energy_model(model)?,
            )
            .map_err(map_error)?;
            ("modified", value(result)?)
        }
        "gquad" => {
            require_rna(model, "gquad")?;
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let result = fold_gquad_with_model(
                &input.sequence,
                model.min_loop,
                model.mea_gamma,
                energy_model(model)?,
            )
            .map_err(map_error)?;
            ("gquad", value(result)?)
        }
        "pseudoknot" => {
            let input: SequenceInput = from_value(request.input.clone(), "input")?;
            let options: PseudoknotApiOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = predict_pseudoknot_with_model(
                &input.sequence,
                model.min_loop,
                &energy,
                &options.decoding,
            )
            .map_err(map_error)?;
            ("pseudoknot", value(result)?)
        }
        "evaluate-pseudoknot" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let options: PseudoknotApiOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = evaluate_pseudoknot_structure_with_model(
                &input.sequence,
                &input.structure,
                model.min_loop,
                &energy,
                &options.decoding,
            )
            .map_err(map_error)?;
            ("pseudoknot-energy", value(result)?)
        }
        "conditional-density2" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalDensity2ApiOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = conditional_density2(
                &input.sequence,
                &input.structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("conditional-density2", value(result)?)
        }
        "conditional-density2-sample" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalSamplingOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = sample_conditional_density2_with_constraints(
                &input.sequence,
                &input.structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
                options.count,
                options.seed,
                options.unique,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("conditional-density2-samples", value(result)?)
        }
        "conditional-density2-suboptimal" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalSuboptimalOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = suboptimal_conditional_density2_with_constraints(
                &input.sequence,
                &input.structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
                options.energy_band_kcal_mol,
                options.limit,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("conditional-density2-suboptimal", value(result)?)
        }
        "conditional-density2-oracle" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalDensity2ApiOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = conditional_density2_ensemble(
                &input.sequence,
                &input.structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
            )
            .map_err(map_error)?;
            ("conditional-density2-oracle", value(result)?)
        }
        "evaluate-conditional-density2" => {
            let input: ConditionalStructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalDensity2ApiOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = evaluate_conditional_density2(
                &input.sequence,
                &input.seed_structure,
                &input.added_structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
                &request.constraints,
            )
            .map_err(map_error)?;
            ("conditional-density2-energy", value(result)?)
        }
        "evaluate-conditional-density2-oracle" => {
            let input: ConditionalStructureInput = from_value(request.input.clone(), "input")?;
            let mut options: ConditionalDensity2ApiOptions = options(request.options.clone())?;
            configure_conditional_model(&mut options.ensemble, model);
            let result = evaluate_conditional_density2_structure(
                &input.sequence,
                &input.seed_structure,
                &input.added_structure,
                model.temperature_celsius,
                model.min_loop,
                model.dangles,
                model.salt_molar,
                &options.ensemble,
            )
            .map_err(map_error)?;
            ("conditional-density2-oracle-energy", value(result)?)
        }
        "comparative" => {
            let input: AlignmentInput = from_value(request.input.clone(), "input")?;
            let options: ComparativeApiOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = comparative_fold_with_model(
                &input.alignment,
                model.min_loop,
                model.mea_gamma,
                &energy,
                &options.folding,
            )
            .map_err(map_error)?;
            ("comparative", value(result)?)
        }
        "landscape" => {
            let input: LandscapeInput = from_value(request.input.clone(), "input")?;
            let sequence = ribon_core::normalize_sequence(&input.sequence).map_err(map_error)?;
            let constraints = ConstraintModel::compile(sequence.len(), &request.constraints)
                .map_err(map_error)?;
            let energy = energy_model(model)?;
            let result = exact_landscape_path(
                &sequence,
                &input.start_structure,
                &input.target_structure,
                model.min_loop,
                &energy,
                &constraints,
            )
            .map_err(map_error)?;
            ("landscape", value(result)?)
        }
        "inverse-design" => {
            let input: InverseDesignInput = from_value(request.input.clone(), "input")?;
            let options: InverseDesignOptions = options(request.options.clone())?;
            let energy = energy_model(model)?;
            let result = inverse_fold_exact(
                &input.target_structure,
                &input.template,
                model.min_loop,
                &energy,
                &request.constraints,
                &options,
            )
            .map_err(map_error)?;
            ("inverse-design", value(result)?)
        }
        "ligand" => {
            let input: LigandInput = from_value(request.input.clone(), "input")?;
            let sequence = ribon_core::normalize_sequence(&input.sequence).map_err(map_error)?;
            let constraints = ConstraintModel::compile(sequence.len(), &request.constraints)
                .map_err(map_error)?;
            let energy = energy_model(model)?;
            let result = ligand_ensemble_exact(
                &sequence,
                &input.motifs,
                model.min_loop,
                model.mea_gamma,
                &energy,
                &constraints,
            )
            .map_err(map_error)?;
            ("ligand", value(result)?)
        }
        "fatgraph-topology" => {
            let input: StructureInput = from_value(request.input.clone(), "input")?;
            let result = fatgraph_topology(&input.sequence, &input.structure).map_err(map_error)?;
            let mut data = value(result)?;
            if let Some(object) = data.as_object_mut() {
                object.insert("sequence".into(), Value::String(input.sequence));
                object.insert("structure".into(), Value::String(input.structure));
            }
            ("fatgraph-topology", data)
        }
        "parameters" => (
            "parameters",
            json!({
                "active": if selected_family(model) == NucleicAcid::Dna {
                    ribon_core::dna_parameter_manifest()
                } else {
                    parameter_manifest()
                },
                "custom_overlay": model.parameter_overrides,
                "catalog": parameter_catalog(),
            }),
        ),
        "capabilities" => (
            "capabilities",
            json!({
                "operations": [
                    "analyze", "fold", "evaluate", "validate", "layout", "sample",
                    "accessibility", "suboptimal", "duplex", "cofold", "local",
                    "circular", "modified", "gquad", "pseudoknot", "evaluate-pseudoknot",
                    "comparative", "parameters", "capabilities"
                ],
                "advanced_operations": [
                    "conditional-density2", "conditional-density2-sample",
                    "conditional-density2-suboptimal", "evaluate-conditional-density2",
                    "landscape", "inverse-design", "ligand", "fatgraph-topology"
                ],
                "layouts": ["naview", "simple", "circular", "turtle", "puzzler", "linear"],
                "thermodynamic_models": [
                    ribon_core::MODEL_ID,
                    ribon_core::DNA_MODEL_ID, ribon_core::CUSTOM_MODEL_ID,
                ],
                "decoders": [
                    "mfe", "centroid", "mea", "probknot-mutual-maximum", "thermodynamic-h-type-seeded-hybrid",
                    "arbitrary-topology-matching-centroid", "arbitrary-topology-matching-mea",
                    "multi-component-h-type-centroid", "multi-component-h-type-mea",
                    "exact-arbitrary-topology-mfe", "exact-arbitrary-topology-centroid",
                    "exact-arbitrary-topology-mea"
                    , "conditional-density2-mfe", "conditional-density2-centroid",
                    "conditional-density2-mea"
                    , "conditional-density2-boltzmann-sampling",
                    "conditional-density2-k-best"
                ],
                "model": {
                    "id": ribon_core::MODEL_ID,
                    "parameter_bundle_sha256": ribon_core::PARAMETER_BUNDLE_SHA256,
                    "mfe_dangles": [0, 1, 2, 3],
                    "ensemble_dangles": [0, 1, 2, 3],
                    "circular_dangles": [0, 1, 2, 3],
                },
                "limits": {
                    "single_sequence_nt": MAX_STANDARD_SEQUENCE_NT,
                    "layout_nt": MAX_LAYOUT_NT,
                    "local_sequence_nt": MAX_LOCAL_SEQUENCE_NT,
                    "local_window_nt": MAX_LOCAL_WINDOW_NT,
                    "duplex_total_nt": MAX_DUPLEX_TOTAL_NT,
                    "duplex_length_product": MAX_DUPLEX_LENGTH_PRODUCT,
                    "alignment_rows": MAX_ALIGNMENT_ROWS,
                    "accessibility_windows": MAX_ACCESSIBILITY_WINDOWS,
                    "sampling_draws": MAX_SAMPLING_DRAWS,
                    "suboptimal_requested_count": MAX_SUBOPTIMAL_STRUCTURES,
                    "odd_dangle_exact_nt": MAX_ODD_DANGLE_EXACT_NT,
                    "arbitrary_matching_exact_nt": MAX_EXACT_ARBITRARY_NT,
                    "landscape_nt": MAX_LANDSCAPE_NT,
                    "inverse_design_candidate_upper_bound": MAX_INVERSE_CANDIDATES,
                    "ligand_nt": MAX_LIGAND_NT,
                    "override": "execution.allow_expensive=true",
                },
                "boundaries": {
                    "global_ensemble": "pseudoknot-free",
                    "odd_dangle_ensemble": "complete planar fixed-structure enumeration with exclusive single dangles/coaxial stacks",
                    "cofold": "complete cut-point noncrossing state space; polynomial for dangles 0/2 and exhaustive for 1/3",
                    "pseudoknot_ensemble": "exact weighted-interval ensemble of all span-disjoint H-type component sets plus open state; optional exact exponential enumeration of every canonical arbitrary-topology matching",
                    "conditional_density2": "complete unambiguous G-conditioned density-2 interval hypergraph for dangles 0/2; O(n^3) time, O(n^2) space, log-domain inside/outside, MFE, base-pair marginals, centroid, and MEA",
                    "conditional_density2_odd_dangles": "exact exhaustive fixed-structure dispatcher for nonlocal exclusive single-dangle/coaxial models 1/3; no approximation or state cap",
                    "gquad_ensemble": "integrated multiple/nested exterior and generalized-multiloop G-quadruplex grammar",
                    "modified_bases": "experimental sparse m6A/pseudouridine/inosine/7-deazaadenosine/purine deltas from an explicit canonical reference state, model-based dihydrouridine correction, and explicit position pseudo-energies",
                    "comparative": "alignment-row Turner MFE/McCaskill grammar with column covariation",
                    "landscape": "complete constraint-filtered planar state graph with exact single-pair-move minimax path; exponential without beam or state cap",
                    "inverse_design": "complete IUPAC sequence-space search maximizing exact constrained-partition target probability; output count never truncates search",
                    "ligand": "complete planar structures times every non-overlapping ligand-site microstate with standard-state concentration correction",
                },
                "guarantees": {
                    "implicit_state_caps": false,
                    "internal_loop_bound": "complete geometric range by default",
                    "pseudoknot_state_count": "exact decimal string plus runtime-usize saturated numeric field",
                    "planar_puzzler_intersections": 0,
                    "extended_dot_bracket_coloring": "valid greedy coloring with exact 30-colorability fallback at the standard alphabet boundary",
                },
                "native_vector_renderer": true,
            }),
        ),
        other => {
            return Err(ApiError::operation(format!(
                "unknown operation {other:?}; call capabilities for the supported set"
            )))
        }
    };
    Ok(Payload {
        kind: kind.into(),
        data,
    })
}

fn encode(response: &Response) -> Vec<u8> {
    serde_json::to_vec(response).unwrap_or_else(|_| {
        br#"{"schema_version": 1,"ok":false,"error":{"code":"internal_error","message":"response serialization failed"}}"#.to_vec()
    })
}

pub fn execute_bytes(input: &[u8]) -> Vec<u8> {
    let request = match serde_json::from_slice::<Request>(input) {
        Ok(request) => request,
        Err(error) => {
            return encode(&Response {
                schema_version: SCHEMA_VERSION,
                id: None,
                ok: false,
                engine: Engine::current(),
                operation: None,
                model: None,
                constraints: None,
                execution: None,
                result: None,
                error: Some(ErrorBody {
                    code: "invalid_request".into(),
                    message: format!("invalid request JSON: {error}"),
                }),
            })
        }
    };
    let operation = request.operation.clone();
    let model = request.model.clone();
    let constraints = request.constraints.clone();
    let execution = request.execution.clone();
    let outcome = if request.schema_version != SCHEMA_VERSION {
        Err(ApiError::request(format!(
            "unsupported schema_version {}; expected {SCHEMA_VERSION}",
            request.schema_version
        )))
    } else {
        validate_model(&request.model)
            .and_then(|_| enforce_resource_policy(&request))
            .and_then(|_| dispatch(&request))
    };
    match outcome {
        Ok(result) => encode(&Response {
            schema_version: SCHEMA_VERSION,
            id: request.id,
            ok: true,
            engine: Engine::current(),
            operation: Some(operation),
            model: Some(model),
            constraints: Some(constraints),
            execution: Some(execution),
            result: Some(result),
            error: None,
        }),
        Err(error) => encode(&Response {
            schema_version: SCHEMA_VERSION,
            id: request.id,
            ok: false,
            engine: Engine::current(),
            operation: Some(operation),
            model: Some(model),
            constraints: Some(constraints),
            execution: Some(execution),
            result: None,
            error: Some(ErrorBody {
                code: error.code.into(),
                message: error.message,
            }),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(request: Value) -> Value {
        serde_json::from_slice(&execute_bytes(&serde_json::to_vec(&request).unwrap())).unwrap()
    }

    #[test]
    fn capabilities_has_a_stable_envelope() {
        let response = call(json!({
            "schema_version": 1,
            "operation": "capabilities",
            "input": {}
        }));
        assert_eq!(response["ok"], true);
        assert_eq!(response["engine"]["backend"], "pure-rust");
        assert_eq!(
            response["engine"]["parameter_bundle_sha256"],
            ribon_core::PARAMETER_BUNDLE_SHA256
        );
        assert_eq!(response["result"]["kind"], "capabilities");
        let capabilities = &response["result"]["data"];
        assert_eq!(capabilities["operations"].as_array().unwrap().len(), 19);
        assert_eq!(
            capabilities["advanced_operations"]
                .as_array()
                .unwrap()
                .len(),
            8
        );
        let advertised = capabilities["operations"]
            .as_array()
            .unwrap()
            .iter()
            .chain(capabilities["advanced_operations"].as_array().unwrap());
        assert!(!advertised
            .clone()
            .any(|name| name.as_str().unwrap().contains("oracle")));
    }

    #[test]
    fn resource_policy_rejects_accidental_document_time_explosions() {
        let long = "A".repeat(MAX_STANDARD_SEQUENCE_NT + 1);
        let sequence = "GGGAAACCC";
        let cases = [
            json!({
                "schema_version": 1,
                "operation": "analyze",
                "input": {"sequence": long}
            }),
            json!({
                "schema_version": 1,
                "operation": "sample",
                "input": {"sequence": sequence},
                "options": {"count": MAX_SAMPLING_DRAWS + 1, "seed": 0}
            }),
            json!({
                "schema_version": 1,
                "operation": "landscape",
                "input": {
                    "sequence": "GGGAAACCCGGGAAA",
                    "start_structure": "...............",
                    "target_structure": "..............."
                }
            }),
        ];
        for request in cases {
            let response = call(request);
            assert_eq!(response["ok"], false, "{response:#}");
            assert_eq!(response["error"]["code"], "resource_limit", "{response:#}");
            assert!(response["error"]["message"]
                .as_str()
                .unwrap()
                .contains("execution.allow_expensive=true"));
        }
    }

    #[test]
    fn explicit_expensive_opt_in_changes_only_the_resource_policy() {
        let response = call(json!({
            "schema_version": 1,
            "operation": "sample",
            "input": {"sequence": "GGGAAACCC"},
            "options": {"count": MAX_SAMPLING_DRAWS + 1, "seed": 7},
            "execution": {"allow-expensive": true}
        }));
        assert_eq!(response["ok"], true, "{response:#}");
        assert_eq!(response["result"]["kind"], "samples");
        assert_eq!(
            response["result"]["data"]["returned"],
            MAX_SAMPLING_DRAWS + 1
        );
        assert_eq!(response["execution"]["allow-expensive"], true);
        assert_eq!(response["model"]["id"], ribon_core::MODEL_ID);
    }

    #[test]
    fn analysis_result_flows_to_layout() {
        let analysis = call(json!({
            "schema_version": 1,
            "id": "hairpin",
            "operation": "analyze",
            "input": {"sequence": "GGGAAACCC"}
        }));
        assert_eq!(analysis["ok"], true, "{analysis:#}");
        let sequence = analysis["result"]["data"]["sequence"].clone();
        let structure = analysis["result"]["data"]["mfe_structure"].clone();
        let layout = call(json!({
            "schema_version": 1,
            "operation": "layout",
            "input": {"sequence": sequence, "structure": structure},
            "options": {"method": "naview"}
        }));
        assert_eq!(layout["ok"], true, "{layout:#}");
        assert_eq!(
            layout["result"]["data"]["points"].as_array().unwrap().len(),
            9
        );
    }

    #[test]
    fn exact_landscape_flows_through_the_stable_protocol() {
        let response = call(json!({
            "schema_version": 1,
            "operation": "landscape",
            "input": {
                "sequence": "GGGAAACCC",
                "start_structure": ".........",
                "target_structure": "(((...)))"
            },
            "model": {"dangles": 2}
        }));
        assert_eq!(response["ok"], true, "{response:#}");
        assert_eq!(response["result"]["kind"], "landscape");
        assert_eq!(response["result"]["data"]["state_space_complete"], true);
        assert_eq!(
            response["result"]["data"]["path"][0]["structure"],
            "........."
        );
    }

    #[test]
    fn exact_inverse_design_flows_through_the_stable_protocol() {
        let response = call(json!({
            "schema_version": 1,
            "operation": "inverse-design",
            "input": {"target_structure": "(...)", "template": "NNNNN"},
            "model": {"dangles": 0},
            "options": {"return-count": 3}
        }));
        assert_eq!(response["ok"], true, "{response:#}");
        assert_eq!(response["result"]["kind"], "inverse-design");
        assert_eq!(response["result"]["data"]["candidate_sequence_count"], 384);
        assert_eq!(response["result"]["data"]["search_complete"], true);
        assert_eq!(response["result"]["data"]["returned"], 3);
    }

    #[test]
    fn exact_ligand_ensemble_flows_through_the_stable_protocol() {
        let response = call(json!({
            "schema_version": 1,
            "operation": "ligand",
            "input": {
                "sequence": "GGGAAACCC",
                "motifs": [{
                    "id": "aptamer",
                    "start": 1,
                    "sequence": "GGGAAACCC",
                    "structure": "(((...)))",
                    "standard-binding-energy-kcal-mol": -3.0,
                    "concentration-molar": 1.0
                }]
            },
            "model": {"dangles": 0}
        }));
        assert_eq!(response["ok"], true, "{response:#}");
        assert_eq!(response["result"]["kind"], "ligand");
        assert_eq!(response["result"]["data"]["state_space_complete"], true);
        assert!(
            response["result"]["data"]["motifs"][0]["occupancy_probability"]
                .as_f64()
                .is_some_and(|probability| probability > 0.0 && probability < 1.0)
        );
    }

    #[test]
    fn dna_model_flows_through_analysis_evaluation_and_conditional_pf() {
        let model = json!({
            "id": ribon_core::DNA_MODEL_ID,
            "temperature_celsius": 37.0,
            "min_loop": 3,
            "dangles": 0,
            "salt_molar": 1.021,
            "mea_gamma": 1.0
        });
        let analysis = call(json!({
            "schema_version": 1,
            "operation": "analyze",
            "input": {"sequence": "GGGTTTCCC"},
            "model": model.clone()
        }));
        assert_eq!(analysis["ok"], true, "{analysis:#}");
        assert_eq!(analysis["result"]["data"]["mfe_structure"], "(((...)))");
        assert!(
            (analysis["result"]["data"]["mfe_energy_kcal_mol"]
                .as_f64()
                .unwrap()
                + 0.2)
                .abs()
                < 1.0e-12
        );
        assert!(analysis["result"]["data"]["model"]["parameter_set"]
            .as_str()
            .unwrap()
            .contains("DNA"));

        let conditional = call(json!({
            "schema_version": 1,
            "operation": "conditional-density2",
            "input": {"sequence": "GGGTTTCCC", "structure": "........."},
            "model": model
        }));
        assert_eq!(conditional["ok"], true, "{conditional:#}");
        assert_eq!(conditional["result"]["data"]["mfe_structure"], "[[[...]]]");
        assert!(
            (conditional["result"]["data"]["mfe_energy_kcal_mol"]
                .as_f64()
                .unwrap()
                + 0.2)
                .abs()
                < 1.0e-12
        );
        assert!(conditional["result"]["data"]["model"]
            .as_str()
            .unwrap()
            .contains("DNA"));
    }

    #[test]
    fn normalized_custom_parameter_overlay_reaches_all_energy_consumers() {
        let custom = json!({
            "id": ribon_core::CUSTOM_MODEL_ID,
            "parameter_base": "rna",
            "parameter_overrides": {
                "schema_version": 1,
                "name": "terminal-au-test",
                "fingerprint_sha256": "abababababababababababababababababababababababababababababababab",
                "misc": [410, 360, 500, 370]
            },
            "temperature_celsius": 37.0,
            "min_loop": 3,
            "dangles": 0,
            "salt_molar": 1.021,
            "mea_gamma": 1.0
        });
        let evaluated = call(json!({
            "schema_version": 1,
            "operation": "evaluate",
            "input": {"sequence": "AAAAU", "structure": "(...)"},
            "model": custom.clone()
        }));
        let ordinary = call(json!({
            "schema_version": 1,
            "operation": "evaluate",
            "input": {"sequence": "AAAAU", "structure": "(...)"},
            "model": {"dangles": 0}
        }));
        assert_eq!(evaluated["ok"], true, "{evaluated:#}");
        let delta = evaluated["result"]["data"]["total_kcal_mol"]
            .as_f64()
            .unwrap()
            - ordinary["result"]["data"]["total_kcal_mol"]
                .as_f64()
                .unwrap();
        assert!((delta - 9.0).abs() < 1.0e-12);

        let conditional = call(json!({
            "schema_version": 1,
            "operation": "conditional-density2",
            "input": {"sequence": "AAAAU", "structure": "....."},
            "model": custom
        }));
        assert_eq!(conditional["ok"], true, "{conditional:#}");
        assert!(conditional["result"]["data"]["model"]
            .as_str()
            .unwrap()
            .contains("custom"));

        let mut duplex_profile = custom;
        duplex_profile["parameter_overrides"]["stack_37"] = json!(vec![0; 49]);
        duplex_profile["parameter_overrides"]["stack_dh"] = json!(vec![0; 49]);
        let custom_duplex = call(json!({
            "schema_version": 1,
            "operation": "duplex",
            "input": {"sequence_a": "GGG", "sequence_b": "CCC"},
            "model": duplex_profile
        }));
        let ordinary_duplex = call(json!({
            "schema_version": 1,
            "operation": "duplex",
            "input": {"sequence_a": "GGG", "sequence_b": "CCC"},
            "model": {"dangles": 0}
        }));
        assert_eq!(custom_duplex["ok"], true, "{custom_duplex:#}");
        assert_ne!(
            custom_duplex["result"]["data"]["mfe_energy_kcal_mol"],
            ordinary_duplex["result"]["data"]["mfe_energy_kcal_mol"]
        );
    }

    #[test]
    fn all_parameter_families_have_no_silent_general_operation_fallbacks() {
        let dna = json!({
            "id": ribon_core::DNA_MODEL_ID,
            "parameter_base": "dna",
            "temperature_celsius": 37.0,
            "min_loop": 3,
            "dangles": 0,
            "salt_molar": 1.021,
            "mea_gamma": 1.0
        });
        let custom = json!({
            "id": ribon_core::CUSTOM_MODEL_ID,
            "parameter_base": "rna",
            "parameter_overrides": {
                "schema_version": 1,
                "name": "operation-matrix",
                "fingerprint_sha256": "efefefefefefefefefefefefefefefefefefefefefefefefefefefefefefefef",
                "misc": [410, 360, 500, 370]
            },
            "temperature_celsius": 37.0,
            "min_loop": 3,
            "dangles": 0,
            "salt_molar": 1.021,
            "mea_gamma": 1.0
        });
        let requests = [
            ("analyze", json!({"sequence": "GGGAAACCC"}), json!({})),
            ("fold", json!({"sequence": "GGGAAACCC"}), json!({})),
            (
                "evaluate",
                json!({"sequence": "GGGAAACCC", "structure": "(((...)))"}),
                json!({}),
            ),
            (
                "sample",
                json!({"sequence": "GGGAAACCC"}),
                json!({"count": 2, "seed": 7}),
            ),
            (
                "accessibility",
                json!({"sequence": "GGGAAACCC"}),
                json!({"windows": [{"from": 1, "to": 1}]}),
            ),
            (
                "suboptimal",
                json!({"sequence": "GGGAAACCC"}),
                json!({"energy_band_kcal_mol": 2.0, "limit": 3}),
            ),
            (
                "duplex",
                json!({"sequence_a": "GGG", "sequence_b": "CCC"}),
                json!({}),
            ),
            (
                "cofold",
                json!({"sequence_a": "GGG", "sequence_b": "CCC"}),
                json!({}),
            ),
            (
                "local",
                json!({"sequence": "GGGAAACCC"}),
                json!({"window_size": 9, "max_pair_span": 8, "max_unpaired": 3}),
            ),
            ("circular", json!({"sequence": "GGGAAACCC"}), json!({})),
            ("pseudoknot", json!({"sequence": "GGGAAACCC"}), json!({})),
            (
                "evaluate-pseudoknot",
                json!({"sequence": "GGGAAACCC", "structure": "(((...)))"}),
                json!({}),
            ),
            (
                "conditional-density2",
                json!({"sequence": "GGGAAACCC", "structure": "........."}),
                json!({}),
            ),
            (
                "conditional-density2-oracle",
                json!({"sequence": "GGGAAACCC", "structure": "........."}),
                json!({}),
            ),
            (
                "conditional-density2-sample",
                json!({"sequence": "GGGAAACCC", "structure": "........."}),
                json!({"count": 2, "seed": 11}),
            ),
            (
                "conditional-density2-suboptimal",
                json!({"sequence": "GGGAAACCC", "structure": "........."}),
                json!({"energy_band_kcal_mol": 2.0, "limit": 3}),
            ),
            (
                "evaluate-conditional-density2",
                json!({
                    "sequence": "GGGAAACCC",
                    "seed_structure": ".........",
                    "added_structure": "[[[...]]]"
                }),
                json!({}),
            ),
            (
                "evaluate-conditional-density2-oracle",
                json!({
                    "sequence": "GGGAAACCC",
                    "seed_structure": ".........",
                    "added_structure": "[[[...]]]"
                }),
                json!({}),
            ),
            (
                "comparative",
                json!({"alignment": ["GGGAAACCC", "GGGAAACCC"]}),
                json!({}),
            ),
            (
                "landscape",
                json!({
                    "sequence": "GGGAAACCC",
                    "start_structure": ".........",
                    "target_structure": "(((...)))"
                }),
                json!({}),
            ),
            (
                "inverse-design",
                json!({"target_structure": ".........", "template": "GGGAAACCC"}),
                json!({"return-count": 1}),
            ),
            (
                "ligand",
                json!({
                    "sequence": "GGGAAACCC",
                    "motifs": [{
                        "id": "unstructured",
                        "start": 1,
                        "sequence": "GGGAAACCC",
                        "structure": ".........",
                        "standard-binding-energy-kcal-mol": -1.0,
                        "concentration-molar": 0.001
                    }]
                }),
                json!({}),
            ),
        ];
        for (name, model) in [("dna", dna.clone()), ("custom", custom.clone())] {
            for (operation, input, options) in &requests {
                let response = call(json!({
                    "schema_version": 1,
                    "operation": operation,
                    "input": input,
                    "model": model,
                    "options": options
                }));
                assert_eq!(
                    response["ok"], true,
                    "{name} {operation} unexpectedly failed: {response:#}"
                );
            }
        }

        for operation in ["modified", "gquad"] {
            let response = call(json!({
                "schema_version": 1,
                "operation": operation,
                "input": {"sequence": "GGGAAACCC"},
                "model": dna
            }));
            assert_eq!(response["ok"], false, "{operation}: {response:#}");
            assert!(response["error"]["message"]
                .as_str()
                .unwrap()
                .contains("RNA-chemistry-specific"));
        }
        for (operation, sequence) in [("modified", "ACGAAACGU"), ("gquad", "GGAGGAGGAGG")] {
            let response = call(json!({
                "schema_version": 1,
                "operation": operation,
                "input": {"sequence": sequence},
                "model": custom
            }));
            assert_eq!(response["ok"], true, "{operation}: {response:#}");
            assert!(response["result"]["data"]["model"]
                .as_str()
                .unwrap()
                .contains("custom"));
        }
    }

    #[test]
    fn conditional_density2_result_flows_to_layout() {
        let ensemble = call(json!({
            "schema_version": 1,
            "operation": "conditional-density2",
            "input": {
                "sequence": "GCGAAACGCU",
                "structure": "(........)"
            },
            "model": {"min_loop": 0, "dangles": 0}
        }));
        assert_eq!(ensemble["ok"], true, "{ensemble:#}");
        assert_eq!(ensemble["result"]["data"]["state_space_complete"], true);
        let layout = call(json!({
            "schema_version": 1,
            "operation": "layout",
            "input": {
                "sequence": ensemble["result"]["data"]["sequence"].clone(),
                "structure": ensemble["result"]["data"]["mfe_structure"].clone()
            }
        }));
        assert_eq!(layout["ok"], true, "{layout:#}");
        assert_eq!(
            layout["result"]["data"]["points"].as_array().unwrap().len(),
            10
        );
    }

    #[test]
    fn conditional_density2_dispatches_odd_dangles_to_the_exact_engine() {
        for dangles in [1, 3] {
            let ensemble = call(json!({
                "schema_version": 1,
                "operation": "conditional-density2",
                "input": {
                    "sequence": "GCGAAACGCU",
                    "structure": "(........)"
                },
                "model": {"min_loop": 0, "dangles": dangles}
            }));
            assert_eq!(ensemble["ok"], true, "{ensemble:#}");
            assert_eq!(ensemble["result"]["data"]["time_complexity"], "exponential");
            assert_eq!(ensemble["result"]["data"]["state_space_complete"], true);
            assert!(ensemble["result"]["data"]["algorithm"]
                .as_str()
                .unwrap()
                .contains("exhaustive"));
            let evaluated = call(json!({
                "schema_version": 1,
                "operation": "evaluate-conditional-density2",
                "input": {
                    "sequence": "GAGAACAACU",
                    "seed_structure": "(....)....",
                    "added_structure": "..(.....)."
                },
                "model": {"min_loop": 0, "dangles": dangles}
            }));
            assert_eq!(evaluated["ok"], true, "{evaluated:#}");
            assert_eq!(evaluated["result"]["data"]["derivation_unique"], true);
            assert!(evaluated["result"]["data"]["energy_kcal_mol"]
                .as_f64()
                .unwrap()
                .is_finite());
        }
    }

    #[test]
    fn conditional_sampling_kbest_constraints_and_topology_share_the_stable_protocol() {
        let base = json!({
            "schema_version": 1,
            "input": {
                "sequence": "GCGCGCGCGCGC",
                "structure": "(......)...."
            },
            "model": {"min_loop": 3, "dangles": 0},
            "constraints": {
                "force-unpaired": [2],
                "probing": {
                    "kind": "shape",
                    "method": "zarringhalam",
                    "reactivities": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, null]
                }
            }
        });
        let mut sample_request = base.clone();
        sample_request["operation"] = json!("conditional-density2-sample");
        sample_request["options"] = json!({"count": 16, "seed": 9, "unique": false});
        let sampled = call(sample_request);
        assert_eq!(sampled["ok"], true, "{sampled:#}");
        assert_eq!(sampled["result"]["data"]["returned"], 16);
        assert_eq!(
            sampled["result"]["data"]["constraints"]["probing_method"],
            "zarringhalam"
        );

        let mut suboptimal_request = base;
        suboptimal_request["operation"] = json!("conditional-density2-suboptimal");
        suboptimal_request["options"] = json!({"energy_band_kcal_mol": 20.0, "limit": 8});
        let suboptimal = call(suboptimal_request);
        assert_eq!(suboptimal["ok"], true, "{suboptimal:#}");
        assert!(suboptimal["result"]["data"]["structures"]
            .as_array()
            .is_some_and(|structures| !structures.is_empty()));

        let topology = call(json!({
            "schema_version": 1,
            "operation": "fatgraph-topology",
            "input": {"sequence": "GCGCGCG", "structure": "(.[.).]"}
        }));
        assert_eq!(topology["ok"], true, "{topology:#}");
        assert_eq!(topology["result"]["data"]["genus"], 1);
    }

    #[test]
    fn unknown_fields_and_versions_are_rejected() {
        let response = call(json!({
            "schema_version": 99,
            "operation": "capabilities",
            "input": {}
        }));
        assert_eq!(response["ok"], false);
        assert_eq!(response["error"]["code"], "invalid_request");
    }
}
