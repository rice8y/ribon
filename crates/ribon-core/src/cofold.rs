//! Two-strand association analysis and equilibrium concentrations.
//!
//! The thermodynamic cycle combines the exact single-strand Turner ensemble
//! implemented by this crate with the connected intermolecular duplex
//! ensemble.  The two macrostates are deliberately explicit: an unbound state
//! containing independently folded monomers and a bound state containing one
//! connected intermolecular interaction.  This makes standard-state handling,
//! pair-probability mixing, and concentration calculations reproducible.

use crate::cofold_standard::{cofold_standard, StandardCofold};
use crate::energy::EnergyModel;
use crate::partition::PairProbability;
use crate::structure::{normalize_sequence, RnaError};
use crate::{analyze_with_model, AnalysisResult, ConstraintConfig, NucleicAcid};
use serde::Serialize;

const GAS_CONSTANT_KCAL: f64 = 0.001_987_17;
const STANDARD_CONCENTRATION_MOLAR: f64 = 1.0;

#[derive(Clone, Debug, Serialize)]
pub struct CofoldMonomer {
    pub sequence: String,
    pub mfe_structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub ensemble_free_energy_kcal_mol: f64,
    pub log_partition_function: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct CofoldComplex {
    pub sequence: String,
    pub structure: String,
    pub mfe_energy_kcal_mol: f64,
    pub unbound_mfe_energy_kcal_mol: f64,
    pub bound_mfe_energy_kcal_mol: Option<f64>,
    pub ensemble_free_energy_kcal_mol: f64,
    pub association_free_energy_kcal_mol: Option<f64>,
    pub log_partition_function: f64,
    pub bound_probability: f64,
    pub pair_probabilities: Vec<PairProbability>,
    pub unpaired_probabilities: Vec<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct EquilibriumConcentrations {
    pub total_a_molar: f64,
    pub total_b_molar: f64,
    pub free_a_molar: f64,
    pub free_b_molar: f64,
    pub aa_molar: f64,
    pub ab_molar: f64,
    pub bb_molar: f64,
    pub log_k_aa_per_molar: Option<f64>,
    pub log_k_ab_per_molar: Option<f64>,
    pub log_k_bb_per_molar: Option<f64>,
    pub mass_balance_a_error_molar: f64,
    pub mass_balance_b_error_molar: f64,
    pub iterations: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct CofoldResult {
    pub sequence_a: String,
    pub sequence_b: String,
    pub monomer_a: CofoldMonomer,
    pub monomer_b: CofoldMonomer,
    pub complex_ab: CofoldComplex,
    pub equilibrium: Option<EquilibriumConcentrations>,
    pub temperature_celsius: f64,
    pub salt_molar: f64,
    pub standard_concentration_molar: f64,
    pub model: &'static str,
}

/// Analyze A, B, and their associated complex.  Optional total concentrations
/// activate the five-species AA/A/AB/B/BB equilibrium calculation.
#[allow(clippy::too_many_arguments)]
pub fn cofold(
    sequence_a: &str,
    sequence_b: &str,
    temperature_celsius: f64,
    min_loop: usize,
    gamma: f64,
    dangles: u8,
    salt_molar: f64,
    concentrations_molar: Option<(f64, f64)>,
) -> Result<CofoldResult, RnaError> {
    let model = EnergyModel::with_dangles_and_salt(temperature_celsius, dangles, salt_molar)?;
    cofold_with_model(
        sequence_a,
        sequence_b,
        min_loop,
        gamma,
        model,
        concentrations_molar,
    )
}

/// Analyze two strands using an explicitly selected RNA or DNA parameter
/// family across monomer and complex states.
pub fn cofold_with_model(
    sequence_a: &str,
    sequence_b: &str,
    min_loop: usize,
    gamma: f64,
    model: EnergyModel,
    concentrations_molar: Option<(f64, f64)>,
) -> Result<CofoldResult, RnaError> {
    let a = normalize_sequence(sequence_a)?;
    let b = normalize_sequence(sequence_b)?;
    let temperature_celsius = model.temperature_celsius();
    let dangles = model.dangles();
    let salt_molar = model.salt_molar();
    let analysis_a = analyze_with_model(
        a.clone(),
        min_loop,
        gamma,
        model.clone(),
        &ConstraintConfig::default(),
    )?;
    let analysis_b = analyze_with_model(
        b.clone(),
        min_loop,
        gamma,
        model.clone(),
        &ConstraintConfig::default(),
    )?;
    let full_ab = cofold_standard(&a, &b, temperature_celsius, min_loop, &model)?;
    let complex_ab = full_complex(&analysis_a, &analysis_b, &full_ab);

    let equilibrium = if let Some((total_a, total_b)) = concentrations_molar {
        validate_concentration(total_a, "total concentration A")?;
        validate_concentration(total_b, "total concentration B")?;
        let full_aa = cofold_standard(&a, &a, temperature_celsius, min_loop, &model)?;
        let full_bb = cofold_standard(&b, &b, temperature_celsius, min_loop, &model)?;
        let (log_kaa, log_kab, log_kbb) = (
            full_association_log_k(&analysis_a, &analysis_a, &full_aa, true),
            full_association_log_k(&analysis_a, &analysis_b, &full_ab, false),
            full_association_log_k(&analysis_b, &analysis_b, &full_bb, true),
        );
        Some(solve_equilibrium(
            total_a, total_b, log_kaa, log_kab, log_kbb,
        )?)
    } else {
        None
    };

    Ok(CofoldResult {
        sequence_a: a,
        sequence_b: b,
        monomer_a: monomer(&analysis_a),
        monomer_b: monomer(&analysis_b),
        complex_ab,
        equilibrium,
        temperature_celsius,
        salt_molar,
        standard_concentration_molar: STANDARD_CONCENTRATION_MOLAR,
        model: if model.parameter_profile_name().is_some() {
            if dangles % 2 == 0 {
                "full polynomial cut-point custom nearest-neighbor cofold ensemble"
            } else {
                "full exhaustive cut-point custom fixed-structure cofold ensemble"
            }
        } else {
            match (model.nucleic_acid(), dangles % 2) {
                (NucleicAcid::Rna, 0) => "full polynomial cut-point RNAstructure 6.6 RNA cofold ensemble with intra- and intermolecular branches",
                (NucleicAcid::Rna, _) => "full exhaustive cut-point fixed-structure RNAstructure 6.6 RNA cofold ensemble",
                (NucleicAcid::Dna, 0) => "full polynomial cut-point RNAstructure 6.6 DNA cofold ensemble with intra- and intermolecular branches",
                (NucleicAcid::Dna, _) => "full exhaustive cut-point fixed-structure RNAstructure 6.6 DNA cofold ensemble",
            }
        },
    })
}

fn monomer(result: &AnalysisResult) -> CofoldMonomer {
    CofoldMonomer {
        sequence: result.sequence.clone(),
        mfe_structure: result.mfe_structure.clone(),
        mfe_energy_kcal_mol: result.mfe_energy_kcal_mol,
        ensemble_free_energy_kcal_mol: result.ensemble_free_energy_kcal_mol,
        log_partition_function: result.log_partition_function,
    }
}

fn full_complex(a: &AnalysisResult, b: &AnalysisResult, full: &StandardCofold) -> CofoldComplex {
    let log_unbound = full.log_unbound_partition_function;
    let p_bound = if full.log_bound_partition_function.is_finite() {
        (full.log_bound_partition_function - full.log_partition_function)
            .exp()
            .clamp(0.0, 1.0)
    } else {
        0.0
    };
    let unbound_mfe = a.mfe_energy_kcal_mol + b.mfe_energy_kcal_mol;
    let rt = GAS_CONSTANT_KCAL * (a.temperature_celsius + 273.15);
    CofoldComplex {
        sequence: format!("{}&{}", a.sequence, b.sequence),
        structure: full.structure.clone(),
        mfe_energy_kcal_mol: full.mfe_energy_kcal_mol,
        unbound_mfe_energy_kcal_mol: unbound_mfe,
        bound_mfe_energy_kcal_mol: full.bound_mfe_energy_kcal_mol,
        ensemble_free_energy_kcal_mol: -rt * full.log_partition_function,
        association_free_energy_kcal_mol: full
            .log_bound_partition_function
            .is_finite()
            .then_some(-rt * (full.log_bound_partition_function - log_unbound)),
        log_partition_function: full.log_partition_function,
        bound_probability: p_bound,
        pair_probabilities: full.pair_probabilities.clone(),
        unpaired_probabilities: full.unpaired_probabilities.clone(),
    }
}

fn full_association_log_k(
    a: &AnalysisResult,
    b: &AnalysisResult,
    full: &StandardCofold,
    homodimer: bool,
) -> Option<f64> {
    full.log_bound_partition_function.is_finite().then(|| {
        debug_assert!(
            (full.log_unbound_partition_function
                - a.log_partition_function
                - b.log_partition_function)
                .abs()
                < 1.0e-8
        );
        let symmetry = if homodimer { 2.0f64.ln() } else { 0.0 };
        full.log_bound_partition_function
            - full.log_unbound_partition_function
            - symmetry
            - STANDARD_CONCENTRATION_MOLAR.ln()
    })
}

fn solve_equilibrium(
    total_a: f64,
    total_b: f64,
    log_k_aa: Option<f64>,
    log_k_ab: Option<f64>,
    log_k_bb: Option<f64>,
) -> Result<EquilibriumConcentrations, RnaError> {
    let kaa = equilibrium_constant(log_k_aa)?;
    let kab = equilibrium_constant(log_k_ab)?;
    let kbb = equilibrium_constant(log_k_bb)?;
    let mut a = total_a / (1.0 + (2.0 * kaa * total_a + kab * total_b).sqrt());
    let mut b = total_b / (1.0 + (2.0 * kbb * total_b + kab * total_a).sqrt());
    let scale = total_a.max(total_b).max(1.0e-30);
    let tolerance = scale * 2.0e-13;
    let mut iterations = 0usize;
    loop {
        iterations = iterations.saturating_add(1);
        let f_a = a + 2.0 * kaa * a * a + kab * a * b - total_a;
        let f_b = b + 2.0 * kbb * b * b + kab * a * b - total_b;
        let residual = f_a.abs().max(f_b.abs());
        if residual <= tolerance {
            break;
        }
        let j_aa = 1.0 + 4.0 * kaa * a + kab * b;
        let j_ab = kab * a;
        let j_ba = kab * b;
        let j_bb = 1.0 + 4.0 * kbb * b + kab * a;
        let determinant = j_aa * j_bb - j_ab * j_ba;
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(RnaError::Numerical(
                "dimer equilibrium Jacobian became singular".into(),
            ));
        }
        let delta_a = (-f_a * j_bb + j_ab * f_b) / determinant;
        let delta_b = (-j_aa * f_b + j_ba * f_a) / determinant;
        let mut damping = 1.0;
        let (next_a, next_b) = loop {
            let candidate_a = a + damping * delta_a;
            let candidate_b = b + damping * delta_b;
            let candidate_f_a = candidate_a
                + 2.0 * kaa * candidate_a * candidate_a
                + kab * candidate_a * candidate_b
                - total_a;
            let candidate_f_b = candidate_b
                + 2.0 * kbb * candidate_b * candidate_b
                + kab * candidate_a * candidate_b
                - total_b;
            if (0.0..=total_a).contains(&candidate_a)
                && (0.0..=total_b).contains(&candidate_b)
                && candidate_f_a.abs().max(candidate_f_b.abs()) < residual
            {
                break (candidate_a, candidate_b);
            }
            damping *= 0.5;
            if damping < 1.0e-12 {
                return Err(RnaError::Numerical(
                    "dimer equilibrium Newton line search stagnated".into(),
                ));
            }
        };
        a = next_a;
        b = next_b;
    }
    let aa = kaa * a * a;
    let ab = kab * a * b;
    let bb = kbb * b * b;
    let error_a = a + 2.0 * aa + ab - total_a;
    let error_b = b + 2.0 * bb + ab - total_b;
    if error_a.abs().max(error_b.abs()) > tolerance * 10.0 {
        return Err(RnaError::Numerical(format!(
            "dimer equilibrium did not converge (mass errors {error_a:e}, {error_b:e})"
        )));
    }
    Ok(EquilibriumConcentrations {
        total_a_molar: total_a,
        total_b_molar: total_b,
        free_a_molar: a,
        free_b_molar: b,
        aa_molar: aa,
        ab_molar: ab,
        bb_molar: bb,
        log_k_aa_per_molar: log_k_aa,
        log_k_ab_per_molar: log_k_ab,
        log_k_bb_per_molar: log_k_bb,
        mass_balance_a_error_molar: error_a,
        mass_balance_b_error_molar: error_b,
        iterations,
    })
}

fn equilibrium_constant(log_k: Option<f64>) -> Result<f64, RnaError> {
    let Some(log_k) = log_k else {
        return Ok(0.0);
    };
    let value = log_k.exp();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RnaError::Numerical(format!(
            "equilibrium constant exp({log_k}) exceeds the f64 output domain"
        )))
    }
}

fn validate_concentration(value: f64, name: &str) -> Result<(), RnaError> {
    if !value.is_finite() || value < 0.0 {
        return Err(RnaError::InvalidOption(format!(
            "{name} must be finite and non-negative"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cofold_probabilities_and_mass_balance_are_normalized() {
        let result = cofold(
            "GGGAAACCC",
            "GGGAAACCC",
            37.0,
            3,
            1.0,
            2,
            1.021,
            Some((1.0e-6, 2.0e-6)),
        )
        .unwrap();
        assert!((0.0..=1.0).contains(&result.complex_ab.bound_probability));
        assert!(result
            .complex_ab
            .pair_probabilities
            .iter()
            .all(|entry| (0.0..=1.0).contains(&entry.probability)));
        let equilibrium = result.equilibrium.unwrap();
        assert!(equilibrium.mass_balance_a_error_molar.abs() < 1.0e-16);
        assert!(equilibrium.mass_balance_b_error_molar.abs() < 1.0e-16);
    }

    #[test]
    fn noninteracting_pair_reduces_to_independent_monomers() {
        let result = cofold("AAAA", "AAAA", 37.0, 3, 1.0, 0, 1.021, None).unwrap();
        assert_eq!(result.complex_ab.bound_probability, 0.0);
        assert_eq!(result.complex_ab.structure, "....&....");
        assert!(
            (result.complex_ab.log_partition_function
                - result.monomer_a.log_partition_function
                - result.monomer_b.log_partition_function)
                .abs()
                < 1.0e-12
        );
    }
}
