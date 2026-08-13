use ribon_core::{
    analyze_with_profile, circular_fold, cofold, comparative_fold, fold_gquad, fold_modified,
    local_fold, parameter_manifest, parse_structure, predict_pseudoknot, ComparativeOptions,
    ModifiedBase, ParameterProfile, PseudoknotOptions,
};

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'U',
            b'C' => 'G',
            b'G' => 'C',
            b'U' => 'A',
            _ => 'N',
        })
        .collect()
}

fn matrix(sequence: &str) -> serde_json::Value {
    let fragment = &sequence[..sequence.len().min(48)];
    let short = &fragment[..fragment.len().min(24)];
    let interaction_a = &fragment[..fragment.len().min(12)];
    let interaction_b = reverse_complement(interaction_a);
    let complex = cofold(
        interaction_a,
        &interaction_b,
        37.0,
        3,
        1.0,
        2,
        1.021,
        Some((1.0e-6, 2.0e-6)),
    )
    .unwrap();
    let modification = ModifiedBase {
        position: 1,
        symbol: "mod".into(),
        canonical_base: fragment.as_bytes()[0] as char,
        kind: None,
        paired_energy_kcal_mol: 0.1,
        unpaired_energy_kcal_mol: 0.0,
        stack_energy_kcal_mol: -0.1,
    };
    let modified = fold_modified(fragment, &[modification], 37.0, 3, 1.0, 2, 1.021).unwrap();
    let gquad = fold_gquad(fragment, 37.0, 3, 1.0, 2, 1.021).unwrap();
    let pseudoknot = predict_pseudoknot(
        fragment,
        37.0,
        3,
        2,
        1.021,
        &PseudoknotOptions {
            threshold: 0.02,
            ..PseudoknotOptions::default()
        },
    )
    .unwrap();
    let comparative = comparative_fold(
        &[
            fragment.to_string(),
            fragment.to_string(),
            fragment.to_string(),
        ],
        37.0,
        3,
        1.0,
        2,
        1.021,
        &ComparativeOptions::default(),
    )
    .unwrap();
    let profiled = analyze_with_profile(fragment, &ParameterProfile::default()).unwrap();
    let equilibrium = complex.equilibrium.unwrap();
    let candidate_probability_sum = gquad
        .candidates
        .iter()
        .map(|candidate| candidate.probability)
        .sum::<f64>();
    serde_json::json!({
        "fragment_length": fragment.len(),
        "short_length": short.len(),
        "cofold_bound_probability": complex.complex_ab.bound_probability,
        "cofold_mass_balance_a_error_molar": equilibrium.mass_balance_a_error_molar,
        "cofold_mass_balance_b_error_molar": equilibrium.mass_balance_b_error_molar,
        "modified_symbol": modified.display_symbols[0],
        "modified_mfe_energy_kcal_mol": modified.analysis.mfe_energy_kcal_mol,
        "gquad_probability": gquad.gquad_probability,
        "gquad_candidate_probability_sum": candidate_probability_sum,
        "gquad_candidate_count": gquad.candidates.len(),
        "pseudoknot_pair_count": pseudoknot.pair_count,
        "pseudoknot_roundtrip_pair_count": parse_structure(fragment, &pseudoknot.structure).unwrap().pairs.len(),
        "pseudoknot_crossing_count": pseudoknot.crossing_count,
        "comparative_sequence_count": comparative.sequence_count,
        "comparative_consensus_matches_input": comparative.consensus_sequence == fragment,
        "profile_mfe_matches_default": profiled.mfe_structure == comparative.analysis.mfe_structure,
        "parameter_bundle_sha256": parameter_manifest().parameter_bundle_sha256,
    })
}

fn main() {
    let arguments = std::env::args().collect::<Vec<_>>();
    if arguments.len() < 3 {
        eprintln!(
            "usage: extended-json circular|gquad|local|matrix SEQUENCE [WINDOW SPAN UNPAIRED]"
        );
        std::process::exit(2);
    }
    let value = match arguments[1].as_str() {
        "circular" => {
            serde_json::to_value(circular_fold(&arguments[2], 37.0, 3, 1.0, 0, 1.021).unwrap())
                .unwrap()
        }
        "gquad" => {
            let temperature = arguments
                .get(3)
                .map(|value| value.parse().unwrap())
                .unwrap_or(37.0);
            serde_json::to_value(fold_gquad(&arguments[2], temperature, 3, 1.0, 0, 1.021).unwrap())
                .unwrap()
        }
        "local" if arguments.len() == 6 => serde_json::to_value(
            local_fold(
                &arguments[2],
                37.0,
                3,
                0,
                1.021,
                arguments[3].parse().unwrap(),
                arguments[4].parse().unwrap(),
                arguments[5].parse().unwrap(),
            )
            .unwrap(),
        )
        .unwrap(),
        "matrix" => matrix(&arguments[2]),
        _ => {
            eprintln!("invalid mode or arguments");
            std::process::exit(2);
        }
    };
    println!("{}", serde_json::to_string(&value).unwrap());
}
