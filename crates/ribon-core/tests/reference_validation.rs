use ribon_core::{
    analyze_with_options, duplex, fold_sequence_with_options, layout_structure, parse_structure,
    ConstraintConfig, LayoutKind,
};
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize)]
struct LayoutCorpus {
    source: String,
    commit: String,
    cases: Vec<LayoutCase>,
}

#[derive(Deserialize)]
struct LayoutCase {
    sequence: String,
    structure: String,
    points: Vec<[f64; 2]>,
}

#[derive(Deserialize)]
struct DuplexCorpus {
    source: String,
    commit: String,
    cases: Vec<DuplexCase>,
}

#[derive(Deserialize)]
struct DuplexCase {
    sequence_a: String,
    sequence_b: String,
    salt_molar: f64,
    full_structure: String,
    mfe_energy: f64,
}

#[derive(Deserialize)]
struct RnaStructureReport {
    schema_version: u32,
    reference: String,
    reference_archive_sha256: String,
    cases: Vec<RnaStructureCase>,
}

#[derive(Deserialize)]
struct RnaStructureCase {
    accession: String,
    rnastructure: RnaStructureValues,
    rnastructure_significant_pair_probabilities: Vec<ReferencePairProbability>,
}

#[derive(Deserialize)]
struct RnaStructureValues {
    sequence: String,
    mfe_structure: String,
    mfe_energy_kcal_mol: f64,
    ensemble_free_energy_kcal_mol: f64,
    centroid_structure: String,
    mea_structure: String,
}

#[derive(Deserialize)]
struct ReferencePairProbability {
    i: usize,
    j: usize,
    probability: f64,
}

fn normalized_distance_signature(points: &[[f64; 2]]) -> Vec<f64> {
    let backbone: f64 = points
        .windows(2)
        .map(|edge| {
            let dx = edge[1][0] - edge[0][0];
            let dy = edge[1][1] - edge[0][1];
            (dx * dx + dy * dy).sqrt()
        })
        .sum::<f64>()
        / (points.len() - 1) as f64;
    let mut signature = Vec::new();
    for i in 0..points.len() {
        for j in (i + 1)..points.len() {
            let dx = points[j][0] - points[i][0];
            let dy = points[j][1] - points[i][1];
            signature.push((dx * dx + dy * dy).sqrt() / backbone);
        }
    }
    signature
}

#[test]
fn naview_geometry_matches_independent_reference_coordinates() {
    let corpus: LayoutCorpus =
        serde_json::from_str(include_str!("fixtures/vienna_naview_2_7_2.json")).unwrap();
    assert_eq!(corpus.source, "ViennaRNA 2.7.2, vrna_plot_coords_naview_pt");
    assert_eq!(corpus.commit, "1ffec79f5e258896160f7362ced8263450f371dc");
    assert_eq!(corpus.cases.len(), 13);

    for case in corpus.cases {
        let layout = layout_structure(&case.sequence, &case.structure, LayoutKind::Naview).unwrap();
        let ours: Vec<[f64; 2]> = layout
            .points
            .iter()
            .map(|point| [point.x * layout.aspect_ratio, point.y])
            .collect();
        let ours = normalized_distance_signature(&ours);
        let reference = normalized_distance_signature(&case.points);
        let rms = (ours
            .iter()
            .zip(reference.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            / ours.len() as f64)
            .sqrt();
        assert!(rms < 1.0e-6, "NAView distance-signature RMS={rms}");
    }
}

#[test]
fn simple_layout_matches_rnaplot_reference_coordinates() {
    let corpus: LayoutCorpus =
        serde_json::from_str(include_str!("fixtures/vienna_simple_2_7_2.json")).unwrap();
    for case in corpus.cases {
        let layout = layout_structure(&case.sequence, &case.structure, LayoutKind::Simple).unwrap();
        let ours: Vec<[f64; 2]> = layout
            .points
            .iter()
            .map(|point| [point.x * layout.aspect_ratio, point.y])
            .collect();
        let ours = normalized_distance_signature(&ours);
        let reference = normalized_distance_signature(&case.points);
        let rms = (ours
            .iter()
            .zip(reference.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            / ours.len() as f64)
            .sqrt();
        assert!(rms < 2.0e-6, "simple-layout distance-signature RMS={rms}");
    }
}

fn pair_set(sequence: &str, structure: &str) -> HashSet<(usize, usize)> {
    parse_structure(sequence, structure)
        .unwrap()
        .pairs
        .into_iter()
        .map(|pair| (pair.i, pair.j))
        .collect()
}

fn structure_scores(sequence: &str, ours: &str, reference: &str) -> (f64, f64) {
    let ours = pair_set(sequence, ours);
    let reference = pair_set(sequence, reference);
    let hits = ours.intersection(&reference).count();
    let sensitivity = if reference.is_empty() {
        1.0
    } else {
        hits as f64 / reference.len() as f64
    };
    let precision = if ours.is_empty() {
        1.0
    } else {
        hits as f64 / ours.len() as f64
    };
    (sensitivity, precision)
}

#[test]
fn twenty_four_rfam_cases_track_rnastructure_6_6_values() {
    let report: RnaStructureReport = serde_json::from_str(include_str!(
        "../../../tests/reports/rnastructure_6_6_validation.json"
    ))
    .unwrap();
    assert_eq!(report.schema_version, 1);
    assert!(report.reference.contains("RNAstructure 6.6"));
    assert_eq!(
        report.reference_archive_sha256,
        "8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8"
    );
    assert_eq!(report.cases.len(), 24);

    let mut mfe_energy_error = 0.0;
    let mut ensemble_energy_error = 0.0;
    let mut maximum_mfe_energy_error = 0.0f64;
    let mut mfe_sensitivity = 0.0;
    let mut mfe_precision = 0.0;
    let mut centroid_sensitivity = 0.0;
    let mut centroid_precision = 0.0;
    let mut mea_sensitivity = 0.0;
    let mut mea_precision = 0.0;
    let mut probability_error = 0.0;
    let mut probability_count = 0usize;

    for case in report.cases {
        let reference = case.rnastructure;
        let analysis = analyze_with_options(
            &reference.sequence,
            37.0,
            3,
            1.0,
            2,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap_or_else(|error| panic!("{}: {error}", case.accession));
        let mfe = fold_sequence_with_options(
            &reference.sequence,
            37.0,
            3,
            3,
            1.021,
            &ConstraintConfig::default(),
        )
        .unwrap();
        let error = (mfe.energy_kcal_mol - reference.mfe_energy_kcal_mol).abs();
        mfe_energy_error += error;
        maximum_mfe_energy_error = maximum_mfe_energy_error.max(error);
        ensemble_energy_error += (analysis.ensemble_free_energy_kcal_mol
            - reference.ensemble_free_energy_kcal_mol)
            .abs();
        let (sensitivity, precision) = structure_scores(
            &reference.sequence,
            &mfe.structure,
            &reference.mfe_structure,
        );
        mfe_sensitivity += sensitivity;
        mfe_precision += precision;
        let (sensitivity, precision) = structure_scores(
            &reference.sequence,
            &analysis.centroid_structure,
            &reference.centroid_structure,
        );
        centroid_sensitivity += sensitivity;
        centroid_precision += precision;
        let (sensitivity, precision) = structure_scores(
            &reference.sequence,
            &analysis.mea_structure,
            &reference.mea_structure,
        );
        mea_sensitivity += sensitivity;
        mea_precision += precision;
        for pair in case.rnastructure_significant_pair_probabilities {
            let observed = analysis
                .pair_probabilities
                .iter()
                .find(|entry| entry.i == pair.i && entry.j == pair.j)
                .map_or(0.0, |entry| entry.probability);
            probability_error += (observed - pair.probability).abs();
            probability_count += 1;
        }
    }
    let count = 24.0;
    assert!(mfe_energy_error / count <= 1.1);
    assert!(maximum_mfe_energy_error <= 3.0);
    assert!(ensemble_energy_error / count <= 1.6);
    assert!(mfe_sensitivity / count >= 0.79);
    assert!(mfe_precision / count >= 0.81);
    assert!(centroid_sensitivity / count >= 0.88);
    assert!(centroid_precision / count >= 0.90);
    assert!(mea_sensitivity / count >= 0.87);
    assert!(mea_precision / count >= 0.88);
    assert!(probability_count > 1_000);
    let significant_probability_mae = probability_error / probability_count as f64;
    assert!(
        significant_probability_mae <= 0.08,
        "significant pair-probability MAE={significant_probability_mae} over {probability_count} pairs"
    );
}

#[test]
fn connected_duplex_mfe_retains_the_published_reference_fixture() {
    let corpus: DuplexCorpus =
        serde_json::from_str(include_str!("fixtures/vienna_duplex_2_7_2.json")).unwrap();
    assert_eq!(corpus.source, "ViennaRNA 2.7.2 RNAduplex, dangles=0");
    assert_eq!(corpus.commit, "1ffec79f5e258896160f7362ced8263450f371dc");
    for case in corpus.cases {
        let ours = duplex(&case.sequence_a, &case.sequence_b, 37.0, case.salt_molar).unwrap();
        assert_eq!(ours.structure, case.full_structure);
        assert!((ours.mfe_energy_kcal_mol - case.mfe_energy).abs() < 1.0e-9);
        assert!(ours.log_bound_partition_function.is_finite());
        assert!((0.0..=1.0).contains(&ours.standard_state_bound_probability));
    }
}
