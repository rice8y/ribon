use ribon_core::{
    accessibility, analyze_with_options, analyze_with_profile, circular_fold, cofold,
    comparative_fold, duplex, ensemble_defect, evaluate_pseudoknot_structure, exact_landscape_path,
    fold_gquad, fold_mfe_with_constraints, fold_modified, fold_sequence_with_options,
    inverse_fold_exact, layout_structure, ligand_ensemble_exact, local_fold, parse_structure,
    partition_with_constraints, predict_pseudoknot, sample_boltzmann, suboptimal_structures,
    AccessibilityWindow, ComparativeOptions, ConstraintConfig, ConstraintModel, EnergyModel,
    InverseDesignOptions, LayoutKind, LigandMotif, ModifiedBase, ModifiedBaseKind,
    ParameterProfile, PseudoknotOptions,
};
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize)]
struct Corpus {
    source: String,
    license: String,
    cases: Vec<RealRna>,
}

#[derive(Deserialize)]
struct RealRna {
    accession: String,
    family_id: String,
    sequence_id: String,
    sequence: String,
    reference_structure: String,
    structure_source: String,
    length: usize,
}

fn corpus() -> Corpus {
    serde_json::from_str(include_str!("../../../tests/data/rfam_real_24.json")).unwrap()
}

#[test]
fn rfam_corpus_is_diverse_provenanced_and_structurally_valid() {
    let corpus = corpus();
    assert_eq!(corpus.source, "Rfam SEED alignments");
    assert_eq!(corpus.license, "CC0-1.0");
    assert_eq!(corpus.cases.len(), 24);
    assert_eq!(
        corpus
            .cases
            .iter()
            .map(|case| case.accession.as_str())
            .collect::<HashSet<_>>()
            .len(),
        24
    );
    assert!(
        corpus
            .cases
            .iter()
            .filter(|case| case.structure_source.starts_with("Published"))
            .count()
            >= 20
    );
    for case in corpus.cases {
        assert!(!case.family_id.is_empty());
        assert!(!case.sequence_id.is_empty());
        assert_eq!(case.sequence.len(), case.length, "{}", case.accession);
        assert_eq!(
            case.reference_structure.len(),
            case.length,
            "{}",
            case.accession
        );
        let parsed = parse_structure(&case.sequence, &case.reference_structure).unwrap();
        assert_eq!(parsed.length, case.length);
    }
}

#[test]
fn all_vector_layouts_are_finite_on_all_24_real_rnas() {
    for case in corpus().cases {
        for method in [
            LayoutKind::Simple,
            LayoutKind::Naview,
            LayoutKind::Circular,
            LayoutKind::Turtle,
            LayoutKind::Puzzler,
            LayoutKind::Linear,
        ] {
            let result = layout_structure(&case.sequence, &case.reference_structure, method)
                .unwrap_or_else(|error| panic!("{} {method:?}: {error}", case.accession));
            assert_eq!(result.points.len(), case.length);
            assert!(result.points.iter().all(|point| {
                point.x.is_finite()
                    && point.y.is_finite()
                    && (0.0..=1.0).contains(&point.x)
                    && (0.0..=1.0).contains(&point.y)
            }));
            if method == LayoutKind::Puzzler
                && parse_structure(&case.sequence, &case.reference_structure)
                    .unwrap()
                    .pairs
                    .iter()
                    .all(|pair| pair.level == 0)
            {
                assert_eq!(
                    straight_segment_crossings(
                        &result.points,
                        &parse_structure(&case.sequence, &case.reference_structure).unwrap(),
                    ),
                    0,
                    "{} Puzzler intersections",
                    case.accession,
                );
            }
        }
    }
}

fn straight_segment_crossings(
    points: &[ribon_core::layout::Point],
    parsed: &ribon_core::ParsedStructure,
) -> usize {
    let mut edges = (0..points.len().saturating_sub(1))
        .filter(|index| !parsed.strand_breaks.contains(&(index + 1)))
        .map(|index| (index, index + 1))
        .collect::<Vec<_>>();
    edges.extend(parsed.pairs.iter().map(|pair| (pair.i - 1, pair.j - 1)));
    let orientation = |a: usize, b: usize, c: usize| {
        (points[b].x - points[a].x) * (points[c].y - points[a].y)
            - (points[b].y - points[a].y) * (points[c].x - points[a].x)
    };
    let mut count = 0;
    for (index, &(a, b)) in edges.iter().enumerate() {
        for &(c, d) in &edges[index + 1..] {
            if [a, b].contains(&c) || [a, b].contains(&d) {
                continue;
            }
            if orientation(a, b, c) * orientation(a, b, d) < -1.0e-12
                && orientation(c, d, a) * orientation(c, d, b) < -1.0e-12
            {
                count += 1;
            }
        }
    }
    count
}

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

#[test]
fn every_analysis_feature_runs_on_all_24_real_rna_families() {
    // Cubic and k-best operations use a real 48-nt prefix from every family;
    // the separate full-length differential test exercises all 82-480 nt
    // records. This matrix detects API integration and numerical-domain bugs
    // without turning the default test suite into a multi-minute benchmark.
    for case in corpus().cases {
        let fragment = &case.sequence[..case.sequence.len().min(48)];
        let config = ConstraintConfig::default();
        let model = EnergyModel::with_dangles_and_salt(37.0, 2, 1.021).unwrap();
        let constraints = ConstraintModel::compile(fragment.len(), &config).unwrap();

        let analysis = analyze_with_options(fragment, 37.0, 3, 1.0, 2, 1.021, &config)
            .unwrap_or_else(|error| panic!("{} prediction: {error}", case.accession));
        assert!(analysis.mfe_energy_kcal_mol.is_finite());
        assert!(analysis.ensemble_free_energy_kcal_mol.is_finite());
        for structure in [
            &analysis.mfe_structure,
            &analysis.centroid_structure,
            &analysis.mea_structure,
        ] {
            assert_eq!(
                parse_structure(fragment, structure).unwrap().length,
                fragment.len()
            );
        }

        for dangles in 0..=3 {
            let fold = fold_sequence_with_options(fragment, 37.0, 3, dangles, 1.021, &config)
                .unwrap_or_else(|error| panic!("{} dangles={dangles}: {error}", case.accession));
            let evaluated = EnergyModel::with_dangles_and_salt(37.0, dangles, 1.021)
                .unwrap()
                .evaluate(fragment, &fold.structure)
                .unwrap();
            assert!(evaluated.total_kcal_mol.is_finite());
        }

        let partition = partition_with_constraints(fragment, 37.0, 3, &model, &constraints)
            .unwrap_or_else(|error| panic!("{} partition: {error}", case.accession));
        let defect = ensemble_defect(fragment, &analysis.mfe_structure, &partition).unwrap();
        assert!((0.0..=1.0).contains(&defect.normalized_ensemble_defect));

        let samples = sample_boltzmann(
            fragment,
            37.0,
            3,
            &model,
            &constraints,
            3,
            20_260_810,
            false,
        )
        .unwrap_or_else(|error| panic!("{} sampling: {error}", case.accession));
        assert_eq!(samples.samples.len(), 3);

        let alternatives = suboptimal_structures(fragment, 3, &model, &constraints, 3.0, 4)
            .unwrap_or_else(|error| panic!("{} suboptimal: {error}", case.accession));
        assert!(!alternatives.structures.is_empty());

        let opening = accessibility(
            fragment,
            37.0,
            3,
            2,
            1.021,
            &config,
            &[
                AccessibilityWindow { from: 1, to: 1 },
                AccessibilityWindow {
                    from: fragment.len() / 2,
                    to: fragment.len() / 2 + 2,
                },
            ],
        )
        .unwrap_or_else(|error| panic!("{} accessibility: {error}", case.accession));
        assert!(opening
            .windows
            .iter()
            .all(|window| (0.0..=1.0).contains(&window.probability_unpaired)));

        let interaction_a = &fragment[..12];
        let interaction_b = reverse_complement(interaction_a);
        let interaction = duplex(interaction_a, &interaction_b, 37.0, 1.021)
            .unwrap_or_else(|error| panic!("{} duplex: {error}", case.accession));
        assert!(interaction.mfe_energy_kcal_mol.is_finite());
        assert_eq!(interaction.structure.len(), 25);

        // v1+ extended-analysis matrix: every feature is exercised on a real
        // fragment from every one of the 24 independently accessioned Rfam
        // families, with invariants that are independent of visual QA.
        let local_fragment = &fragment[..fragment.len().min(24)];
        let local = local_fold(
            local_fragment,
            37.0,
            3,
            2,
            1.021,
            local_fragment.len().min(18),
            local_fragment.len().min(12).saturating_sub(1).max(1),
            1,
        )
        .unwrap_or_else(|error| panic!("{} local fold: {error}", case.accession));
        assert!(local
            .pair_probabilities
            .iter()
            .all(|entry| (0.0..=1.0).contains(&entry.probability)));

        let circular_fragment = &fragment[..fragment.len().min(24)];
        let circular = circular_fold(circular_fragment, 37.0, 3, 1.0, 0, 1.021)
            .unwrap_or_else(|error| panic!("{} circular: {error}", case.accession));
        assert!(circular.mfe_energy_kcal_mol.is_finite());
        assert!(circular.log_partition_function.is_finite());
        assert_eq!(
            circular.unpaired_probabilities.len(),
            circular_fragment.len()
        );

        let (modified_index, modified_kind, modified_symbol) = fragment
            .bytes()
            .position(|base| base == b'A')
            .map(|index| (index, Some(ModifiedBaseKind::M6a), "m6A"))
            .or_else(|| {
                fragment
                    .bytes()
                    .position(|base| base == b'U')
                    .map(|index| (index, Some(ModifiedBaseKind::Pseudouridine), "Psi"))
            })
            .unwrap_or((0, None, "mod"));
        let modification = ModifiedBase {
            position: modified_index + 1,
            symbol: modified_symbol.into(),
            canonical_base: fragment.as_bytes()[modified_index] as char,
            kind: modified_kind,
            paired_energy_kcal_mol: 0.1,
            unpaired_energy_kcal_mol: 0.0,
            stack_energy_kcal_mol: -0.1,
        };
        let modified = fold_modified(fragment, &[modification], 37.0, 3, 1.0, 2, 1.021)
            .unwrap_or_else(|error| panic!("{} modified: {error}", case.accession));
        assert_eq!(modified.display_symbols[modified_index], modified_symbol);
        assert!(modified.analysis.mfe_energy_kcal_mol.is_finite());

        let gquad = fold_gquad(fragment, 37.0, 3, 1.0, 2, 1.021)
            .unwrap_or_else(|error| panic!("{} gquad: {error}", case.accession));
        assert!((0.0..=1.0).contains(&gquad.gquad_probability));
        assert!(gquad.log_partition_function.is_finite());

        let pseudoknot = predict_pseudoknot(
            fragment,
            37.0,
            3,
            2,
            1.021,
            &PseudoknotOptions {
                threshold: 0.02,
                // This cross-operation smoke matrix checks API integration;
                // the separate 24-case ShPK validation exercises the default
                // exhaustive H-type state space.
                max_components: Some(4),
                max_ensemble_states: Some(512),
                ..PseudoknotOptions::default()
            },
        )
        .unwrap_or_else(|error| panic!("{} pseudoknot: {error}", case.accession));
        assert_eq!(
            parse_structure(fragment, &pseudoknot.structure)
                .unwrap()
                .pairs
                .len(),
            pseudoknot.pair_count
        );
        let evaluated_pseudoknot = evaluate_pseudoknot_structure(
            fragment,
            &pseudoknot.structure,
            37.0,
            3,
            2,
            1.021,
            &PseudoknotOptions::default(),
        )
        .unwrap_or_else(|error| panic!("{} pseudoknot evaluation: {error}", case.accession));
        assert!(evaluated_pseudoknot.energy.total_kcal_mol.is_finite());
        let arbitrary_fragment = fragment.chars().take(10).collect::<String>();
        let exact_arbitrary = predict_pseudoknot(
            &arbitrary_fragment,
            37.0,
            3,
            0,
            1.021,
            &PseudoknotOptions {
                exact_arbitrary_ensemble: true,
                ..PseudoknotOptions::default()
            },
        )
        .unwrap_or_else(|error| panic!("{} exact arbitrary pseudoknot: {error}", case.accession))
        .exact_arbitrary_ensemble
        .expect("exact arbitrary ensemble was requested");
        assert!(exact_arbitrary.state_space_complete);
        assert!(exact_arbitrary.log_partition_function.is_finite());
        for position in 1..=arbitrary_fragment.len() {
            let paired = exact_arbitrary
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == position || pair.j == position)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!(
                (paired + exact_arbitrary.unpaired_probabilities[position - 1] - 1.0).abs()
                    < 1.0e-10
            );
        }

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
        .unwrap_or_else(|error| panic!("{} comparative: {error}", case.accession));
        assert_eq!(comparative.consensus_sequence, fragment);
        assert!(comparative.consensus_energy_kcal_mol.is_finite());
        assert!(
            (comparative.analysis.mfe_energy_kcal_mol - analysis.mfe_energy_kcal_mol).abs()
                < 1.0e-10,
            "{} identical alignment MFE reduction drift: {} versus {}",
            case.accession,
            comparative.analysis.mfe_energy_kcal_mol,
            analysis.mfe_energy_kcal_mol,
        );
        assert!(
            (comparative.analysis.log_partition_function - analysis.log_partition_function).abs()
                < 1.0e-10,
            "{} identical alignment PF reduction drift",
            case.accession
        );
        assert_eq!(
            comparative.analysis.pair_probabilities.len(),
            analysis.pair_probabilities.len(),
            "{} identical alignment pair table length drift",
            case.accession
        );
        for (actual, expected) in comparative
            .analysis
            .pair_probabilities
            .iter()
            .zip(&analysis.pair_probabilities)
        {
            assert_eq!((actual.i, actual.j), (expected.i, expected.j));
            assert!(
                (actual.probability - expected.probability).abs() < 1.0e-10,
                "{} identical alignment pair probability drift at {}-{}",
                case.accession,
                actual.i,
                actual.j
            );
        }

        // Exact exponential features run on a real 10-nt prefix from every
        // family. Completeness is the contract; no sequence/structure beam is
        // introduced merely to keep this integration matrix fast.
        let exact_fragment = fragment.chars().take(10).collect::<String>();
        let exact_model = EnergyModel::with_dangles_and_salt(37.0, 0, 1.021).unwrap();
        let exact_constraints =
            ConstraintModel::compile(exact_fragment.len(), &ConstraintConfig::default()).unwrap();
        let exact_fold =
            fold_mfe_with_constraints(&exact_fragment, 3, &exact_model, &exact_constraints)
                .unwrap_or_else(|error| panic!("{} exact feature MFE: {error}", case.accession));
        let open_structure = ".".repeat(exact_fragment.len());
        let landscape = exact_landscape_path(
            &exact_fragment,
            &open_structure,
            &exact_fold.structure,
            3,
            &exact_model,
            &exact_constraints,
        )
        .unwrap_or_else(|error| panic!("{} landscape: {error}", case.accession));
        assert!(landscape.state_space_complete);
        assert_eq!(landscape.path.first().unwrap().structure, open_structure);
        assert_eq!(
            landscape.path.last().unwrap().structure,
            exact_fold.structure
        );
        assert!(
            (landscape
                .path
                .iter()
                .map(|step| step.energy_kcal_mol)
                .fold(f64::NEG_INFINITY, f64::max)
                - landscape.saddle_energy_kcal_mol)
                .abs()
                < 1.0e-12
        );

        let design = inverse_fold_exact(
            &open_structure,
            &exact_fragment,
            3,
            &exact_model,
            &ConstraintConfig::default(),
            &InverseDesignOptions {
                return_count: 1,
                ..InverseDesignOptions::default()
            },
        )
        .unwrap_or_else(|error| panic!("{} inverse design: {error}", case.accession));
        assert!(design.search_complete);
        assert_eq!(design.candidate_sequence_count, 1);
        assert_eq!(design.evaluated_sequence_count, 1);
        assert_eq!(design.candidates[0].sequence, exact_fragment);

        let ligand = ligand_ensemble_exact(
            &exact_fragment,
            &[LigandMotif {
                id: format!("{}-unstructured-domain", case.accession),
                start: 1,
                sequence: exact_fragment.clone(),
                structure: open_structure,
                standard_binding_energy_kcal_mol: -1.0,
                concentration_molar: 1.0e-3,
            }],
            3,
            1.0,
            &exact_model,
            &exact_constraints,
        )
        .unwrap_or_else(|error| panic!("{} ligand ensemble: {error}", case.accession));
        assert!(ligand.state_space_complete);
        assert!((0.0..=1.0).contains(&ligand.motifs[0].occupancy_probability));
        for position in 1..=exact_fragment.len() {
            let paired = ligand
                .pair_probabilities
                .iter()
                .filter(|pair| pair.i == position || pair.j == position)
                .map(|pair| pair.probability)
                .sum::<f64>();
            assert!(
                (paired + ligand.unpaired_probabilities[position - 1] - 1.0).abs() < 1.0e-10,
                "{} ligand probability normalization at {position}",
                case.accession
            );
        }

        let profile = ParameterProfile::default();
        let profiled = analyze_with_profile(fragment, &profile)
            .unwrap_or_else(|error| panic!("{} profile: {error}", case.accession));
        assert_eq!(profiled.mfe_structure, analysis.mfe_structure);

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
        .unwrap_or_else(|error| panic!("{} cofold: {error}", case.accession));
        let equilibrium = complex.equilibrium.unwrap();
        assert!(equilibrium.mass_balance_a_error_molar.abs() < 1.0e-15);
        assert!(equilibrium.mass_balance_b_error_molar.abs() < 1.0e-15);
    }
}
