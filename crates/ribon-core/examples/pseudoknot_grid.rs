use ribon_core::{parse_structure, predict_pseudoknot, PseudoknotOptions};
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Deserialize)]
struct Corpus {
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    id: String,
    sequence: String,
    structure: String,
}

fn pairs(sequence: &str, structure: &str) -> HashSet<(usize, usize)> {
    parse_structure(sequence, structure)
        .expect("valid corpus structure")
        .pairs
        .into_iter()
        .map(|pair| (pair.i, pair.j))
        .collect()
}

fn main() {
    let corpus: Corpus =
        serde_json::from_str(include_str!("../../../tests/data/pseudoknot_real_24.json"))
            .expect("valid corpus");
    for threshold in [0.0, 0.01, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2] {
        for iterations in [1, 2, 3, 5] {
            for min_helix in [1, 2] {
                let options = PseudoknotOptions {
                    threshold,
                    iterations,
                    min_helix,
                    ..PseudoknotOptions::default()
                };
                let mut sensitivity = 0.0;
                let mut precision = 0.0;
                let mut crossing = 0;
                for case in &corpus.cases {
                    let prediction =
                        predict_pseudoknot(&case.sequence, 37.0, 3, 2, 1.021, &options)
                            .expect("prediction");
                    let expected = pairs(&case.sequence, &case.structure);
                    let predicted = pairs(&case.sequence, &prediction.structure);
                    let hits = expected.intersection(&predicted).count();
                    sensitivity += hits as f64 / expected.len() as f64;
                    precision += if predicted.is_empty() {
                        0.0
                    } else {
                        hits as f64 / predicted.len() as f64
                    };
                    crossing += usize::from(prediction.crossing_count > 0);
                }
                let n = corpus.cases.len() as f64;
                println!(
                    "threshold={threshold:.2} iterations={iterations} min_helix={min_helix} sensitivity={:.6} precision={:.6} crossing={crossing}",
                    sensitivity / n,
                    precision / n
                );
            }
        }
    }
    for evidence_weight_kcal_mol in [0.0, 0.35, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0] {
        let options = PseudoknotOptions {
            threshold: 0.2,
            iterations: 1,
            min_helix: 2,
            evidence_weight_kcal_mol,
            ..PseudoknotOptions::default()
        };
        let mut sensitivity = 0.0;
        let mut precision = 0.0;
        for case in &corpus.cases {
            let prediction =
                predict_pseudoknot(&case.sequence, 37.0, 3, 2, 1.021, &options).unwrap();
            let expected = pairs(&case.sequence, &case.structure);
            let predicted = pairs(&case.sequence, &prediction.structure);
            let hits = expected.intersection(&predicted).count();
            sensitivity += hits as f64 / expected.len() as f64;
            precision += hits as f64 / predicted.len() as f64;
        }
        println!(
            "evidence={evidence_weight_kcal_mol:.2} sensitivity={:.6} precision={:.6}",
            sensitivity / corpus.cases.len() as f64,
            precision / corpus.cases.len() as f64
        );
    }
    for unpaired_kcal_mol in [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5] {
        let options = PseudoknotOptions {
            threshold: 0.2,
            iterations: 1,
            min_helix: 2,
            evidence_weight_kcal_mol: 0.0,
            unpaired_kcal_mol,
            ..PseudoknotOptions::default()
        };
        let mut sensitivity = 0.0;
        let mut precision = 0.0;
        for case in &corpus.cases {
            let prediction =
                predict_pseudoknot(&case.sequence, 37.0, 3, 2, 1.021, &options).unwrap();
            let expected = pairs(&case.sequence, &case.structure);
            let predicted = pairs(&case.sequence, &prediction.structure);
            let hits = expected.intersection(&predicted).count();
            sensitivity += hits as f64 / expected.len() as f64;
            precision += hits as f64 / predicted.len() as f64;
        }
        println!(
            "unpaired={unpaired_kcal_mol:.2} sensitivity={:.6} precision={:.6}",
            sensitivity / corpus.cases.len() as f64,
            precision / corpus.cases.len() as f64
        );
    }
    let options = PseudoknotOptions {
        threshold: 0.2,
        iterations: 1,
        min_helix: 2,
        evidence_weight_kcal_mol: 0.0,
        unpaired_kcal_mol: 0.05,
        ..PseudoknotOptions::default()
    };
    for case in &corpus.cases {
        let prediction = predict_pseudoknot(&case.sequence, 37.0, 3, 2, 1.021, &options).unwrap();
        let expected = pairs(&case.sequence, &case.structure);
        let predicted = pairs(&case.sequence, &prediction.structure);
        let hits = expected.intersection(&predicted).count();
        println!(
            "case={} sensitivity={:.3} precision={:.3} core={:?} score={:?} expected={} predicted={}",
            case.id,
            hits as f64 / expected.len() as f64,
            hits as f64 / predicted.len() as f64,
            prediction
                .thermodynamic_core_pairs
                .iter()
                .map(|pair| (pair.i, pair.j))
                .collect::<Vec<_>>(),
            prediction.thermodynamic_core_score_kcal_mol,
            case.structure,
            prediction.structure,
        );
    }
}
