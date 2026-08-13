use ribon_core::{fold_mfe, EnergyModel};

fn main() {
    let mut arguments = std::env::args().skip(1);
    let sequence = arguments
        .next()
        .expect("usage: fold SEQUENCE [DANGLES] [STRUCTURE]");
    let dangles = arguments
        .next()
        .map(|value| value.parse::<u8>().expect("DANGLES must be 0, 1, 2, or 3"))
        .unwrap_or(2);
    let model =
        EnergyModel::with_dangles_and_salt(37.0, dangles, 1.021).expect("invalid dangle model");
    if let Some(structure) = arguments.next() {
        let result = model
            .evaluate(&sequence, &structure)
            .expect("structure evaluation failed");
        println!("{}", serde_json::to_string_pretty(&result).unwrap());
        return;
    }
    let result = fold_mfe(&sequence, 3, &model).expect("MFE folding failed");
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}
