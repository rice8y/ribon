use ribon_core::{suboptimal_structures, ConstraintConfig, ConstraintModel, EnergyModel};

fn main() {
    let mut args = std::env::args().skip(1);
    let sequence = args
        .next()
        .expect("usage: suboptimal-json SEQUENCE BAND [DANGLES] [LIMIT]");
    let band = args
        .next()
        .expect("missing BAND")
        .parse::<f64>()
        .expect("BAND must be a number");
    let dangles = args
        .next()
        .map(|value| value.parse::<u8>().expect("DANGLES must be an integer"))
        .unwrap_or(2);
    let limit = args
        .next()
        .map(|value| value.parse::<usize>().expect("LIMIT must be an integer"))
        .unwrap_or(50);
    let model = EnergyModel::with_dangles(37.0, dangles).expect("invalid model");
    let constraints = ConstraintModel::compile(sequence.len(), &ConstraintConfig::default())
        .expect("invalid constraints");
    let result = suboptimal_structures(&sequence, 3, &model, &constraints, band, limit)
        .expect("suboptimal prediction failed");
    println!("{}", serde_json::to_string(&result).unwrap());
}
