use ribon_core::{analyze_with_constraints, fold_sequence_with_constraints, ConstraintConfig};

fn main() {
    let mut arguments = std::env::args().skip(1);
    let sequence = arguments
        .next()
        .expect("usage: constrained SEQUENCE DANGLES CONSTRAINTS_JSON [fold]");
    let dangles = arguments
        .next()
        .expect("missing DANGLES")
        .parse::<u8>()
        .expect("DANGLES must be 0, 1, 2, or 3");
    let config: ConstraintConfig =
        serde_json::from_str(&arguments.next().expect("missing CONSTRAINTS_JSON"))
            .expect("invalid CONSTRAINTS_JSON");
    let output = if arguments.next().as_deref() == Some("fold") {
        serde_json::to_value(
            fold_sequence_with_constraints(&sequence, 37.0, 3, dangles, &config)
                .expect("constrained MFE failed"),
        )
        .unwrap()
    } else {
        serde_json::to_value(
            analyze_with_constraints(&sequence, 37.0, 3, 1.0, dangles, &config)
                .expect("constrained analysis failed"),
        )
        .unwrap()
    };
    println!("{}", serde_json::to_string(&output).unwrap());
}
