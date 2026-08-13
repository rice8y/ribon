use ribon_core::{analyze_with_options, ConstraintConfig};

fn main() {
    let sequence = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "GGGAAACCC".to_string());
    let dangles = std::env::args()
        .nth(2)
        .map(|value| value.parse::<u8>().expect("DANGLES must be 0 or 2"))
        .unwrap_or(2);
    let salt = std::env::args()
        .nth(3)
        .map(|value| value.parse::<f64>().expect("SALT must be molarity"))
        .unwrap_or(1.021);
    match analyze_with_options(
        &sequence,
        37.0,
        3,
        1.0,
        dangles,
        salt,
        &ConstraintConfig::default(),
    ) {
        Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
        Err(error) => {
            eprintln!("ribon: {error}");
            std::process::exit(2);
        }
    }
}
