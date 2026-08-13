use ribon_core::duplex;

fn main() {
    let mut args = std::env::args().skip(1);
    let sequence_a = args.next().expect("usage: duplex-json A B [SALT]");
    let sequence_b = args.next().expect("usage: duplex-json A B [SALT]");
    let salt = args
        .next()
        .map(|value| value.parse::<f64>().expect("SALT must be a number"))
        .unwrap_or(1.021);
    let result = duplex(&sequence_a, &sequence_b, 37.0, salt).expect("duplex prediction failed");
    println!("{}", serde_json::to_string(&result).unwrap());
}
