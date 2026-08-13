use ribon_core::{layout_structure, LayoutKind};
use std::env;
use std::process::ExitCode;

fn main() -> ExitCode {
    let mut arguments = env::args().skip(1);
    let Some(sequence) = arguments.next() else {
        eprintln!("usage: cargo run -p ribon-core --example layout -- SEQUENCE STRUCTURE [METHOD]");
        return ExitCode::from(2);
    };
    let Some(structure) = arguments.next() else {
        eprintln!("missing dot-bracket structure");
        return ExitCode::from(2);
    };
    let method = arguments.next().unwrap_or_else(|| "naview".into());
    let kind = match method.parse::<LayoutKind>() {
        Ok(kind) => kind,
        Err(error) => {
            eprintln!("{error}");
            return ExitCode::from(2);
        }
    };
    match layout_structure(&sequence, &structure, kind) {
        Ok(layout) => {
            println!(
                "{}",
                serde_json::to_string(&layout).expect("layout serializes")
            );
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(1)
        }
    }
}
