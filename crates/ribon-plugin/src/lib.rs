//! Minimal Typst/WASM boundary for the Ribon request protocol.

use wasm_minimal_protocol::*;

mod protocol;

initiate_protocol!();

/// Execute one stable `ribon.analysis/1` JSON request.
///
/// The function always returns a JSON response envelope. Validation and
/// analysis failures therefore remain structured data instead of crossing the
/// WASM ABI as implementation-specific traps.
#[wasm_func]
pub fn run(request: &[u8]) -> Vec<u8> {
    protocol::execute_bytes(request)
}
