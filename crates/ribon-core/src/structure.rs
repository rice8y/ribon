use serde::Serialize;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Pair {
    /// One-based nucleotide index.
    pub i: usize,
    /// One-based nucleotide index.
    pub j: usize,
    /// Bracket alphabet level: 0=(), 1=[], 2={}, 3=<>, 4+=A/a ... Z/z.
    pub level: usize,
    pub canonical: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct ParsedStructure {
    pub sequence: String,
    pub structure: String,
    pub length: usize,
    pub pairs: Vec<Pair>,
    /// One-based positions after which a strand break occurred.
    pub strand_breaks: Vec<usize>,
    #[serde(skip)]
    pub partner: Vec<Option<usize>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RnaError {
    EmptySequence,
    InvalidSequence { position: usize, symbol: char },
    LengthMismatch { sequence: usize, structure: usize },
    InvalidStructure { position: usize, symbol: char },
    UnmatchedOpening { position: usize, symbol: char },
    UnmatchedClosing { position: usize, symbol: char },
    MultiplePartners { position: usize },
    PseudoknotUnsupported(&'static str),
    InvalidOption(String),
    Numerical(String),
}

impl Display for RnaError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequence => write!(f, "RNA sequence is empty"),
            Self::InvalidSequence { position, symbol } => write!(
                f,
                "invalid RNA symbol {symbol:?} at sequence position {position}"
            ),
            Self::LengthMismatch {
                sequence,
                structure,
            } => write!(
                f,
                "sequence/structure length mismatch ({sequence} nt versus {structure} symbols)"
            ),
            Self::InvalidStructure { position, symbol } => write!(
                f,
                "invalid dot-bracket symbol {symbol:?} at structure position {position}"
            ),
            Self::UnmatchedOpening { position, symbol } => write!(
                f,
                "unmatched opening bracket {symbol:?} at structure position {position}"
            ),
            Self::UnmatchedClosing { position, symbol } => write!(
                f,
                "unmatched closing bracket {symbol:?} at structure position {position}"
            ),
            Self::MultiplePartners { position } => {
                write!(f, "nucleotide {position} is paired more than once")
            }
            Self::PseudoknotUnsupported(operation) => {
                write!(f, "{operation} requires a pseudoknot-free structure")
            }
            Self::InvalidOption(message) | Self::Numerical(message) => f.write_str(message),
        }
    }
}

impl Error for RnaError {}

fn normalize_sequence_with_breaks(input: &str) -> Result<(String, Vec<usize>), RnaError> {
    let mut out = String::new();
    let mut strand_breaks = Vec::new();
    let mut position = 0usize;
    for raw in input.chars() {
        if raw.is_whitespace() {
            continue;
        }
        if raw == '&' {
            if out.is_empty() || strand_breaks.last() == Some(&out.len()) {
                return Err(RnaError::InvalidOption(
                    "strand separators must occur once between non-empty strands".into(),
                ));
            }
            strand_breaks.push(out.len());
            continue;
        }
        position += 1;
        let c = raw.to_ascii_uppercase();
        // IUPAC ambiguity codes are accepted but only A/C/G/U/T participate
        // in the built-in energy models. T uses the U-shaped thermodynamic
        // lookup index while remaining T in public results and layouts.
        if !matches!(
            c,
            'A' | 'C'
                | 'G'
                | 'U'
                | 'T'
                | 'R'
                | 'Y'
                | 'S'
                | 'W'
                | 'K'
                | 'M'
                | 'B'
                | 'D'
                | 'H'
                | 'V'
                | 'N'
        ) {
            return Err(RnaError::InvalidSequence {
                position,
                symbol: raw,
            });
        }
        out.push(c);
    }
    if out.is_empty() {
        return Err(RnaError::EmptySequence);
    }
    if strand_breaks.last() == Some(&out.len()) {
        return Err(RnaError::InvalidOption(
            "strand separators must occur once between non-empty strands".into(),
        ));
    }
    Ok((out, strand_breaks))
}

/// Normalize whitespace and case while preserving DNA thymine. Strand
/// separators are removed; structure parsing uses an internal variant that
/// also retains their indices.
pub fn normalize_sequence(input: &str) -> Result<String, RnaError> {
    normalize_sequence_with_breaks(input).map(|(sequence, _)| sequence)
}

fn bracket_kind(c: char) -> Option<(usize, bool)> {
    match c {
        '(' => Some((0, true)),
        ')' => Some((0, false)),
        '[' => Some((1, true)),
        ']' => Some((1, false)),
        '{' => Some((2, true)),
        '}' => Some((2, false)),
        '<' => Some((3, true)),
        '>' => Some((3, false)),
        'A'..='Z' => Some((4 + (c as usize - 'A' as usize), true)),
        'a'..='z' => Some((4 + (c as usize - 'a' as usize), false)),
        _ => None,
    }
}

fn pair_is_canonical(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'U')
            | (b'U', b'A')
            | (b'A', b'T')
            | (b'T', b'A')
            | (b'C', b'G')
            | (b'G', b'C')
            | (b'G', b'U')
            | (b'U', b'G')
            | (b'G', b'T')
            | (b'T', b'G')
    )
}

pub fn parse_structure(sequence: &str, structure: &str) -> Result<ParsedStructure, RnaError> {
    let (sequence, sequence_breaks) = normalize_sequence_with_breaks(sequence)?;
    let mut symbols = Vec::new();
    let mut strand_breaks = Vec::new();
    for raw in structure.chars() {
        if raw.is_whitespace() {
            continue;
        }
        if raw == '&' {
            if symbols.is_empty() || strand_breaks.last() == Some(&symbols.len()) {
                return Err(RnaError::InvalidOption(
                    "structure strand separators must occur once between non-empty strands".into(),
                ));
            }
            strand_breaks.push(symbols.len());
        } else {
            symbols.push(raw);
        }
    }
    if symbols.len() != sequence.len() {
        return Err(RnaError::LengthMismatch {
            sequence: sequence.len(),
            structure: symbols.len(),
        });
    }
    if strand_breaks.last() == Some(&symbols.len()) {
        return Err(RnaError::InvalidOption(
            "structure strand separators must occur once between non-empty strands".into(),
        ));
    }
    if !sequence_breaks.is_empty() && !strand_breaks.is_empty() && sequence_breaks != strand_breaks
    {
        return Err(RnaError::InvalidOption(format!(
            "sequence and structure strand breaks differ ({sequence_breaks:?} versus {strand_breaks:?})"
        )));
    }
    if strand_breaks.is_empty() {
        strand_breaks = sequence_breaks;
    }

    let mut stacks: Vec<Vec<usize>> = vec![Vec::new(); 30];
    let mut partner = vec![None; symbols.len()];
    let mut pairs = Vec::new();
    let bases = sequence.as_bytes();

    for (index, &symbol) in symbols.iter().enumerate() {
        // `+` and `~` are ViennaRNA's G-quadruplex run/linker markers. Their
        // topology is carried by the analysis result; for ordinary pair
        // parsing both occupy unpaired positions so vector layouts render them.
        if matches!(symbol, '.' | '_' | '-' | ',' | ':' | '+' | '~') {
            continue;
        }
        let Some((level, opening)) = bracket_kind(symbol) else {
            return Err(RnaError::InvalidStructure {
                position: index + 1,
                symbol,
            });
        };
        if opening {
            stacks[level].push(index);
        } else {
            let Some(left) = stacks[level].pop() else {
                return Err(RnaError::UnmatchedClosing {
                    position: index + 1,
                    symbol,
                });
            };
            if partner[left].is_some() || partner[index].is_some() {
                return Err(RnaError::MultiplePartners {
                    position: index + 1,
                });
            }
            partner[left] = Some(index);
            partner[index] = Some(left);
            pairs.push(Pair {
                i: left + 1,
                j: index + 1,
                level,
                canonical: pair_is_canonical(bases[left], bases[index]),
            });
        }
    }

    for (level, stack) in stacks.iter().enumerate() {
        if let Some(&position) = stack.last() {
            let symbol = if level == 0 {
                '('
            } else if level == 1 {
                '['
            } else if level == 2 {
                '{'
            } else if level == 3 {
                '<'
            } else {
                char::from_u32(('A' as u32) + (level - 4) as u32).unwrap_or('A')
            };
            return Err(RnaError::UnmatchedOpening {
                position: position + 1,
                symbol,
            });
        }
    }

    pairs.sort_by_key(|pair| (pair.i, pair.j));
    Ok(ParsedStructure {
        length: sequence.len(),
        sequence,
        structure: symbols.into_iter().collect(),
        pairs,
        strand_breaks,
        partner,
    })
}

pub fn is_pseudoknotted(pairs: &[Pair]) -> bool {
    for (index, a) in pairs.iter().enumerate() {
        for b in &pairs[index + 1..] {
            if (a.i < b.i && b.i < a.j && a.j < b.j) || (b.i < a.i && a.i < b.j && b.j < a.j) {
                return true;
            }
        }
    }
    false
}

pub fn pairs_to_dot_bracket(length: usize, pairs: &[(usize, usize)]) -> String {
    let mut chars = vec!['.'; length];
    for &(i, j) in pairs {
        if i < j && j < length {
            chars[i] = '(';
            chars[j] = ')';
        }
    }
    chars.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_extended_dot_bracket_and_preserves_dna_thymine() {
        let parsed = parse_structure("atgcaa", "([.)].").unwrap();
        assert_eq!(parsed.sequence, "ATGCAA");
        assert_eq!(parsed.pairs.len(), 2);
        assert!(is_pseudoknotted(&parsed.pairs));
    }

    #[test]
    fn rejects_mismatch() {
        let error = parse_structure("ACGU", "...").unwrap_err();
        assert!(matches!(error, RnaError::LengthMismatch { .. }));
    }

    #[test]
    fn retains_and_cross_checks_strand_breaks() {
        let parsed = parse_structure("GGG&CCC", "(((&)))").unwrap();
        assert_eq!(parsed.strand_breaks, vec![3]);
        assert!(parse_structure("GG&GCCC", "(((&)))").is_err());
        assert!(parse_structure("&GGGCCC", "...&...").is_err());
    }
}
