//! Sparse nearest-neighbor parameters for common naturally modified bases.
//!
//! The tables are transcribed from the primary thermodynamic studies listed in
//! [`ModifiedBaseKind::source_url`].  A key stores a nearest neighbor as
//! `5'-xy-3' / 3'-zw-5'`, i.e. `[x, y, z, w]`.  Values are kcal/mol.

use serde::{Deserialize, Serialize};

const T37_KELVIN: f64 = 310.15;

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, Hash)]
#[serde(rename_all = "kebab-case")]
pub enum ModifiedBaseKind {
    #[serde(alias = "m6A", alias = "m6a")]
    M6a,
    #[serde(alias = "psi", alias = "Psi", alias = "pseudouridine")]
    Pseudouridine,
    #[serde(alias = "I", alias = "inosine")]
    Inosine,
    #[serde(alias = "7DA", alias = "7da", alias = "7-deazaadenosine")]
    SevenDeazaadenosine,
    #[serde(alias = "P", alias = "nebularine")]
    Purine,
    #[serde(alias = "D", alias = "dihydrouridine")]
    Dihydrouridine,
}

impl ModifiedBaseKind {
    pub(crate) fn precursor(self) -> u8 {
        match self {
            Self::M6a | Self::Inosine | Self::SevenDeazaadenosine | Self::Purine => b'A',
            Self::Pseudouridine | Self::Dihydrouridine => b'U',
        }
    }

    pub(crate) fn folding_base(self) -> u8 {
        match self {
            // Inosine pairs with C and U in the supported sparse parameter set.
            Self::Inosine => b'G',
            _ => self.precursor(),
        }
    }

    fn code(self) -> u8 {
        match self {
            Self::M6a => b'M',
            Self::Pseudouridine => b'Y',
            Self::Inosine => b'I',
            Self::SevenDeazaadenosine => b'7',
            Self::Purine => b'P',
            Self::Dihydrouridine => b'D',
        }
    }

    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::M6a => "N6-methyladenosine",
            Self::Pseudouridine => "pseudouridine",
            Self::Inosine => "inosine",
            Self::SevenDeazaadenosine => "7-deazaadenosine",
            Self::Purine => "purine/nebularine",
            Self::Dihydrouridine => "dihydrouridine",
        }
    }

    pub(crate) fn source_url(self) -> &'static str {
        match self {
            Self::M6a => "https://doi.org/10.1038/s41467-022-28817-4",
            Self::Pseudouridine => "https://doi.org/10.1261/rna.039610.113",
            Self::Inosine => "https://doi.org/10.1093/nar/gky907",
            Self::SevenDeazaadenosine => "https://doi.org/10.1261/rna.055277.115",
            Self::Purine => "https://doi.org/10.1093/nar/gkw830",
            Self::Dihydrouridine => "https://doi.org/10.1093/bioinformatics/btad696",
        }
    }

    pub(crate) fn calibration(self) -> &'static str {
        match self {
            Self::M6a => {
                "experimental nearest-neighbor delta-G at 37 C; experimental terminal m6A-U term"
            }
            Self::Pseudouridine => "experimental nearest-neighbor delta-H and delta-G at 37 C",
            Self::Inosine => {
                "experimental I-C delta-H/delta-G plus experimental I-U delta-G at 37 C"
            }
            Self::SevenDeazaadenosine => {
                "experimental nearest-neighbor delta-H and delta-G at 37 C"
            }
            Self::Purine => "experimental nearest-neighbor delta-H and delta-G at 37 C",
            Self::Dihydrouridine => {
                "published model-based +1.5 kcal/mol correction per affected stack"
            }
        }
    }
}

#[derive(Clone, Copy)]
struct Parameter {
    key: [u8; 4],
    dg37: f64,
    dh: Option<f64>,
}

impl Parameter {
    const fn g(key: &[u8; 4], dg37: f64) -> Self {
        Self {
            key: *key,
            dg37,
            dh: None,
        }
    }

    const fn gh(key: &[u8; 4], dg37: f64, dh: f64) -> Self {
        Self {
            key: *key,
            dg37,
            dh: Some(dh),
        }
    }

    fn free_energy(self, temperature_celsius: f64) -> f64 {
        if (temperature_celsius - 37.0).abs() < 1.0e-12 {
            return self.dg37;
        }
        let Some(dh) = self.dh else {
            return self.dg37;
        };
        let entropy = (dh - self.dg37) / T37_KELVIN;
        dh - (temperature_celsius + 273.15) * entropy
    }
}

const M6A: &[Parameter] = &[
    Parameter::g(b"MCUG", -1.79),
    Parameter::g(b"UCMG", -1.72),
    Parameter::g(b"MGUC", -1.56),
    Parameter::g(b"UGMC", -1.24),
    Parameter::g(b"MUUA", -1.10),
    Parameter::g(b"MAUU", -0.92),
    Parameter::g(b"UUMA", -0.83),
    Parameter::g(b"UAMU", -0.73),
    Parameter::g(b"MUUG", -0.69),
    Parameter::g(b"MUUM", -0.46),
    Parameter::g(b"UGMU", -0.32),
    Parameter::g(b"UUMG", -0.32),
    Parameter::g(b"MMUU", -0.21),
    Parameter::g(b"MGUU", -0.03),
    Parameter::g(b"UMMU", 1.45),
];

const PSI: &[Parameter] = &[
    Parameter::gh(b"AYUA", -2.80, -22.08),
    Parameter::gh(b"CYGA", -2.77, -16.23),
    Parameter::gh(b"GYCA", -3.29, -24.07),
    Parameter::gh(b"UYAA", -1.62, -20.81),
];

// Rust byte-string literals cannot express whitespace-free visual grouping;
// the second orientation of the pseudouridine table is kept separately.
const PSI_REVERSE: &[Parameter] = &[
    Parameter::gh(b"YAAU", -2.10, -12.47),
    Parameter::gh(b"YCAG", -2.49, -17.29),
    Parameter::gh(b"YGAC", -2.20, -11.19),
    Parameter::gh(b"YUAA", -2.74, -26.94),
];

const SEVEN_DA: &[Parameter] = &[
    Parameter::gh(b"A7UU", -0.59, -8.4),
    Parameter::gh(b"C7GU", -1.81, -11.8),
    Parameter::gh(b"G7CU", -1.66, -10.8),
    Parameter::gh(b"U7AU", -1.07, -9.4),
    Parameter::gh(b"7AUU", -0.68, -9.9),
    Parameter::gh(b"7CUG", -2.10, -14.8),
    Parameter::gh(b"7GUC", -1.98, -15.1),
    Parameter::gh(b"7UUA", -1.46, -13.9),
];

const INOSINE_C: &[Parameter] = &[
    Parameter::gh(b"IGCC", -2.23, -14.5),
    Parameter::gh(b"ICCG", -1.89, -10.6),
    Parameter::gh(b"IACU", -1.18, -15.3),
    Parameter::gh(b"IUCA", -1.02, -7.7),
    Parameter::gh(b"GICC", -2.62, -16.8),
    Parameter::gh(b"CIGC", -1.86, -12.7),
    Parameter::gh(b"AIUC", -1.57, -14.2),
    Parameter::gh(b"UIAC", -0.96, -11.8),
];

const INOSINE_U: &[Parameter] = &[
    Parameter::g(b"AIUU", -0.41),
    Parameter::g(b"CIGU", -0.77),
    Parameter::g(b"GICU", -1.34),
    Parameter::g(b"UIAU", 0.37),
    Parameter::g(b"IAUU", 0.43),
    Parameter::g(b"ICUG", -1.03),
    Parameter::g(b"IGUC", -1.22),
    Parameter::g(b"IUUA", -0.50),
];

const PURINE: &[Parameter] = &[
    Parameter::gh(b"APUU", 0.43, -14.0),
    Parameter::gh(b"CPGU", -0.76, -12.4),
    Parameter::gh(b"GPCU", -1.10, -14.2),
    Parameter::gh(b"UPAU", 0.33, -8.7),
    Parameter::gh(b"PAUU", -0.68, -10.4),
    Parameter::gh(b"PCUG", -1.98, -15.7),
    Parameter::gh(b"PGUC", -1.88, -14.5),
    Parameter::gh(b"PUUA", -0.32, -11.9),
];

fn reverse_key(key: [u8; 4]) -> [u8; 4] {
    [key[3], key[2], key[1], key[0]]
}

fn lookup_table(table: &[Parameter], key: [u8; 4], temperature_celsius: f64) -> Option<f64> {
    table
        .iter()
        .copied()
        .find(|parameter| parameter.key == key || parameter.key == reverse_key(key))
        .map(|parameter| parameter.free_energy(temperature_celsius))
}

pub(crate) fn stack_energy(
    kind: ModifiedBaseKind,
    canonical_key: [u8; 4],
    modified_offsets: &[usize],
    temperature_celsius: f64,
) -> Option<f64> {
    if kind == ModifiedBaseKind::Dihydrouridine {
        return None;
    }
    let mut key = canonical_key;
    for &offset in modified_offsets {
        key[offset] = kind.code();
    }
    let tables: &[&[Parameter]] = match kind {
        ModifiedBaseKind::M6a => &[M6A],
        ModifiedBaseKind::Pseudouridine => &[PSI, PSI_REVERSE],
        ModifiedBaseKind::Inosine => &[INOSINE_C, INOSINE_U],
        ModifiedBaseKind::SevenDeazaadenosine => &[SEVEN_DA],
        ModifiedBaseKind::Purine => &[PURINE],
        ModifiedBaseKind::Dihydrouridine => &[],
    };
    tables
        .iter()
        .find_map(|table| lookup_table(table, key, temperature_celsius))
}

pub(crate) fn terminal_energy(
    kind: ModifiedBaseKind,
    other_base: u8,
    temperature_celsius: f64,
) -> Option<f64> {
    let parameter = match (kind, other_base) {
        (ModifiedBaseKind::M6a, b'U') => Parameter::g(b"MUAU", 0.13),
        (ModifiedBaseKind::Pseudouridine, b'A') => Parameter::gh(b"YAAU", 0.31, -2.04),
        (ModifiedBaseKind::Inosine, b'C') => Parameter::gh(b"ICCG", -0.08, 2.0),
        (ModifiedBaseKind::SevenDeazaadenosine, b'U') => Parameter::gh(b"7UAU", 0.31, 9.3),
        (ModifiedBaseKind::Purine, b'U') => Parameter::gh(b"PUAU", 0.86, 2.3),
        _ => return None,
    };
    Some(parameter.free_energy(temperature_celsius))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_values_are_exact_at_37_celsius() {
        assert_eq!(
            stack_energy(ModifiedBaseKind::M6a, *b"ACUG", &[0], 37.0),
            Some(-1.79)
        );
        assert_eq!(
            stack_energy(ModifiedBaseKind::Pseudouridine, *b"AUUA", &[1], 37.0),
            Some(-2.80)
        );
        assert_eq!(
            stack_energy(ModifiedBaseKind::Inosine, *b"AGCC", &[0], 37.0),
            Some(-2.23)
        );
    }

    #[test]
    fn enthalpy_parameters_scale_with_temperature() {
        let at_37 = stack_energy(ModifiedBaseKind::Pseudouridine, *b"AUUA", &[1], 37.0).unwrap();
        let at_10 = stack_energy(ModifiedBaseKind::Pseudouridine, *b"AUUA", &[1], 10.0).unwrap();
        assert_eq!(at_37, -2.80);
        assert!(at_10 < at_37);
    }

    #[test]
    fn every_calibrated_table_family_reproduces_a_published_value() {
        let cases = [
            (ModifiedBaseKind::Pseudouridine, *b"UUAA", vec![0], -2.74),
            (ModifiedBaseKind::Inosine, *b"AGUU", vec![1], -0.41),
            (
                ModifiedBaseKind::SevenDeazaadenosine,
                *b"AAUU",
                vec![1],
                -0.59,
            ),
            (ModifiedBaseKind::Purine, *b"AAUU", vec![1], 0.43),
        ];
        for (kind, key, offsets, expected) in cases {
            let actual = stack_energy(kind, key, &offsets, 37.0).unwrap();
            assert!((actual - expected).abs() < 1e-12, "{kind:?}: {actual}");
        }
        assert_eq!(
            terminal_energy(ModifiedBaseKind::M6a, b'U', 37.0),
            Some(0.13)
        );
        assert_eq!(
            terminal_energy(ModifiedBaseKind::Pseudouridine, b'A', 37.0),
            Some(0.31)
        );
        assert_eq!(
            terminal_energy(ModifiedBaseKind::Inosine, b'C', 37.0),
            Some(-0.08)
        );
        assert_eq!(
            terminal_energy(ModifiedBaseKind::SevenDeazaadenosine, b'U', 37.0),
            Some(0.31)
        );
        assert_eq!(
            terminal_energy(ModifiedBaseKind::Purine, b'U', 37.0),
            Some(0.86)
        );
    }

    #[test]
    fn provenance_urls_identify_the_parameter_sources() {
        let expected = [
            (
                ModifiedBaseKind::M6a,
                "https://doi.org/10.1038/s41467-022-28817-4",
            ),
            (
                ModifiedBaseKind::Pseudouridine,
                "https://doi.org/10.1261/rna.039610.113",
            ),
            (
                ModifiedBaseKind::Inosine,
                "https://doi.org/10.1093/nar/gky907",
            ),
            (
                ModifiedBaseKind::SevenDeazaadenosine,
                "https://doi.org/10.1261/rna.055277.115",
            ),
            (
                ModifiedBaseKind::Purine,
                "https://doi.org/10.1093/nar/gkw830",
            ),
            (
                ModifiedBaseKind::Dihydrouridine,
                "https://doi.org/10.1093/bioinformatics/btad696",
            ),
        ];
        for (kind, source) in expected {
            assert_eq!(kind.source_url(), source, "{kind:?}");
        }
    }
}
