use std::env;
use std::fs;
use std::path::{Path, PathBuf};

const INF: i32 = 10_000_000;
const DATA_DIR: &str = "data/rnastructure-6.6";
const STACK_PAIR_ORDER: [usize; 6] = [4, 0, 1, 2, 5, 3];
const INTERNAL_PAIR_ORDER: [usize; 6] = [4, 0, 1, 5, 2, 3];

fn main() {
    let data = Path::new(DATA_DIR);
    println!("cargo:rerun-if-changed={}", data.display());

    generate_parameter_module(
        data,
        "rna",
        50,
        370,
        terminal_pair_penalty(data, "rna.helix_ends.dg"),
        terminal_pair_penalty(data, "rna.helix_ends.dh"),
        "Turner 2004 RNA",
        "turner2004_generated.rs",
    );
    generate_parameter_module(
        data,
        "dna",
        0,
        named_scalar_centi(data, "dna.miscloop.dh", "terminal AU penalty"),
        terminal_pair_penalty(data, "dna.helix_ends.dg"),
        named_scalar_centi(data, "dna.miscloop.dh", "terminal AU penalty"),
        "Mathews 2004 DNA",
        "mathews2004_dna_generated.rs",
    );
}

#[allow(clippy::too_many_arguments)]
fn generate_parameter_module(
    data: &Path,
    prefix: &str,
    terminal_37: i32,
    terminal_dh: i32,
    terminal_pair_37: i32,
    terminal_pair_dh: i32,
    description: &str,
    output_name: &str,
) {
    let file = |stem: &str, suffix: &str| format!("{prefix}.{stem}.{suffix}");
    let stack_37 = stack_table(data, &file("stack", "dg"));
    let stack_dh = stack_table(data, &file("stack", "dh"));
    let mismatch_h_37 = mismatch_table(data, &file("tstackh", "dg"), terminal_37);
    let mismatch_h_dh = mismatch_table(data, &file("tstackh", "dh"), terminal_dh);
    let mismatch_i_37 = mismatch_table(data, &file("tstacki", "dg"), terminal_37);
    let mismatch_i_dh = mismatch_table(data, &file("tstacki", "dh"), terminal_dh);
    let mismatch_1n_37 = mismatch_table(data, &file("tstacki1n", "dg"), terminal_37);
    let mismatch_1n_dh = mismatch_table(data, &file("tstacki1n", "dh"), terminal_dh);
    let mismatch_23_37 = mismatch_table(data, &file("tstacki23", "dg"), terminal_37);
    let mismatch_23_dh = mismatch_table(data, &file("tstacki23", "dh"), terminal_dh);
    let mismatch_m_37 = terminal_mismatch_table(data, &file("tstackm", "dg"));
    let mismatch_m_dh = terminal_mismatch_table(data, &file("tstackm", "dh"));
    let mismatch_ext_37 = terminal_mismatch_table(data, &file("tstack", "dg"));
    let mismatch_ext_dh = terminal_mismatch_table(data, &file("tstack", "dh"));
    let (dangle5_37, dangle3_37) = dangle_tables(data, &file("dangle", "dg"));
    let (dangle5_dh, dangle3_dh) = dangle_tables(data, &file("dangle", "dh"));
    let int11_37 = int11_table(data, &file("int11", "dg"), terminal_37);
    let int11_dh = int11_table(data, &file("int11", "dh"), terminal_dh);
    let int21_37 = int21_table(data, &file("int21", "dg"), terminal_37);
    let int21_dh = int21_table(data, &file("int21", "dh"), terminal_dh);
    let int22_37 = int22_table(data, &file("int22", "dg"), terminal_37);
    let int22_dh = int22_table(data, &file("int22", "dh"), terminal_dh);
    let (internal_37, bulge_37, hairpin_37) = loop_tables(data, &file("loop", "dg"));
    let (internal_dh, bulge_dh, hairpin_dh) = loop_tables(data, &file("loop", "dh"));
    let misc_37 = miscellaneous_rows(data, &file("miscloop", "dg"));
    let misc_dh = miscellaneous_rows(data, &file("miscloop", "dh"));
    let triloops = special_loops(
        data,
        &file("triloop", "dg"),
        &file("triloop", "dh"),
        terminal_37,
        terminal_dh,
    );
    let tetraloops = special_loops(
        data,
        &file("tloop", "dg"),
        &file("tloop", "dh"),
        terminal_37,
        terminal_dh,
    );
    let hexaloops = special_loops(
        data,
        &file("hexaloop", "dg"),
        &file("hexaloop", "dh"),
        terminal_37,
        terminal_dh,
    );

    let mut generated = format!(
        "// Generated from RNAstructure 6.6 GPL-2.0-only {description} data tables.\n\
         // Source distribution SHA-256: 8a2904c4b9e16854a2aac3c6f3e510c844685f8cf330601e986d12f7d97dadc8.\n\
         // Integer values are 0.01 kcal/mol.\n"
    );
    for (name, values) in [
        ("STACK_37", stack_37),
        ("STACK_DH", stack_dh),
        ("MISMATCH_H_37", mismatch_h_37),
        ("MISMATCH_H_DH", mismatch_h_dh),
        ("MISMATCH_I_37", mismatch_i_37),
        ("MISMATCH_I_DH", mismatch_i_dh),
        ("MISMATCH_1N_37", mismatch_1n_37),
        ("MISMATCH_1N_DH", mismatch_1n_dh),
        ("MISMATCH_23_37", mismatch_23_37),
        ("MISMATCH_23_DH", mismatch_23_dh),
        ("MISMATCH_M_37", mismatch_m_37),
        ("MISMATCH_M_DH", mismatch_m_dh),
        ("MISMATCH_EXT_37", mismatch_ext_37),
        ("MISMATCH_EXT_DH", mismatch_ext_dh),
        ("DANGLE5_37", dangle5_37),
        ("DANGLE5_DH", dangle5_dh),
        ("DANGLE3_37", dangle3_37),
        ("DANGLE3_DH", dangle3_dh),
        ("INT11_37", int11_37),
        ("INT11_DH", int11_dh),
        ("INT21_37", int21_37),
        ("INT21_DH", int21_dh),
        ("INT22_37", int22_37),
        ("INT22_DH", int22_dh),
        ("HAIRPIN_37", hairpin_37),
        ("HAIRPIN_DH", hairpin_dh),
        ("BULGE_37", bulge_37),
        ("BULGE_DH", bulge_dh),
        ("INTERNAL_37", internal_37),
        ("INTERNAL_DH", internal_dh),
    ] {
        write_array(&mut generated, name, &values);
    }
    // Compact runtime order: `(g37, enthalpy)` for multiloop unpaired,
    // closing, and branch terms; then Ninio slope/enthalpy/maximum; then
    // duplex initiation and weak terminal-pair penalty.
    write_array(
        &mut generated,
        "ML_PARAMS",
        &[
            misc_37.multiloop[1],
            misc_dh.multiloop[1],
            misc_37.multiloop[0],
            misc_dh.multiloop[0],
            misc_37.multiloop[2],
            misc_dh.multiloop[2],
        ],
    );
    write_array(
        &mut generated,
        "NINIO",
        &[
            misc_37.ninio_slope,
            misc_dh.ninio_slope,
            misc_37.ninio_maximum,
        ],
    );
    write_array(
        &mut generated,
        "MISC",
        &[
            misc_37.duplex_initiation,
            misc_dh.duplex_initiation,
            terminal_pair_37,
            terminal_pair_dh,
        ],
    );
    generated.push_str(&format!("pub const LXC_37: f64 = {:?};\n", misc_37.lxc));
    write_special(&mut generated, "TRILOOPS", &triloops);
    write_special(&mut generated, "TETRALOOPS", &tetraloops);
    write_special(&mut generated, "HEXALOOPS", &hexaloops);

    let output = PathBuf::from(env::var_os("OUT_DIR").expect("OUT_DIR")).join(output_name);
    fs::write(output, generated).expect("write generated thermodynamic constants");
}

fn read(data: &Path, name: &str) -> String {
    fs::read_to_string(data.join(name)).unwrap_or_else(|error| panic!("read {name}: {error}"))
}

fn energy(token: &str) -> Option<i32> {
    if token == "." {
        Some(INF)
    } else {
        token
            .parse::<f64>()
            .ok()
            .map(|value| (value * 100.0).round() as i32)
    }
}

fn base(token: &str) -> Option<usize> {
    match token {
        "A" => Some(0),
        "C" => Some(1),
        "G" => Some(2),
        "U" | "T" => Some(3),
        _ => None,
    }
}

fn pair(top: usize, bottom: usize) -> Option<usize> {
    match (top, bottom) {
        (1, 2) => Some(0),
        (2, 1) => Some(1),
        (2, 3) => Some(2),
        (3, 2) => Some(3),
        (0, 3) => Some(4),
        (3, 0) => Some(5),
        _ => None,
    }
}

fn reverse_pair(value: usize) -> usize {
    [1, 0, 3, 2, 5, 4][value]
}

fn weak_penalty(pair: usize, penalty: i32) -> i32 {
    if pair >= 2 {
        penalty
    } else {
        0
    }
}

fn matrix4_rows(input: &str) -> Vec<[i32; 4]> {
    input
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            if tokens.len() != 5 || base(tokens[0]).is_none() {
                return None;
            }
            Some([
                energy(tokens[1])?,
                energy(tokens[2])?,
                energy(tokens[3])?,
                energy(tokens[4])?,
            ])
        })
        .collect()
}

fn matrix4(data: &Path, name: &str) -> Vec<[i32; 16]> {
    let rows = matrix4_rows(&read(data, name));
    assert_eq!(rows.len() % 4, 0, "incomplete matrix in {name}");
    rows.chunks_exact(4)
        .map(|chunk| {
            let mut matrix = [INF; 16];
            for row in 0..4 {
                matrix[row * 4..row * 4 + 4].copy_from_slice(&chunk[row]);
            }
            matrix
        })
        .collect()
}

fn stack_table(data: &Path, name: &str) -> Vec<i32> {
    let matrices = matrix4(data, name);
    assert_eq!(matrices.len(), 6, "unexpected stack matrix count");
    let mut output = vec![INF; 49];
    for (block, matrix) in matrices.iter().enumerate() {
        let outer = STACK_PAIR_ORDER[block];
        for top in 0..4 {
            for bottom in 0..4 {
                if let Some(inner) = pair(bottom, top) {
                    output[outer * 7 + inner] = matrix[top * 4 + bottom];
                }
            }
        }
    }
    output
}

fn mismatch_table(data: &Path, name: &str, terminal_penalty: i32) -> Vec<i32> {
    let matrices = matrix4(data, name);
    assert_eq!(matrices.len(), 6, "unexpected mismatch matrix count");
    let mut output = vec![0; 175];
    for (block, matrix) in matrices.iter().enumerate() {
        let outer = STACK_PAIR_ORDER[block];
        let terminal = weak_penalty(outer, terminal_penalty);
        for left in 0..4 {
            for right in 0..4 {
                let value = matrix[left * 4 + right];
                output[((outer * 5 + left + 1) * 5) + right + 1] =
                    if value >= INF { INF } else { value + terminal };
            }
        }
    }
    output
}

fn terminal_mismatch_table(data: &Path, name: &str) -> Vec<i32> {
    let matrices = matrix4(data, name);
    assert_eq!(
        matrices.len(),
        6,
        "unexpected terminal-mismatch matrix count"
    );
    let mut output = vec![0; 175];
    for (block, matrix) in matrices.iter().enumerate() {
        // In the RNAstructure terminal-stack files the closing pair and the
        // two flanking nucleotides are shown from the loop-facing side.  The
        // runtime indexes them from the helix-facing side, which reverses the
        // pair and transposes the 4x4 block.  Unlike hairpin and internal-loop
        // terminal mismatches, exterior/multiloop stem scoring adds the weak
        // terminal-pair penalty separately.
        let pair = reverse_pair(STACK_PAIR_ORDER[block]);
        for left in 0..4 {
            for right in 0..4 {
                output[((pair * 5 + left + 1) * 5) + right + 1] = matrix[right * 4 + left];
            }
        }
    }
    output
}

fn dangle_tables(data: &Path, name: &str) -> (Vec<i32>, Vec<i32>) {
    let rows = read(data, name)
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            if tokens.len() != 4 {
                return None;
            }
            Some([
                energy(tokens[0])?,
                energy(tokens[1])?,
                energy(tokens[2])?,
                energy(tokens[3])?,
            ])
        })
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 32, "unexpected dangle row count in {name}");
    let mut dangle5 = vec![0; 35];
    let mut dangle3 = vec![0; 35];
    for section in 0..2 {
        for top in 0..4 {
            for bottom in 0..4 {
                let Some(pair) = pair(bottom, top) else {
                    continue;
                };
                let values = rows[section * 16 + top * 4 + bottom];
                let target = if section == 0 {
                    &mut dangle3
                } else {
                    &mut dangle5
                };
                for nucleotide in 0..4 {
                    target[pair * 5 + nucleotide + 1] = values[nucleotide];
                }
            }
        }
    }
    (dangle5, dangle3)
}

fn int11_table(data: &Path, name: &str, terminal_penalty: i32) -> Vec<i32> {
    let matrices = matrix4(data, name);
    assert_eq!(matrices.len(), 36, "unexpected int11 matrix count");
    let mut output = vec![INF; 1225];
    for (block, matrix) in matrices.iter().enumerate() {
        let outer = INTERNAL_PAIR_ORDER[block / 6];
        let inner = reverse_pair(INTERNAL_PAIR_ORDER[block % 6]);
        let terminal =
            weak_penalty(outer, terminal_penalty) + weak_penalty(inner, terminal_penalty);
        for left in 0..4 {
            for right in 0..4 {
                let value = matrix[left * 4 + right];
                let index = ((outer * 7 + inner) * 5 + left + 1) * 5 + right + 1;
                output[index] = if value >= INF { INF } else { value + terminal };
            }
        }
    }
    output
}

fn int21_table(data: &Path, name: &str, terminal_penalty: i32) -> Vec<i32> {
    let matrices = matrix4(data, name);
    assert_eq!(matrices.len(), 144, "unexpected int21 matrix count");
    let mut output = vec![INF; 6125];
    for (block, matrix) in matrices.iter().enumerate() {
        let outer = INTERNAL_PAIR_ORDER[block / 24];
        let inner = reverse_pair(INTERNAL_PAIR_ORDER[(block / 4) % 6]);
        let third = block % 4 + 1;
        let terminal =
            weak_penalty(outer, terminal_penalty) + weak_penalty(inner, terminal_penalty);
        for first in 0..4 {
            for second in 0..4 {
                let value = matrix[first * 4 + second];
                // RNAstructure groups the 1x2-loop tables by the nucleotide
                // adjacent to the inner pair on the longer strand.  The
                // matrix columns are the nucleotide adjacent to the outer
                // pair, so the latter two runtime axes are the reverse of
                // their order in the source file.
                let index = (((outer * 7 + inner) * 5 + first + 1) * 5 + third) * 5 + second + 1;
                output[index] = if value >= INF { INF } else { value + terminal };
            }
        }
    }
    output
}

fn int22_rows(input: &str) -> Vec<([usize; 2], [i32; 16])> {
    input
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            if tokens.len() != 17 || tokens[0].len() != 2 {
                return None;
            }
            let label = tokens[0].as_bytes();
            let row = [
                base(std::str::from_utf8(&label[0..1]).ok()?)?,
                base(std::str::from_utf8(&label[1..2]).ok()?)?,
            ];
            let mut values = [INF; 16];
            for (index, token) in tokens[1..].iter().enumerate() {
                values[index] = energy(token)?;
            }
            Some((row, values))
        })
        .collect()
}

fn int22_table(data: &Path, name: &str, terminal_penalty: i32) -> Vec<i32> {
    let rows = int22_rows(&read(data, name));
    assert_eq!(rows.len(), 576, "unexpected int22 row count");
    let mut output = vec![INF; 9216];
    for block in 0..36 {
        let outer = INTERNAL_PAIR_ORDER[block / 6];
        let inner = reverse_pair(INTERNAL_PAIR_ORDER[block % 6]);
        let terminal =
            weak_penalty(outer, terminal_penalty) + weak_penalty(inner, terminal_penalty);
        for row_index in 0..16 {
            let (row, values) = &rows[block * 16 + row_index];
            for (column, value) in values.iter().copied().enumerate() {
                let second = column / 4;
                let third = column % 4;
                let index =
                    (((((outer * 6 + inner) * 4 + row[0]) * 4 + second) * 4 + third) * 4) + row[1];
                output[index] = if value >= INF { INF } else { value + terminal };
            }
        }
    }
    output
}

fn loop_tables(data: &Path, name: &str) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    let mut internal = vec![INF; 31];
    let mut bulge = vec![INF; 31];
    let mut hairpin = vec![INF; 31];
    for line in read(data, name).lines() {
        if line.trim_start().starts_with('#') {
            continue;
        }
        let tokens = line.split_whitespace().collect::<Vec<_>>();
        if tokens.len() != 4 {
            continue;
        }
        let Ok(size) = tokens[0].parse::<usize>() else {
            continue;
        };
        if size > 30 {
            continue;
        }
        internal[size] = energy(tokens[1]).expect("internal loop value");
        bulge[size] = energy(tokens[2]).expect("bulge loop value");
        hairpin[size] = energy(tokens[3]).expect("hairpin loop value");
    }
    (internal, bulge, hairpin)
}

struct MiscellaneousParameters {
    lxc: f64,
    ninio_maximum: i32,
    ninio_slope: i32,
    multiloop: [i32; 3],
    duplex_initiation: i32,
}

fn numeric_rows(data: &Path, name: &str) -> Vec<Vec<f64>> {
    read(data, name)
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let values = line
                .split_whitespace()
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
                .ok()?;
            (!values.is_empty()).then_some(values)
        })
        .collect()
}

fn centi(value: f64) -> i32 {
    (value * 100.0).round() as i32
}

fn miscellaneous_rows(data: &Path, name: &str) -> MiscellaneousParameters {
    let rows = numeric_rows(data, name);
    assert!(
        rows.len() >= 13,
        "unexpected miscellaneous row count in {name}"
    );
    assert_eq!(rows[2].len(), 4, "unexpected Ninio slope count");
    assert!(
        rows[2]
            .iter()
            .all(|value| (*value - rows[2][0]).abs() < 1.0e-12),
        "runtime requires a position-independent Ninio slope"
    );
    assert_eq!(rows[3].len(), 3, "unexpected multiloop term count");
    assert_eq!(rows[4].len(), 3, "unexpected efn2 multiloop term count");
    MiscellaneousParameters {
        lxc: rows[0][0],
        ninio_maximum: centi(rows[1][0]),
        ninio_slope: centi(rows[2][0]),
        multiloop: [centi(rows[3][0]), centi(rows[3][1]), centi(rows[3][2])],
        duplex_initiation: named_scalar_centi(data, name, "Intermolecular initiation"),
    }
}

/// Read the first scalar value following a named comment section.  DNA and
/// RNA `miscloop` files do not have identical optional sections, so indexing
/// the numeric rows would silently select the wrong quantity for DNA.
fn named_scalar_centi(data: &Path, name: &str, section: &str) -> i32 {
    let input = read(data, name);
    let mut found = false;
    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with('#') {
            let heading = trimmed.trim_start_matches('#').trim().to_ascii_lowercase();
            if heading.contains(&section.to_ascii_lowercase()) {
                found = true;
            }
            continue;
        }
        if found {
            if let Some(value) = trimmed
                .split_whitespace()
                .next()
                .and_then(|token| token.parse::<f64>().ok())
            {
                return centi(value);
            }
        }
    }
    panic!("missing {section:?} scalar in {name}")
}

fn terminal_pair_penalty(data: &Path, name: &str) -> i32 {
    let matrices = matrix4(data, name);
    assert_eq!(matrices.len(), 6, "unexpected terminal-pair table count");
    let penalty = matrices[0][0];
    for (block, matrix) in matrices.iter().enumerate() {
        let expected = if STACK_PAIR_ORDER[block] >= 2 {
            penalty
        } else {
            0
        };
        assert!(
            matrix.iter().all(|&value| value == expected),
            "nonuniform terminal-pair table is not supported"
        );
    }
    penalty
}

fn special_rows(data: &Path, name: &str) -> Vec<(String, i32)> {
    read(data, name)
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let tokens = line.split_whitespace().collect::<Vec<_>>();
            if tokens.len() != 2 {
                return None;
            }
            Some((tokens[0].to_string(), energy(tokens[1])?))
        })
        .collect()
}

fn special_loops(
    data: &Path,
    free_energy: &str,
    enthalpy: &str,
    terminal_penalty_37: i32,
    terminal_penalty_dh: i32,
) -> Vec<(String, i32, i32)> {
    let dg = special_rows(data, free_energy);
    let dh = special_rows(data, enthalpy);
    assert_eq!(dg.len(), dh.len(), "special-loop row count mismatch");
    dg.into_iter()
        .zip(dh)
        .map(|((sequence, dg), (enthalpy_sequence, dh))| {
            assert_eq!(
                sequence, enthalpy_sequence,
                "special-loop sequence mismatch"
            );
            let bases = sequence.as_bytes();
            let pair = pair(
                base(std::str::from_utf8(&bases[0..1]).expect("ASCII base")).expect("base"),
                base(std::str::from_utf8(&bases[bases.len() - 1..]).expect("ASCII closing base"))
                    .expect("closing base"),
            )
            .expect("canonical special-loop closing pair");
            (
                sequence,
                dg + weak_penalty(pair, terminal_penalty_37),
                dh + weak_penalty(pair, terminal_penalty_dh),
            )
        })
        .collect()
}

fn write_array(output: &mut String, name: &str, values: &[i32]) {
    output.push_str(&format!(
        "pub static {name}: [i32; {}] = {values:?};\n",
        values.len()
    ));
}

fn write_special(output: &mut String, name: &str, values: &[(String, i32, i32)]) {
    output.push_str(&format!("pub static {name}: &[(&str, i32, i32)] = &[\n"));
    for (sequence, dg, dh) in values {
        output.push_str(&format!("    ({sequence:?}, {dg}, {dh}),\n"));
    }
    output.push_str("];\n");
}
