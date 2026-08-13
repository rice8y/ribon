use crate::structure::{parse_structure, Pair, ParsedStructure, RnaError};
use serde::Serialize;
use std::f64::consts::{FRAC_PI_2, PI};
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LayoutKind {
    Naview,
    Simple,
    Circular,
    Turtle,
    Puzzler,
    Linear,
}

impl FromStr for LayoutKind {
    type Err = RnaError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "naview" | "radial" => Ok(Self::Naview),
            "simple" | "rna-plot" | "rnaplot" => Ok(Self::Simple),
            "circular" | "circle" => Ok(Self::Circular),
            "turtle" | "rnaturtle" => Ok(Self::Turtle),
            "puzzler" => Ok(Self::Puzzler),
            "linear" | "line" => Ok(Self::Linear),
            _ => Err(RnaError::InvalidOption(format!(
                "unknown layout {value:?}; expected naview, simple, circular, turtle, puzzler, or linear"
            ))),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LayoutPair {
    pub i: usize,
    pub j: usize,
    pub level: usize,
    pub canonical: bool,
    pub interstrand: bool,
}

/// Circular backbone segment emitted by RNApuzzler/RNAturtle. Coordinates are
/// normalized with the nucleotide positions; radii are separated by axis so
/// Typst can retain the reference arc under non-square figure dimensions.
#[derive(Clone, Debug, Serialize)]
pub struct LayoutArc {
    pub from: usize,
    pub to: usize,
    pub center_x: f64,
    pub center_y: f64,
    pub radius_x: f64,
    pub radius_y: f64,
    pub start_degrees: f64,
    pub end_degrees: f64,
    pub clockwise: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct LayoutResult {
    pub sequence: String,
    pub structure: String,
    pub method: &'static str,
    pub algorithm: &'static str,
    pub points: Vec<Point>,
    pub backbone_arcs: Vec<LayoutArc>,
    pub pairs: Vec<LayoutPair>,
    pub strand_breaks: Vec<usize>,
    pub aspect_ratio: f64,
    pub crossings: usize,
}

pub fn layout_structure(
    sequence: &str,
    structure: &str,
    kind: LayoutKind,
) -> Result<LayoutResult, RnaError> {
    let parsed = parse_structure(sequence, structure)?;
    let (raw, raw_arcs) = match kind {
        LayoutKind::Simple => (simple_coordinates(&parsed), Vec::new()),
        LayoutKind::Circular => (circular_coordinates(parsed.length), Vec::new()),
        LayoutKind::Linear => (linear_coordinates(&parsed), Vec::new()),
        LayoutKind::Naview => (naview_coordinates(&parsed), Vec::new()),
        LayoutKind::Turtle => {
            let (points, arcs) = crate::turtle::coordinates_and_arcs(&parsed);
            (points, raw_turtle_arcs(arcs))
        }
        LayoutKind::Puzzler => puzzler_coordinates(&parsed),
    };
    let (points, backbone_arcs, aspect_ratio) = normalize_points_and_arcs(raw, raw_arcs);
    let crossings = count_segment_crossings(&points, &parsed.pairs, &parsed.strand_breaks);
    let (method, algorithm) = match kind {
        LayoutKind::Simple => ("simple", "regular-bond radial layout"),
        LayoutKind::Circular => ("circular", "equidistant circular"),
        LayoutKind::Linear if parsed.strand_breaks.is_empty() => ("linear", "linear arc diagram"),
        LayoutKind::Linear => ("linear", "strand-aware antiparallel linear diagram"),
        LayoutKind::Turtle => ("turtle", "RNAturtle affine loop geometry"),
        LayoutKind::Puzzler => (
            "puzzler",
            "RNAturtle geometry with deterministic bounded collision reduction",
        ),
        LayoutKind::Naview => (
            "naview",
            "independently authored classic NAView modified-radial geometry",
        ),
    };
    Ok(LayoutResult {
        sequence: parsed.sequence,
        structure: parsed.structure,
        method,
        algorithm,
        points,
        backbone_arcs,
        pairs: parsed
            .pairs
            .iter()
            .map(|p| LayoutPair {
                i: p.i,
                j: p.j,
                level: p.level,
                canonical: p.canonical,
                interstrand: strand_index(p.i, &parsed.strand_breaks)
                    != strand_index(p.j, &parsed.strand_breaks),
            })
            .collect(),
        strand_breaks: parsed.strand_breaks,
        aspect_ratio,
        crossings,
    })
}

#[derive(Clone, Copy)]
struct RawBackboneArc {
    segment: usize,
    center_x: f64,
    center_y: f64,
    radius: f64,
    start_degrees: f64,
    end_degrees: f64,
    clockwise: bool,
}

fn puzzler_coordinates(parsed: &ParsedStructure) -> (Vec<Point>, Vec<RawBackboneArc>) {
    let (points, arcs) = crate::puzzler::coordinates_and_arcs(parsed);
    (points, raw_turtle_arcs(arcs))
}

fn raw_turtle_arcs(arcs: Vec<crate::turtle::TurtleArc>) -> Vec<RawBackboneArc> {
    arcs.into_iter()
        .map(|arc| RawBackboneArc {
            segment: arc.segment,
            center_x: arc.center_x,
            center_y: arc.center_y,
            radius: arc.radius,
            start_degrees: arc.start_degrees,
            end_degrees: arc.end_degrees,
            clockwise: arc.clockwise,
        })
        .collect()
}

fn circular_coordinates(n: usize) -> Vec<Point> {
    let step = 2.0 * PI / n as f64;
    (0..n)
        .map(|i| Point {
            x: (i as f64 * step - FRAC_PI_2).cos(),
            y: (i as f64 * step - FRAC_PI_2).sin(),
        })
        .collect()
}

fn strand_index(position: usize, strand_breaks: &[usize]) -> usize {
    strand_breaks.partition_point(|&strand_end| strand_end < position)
}

fn strand_ranges(length: usize, strand_breaks: &[usize]) -> Vec<(usize, usize)> {
    let mut ranges = Vec::with_capacity(strand_breaks.len() + 1);
    let mut start = 0usize;
    for &strand_end in strand_breaks {
        ranges.push((start, strand_end));
        start = strand_end;
    }
    ranges.push((start, length));
    ranges
}

fn linear_coordinates(parsed: &ParsedStructure) -> Vec<Point> {
    if parsed.strand_breaks.is_empty() {
        return (0..parsed.length)
            .map(|i| Point {
                x: i as f64,
                y: 0.0,
            })
            .collect();
    }

    let ranges = strand_ranges(parsed.length, &parsed.strand_breaks);
    let maximum_length = ranges
        .iter()
        .map(|(start, end)| end - start)
        .max()
        .unwrap_or(1);
    let mut points = vec![Point::default(); parsed.length];
    for (strand, &(start, end)) in ranges.iter().enumerate() {
        let length = end - start;
        let offset = (maximum_length - length) as f64 / 2.0;
        for local in 0..length {
            let x = if strand % 2 == 0 {
                offset + local as f64
            } else {
                offset + (length - 1 - local) as f64
            };
            points[start + local] = Point {
                x,
                y: strand as f64,
            };
        }
    }
    points
}

/// Reimplementation of the classic RNAplot simple-radial polygon rule. It is
/// kept as a separate algorithm because reference RNAplot exposes it as layout
/// type 0 and its regular bond length is valuable for regression comparisons.
fn simple_coordinates(parsed: &ParsedStructure) -> Vec<Point> {
    let n = parsed.length;
    let mut table = vec![0usize; n + 1];
    table[0] = n;
    for pair in &parsed.pairs {
        // The classic radial algorithm is a planar secondary-structure
        // algorithm. Keep the first bracket level as the planar scaffold and
        // render higher levels as additional pair edges.
        if pair.level == 0 {
            table[pair.i] = pair.j;
            table[pair.j] = pair.i;
        }
    }
    let mut angles = vec![0.0f64; n + 6];
    radial_loop_angles(&table, 0, n, &mut angles);
    let mut points = vec![Point { x: 0.0, y: 0.0 }; n];
    let mut alpha = 0.0f64;
    for index in 1..n {
        points[index] = Point {
            x: points[index - 1].x + alpha.cos(),
            y: points[index - 1].y + alpha.sin(),
        };
        alpha += PI - angles[index + 1];
    }
    points
}

fn radial_loop_angles(table: &[usize], mut i: usize, mut j: usize, angles: &mut [f64]) {
    let mut vertices = 2usize;
    let mut boundaries = Vec::new();
    let previous = i as isize - 1;
    j += 1;
    while i != j {
        let mate = table[i];
        if mate == 0 || i == 0 {
            i += 1;
            vertices += 1;
        } else {
            vertices += 2;
            let (start_left, start_right) = (i, mate);
            boundaries.push(start_left);
            boundaries.push(start_right);
            i = mate + 1;
            let (mut left, mut right) = (start_left, start_right);
            let mut ladder = 0usize;
            loop {
                left += 1;
                right = right.saturating_sub(1);
                ladder += 1;
                if left >= table.len() || table[left] != right || table[left] <= left {
                    break;
                }
            }
            if ladder >= 2 {
                let fill = ladder - 2;
                angles[start_left + 1 + fill] += FRAC_PI_2;
                angles[start_right - 1 - fill] += FRAC_PI_2;
                angles[start_left] += FRAC_PI_2;
                angles[start_right] += FRAC_PI_2;
                if ladder > 2 {
                    for offset in 1..=fill {
                        angles[start_left + offset] = PI;
                        angles[start_right - offset] = PI;
                    }
                }
            }
            if left <= right {
                radial_loop_angles(table, left, right, angles);
            }
        }
    }
    let polygon_angle = PI * (vertices as f64 - 2.0) / vertices as f64;
    boundaries.push(j);
    let mut begin = previous.max(0) as usize;
    let mut cursor = 0usize;
    while cursor < boundaries.len() {
        let end = boundaries[cursor];
        for position in begin..=end.min(angles.len() - 1) {
            angles[position] += polygon_angle;
        }
        cursor += 1;
        if cursor >= boundaries.len() {
            break;
        }
        begin = boundaries[cursor];
        cursor += 1;
    }
}

fn naview_coordinates(parsed: &ParsedStructure) -> Vec<Point> {
    crate::naview::coordinates(parsed)
}

fn normalize_points_and_arcs(
    mut points: Vec<Point>,
    arcs: Vec<RawBackboneArc>,
) -> (Vec<Point>, Vec<LayoutArc>, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for point in &points {
        min_x = min_x.min(point.x);
        max_x = max_x.max(point.x);
        min_y = min_y.min(point.y);
        max_y = max_y.max(point.y);
    }
    let raw_width = max_x - min_x;
    let raw_height = max_y - min_y;
    let width = raw_width.max(1.0);
    let height = raw_height.max(1.0);
    for point in &mut points {
        point.x = (point.x - min_x) / width;
        point.y = if raw_height.abs() < 1.0e-12 {
            0.5
        } else {
            (point.y - min_y) / height
        };
    }
    let arcs = arcs
        .into_iter()
        .filter(|arc| arc.segment + 1 < points.len())
        .map(|arc| LayoutArc {
            from: arc.segment + 1,
            to: arc.segment + 2,
            center_x: (arc.center_x - min_x) / width,
            center_y: if raw_height.abs() < 1.0e-12 {
                0.5
            } else {
                (arc.center_y - min_y) / height
            },
            radius_x: arc.radius / width,
            radius_y: arc.radius / height,
            start_degrees: arc.start_degrees,
            end_degrees: arc.end_degrees,
            clockwise: arc.clockwise,
        })
        .collect();
    (points, arcs, (width / height).clamp(0.08, 20.0))
}

fn orientation(a: Point, b: Point, c: Point) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

fn proper_intersection(a: Point, b: Point, c: Point, d: Point) -> bool {
    let ab_c = orientation(a, b, c);
    let ab_d = orientation(a, b, d);
    let cd_a = orientation(c, d, a);
    let cd_b = orientation(c, d, b);
    ab_c * ab_d < -1.0e-10 && cd_a * cd_b < -1.0e-10
}

pub(crate) fn count_segment_crossings(
    points: &[Point],
    pairs: &[Pair],
    strand_breaks: &[usize],
) -> usize {
    let mut segments = Vec::new();
    for i in 0..points.len().saturating_sub(1) {
        if !strand_breaks.contains(&(i + 1)) {
            segments.push((i, i + 1));
        }
    }
    segments.extend(pairs.iter().map(|pair| (pair.i - 1, pair.j - 1)));
    let mut count = 0usize;
    for index in 0..segments.len() {
        let (a, b) = segments[index];
        for &(c, d) in &segments[index + 1..] {
            if a == c || a == d || b == c || b == d {
                continue;
            }
            if proper_intersection(points[a], points[b], points[c], points[d]) {
                count += 1;
            }
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_layouts_return_finite_normalized_coordinates() {
        for kind in [
            LayoutKind::Simple,
            LayoutKind::Naview,
            LayoutKind::Circular,
            LayoutKind::Linear,
            LayoutKind::Turtle,
            LayoutKind::Puzzler,
        ] {
            let result = layout_structure("GGGAAACCC", "(((...)))", kind).unwrap();
            assert_eq!(result.points.len(), 9);
            assert!(result.points.iter().all(|p| {
                p.x.is_finite()
                    && p.y.is_finite()
                    && (0.0..=1.0).contains(&p.x)
                    && (0.0..=1.0).contains(&p.y)
            }));
        }
    }

    #[test]
    fn circular_layout_is_square() {
        let result = layout_structure("ACGUACGU", "........", LayoutKind::Circular).unwrap();
        assert!((result.aspect_ratio - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn layout_preserves_dna_thymine() {
        let result = layout_structure("atg&t", "...&.", LayoutKind::Linear).unwrap();
        assert_eq!(result.sequence, "ATGT");
        assert_eq!(result.strand_breaks, vec![3]);
    }

    #[test]
    fn multi_strand_linear_layout_uses_antiparallel_rows_and_rungs() {
        let result = layout_structure("GGGG&CCCC", "((((&))))", LayoutKind::Linear).unwrap();
        assert_eq!(result.algorithm, "strand-aware antiparallel linear diagram");
        assert_eq!(result.strand_breaks, vec![4]);
        assert!((result.aspect_ratio - 3.0).abs() < 1.0e-12);

        let top = &result.points[..4];
        let bottom = &result.points[4..];
        assert!(top.windows(2).all(|window| window[0].x < window[1].x));
        assert!(bottom.windows(2).all(|window| window[0].x > window[1].x));
        assert!(top.iter().all(|point| (point.y - 0.0).abs() < 1.0e-12));
        assert!(bottom.iter().all(|point| (point.y - 1.0).abs() < 1.0e-12));
        assert!(result.pairs.iter().all(|pair| pair.interstrand));
        for pair in &result.pairs {
            let first = result.points[pair.i - 1];
            let second = result.points[pair.j - 1];
            assert!((first.x - second.x).abs() < 1.0e-12);
        }
        assert_eq!(result.crossings, 0);
    }

    #[test]
    fn uneven_multi_strand_linear_layout_is_centered_and_finite() {
        let result = layout_structure("AAAA&CC&GGG", "....&..&...", LayoutKind::Linear).unwrap();
        assert_eq!(result.strand_breaks, vec![4, 6]);
        assert!(result.points.iter().all(|point| {
            point.x.is_finite()
                && point.y.is_finite()
                && (0.0..=1.0).contains(&point.x)
                && (0.0..=1.0).contains(&point.y)
        }));
        assert!((result.points[4].x - 2.0 / 3.0).abs() < 1.0e-12);
        assert!((result.points[5].x - 1.0 / 3.0).abs() < 1.0e-12);
    }

    #[test]
    fn multi_strand_linear_pairs_distinguish_arcs_from_rungs() {
        let result = layout_structure("NNNN&NNNN", "([.)&..].", LayoutKind::Linear).unwrap();
        assert_eq!(result.pairs.len(), 2);
        assert!(
            !result
                .pairs
                .iter()
                .find(|pair| pair.i == 1 && pair.j == 4)
                .unwrap()
                .interstrand
        );
        assert!(
            result
                .pairs
                .iter()
                .find(|pair| pair.i == 2 && pair.j == 7)
                .unwrap()
                .interstrand
        );
    }

    #[test]
    fn naview_uses_the_round_bracket_level_as_its_planar_scaffold() {
        let pseudoknot =
            layout_structure("NNNNNNNNNNNN", "(([[..))..]]", LayoutKind::Naview).unwrap();
        let scaffold =
            layout_structure("NNNNNNNNNNNN", "((....))....", LayoutKind::Naview).unwrap();
        for (actual, expected) in pseudoknot.points.iter().zip(&scaffold.points) {
            assert!((actual.x - expected.x).abs() < 1.0e-12);
            assert!((actual.y - expected.y).abs() < 1.0e-12);
        }
        assert_eq!(pseudoknot.pairs.len(), 4);
        assert!(pseudoknot.crossings > 0);
    }

    #[test]
    fn naview_keeps_a_planar_five_way_loop_crossing_free() {
        let structure = "((.((...)).((...)).((...)).((...)).((...)).))";
        let result =
            layout_structure(&"N".repeat(structure.len()), structure, LayoutKind::Naview).unwrap();
        assert_eq!(result.crossings, 0);
        assert!(result
            .points
            .iter()
            .all(|point| point.x.is_finite() && point.y.is_finite()));
    }

    #[test]
    fn turtle_and_puzzler_emit_geometrically_exact_backbone_arcs() {
        let structure = "((.((...)).((....)).))";
        for kind in [LayoutKind::Turtle, LayoutKind::Puzzler] {
            let result = layout_structure(&"N".repeat(structure.len()), structure, kind).unwrap();
            assert!(!result.backbone_arcs.is_empty());
            for arc in &result.backbone_arcs {
                assert_eq!(arc.to, arc.from + 1);
                assert!(arc.radius_x.is_finite() && arc.radius_x > 0.0);
                assert!(arc.radius_y.is_finite() && arc.radius_y > 0.0);
                for point in [&result.points[arc.from - 1], &result.points[arc.to - 1]] {
                    let ellipse = ((point.x - arc.center_x) / arc.radius_x).powi(2)
                        + ((point.y - arc.center_y) / arc.radius_y).powi(2);
                    assert!((ellipse - 1.0).abs() < 1.0e-8, "ellipse={ellipse}");
                }
            }
        }
    }
}
