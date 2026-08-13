//! Independent Rust implementation of the RNAturtle affine loop layout.
//!
//! The implementation follows the published/configuration geometry used by
//! RNAplot: paired chords are 35 units, backbone chords are 25 units, loops
//! are placed on the smallest circle satisfying those chord lengths, and the
//! resulting turn/distance representation is walked as a turtle path.  No C
//! code is linked into the package.

use crate::layout::Point;
use crate::structure::ParsedStructure;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const PAIRED: f64 = 35.0;
const UNPAIRED: f64 = 25.0;
const EPSILON: f64 = 1.0e-10;

#[derive(Clone, Debug)]
struct ArcConfig {
    angle: f64,
    segments: usize,
}

#[derive(Clone, Debug)]
struct LoopConfig {
    radius: f64,
    arcs: Vec<ArcConfig>,
}

#[derive(Clone, Debug)]
struct BaseInfo {
    angle: f64,
    distance: f64,
    arc_radius: Option<f64>,
    config: Option<LoopConfig>,
}

impl Default for BaseInfo {
    fn default() -> Self {
        Self {
            angle: 0.0,
            distance: UNPAIRED,
            arc_radius: None,
            config: None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct TurtleArc {
    pub segment: usize,
    pub center_x: f64,
    pub center_y: f64,
    pub radius: f64,
    pub start_degrees: f64,
    pub end_degrees: f64,
    pub clockwise: bool,
}

#[cfg(test)]
pub(crate) fn coordinates(parsed: &ParsedStructure) -> Vec<Point> {
    coordinates_and_arcs(parsed).0
}

pub(crate) fn coordinates_and_arcs(parsed: &ParsedStructure) -> (Vec<Point>, Vec<TurtleArc>) {
    let table = pair_table(parsed);
    let n = parsed.length;
    if n <= 1 {
        return (vec![Point { x: 0.0, y: 0.0 }; n], Vec::new());
    }
    let mut info = vec![BaseInfo::default(); n + 1];
    generate_configs(&table, &mut info);
    compute_affine(&table, &mut info);

    let mut points = vec![Point::default(); n];
    let mut heading = 0.0;
    for index in 1..n {
        heading -= info[index + 1].angle;
        points[index] = Point {
            x: points[index - 1].x + info[index].distance * heading.cos(),
            y: points[index - 1].y + info[index].distance * heading.sin(),
        };
    }
    let arcs = refit_arcs(
        &points,
        &(1..n)
            .filter_map(|index| {
                info[index].arc_radius.map(|radius| TurtleArc {
                    segment: index - 1,
                    center_x: 0.0,
                    center_y: 0.0,
                    radius,
                    start_degrees: 0.0,
                    end_degrees: 0.0,
                    clockwise: true,
                })
            })
            .collect::<Vec<_>>(),
    );
    (points, arcs)
}

pub(crate) fn refit_arcs(points: &[Point], templates: &[TurtleArc]) -> Vec<TurtleArc> {
    templates
        .iter()
        .filter_map(|template| {
            let from = *points.get(template.segment)?;
            let to = *points.get(template.segment + 1)?;
            let dx = to.x - from.x;
            let dy = to.y - from.y;
            let chord = dx.hypot(dy);
            if chord < EPSILON {
                return None;
            }
            let radius = template.radius.max(chord * 0.5 + EPSILON);
            let half = chord * 0.5;
            let height = (radius * radius - half * half).max(0.0).sqrt();
            let middle_x = (from.x + to.x) * 0.5;
            let middle_y = (from.y + to.y) * 0.5;
            // RNAturtle's global direction is clockwise. The circle center is
            // therefore on the right side of the directed chord.
            let center_x = middle_x + dy / chord * height;
            let center_y = middle_y - dx / chord * height;
            Some(TurtleArc {
                segment: template.segment,
                center_x,
                center_y,
                radius,
                start_degrees: (from.y - center_y).atan2(from.x - center_x).to_degrees(),
                end_degrees: (to.y - center_y).atan2(to.x - center_x).to_degrees(),
                clockwise: true,
            })
        })
        .collect()
}

fn pair_table(parsed: &ParsedStructure) -> Vec<usize> {
    let mut table = vec![0; parsed.length + 1];
    table[0] = parsed.length;
    for pair in &parsed.pairs {
        if pair.level == 0 {
            table[pair.i] = pair.j;
            table[pair.j] = pair.i;
        }
    }
    table
}

fn generate_configs(table: &[usize], info: &mut [BaseInfo]) {
    let n = table[0];
    let mut i = 1;
    while i < n {
        if table[i] == 0 {
            i += 1;
        } else if table[i] > i {
            generate_stem_config(i, table, info);
            i = table[i] + 1;
        } else {
            i += 1;
        }
    }
}

fn generate_stem_config(mut i: usize, table: &[usize], info: &mut [BaseInfo]) {
    while i + 1 < table.len() && table[i + 1] + 1 == table[i] {
        i += 1;
    }
    generate_loop_config(i, table, info);
}

fn generate_loop_config(start: usize, table: &[usize], info: &mut [BaseInfo]) {
    let end = table[start];
    if end <= start {
        return;
    }
    let mut unpaired = 0usize;
    let mut children = Vec::new();
    let mut i = start + 1;
    while i < end {
        if table[i] == 0 {
            unpaired += 1;
            i += 1;
        } else if table[i] > i {
            children.push(i);
            i = table[i] + 1;
        } else {
            i += 1;
        }
    }

    if children.len() == 1 && unpaired == 1 {
        generate_stem_config(children[0], table, info);
        return;
    }

    let stems = children.len() + 1;
    let backbone_segments = unpaired + stems;
    let radius = solve_radius(PAIRED, UNPAIRED, stems, backbone_segments, TAU);
    let angle_paired = chord_angle(PAIRED, radius);
    let angle_unpaired = chord_angle(UNPAIRED, radius);
    let mut arcs = Vec::with_capacity(stems);
    let mut cursor = start + 1;
    let mut arc_unpaired = 0usize;
    while cursor <= end {
        if table[cursor] == 0 {
            arc_unpaired += 1;
            cursor += 1;
        } else {
            let segments = arc_unpaired + 1;
            arcs.push(ArcConfig {
                angle: angle_paired + segments as f64 * angle_unpaired,
                segments,
            });
            if cursor == end {
                break;
            }
            arc_unpaired = 0;
            cursor = table[cursor] + 1;
        }
    }
    info[start].config = Some(LoopConfig { radius, arcs });
    for child in children {
        generate_stem_config(child, table, info);
    }
}

fn chord_angle(chord: f64, radius: f64) -> f64 {
    2.0 * (chord / (2.0 * radius)).clamp(-1.0, 1.0).asin()
}

fn solve_radius(paired: f64, unpaired: f64, stems: usize, segments: usize, angle: f64) -> f64 {
    let residual = |radius: f64| {
        stems as f64 * (paired / (2.0 * radius)).clamp(-1.0, 1.0).asin()
            + segments as f64 * (unpaired / (2.0 * radius)).clamp(-1.0, 1.0).asin()
            - angle * 0.5
    };
    let mut lower = (paired.max(unpaired) * 0.5) * (1.0 + f64::EPSILON.sqrt());
    let mut upper = lower * 2.0;
    while residual(upper) > 0.0 {
        upper *= 2.0;
    }
    while upper - lower > 1.0e-12 * upper.max(1.0) {
        let midpoint = lower + (upper - lower) * 0.5;
        if residual(midpoint) > 0.0 {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }
    lower + (upper - lower) * 0.5
}

fn compute_affine(table: &[usize], info: &mut [BaseInfo]) {
    let n = table[0];
    info[0].angle = 0.0;
    if n >= 2 {
        info[1].angle = 0.0;
        info[2].angle = 0.0;
    }
    let direction = -1.0;
    let mut current = 1usize;
    let mut exterior_runs = 0usize;
    while current < n {
        if table[current] == 0 {
            current = handle_exterior(table, current, info, direction);
            exterior_runs += 1;
        }
        if current >= n {
            break;
        }
        if current > 1
            && table[current] != 0
            && table[current - 1] != 0
            && table[current] != table[current - 1] + 1
        {
            if current == 1 && exterior_runs == 0 {
                info[0].angle = -FRAC_PI_2;
                info[1].angle = -FRAC_PI_2;
                info[2].angle = -FRAC_PI_2;
            } else if current < n {
                info[current].angle += direction * FRAC_PI_2;
                info[current + 1].distance = UNPAIRED;
                info[current + 1].angle += direction * FRAC_PI_2;
                exterior_runs += 1;
            }
        } else if current == 1 && exterior_runs == 0 {
            info[0].angle = -FRAC_PI_2;
            info[1].angle = -FRAC_PI_2;
            info[2].angle = -FRAC_PI_2;
        }
        handle_stem(table, current, info, direction);
        current = table[current] + 1;
        if current == n {
            current = handle_exterior(table, current, info, direction);
        }
    }
}

fn handle_exterior(
    table: &[usize],
    mut current: usize,
    info: &mut [BaseInfo],
    direction: f64,
) -> usize {
    let n = table[0];
    if current > 1 {
        info[current].angle += direction * FRAC_PI_2;
    }
    while current < n && table[current] == 0 {
        info[current + 1].angle = 0.0;
        current += 1;
    }
    if current < n {
        info[current + 1].angle = direction * FRAC_PI_2;
    }
    current
}

fn handle_stem(table: &[usize], mut i: usize, info: &mut [BaseInfo], direction: f64) {
    if i == 0 || i >= table.len() || table[i] <= i {
        return;
    }
    let end = table[i] + 1;
    i += 1;
    while i < table.len()
        && table[i] > 0
        && (table[i] == end.saturating_sub(1) || table[i] + 1 == table[i - 1])
    {
        if i + 1 < info.len() {
            info[i + 1].angle = 0.0;
        }
        i += 1;
    }
    if i > 0 && table.get(i).copied() != Some(end.saturating_sub(1)) {
        i -= 1;
        handle_loop(i, table, info, direction);
    }
}

fn handle_loop(mut i: usize, table: &[usize], info: &mut [BaseInfo], direction: f64) {
    let start = i;
    let end = table[i];
    if end <= i {
        return;
    }
    let (stems, unpaired) = loop_counts(start, table);
    let child = only_child(start, table);
    if stems == 2 && unpaired == 1 {
        let length = UNPAIRED;
        let alpha = (UNPAIRED / (2.0 * length)).acos();
        if table[i + 1] == 0 {
            info[i + 1].angle += direction * alpha;
            i += 1;
            info[i + 1].angle = -direction * alpha * 2.0;
            i += 1;
            if i + 1 < info.len() {
                info[i + 1].angle = direction * alpha;
            }
            if let Some(child) = child {
                handle_stem(table, child, info, direction);
            }
        } else if let Some(child) = child {
            info[i + 1].angle += 0.0;
            i += 1;
            if i + 1 < info.len() {
                info[i + 1].angle += 0.0;
            }
            if i + 2 < info.len() {
                info[i + 2].angle += 0.0;
            }
            handle_stem(table, child, info, direction);
            i = table[child];
            if i + 1 < info.len() {
                info[i + 1].angle += direction * alpha;
            }
            i += 1;
            if i + 1 < info.len() {
                info[i + 1].angle = -direction * alpha * 2.0;
            }
            i += 1;
            if i + 1 < info.len() {
                info[i + 1].angle = direction * alpha;
            }
        }
        return;
    }

    let Some(config) = info[start].config.clone() else {
        return;
    };
    if config.arcs.is_empty() {
        return;
    }
    let paired_angle = chord_angle(PAIRED, config.radius);
    let mut arc_index = 0usize;
    let mut geometry = arc_geometry(&config, arc_index, paired_angle);
    arc_index += 1;
    info[i + 1].angle += direction * (PI - geometry.1);
    info[i].distance = geometry.0;
    info[i].arc_radius = Some(config.radius);
    i += 1;
    let mut child_stem = false;
    while i < end {
        if table[i] == 0 {
            info[i + 1].angle = -direction * (geometry.2 - PI);
            info[i].distance = geometry.0;
            info[i].arc_radius = Some(config.radius);
            i += 1;
        } else if table[i] > i {
            info[i + 1].angle = direction * (PI - geometry.1);
            child_stem = true;
            handle_stem(table, i, info, direction);
            i = table[i];
        } else {
            if child_stem && arc_index < config.arcs.len() {
                child_stem = false;
                geometry = arc_geometry(&config, arc_index, paired_angle);
                arc_index += 1;
            }
            if i + 1 < info.len() {
                info[i + 1].angle += direction * (PI - geometry.1);
            }
            info[i].distance = geometry.0;
            info[i].arc_radius = Some(config.radius);
            i += 1;
        }
    }
    if i + 1 < info.len() {
        info[i + 1].angle = direction * (PI - geometry.1);
    }
}

fn arc_geometry(config: &LoopConfig, index: usize, paired_angle: f64) -> (f64, f64, f64) {
    let arc = &config.arcs[index.min(config.arcs.len() - 1)];
    let backbone_angle = (arc.angle - paired_angle) / arc.segments.max(1) as f64;
    let distance = (2.0 * config.radius * config.radius * (1.0 - backbone_angle.cos())).sqrt();
    let delta_pair_backbone = 0.5 * (PI + paired_angle + backbone_angle);
    let delta_backbone = PI + backbone_angle;
    (distance, delta_pair_backbone, delta_backbone)
}

fn loop_counts(start: usize, table: &[usize]) -> (usize, usize) {
    let end = table[start];
    let mut stems = 1usize;
    let mut unpaired = 0usize;
    let mut i = start + 1;
    while i < end {
        if table[i] == 0 {
            unpaired += 1;
            i += 1;
        } else if table[i] > i {
            stems += 1;
            i = table[i] + 1;
        } else {
            i += 1;
        }
    }
    (stems, unpaired)
}

fn only_child(start: usize, table: &[usize]) -> Option<usize> {
    let end = table[start];
    let mut found = None;
    let mut i = start + 1;
    while i < end {
        if table[i] > i {
            if found.is_some() {
                return None;
            }
            found = Some(i);
            i = table[i] + 1;
        } else {
            i += 1;
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::parse_structure;

    #[test]
    fn turtle_coordinates_are_finite_for_diverse_loops() {
        for structure in [
            "(((...)))",
            "((..((...))...))",
            "((.((...)).((...)).((...)).))",
            "...(((...)))....(((....)))..",
        ] {
            let parsed = parse_structure(&"N".repeat(structure.len()), structure).unwrap();
            let points = coordinates(&parsed);
            assert_eq!(points.len(), structure.len());
            assert!(points
                .iter()
                .all(|point| point.x.is_finite() && point.y.is_finite()));
        }
    }
}
