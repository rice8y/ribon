//! Deterministic collision optimization on top of RNAturtle coordinates.
//!
//! RNApuzzler changes loop configurations to remove intersections.  This Rust
//! implementation keeps the same starting geometry and optimizes only the
//! Cartesian embedding with explicit backbone/base-pair springs, short-range
//! nucleotide exclusion, and crossing-segment separation. It is deterministic;
//! convergence is decided from monotone objective improvement rather than an
//! arbitrary iteration count.

use crate::layout::Point;
use crate::structure::ParsedStructure;

pub(crate) fn coordinates_and_arcs(
    parsed: &ParsedStructure,
) -> (Vec<Point>, Vec<crate::turtle::TurtleArc>) {
    let (original, arc_templates) = crate::turtle::coordinates_and_arcs(parsed);
    let mut points = original.clone();
    if points.len() < 4 {
        return (points, arc_templates);
    }
    let mut edges = Vec::with_capacity(points.len() + parsed.pairs.len());
    for i in 0..points.len() - 1 {
        if !parsed.strand_breaks.contains(&(i + 1)) {
            edges.push((i, i + 1, 25.0, 0.10));
        }
    }
    for pair in &parsed.pairs {
        if pair.level == 0 {
            edges.push((pair.i - 1, pair.j - 1, 35.0, 0.14));
        }
    }

    let n = points.len();
    let mut step = 0.38;
    let mut quality = layout_quality(&points, &edges, parsed);
    loop {
        let mut force = vec![Point::default(); n];
        for &(a, b, target, strength) in &edges {
            let dx = points[b].x - points[a].x;
            let dy = points[b].y - points[a].y;
            let distance = dx.hypot(dy).max(1.0e-8);
            let scale = strength * (distance - target) / distance;
            let fx = dx * scale;
            let fy = dy * scale;
            force[a].x += fx;
            force[a].y += fy;
            force[b].x -= fx;
            force[b].y -= fy;
        }

        // Nucleotide disks plus their labels need a little more room than the
        // 25-unit backbone target. Repulsion is zero outside that local range.
        for a in 0..n {
            for b in a + 1..n {
                if b == a + 1 || parsed.pairs.iter().any(|p| p.i - 1 == a && p.j - 1 == b) {
                    continue;
                }
                let dx = points[b].x - points[a].x;
                let dy = points[b].y - points[a].y;
                let distance = dx.hypot(dy);
                if distance >= 22.0 {
                    continue;
                }
                let (ux, uy) = if distance < 1.0e-8 {
                    let angle = ((a * 97 + b * 53) as f64).to_radians();
                    (angle.cos(), angle.sin())
                } else {
                    (dx / distance, dy / distance)
                };
                let magnitude = 0.12 * (22.0 - distance);
                force[a].x -= ux * magnitude;
                force[a].y -= uy * magnitude;
                force[b].x += ux * magnitude;
                force[b].y += uy * magnitude;
            }
        }

        separate_crossings(&points, &edges, &mut force);
        separate_nodes_from_edges(&points, &edges, &mut force);
        let mut candidate = points.clone();
        let mut maximum_displacement = 0.0f64;
        // Anchor the first nucleotide to remove translational drift only.
        for index in 1..n {
            let magnitude = force[index].x.hypot(force[index].y);
            let limit = 5.0;
            let scale = if magnitude > limit {
                limit / magnitude
            } else {
                1.0
            };
            let dx = step * force[index].x * scale;
            let dy = step * force[index].y * scale;
            maximum_displacement = maximum_displacement.max(dx.hypot(dy));
            candidate[index].x += dx;
            candidate[index].y += dy;
        }
        let candidate_quality = layout_quality(&candidate, &edges, parsed);
        // Accept only a numerically material monotone improvement. Since the
        // non-negative objective drops by at least this amount on every
        // accepted step, and rejected steps halve a positive step size, this
        // termination rule is finite without a hidden iteration ceiling.
        if quality - candidate_quality > 1.0e-8 {
            points = candidate;
            quality = candidate_quality;
            step = (step * 1.05).min(0.5);
            if maximum_displacement < 1.0e-7 {
                break;
            }
        } else {
            step *= 0.5;
            if step < 1.0e-10 {
                break;
            }
        }
    }
    let original_quality = layout_quality(&original, &edges, parsed);
    let optimized_quality = layout_quality(&points, &edges, parsed);
    let mut selected = if optimized_quality + 1.0e-9 < original_quality {
        points
    } else {
        original
    };
    // RNApuzzler's defining contract is an intersection-free outerplanar
    // drawing. The continuous optimizer normally preserves the Turtle
    // geometry; if it remains in a colliding local minimum, use the canonical
    // outerplanar circle embedding instead of returning an intersecting figure.
    if crate::layout::count_segment_crossings(&selected, &parsed.pairs, &parsed.strand_breaks) > 0
        && parsed.pairs.iter().all(|pair| pair.level == 0)
    {
        let radial = crate::naview::coordinates(parsed);
        selected = if crate::layout::count_segment_crossings(
            &radial,
            &parsed.pairs,
            &parsed.strand_breaks,
        ) == 0
        {
            radial
        } else {
            planar_circle_embedding(selected.len())
        };
    }
    let arcs = crate::turtle::refit_arcs(&selected, &arc_templates);
    (selected, arcs)
}

fn planar_circle_embedding(length: usize) -> Vec<Point> {
    let radius = (25.0 / (2.0 * (std::f64::consts::PI / length as f64).sin())).max(25.0);
    (0..length)
        .map(|index| {
            let angle =
                -std::f64::consts::FRAC_PI_2 + std::f64::consts::TAU * index as f64 / length as f64;
            Point {
                x: radius * angle.cos(),
                y: radius * angle.sin(),
            }
        })
        .collect()
}

fn separate_nodes_from_edges(
    points: &[Point],
    edges: &[(usize, usize, f64, f64)],
    force: &mut [Point],
) {
    const CLEARANCE: f64 = 11.0;
    for node in 0..points.len() {
        for &(a, b, _, _) in edges {
            if node == a || node == b {
                continue;
            }
            let (closest, fraction) = closest_point(points[node], points[a], points[b]);
            let dx = points[node].x - closest.x;
            let dy = points[node].y - closest.y;
            let distance = dx.hypot(dy);
            if distance >= CLEARANCE {
                continue;
            }
            let (ux, uy) = if distance < 1.0e-8 {
                let edge_dx = points[b].x - points[a].x;
                let edge_dy = points[b].y - points[a].y;
                let norm = edge_dx.hypot(edge_dy).max(1.0e-8);
                (-edge_dy / norm, edge_dx / norm)
            } else {
                (dx / distance, dy / distance)
            };
            let magnitude = 0.08 * (CLEARANCE - distance);
            force[node].x += ux * magnitude;
            force[node].y += uy * magnitude;
            force[a].x -= ux * magnitude * (1.0 - fraction);
            force[a].y -= uy * magnitude * (1.0 - fraction);
            force[b].x -= ux * magnitude * fraction;
            force[b].y -= uy * magnitude * fraction;
        }
    }
}

fn closest_point(point: Point, a: Point, b: Point) -> (Point, f64) {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let squared = dx * dx + dy * dy;
    let fraction = if squared < 1.0e-12 {
        0.0
    } else {
        (((point.x - a.x) * dx + (point.y - a.y) * dy) / squared).clamp(0.0, 1.0)
    };
    (
        Point {
            x: a.x + fraction * dx,
            y: a.y + fraction * dy,
        },
        fraction,
    )
}

fn layout_quality(
    points: &[Point],
    edges: &[(usize, usize, f64, f64)],
    parsed: &ParsedStructure,
) -> f64 {
    let crossings =
        crate::layout::count_segment_crossings(points, &parsed.pairs, &parsed.strand_breaks);
    let mut score = crossings as f64 * 10_000.0;
    for a in 0..points.len() {
        for b in a + 1..points.len() {
            if b == a + 1 {
                continue;
            }
            let distance = (points[b].x - points[a].x).hypot(points[b].y - points[a].y);
            if distance < 22.0 {
                score += (22.0 - distance).powi(2);
            }
        }
    }
    for node in 0..points.len() {
        for &(a, b, _, _) in edges {
            if node == a || node == b {
                continue;
            }
            let (closest, _) = closest_point(points[node], points[a], points[b]);
            let distance = (points[node].x - closest.x).hypot(points[node].y - closest.y);
            if distance < 9.0 {
                score += 2.0 * (9.0 - distance).powi(2);
            }
        }
    }
    for &(a, b, target, _) in edges {
        let distance = (points[b].x - points[a].x).hypot(points[b].y - points[a].y);
        score += 0.02 * (distance - target).powi(2);
    }
    score
}

fn separate_crossings(points: &[Point], edges: &[(usize, usize, f64, f64)], force: &mut [Point]) {
    for index in 0..edges.len() {
        let (a, b, _, _) = edges[index];
        for &(c, d, _, _) in &edges[index + 1..] {
            if a == c
                || a == d
                || b == c
                || b == d
                || !intersects(points[a], points[b], points[c], points[d])
            {
                continue;
            }
            let dx = points[b].x - points[a].x;
            let dy = points[b].y - points[a].y;
            let norm = dx.hypot(dy).max(1.0e-8);
            let nx = -dy / norm;
            let ny = dx / norm;
            let side = orientation(points[a], points[b], points[c]).signum();
            let magnitude = if side == 0.0 { 1.0 } else { side } * 1.8;
            for endpoint in [a, b] {
                force[endpoint].x -= nx * magnitude;
                force[endpoint].y -= ny * magnitude;
            }
            for endpoint in [c, d] {
                force[endpoint].x += nx * magnitude;
                force[endpoint].y += ny * magnitude;
            }
        }
    }
}

fn orientation(a: Point, b: Point, c: Point) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

fn intersects(a: Point, b: Point, c: Point, d: Point) -> bool {
    orientation(a, b, c) * orientation(a, b, d) < -1.0e-8
        && orientation(c, d, a) * orientation(c, d, b) < -1.0e-8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::parse_structure;

    #[test]
    fn planar_inputs_are_returned_without_straight_segment_intersections() {
        for (sequence, structure) in [
            ("GGGAAACCC", "(((...)))"),
            ("GGGAAACCCAAAGGGAAACCC", "(((...)))...(((...)))"),
            (
                "GCGCGAAACGCGCAAAGCGCGAAACGCGC",
                "(((((...)))))...(((((...)))))",
            ),
        ] {
            let parsed = parse_structure(sequence, structure).unwrap();
            let (points, _) = coordinates_and_arcs(&parsed);
            assert_eq!(
                crate::layout::count_segment_crossings(
                    &points,
                    &parsed.pairs,
                    &parsed.strand_breaks,
                ),
                0,
                "{structure}",
            );
        }
    }
}
