//! Independently authored Rust implementation of the modified radial layout
//! described by Bruccoleri and Heinrich (1988).
//!
//! The implementation models helices as regions and the intervening loops as
//! a cyclic graph. It does not copy or link the historically distributed
//! NAView C implementation, whose additional copying terms are not suitable
//! for this package.

use crate::layout::Point;
use crate::structure::ParsedStructure;
use std::f64::consts::{FRAC_PI_2, PI, SQRT_2};

const TAU: f64 = 2.0 * PI;
const MIN_ARC_SEPARATION: f64 = 0.5;
const UNSET: f64 = 9_999.0;

#[derive(Clone, Debug)]
struct Nucleotide {
    mate: usize,
    region: Option<usize>,
    point: Point,
}

#[derive(Clone, Copy, Debug)]
struct HelixRegion {
    outer_left: usize,
    inner_left: usize,
    inner_right: usize,
    outer_right: usize,
}

#[derive(Clone, Debug, Default)]
struct LoopNode {
    links: Vec<usize>,
    depth: usize,
    radius: f64,
    center: Point,
}

#[derive(Clone, Debug)]
struct LoopLink {
    neighbor: usize,
    region: usize,
    start: usize,
    end: usize,
    radial: Point,
    angle: f64,
    extruded: bool,
    broken: bool,
}

struct NaviewLayout {
    n: usize,
    bases: Vec<Nucleotide>,
    regions: Vec<HelixRegion>,
    loops: Vec<LoopNode>,
    links: Vec<LoopLink>,
    region_claimed: Vec<bool>,
}

pub(crate) fn coordinates(parsed: &ParsedStructure) -> Vec<Point> {
    if parsed.length <= 2 {
        return (0..parsed.length)
            .map(|index| Point {
                x: index as f64,
                y: 0.0,
            })
            .collect();
    }

    let mut layout = NaviewLayout::new(parsed);
    layout.build_regions();
    let exterior = layout.build_loop(0);
    let root = layout.central_loop(exterior);
    layout.place_loop(root, None);
    layout.bases.iter().skip(1).map(|base| base.point).collect()
}

impl NaviewLayout {
    fn new(parsed: &ParsedStructure) -> Self {
        let n = parsed.length;
        let mut bases = vec![
            Nucleotide {
                mate: 0,
                region: None,
                point: Point { x: UNSET, y: UNSET },
            };
            n + 1
        ];
        let mut pair_count = 0usize;
        for pair in parsed.pairs.iter().filter(|pair| pair.level == 0) {
            bases[pair.i].mate = pair.j;
            bases[pair.j].mate = pair.i;
            pair_count += 1;
        }
        // Classic NAView treats an unpaired chain as a single artificial
        // closing region.  The artificial pair only influences coordinates;
        // it is never returned as a structure edge.
        if pair_count == 0 {
            bases[1].mate = n;
            bases[n].mate = 1;
        }

        Self {
            n,
            bases,
            regions: Vec::new(),
            loops: Vec::new(),
            links: Vec::new(),
            region_claimed: Vec::new(),
        }
    }

    fn build_regions(&mut self) {
        let mut assigned = vec![false; self.n + 1];
        for index in 0..=self.n {
            let mate = self.bases[index].mate;
            if mate == 0 || assigned[index] {
                continue;
            }

            let outer_left = index;
            let outer_right = mate;
            let mut left = index;
            let mut right = mate;
            while left < right && self.bases[left].mate == right {
                assigned[left] = true;
                assigned[right] = true;
                left += 1;
                right = right.saturating_sub(1);
            }
            let region = HelixRegion {
                outer_left,
                inner_left: left - 1,
                inner_right: right + 1,
                outer_right,
            };
            let id = self.regions.len();
            for offset in 0..=(region.inner_left - region.outer_left) {
                self.bases[region.outer_left + offset].region = Some(id);
                self.bases[region.outer_right - offset].region = Some(id);
            }
            self.regions.push(region);
        }
        self.region_claimed = vec![false; self.regions.len()];
    }

    fn build_loop(&mut self, start: usize) -> usize {
        let loop_id = self.loops.len();
        self.loops.push(LoopNode::default());
        let mut cursor = start;
        loop {
            let mate = self.bases[cursor].mate;
            if mate != 0 {
                let region_id = self.bases[cursor]
                    .region
                    .expect("paired bases belong to a helix region");
                if !self.region_claimed[region_id] {
                    self.region_claimed[region_id] = true;
                    let region = self.regions[region_id];
                    let (forward_start, forward_end, reverse_start, reverse_end, child_start) =
                        if cursor == region.outer_left {
                            (
                                region.outer_left,
                                region.outer_right,
                                region.inner_right,
                                region.inner_left,
                                self.cyclic_next(region.inner_left),
                            )
                        } else {
                            debug_assert_eq!(cursor, region.inner_right);
                            (
                                region.inner_right,
                                region.inner_left,
                                region.outer_left,
                                region.outer_right,
                                self.cyclic_next(region.outer_right),
                            )
                        };
                    let child = self.build_loop(child_start);
                    let outward = self.links.len();
                    self.links
                        .push(LoopLink::new(child, region_id, forward_start, forward_end));
                    self.loops[loop_id].links.push(outward);
                    let inward = self.links.len();
                    self.links.push(LoopLink::new(
                        loop_id,
                        region_id,
                        reverse_start,
                        reverse_end,
                    ));
                    self.loops[child].links.push(inward);
                }
                cursor = mate;
            }
            cursor = self.cyclic_next(cursor);
            if cursor == start {
                break;
            }
        }
        loop_id
    }

    fn central_loop(&mut self, fallback: usize) -> usize {
        let mut root = fallback;
        let mut best_degree = 0usize;
        let mut best_depth = 0usize;
        for candidate in 0..self.loops.len() {
            let mut seen = vec![false; self.loops.len()];
            let depth = self.leaf_distance(candidate, None, &mut seen);
            self.loops[candidate].depth = depth;
            let degree = self.loops[candidate].links.len();
            if degree > best_degree || (degree == best_degree && depth > best_depth) {
                root = candidate;
                best_degree = degree;
                best_depth = depth;
            }
        }
        root
    }

    fn leaf_distance(&self, node: usize, parent: Option<usize>, seen: &mut [bool]) -> usize {
        if self.loops[node].links.len() <= 1 {
            return 0;
        }
        if seen[node] {
            return 0;
        }
        seen[node] = true;
        let distance = self.loops[node]
            .links
            .iter()
            .filter_map(|&link_id| {
                let neighbor = self.links[link_id].neighbor;
                (Some(neighbor) != parent).then(|| self.leaf_distance(neighbor, Some(node), seen))
            })
            .min()
            .unwrap_or(0);
        seen[node] = false;
        distance + 1
    }

    fn place_loop(&mut self, loop_id: usize, anchor: Option<usize>) {
        let link_ids = self.loops[loop_id].links.clone();
        if link_ids.is_empty() {
            return;
        }
        let step = TAU / (self.n + 1) as f64;
        for &link_id in &link_ids {
            let start = self.links[link_id].start;
            let end = self.links[link_id].end;
            let start_point = Point {
                x: -(step * start as f64).sin(),
                y: (step * start as f64).cos(),
            };
            let end_point = Point {
                x: -(step * end as f64).sin(),
                y: (step * end as f64).cos(),
            };
            let vector = Point {
                x: end_point.y - start_point.y,
                y: start_point.x - end_point.x,
            };
            let length = vector.x.hypot(vector.y);
            let radial = Point {
                x: vector.x / length,
                y: vector.y / length,
            };
            self.links[link_id].radial = radial;
            self.links[link_id].angle = positive_angle(radial.y.atan2(radial.x));
        }

        let anchor_index = anchor.and_then(|anchor_id| {
            let region = self.links[anchor_id].region;
            link_ids
                .iter()
                .position(|&link_id| self.links[link_id].region == region)
        });

        'radius_retry: loop {
            self.fit_radius(loop_id);
            let radius = self.loops[loop_id].radius;
            let mut center = if let Some(index) = anchor_index {
                let link = &self.links[link_ids[index]];
                let midpoint = midpoint(self.bases[link.start].point, self.bases[link.end].point);
                Point {
                    x: midpoint.x - radius * link.radial.x,
                    y: midpoint.y - radius * link.radial.y,
                }
            } else {
                Point { x: 0.0, y: 0.0 }
            };

            let first_block = self.block_start(&link_ids, anchor_index.unwrap_or(0));
            let mut block_start = first_block;
            loop {
                let mut block_end = block_start;
                let mut includes_anchor = anchor_index == Some(block_end);
                let mut block_size = 1usize;
                while block_size < link_ids.len()
                    && self.links_connected(
                        link_ids[block_end],
                        link_ids[next_index(block_end, link_ids.len())],
                    )
                {
                    block_end = next_index(block_end, link_ids.len());
                    includes_anchor |= anchor_index == Some(block_end);
                    block_size += 1;
                }

                let middle = if let Some(index) = anchor_index.filter(|_| includes_anchor) {
                    index
                } else {
                    advance_index(block_start, (block_size - 1) / 2, link_ids.len())
                };
                self.place_connected_block(
                    &link_ids,
                    block_start,
                    block_end,
                    middle,
                    anchor_index,
                    center,
                    radius,
                );

                let after = next_index(block_end, link_ids.len());
                if block_end != block_start && !(block_start == first_block && after == first_block)
                {
                    self.align_connected_block(
                        &link_ids,
                        block_start,
                        block_end,
                        includes_anchor,
                        &mut center,
                        radius,
                    );
                }
                block_start = after;
                if block_start == first_block {
                    break;
                }
            }

            for index in 0..link_ids.len() {
                let next = next_index(index, link_ids.len());
                let link_id = link_ids[index];
                let next_id = link_ids[next];
                let current = self.links[link_id].clone();
                let following = self.links[next_id].clone();
                let current_vector = subtract(self.bases[current.end].point, center);
                let next_vector = subtract(self.bases[following.start].point, center);
                let current_radius = norm(current_vector);
                let next_radius = norm(next_vector);
                let current_angle = positive_angle(current_vector.y.atan2(current_vector.x));
                let mut next_angle = positive_angle(next_vector.y.atan2(next_vector.x));
                if next_angle < current_angle {
                    next_angle += TAU;
                }
                let drawn_span = next_angle - current_angle;
                let expected_span = loop_span(current.angle, following.angle);
                if (drawn_span - expected_span).abs() > PI
                    && !current.extruded
                    && following.start.wrapping_sub(current.end) != 1
                {
                    self.links[link_id].extruded = true;
                    continue 'radius_retry;
                }

                if current.extruded {
                    self.place_extruded_segment(link_id, next_id);
                } else {
                    let count = self.cyclic_distance(current.end, following.start);
                    if count > 0 {
                        let increment = drawn_span / count as f64;
                        for offset in 1..count {
                            let base = self.cyclic_advance(current.end, offset);
                            let angle = current_angle + offset as f64 * increment;
                            let radius_at = current_radius
                                + (next_radius - current_radius) * (angle - current_angle)
                                    / drawn_span.max(1.0e-12);
                            self.bases[base].point = Point {
                                x: center.x + radius_at * angle.cos(),
                                y: center.y + radius_at * angle.sin(),
                            };
                        }
                    }
                }
            }
            self.loops[loop_id].center = center;
            break;
        }

        for (index, &link_id) in link_ids.iter().enumerate() {
            if anchor_index == Some(index) {
                continue;
            }
            self.extend_helix(link_id);
            let child = self.links[link_id].neighbor;
            self.place_loop(child, Some(link_id));
        }
    }

    fn fit_radius(&mut self, loop_id: usize) {
        let link_ids = self.loops[loop_id].links.clone();
        loop {
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            let mut tightest = f64::INFINITY;
            let mut tightest_link = None;
            for index in 0..link_ids.len() {
                let link_id = link_ids[index];
                let next_id = link_ids[next_index(index, link_ids.len())];
                let link = &self.links[link_id];
                let next = &self.links[next_id];
                let span = loop_span(link.angle, next.angle);
                let available = if link.extruded {
                    if span <= FRAC_PI_2 {
                        2.0
                    } else {
                        1.5
                    }
                } else {
                    self.cyclic_distance(link.end, next.start) as f64
                };
                numerator += span * (1.0 / available + 1.0);
                denominator += span * span / available;
                let per_base_angle = span / available;
                if !link.extruded && available > 1.0 && per_base_angle < tightest {
                    tightest = per_base_angle;
                    tightest_link = Some(link_id);
                }
            }
            let mut radius = numerator / denominator;
            radius = radius.max(SQRT_2 / 2.0);
            if tightest * radius < MIN_ARC_SEPARATION {
                self.links[tightest_link.expect("a tight loop segment exists")].extruded = true;
                continue;
            }
            if self.loops[loop_id].radius == 0.0 {
                self.loops[loop_id].radius = radius;
            }
            break;
        }
    }

    fn block_start(&mut self, link_ids: &[usize], initial: usize) -> usize {
        let mut start = initial;
        for _ in 0..link_ids.len() {
            let previous = prev_index(start, link_ids.len());
            if !self.links_connected(link_ids[previous], link_ids[start]) {
                return start;
            }
            start = previous;
        }

        let (largest_gap_index, _) = link_ids
            .iter()
            .enumerate()
            .map(|(index, &link_id)| {
                let next = link_ids[next_index(index, link_ids.len())];
                (
                    index,
                    loop_span(self.links[link_id].angle, self.links[next].angle),
                )
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .expect("loop has a connection");
        self.links[link_ids[largest_gap_index]].broken = true;
        next_index(largest_gap_index, link_ids.len())
    }

    #[allow(clippy::too_many_arguments)]
    fn place_connected_block(
        &mut self,
        link_ids: &[usize],
        block_start: usize,
        block_end: usize,
        middle: usize,
        anchor_index: Option<usize>,
        center: Point,
        radius: f64,
    ) {
        let mut upward = Some(middle);
        let mut downward = Some(middle);
        let mut direction = 0i8;
        while upward.is_some() || downward.is_some() {
            let index = match direction {
                -1 => upward,
                0 => Some(middle),
                _ => downward,
            };
            if let Some(index) = index {
                if anchor_index != Some(index) {
                    match direction {
                        0 => self.place_radial_pair(link_ids[index], center, radius),
                        -1 => {
                            let next = next_index(index, link_ids.len());
                            self.place_before(link_ids[index], link_ids[next]);
                        }
                        _ => {
                            let previous = prev_index(index, link_ids.len());
                            self.place_after(link_ids[previous], link_ids[index]);
                        }
                    }
                }
            }

            if direction < 0 {
                downward = downward.and_then(|value| {
                    (value != block_end).then(|| next_index(value, link_ids.len()))
                });
                direction = 1;
            } else {
                upward = upward.and_then(|value| {
                    (value != block_start).then(|| prev_index(value, link_ids.len()))
                });
                direction = -1;
            }
        }
    }

    fn place_radial_pair(&mut self, link_id: usize, center: Point, radius: f64) {
        let link = self.links[link_id].clone();
        let half_chord = (1.0 / (2.0 * radius)).asin();
        let start_angle = link.angle - half_chord;
        let end_angle = link.angle + half_chord;
        self.bases[link.start].point = Point {
            x: center.x + radius * start_angle.cos(),
            y: center.y + radius * start_angle.sin(),
        };
        self.bases[link.end].point = Point {
            x: center.x + radius * end_angle.cos(),
            y: center.y + radius * end_angle.sin(),
        };
    }

    fn place_before(&mut self, link_id: usize, next_id: usize) {
        let link = self.links[link_id].clone();
        let next = self.links[next_id].clone();
        let bisector = circular_midpoint(link.angle, next.angle);
        let tangent = Point {
            x: bisector.sin(),
            y: -bisector.cos(),
        };
        let gap = angle_delta(link.angle, next.angle);
        let length = connector_gap_length(link.extruded, gap);
        self.bases[link.end].point = add_scaled(self.bases[next.start].point, tangent, length);
        self.bases[link.start].point = Point {
            x: self.bases[link.end].point.x + link.radial.y,
            y: self.bases[link.end].point.y - link.radial.x,
        };
    }

    fn place_after(&mut self, link_id: usize, next_id: usize) {
        let link = self.links[link_id].clone();
        let next = self.links[next_id].clone();
        let bisector = circular_midpoint(link.angle, next.angle);
        let tangent = Point {
            x: -bisector.sin(),
            y: bisector.cos(),
        };
        let gap = angle_delta(link.angle, next.angle);
        let length = connector_gap_length(link.extruded, gap);
        self.bases[next.start].point = add_scaled(self.bases[link.end].point, tangent, length);
        self.bases[next.end].point = Point {
            x: self.bases[next.start].point.x - next.radial.y,
            y: self.bases[next.start].point.y + next.radial.x,
        };
    }

    #[allow(clippy::too_many_arguments)]
    fn align_connected_block(
        &mut self,
        link_ids: &[usize],
        block_start: usize,
        block_end: usize,
        includes_anchor: bool,
        center: &mut Point,
        radius: f64,
    ) {
        let first = self.links[link_ids[block_start]].clone();
        let last = self.links[link_ids[block_end]].clone();
        let a = self.bases[first.start].point;
        let b = self.bases[last.end].point;
        let chord = subtract(b, a);
        let chord_length = norm(chord);
        if chord_length <= 1.0e-12 {
            return;
        }
        let middle = midpoint(a, b);
        let chord_unit = scale(chord, 1.0 / chord_length);
        let center_vector = subtract(*center, middle);
        let projected = dot(center_vector, chord_unit);
        let mut normal = subtract(scale(chord_unit, projected), center_vector);
        let normal_length = norm(normal);
        if normal_length <= 1.0e-12 {
            normal = Point {
                x: -chord_unit.y,
                y: chord_unit.x,
            };
        } else {
            normal = scale(normal, 1.0 / normal_length);
        }

        let start_angle = polar_angle(subtract(a, *center));
        let mut end_angle = polar_angle(subtract(b, *center));
        if end_angle < start_angle {
            end_angle += TAU;
        }
        let sign = if end_angle - start_angle > PI {
            -1.0
        } else {
            1.0
        };
        let target_middle = add_scaled(*center, normal, sign * radius);
        let displacement = subtract(target_middle, middle);
        if includes_anchor {
            *center = subtract(*center, displacement);
        } else {
            let mut index = block_start;
            loop {
                let link = &self.links[link_ids[index]];
                self.bases[link.start].point = add(self.bases[link.start].point, displacement);
                self.bases[link.end].point = add(self.bases[link.end].point, displacement);
                if index == block_end {
                    break;
                }
                index = next_index(index, link_ids.len());
            }
        }
    }

    fn extend_helix(&mut self, link_id: usize) {
        let link = self.links[link_id].clone();
        let region = self.regions[link.region];
        let (start, end) = if link.start == region.outer_left {
            (region.outer_left, region.inner_left)
        } else {
            (region.inner_right, region.outer_right)
        };
        for offset in 1..=(end - start) {
            let index = start + offset;
            let mate = self.bases[index].mate;
            self.bases[index].point =
                add_scaled(self.bases[link.start].point, link.radial, offset as f64);
            self.bases[mate].point =
                add_scaled(self.bases[link.end].point, link.radial, offset as f64);
        }
    }

    fn place_extruded_segment(&mut self, link_id: usize, next_id: usize) {
        let link = self.links[link_id].clone();
        let next = self.links[next_id].clone();
        let mut start = link.end;
        let mut end = next.start;
        let mut count = self.cyclic_distance(start, end);
        let span = angle_delta(link.angle, next.angle);
        if count == 2 {
            self.place_circle_arc(start, end);
            return;
        }

        let chord = subtract(self.bases[end].point, self.bases[start].point);
        let chord_length = norm(chord);
        if chord_length >= 1.5 && span <= FRAC_PI_2 {
            let direction = scale(chord, 1.0 / chord_length);
            let next_start = self.cyclic_next(start);
            let previous_end = self.cyclic_previous(end);
            self.bases[next_start].point = add_scaled(self.bases[start].point, direction, 0.5);
            self.bases[previous_end].point = add_scaled(self.bases[end].point, direction, -0.5);
            start = next_start;
            end = previous_end;
        }

        loop {
            self.place_circle_arc(start, end);
            let next_start = self.cyclic_next(start);
            let previous_end = self.cyclic_previous(end);
            let first_angle = polar_angle(subtract(
                self.bases[next_start].point,
                self.bases[start].point,
            ));
            let last_angle = polar_angle(subtract(
                self.bases[previous_end].point,
                self.bases[end].point,
            ));
            let first_collision = angle_delta(link.angle, first_angle) > PI;
            let last_collision = angle_delta(last_angle, next.angle) > PI;
            if !(first_collision || last_collision) || count <= 1 {
                break;
            }
            let average = link.angle + span / 2.0;
            let first_direction = average.min(link.angle + 0.5);
            let last_direction = average.max(link.angle + span - 0.5);
            self.bases[next_start].point = add_scaled(
                self.bases[start].point,
                Point {
                    x: first_direction.cos(),
                    y: first_direction.sin(),
                },
                1.0,
            );
            self.bases[previous_end].point = add_scaled(
                self.bases[end].point,
                Point {
                    x: last_direction.cos(),
                    y: last_direction.sin(),
                },
                1.0,
            );
            start = next_start;
            end = previous_end;
            count = count.saturating_sub(2);
        }
    }

    fn place_circle_arc(&mut self, start: usize, end: usize) {
        let chord = subtract(self.bases[end].point, self.bases[start].point);
        let length = norm(chord);
        let steps = self.cyclic_distance(start, end);
        if steps <= 1 || length <= 1.0e-12 {
            return;
        }
        if length >= steps as f64 {
            let direction = scale(chord, 1.0 / length);
            for offset in 1..steps {
                let index = self.cyclic_advance(start, offset);
                self.bases[index].point = add_scaled(
                    self.bases[start].point,
                    direction,
                    offset as f64 / steps as f64,
                );
            }
            return;
        }

        let (height, increment) = arc_center(steps - 1, length);
        let direction = scale(chord, 1.0 / length);
        let chord_middle = add_scaled(self.bases[start].point, direction, length / 2.0);
        let normal = Point {
            x: direction.y,
            y: -direction.x,
        };
        let center = add_scaled(chord_middle, normal, height);
        let radius_vector = subtract(self.bases[start].point, center);
        let radius = norm(radius_vector);
        let initial = radius_vector.y.atan2(radius_vector.x);
        for offset in 1..steps {
            let index = self.cyclic_advance(start, offset);
            let angle = initial + offset as f64 * increment;
            self.bases[index].point = Point {
                x: center.x + radius * angle.cos(),
                y: center.y + radius * angle.sin(),
            };
        }
    }

    fn links_connected(&self, link_id: usize, _next_id: usize) -> bool {
        let link = &self.links[link_id];
        link.extruded || link.end + 1 == self.links[_next_id].start
    }

    fn cyclic_next(&self, index: usize) -> usize {
        if index == self.n {
            0
        } else {
            index + 1
        }
    }

    fn cyclic_previous(&self, index: usize) -> usize {
        if index == 0 {
            self.n
        } else {
            index - 1
        }
    }

    fn cyclic_advance(&self, index: usize, amount: usize) -> usize {
        (index + amount) % (self.n + 1)
    }

    fn cyclic_distance(&self, start: usize, end: usize) -> usize {
        if end >= start {
            end - start
        } else {
            end + self.n + 1 - start
        }
    }
}

impl LoopLink {
    fn new(neighbor: usize, region: usize, start: usize, end: usize) -> Self {
        Self {
            neighbor,
            region,
            start,
            end,
            radial: Point { x: 0.0, y: 0.0 },
            angle: 0.0,
            extruded: false,
            broken: false,
        }
    }
}

fn arc_center(n: usize, chord: f64) -> (f64, f64) {
    let mut upper = (n + 1) as f64 / PI;
    let mut lower = -upper - chord / ((n + 1) as f64 + 1.0e-6 - chord);
    if chord < 1.0 {
        lower = 0.0;
    }
    let mut height = 0.0;
    let theta = loop {
        let previous_height = height;
        height = (upper + lower) / 2.0;
        let radius = height.hypot(chord / 2.0);
        let theta = (1.0 - 0.5 / (radius * radius)).clamp(-1.0, 1.0).acos();
        let phi = (height / radius).clamp(-1.0, 1.0).acos();
        let error = theta * (n + 1) as f64 + 2.0 * phi - TAU;
        if error.abs() <= 1.0e-12 || height == previous_height {
            break theta;
        }
        if error > 0.0 {
            lower = height;
        } else {
            upper = height;
        }
    };
    (height, theta)
}

fn positive_angle(angle: f64) -> f64 {
    angle.rem_euclid(TAU)
}

fn polar_angle(vector: Point) -> f64 {
    positive_angle(vector.y.atan2(vector.x))
}

fn angle_delta(from: f64, to: f64) -> f64 {
    (to - from).rem_euclid(TAU)
}

fn loop_span(from: f64, to: f64) -> f64 {
    let span = to - from;
    if span <= 0.0 {
        span + TAU
    } else {
        span
    }
}

fn circular_midpoint(from: f64, to: f64) -> f64 {
    positive_angle(from + angle_delta(from, to) / 2.0)
}

fn connector_gap_length(extruded: bool, span: f64) -> f64 {
    if !extruded {
        1.0
    } else if span <= FRAC_PI_2 {
        2.0
    } else {
        1.5
    }
}

fn next_index(index: usize, length: usize) -> usize {
    (index + 1) % length
}

fn prev_index(index: usize, length: usize) -> usize {
    if index == 0 {
        length - 1
    } else {
        index - 1
    }
}

fn advance_index(index: usize, amount: usize, length: usize) -> usize {
    (index + amount) % length
}

fn add(a: Point, b: Point) -> Point {
    Point {
        x: a.x + b.x,
        y: a.y + b.y,
    }
}

fn subtract(a: Point, b: Point) -> Point {
    Point {
        x: a.x - b.x,
        y: a.y - b.y,
    }
}

fn scale(point: Point, amount: f64) -> Point {
    Point {
        x: point.x * amount,
        y: point.y * amount,
    }
}

fn add_scaled(point: Point, vector: Point, amount: f64) -> Point {
    Point {
        x: point.x + vector.x * amount,
        y: point.y + vector.y * amount,
    }
}

fn midpoint(a: Point, b: Point) -> Point {
    Point {
        x: (a.x + b.x) / 2.0,
        y: (a.y + b.y) / 2.0,
    }
}

fn dot(a: Point, b: Point) -> f64 {
    a.x * b.x + a.y * b.y
}

fn norm(point: Point) -> f64 {
    point.x.hypot(point.y)
}
