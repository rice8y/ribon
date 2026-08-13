use ribon_core::{conditional_density2_polynomial, ConditionalDensity2Options};
use serde::Serialize;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

struct TrackingAllocator;

static CURRENT_HEAP_BYTES: AtomicUsize = AtomicUsize::new(0);
static PEAK_HEAP_BYTES: AtomicUsize = AtomicUsize::new(0);

fn add_heap_bytes(bytes: usize) {
    let current = CURRENT_HEAP_BYTES.fetch_add(bytes, Ordering::Relaxed) + bytes;
    PEAK_HEAP_BYTES.fetch_max(current, Ordering::Relaxed);
}

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            add_heap_bytes(layout.size());
        }
        pointer
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() {
            add_heap_bytes(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) };
        CURRENT_HEAP_BYTES.fetch_sub(layout.size(), Ordering::Relaxed);
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let new_pointer = unsafe { System.realloc(pointer, layout, new_size) };
        if !new_pointer.is_null() {
            if new_size >= layout.size() {
                add_heap_bytes(new_size - layout.size());
            } else {
                CURRENT_HEAP_BYTES.fetch_sub(layout.size() - new_size, Ordering::Relaxed);
            }
        }
        new_pointer
    }
}

#[global_allocator]
static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator;

#[derive(Serialize)]
struct Row {
    length: usize,
    elapsed_seconds: f64,
    normalized_cubic_seconds: f64,
    peak_heap_bytes: usize,
    pair_probability_entries: usize,
    log_partition_function: f64,
}

fn input(length: usize) -> (String, String) {
    let sequence = (0..length)
        .map(|index| if index % 2 == 0 { 'G' } else { 'C' })
        .collect::<String>();
    let mut structure = vec!['.'; length];
    for start in (0..length.saturating_sub(21)).step_by(24) {
        let end = start + 21;
        structure[start] = '(';
        structure[end] = ')';
        structure[start + 1] = '(';
        structure[end - 1] = ')';
    }
    (sequence, structure.into_iter().collect())
}

fn main() {
    let lengths = std::env::args()
        .skip(1)
        .map(|value| value.parse::<usize>().expect("length must be an integer"))
        .collect::<Vec<_>>();
    let lengths = if lengths.is_empty() {
        vec![120, 180, 240]
    } else {
        lengths
    };
    let options = ConditionalDensity2Options::default();
    let mut rows = Vec::new();
    for length in lengths {
        let (sequence, structure) = input(length);
        let baseline_heap = CURRENT_HEAP_BYTES.load(Ordering::Relaxed);
        PEAK_HEAP_BYTES.store(baseline_heap, Ordering::Relaxed);
        let start = Instant::now();
        let result =
            conditional_density2_polynomial(&sequence, &structure, 37.0, 3, 0, 1.021, &options)
                .expect("benchmark input must be valid");
        let elapsed = start.elapsed().as_secs_f64();
        let peak_heap_bytes = PEAK_HEAP_BYTES
            .load(Ordering::Relaxed)
            .saturating_sub(baseline_heap);
        rows.push(Row {
            length,
            elapsed_seconds: elapsed,
            normalized_cubic_seconds: elapsed / (length as f64).powi(3),
            peak_heap_bytes,
            pair_probability_entries: result.pair_probabilities.len(),
            log_partition_function: result.log_partition_function,
        });
    }
    println!(
        "{}",
        serde_json::to_string(&rows).expect("serialize benchmark")
    );
}
