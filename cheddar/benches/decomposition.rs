use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaChaRng;
use square_decomp::decomposition::{decompose_four_squares, decompose_three_squares_rs};

fn bench_decompose_three_squares(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs: Vec<u64> = (0..size)
        .map(|_| loop {
            let z = rng.gen();
            if z % 4 == 1 {
                return z;
            }
        })
        .collect();
    let mut group = c.benchmark_group("decompose_three_squares");

    group.throughput(Throughput::Elements(size));
    group.bench_function("decompose_three_squares_rs", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(decompose_three_squares_rs(black_box(n)));
            }
        })
    });
}

fn bench_decompose_four_squares(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs: Vec<u64> = (0..size)
        .map(|_| loop {
            let z = rng.gen();
            if z >> 63 == 0u64 || z & 1 == 0u64 {
                return z;
            }
        })
        .collect();
    let mut group = c.benchmark_group("decompose_four_squares");

    group.throughput(Throughput::Elements(size));
    group.bench_function("decompose_four_squares", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(decompose_four_squares(black_box(n)));
            }
        })
    });
}

criterion_group!(
    is_prime,
    bench_decompose_three_squares,
    bench_decompose_four_squares
);
criterion_main!(is_prime);
