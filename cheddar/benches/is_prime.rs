use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaChaRng;
use square_decomp::primality::{is_prime_32, is_prime_32_faster, is_prime_64, is_prime_64_faster};

fn bench_is_prime_32(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs: Vec<u32> = (0..size).map(|_| rng.gen()).collect();
    let mut group = c.benchmark_group("is_prime_32");

    group.throughput(Throughput::Elements(size));
    group.bench_function("is_prime_32", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(is_prime_32(black_box(n)));
            }
        })
    });

    group.bench_function("is_prime_32_faster", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(is_prime_32_faster(black_box(n)));
            }
        })
    });
}

fn bench_is_prime_64(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs: Vec<u64> = (0..size).map(|_| rng.gen()).collect();
    let mut group = c.benchmark_group("is_prime_64");

    group.throughput(Throughput::Elements(size));
    group.bench_function("is_prime_64", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(is_prime_64(black_box(n)));
            }
        })
    });

    group.bench_function("is_prime_64_faster", |b| {
        b.iter(|| {
            for &n in inputs.iter() {
                black_box(is_prime_64_faster(black_box(n)));
            }
        })
    });
}

criterion_group!(is_prime, bench_is_prime_32, bench_is_prime_64);
criterion_main!(is_prime);
