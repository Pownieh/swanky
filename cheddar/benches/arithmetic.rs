use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use rand;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaChaRng;
use square_decomp::arithmetic::{montgomery_mulmod, montgomery_powmod, mulmod, powmod};

fn bench_mulmod_64(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs_mod: Vec<u64> = (0..size).map(|_| rng.gen::<u64>() | 1).collect();
    let inputs_lhs: Vec<u64> = inputs_mod.iter().map(|&m| rng.gen_range(0..m)).collect();
    let inputs_rhs: Vec<u64> = inputs_mod.iter().map(|&m| rng.gen_range(0..m)).collect();
    let mut group = c.benchmark_group("mulmod_64");

    group.throughput(Throughput::Elements(size));
    group.bench_function("mulmod", |b| {
        b.iter(|| {
            for ((&x, &e), &m) in inputs_lhs
                .iter()
                .zip(inputs_rhs.iter())
                .zip(inputs_mod.iter())
            {
                black_box(mulmod(black_box(x), black_box(e), black_box(m)));
            }
        })
    });

    group.bench_function("montgomery_mulmod", |b| {
        b.iter(|| {
            for ((&x, &e), &m) in inputs_lhs
                .iter()
                .zip(inputs_rhs.iter())
                .zip(inputs_mod.iter())
            {
                black_box(montgomery_mulmod(black_box(x), black_box(e), black_box(m)));
            }
        })
    });
}

fn bench_powmod_64(c: &mut Criterion) {
    let size = 1_000;
    let mut rng = ChaChaRng::from_seed([42u8; 32]);
    let inputs_mod: Vec<u64> = (0..size).map(|_| rng.gen::<u64>() | 1).collect();
    let inputs_base: Vec<u64> = inputs_mod.iter().map(|&m| rng.gen_range(0..m)).collect();
    let inputs_exp: Vec<u64> = inputs_mod.iter().map(|&m| rng.gen_range(0..m)).collect();
    let mut group = c.benchmark_group("powmod_64");

    group.throughput(Throughput::Elements(size));
    group.bench_function("powmod", |b| {
        b.iter(|| {
            for ((&x, &e), &m) in inputs_base
                .iter()
                .zip(inputs_exp.iter())
                .zip(inputs_mod.iter())
            {
                black_box(powmod(black_box(x), black_box(e), black_box(m)));
            }
        })
    });

    group.bench_function("montgomery_powmod", |b| {
        b.iter(|| {
            for ((&x, &e), &m) in inputs_base
                .iter()
                .zip(inputs_exp.iter())
                .zip(inputs_mod.iter())
            {
                black_box(montgomery_powmod(black_box(x), black_box(e), black_box(m)));
            }
        })
    });
}

criterion_group!(arithmetic, bench_mulmod_64, bench_powmod_64);
criterion_main!(arithmetic);
