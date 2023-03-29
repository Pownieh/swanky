use ocelot::svole::wykw::{
    LpnParams, LPN_EXTEND_MEDIUM, LPN_EXTEND_SMALL, LPN_SETUP_MEDIUM, LPN_SETUP_SMALL,
};
use rand::distributions::Uniform;
use rand::Rng;
use scuttlebutt::field::F61p;
use scuttlebutt::{
    track_unix_channel_pair, unix_channel_pair, AbstractChannel, AesRng, TrackChannel,
};
use square_decomp::{prover::RangeProver, verifier::RangeVerifier};
use std::time::{Duration, Instant};

use serde::Serialize;
use serde_json;

const N: usize = 2_000;
const STAT_SEC: usize = 40;
const PROD: usize = N * STAT_SEC;
const C_VEC: usize = N * 4;
// todo: C_VEC SIZE = 4 * N
const MODULUS: u64 = (1 << 61) - 1;

fn compute_bound() -> (u64, u64, u64) {
    // bound = sqrt(p) / (39 * sqrt(4m) + 4))
    let infinity_bound = f64::ceil(f64::sqrt((MODULUS as f64 - 1f64) / (2f64 * 4f64)));
    let numerator: f64 = infinity_bound;
    let denominator: f64 = 2f64 * (9.75 * f64::sqrt(4f64 * N as f64) + 1f64);

    // We have to add a mask that essentially adds an extra bit, so we what we have to add is really
    // just the size of the infinity_bound / denominator. Now, the check changes to ensuring that
    // our number is really below twice that.

    let bound = f64::ceil(numerator / denominator) as u64;
    let mask = bound;
    let mask_bit_size = f64::ceil((mask as f64).log2());

    return (bound, mask, mask_bit_size as u64);
}

#[derive(Clone, Debug, Default, Serialize)]
struct StatTuple {
    pub n: usize,
    pub ns_avg: f64,
    pub ns_avg_element: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
struct RuntimeStats {
    proving_ranges: Vec<Duration>,
    total_time: u128,
    computing_four_squares: Vec<Duration>,
    four_squares: StatTuple,
    proving_range: StatTuple,
    kilobytes_sent: f64,
    kilobytes_read: f64,
}

impl RuntimeStats {
    fn init() -> Self {
        Self {
            proving_ranges: Default::default(),
            total_time: 0,
            computing_four_squares: Default::default(),
            four_squares: Default::default(),
            proving_range: Default::default(),
            kilobytes_sent: 0.0,
            kilobytes_read: 0.0,
        }
    }

    fn analyse_times(times: &Vec<Duration>, elements: usize) -> StatTuple {
        let total_time: u128 = times.iter().map(|d| d.as_nanos()).sum();
        let average_iteration: f64 = total_time as f64 / times.len() as f64;
        let average_per_element: f64 = average_iteration / elements as f64;

        StatTuple {
            n: elements,
            ns_avg: average_iteration,
            ns_avg_element: average_per_element,
        }
    }

    pub fn compute_stats(&mut self, iterations: usize, elements: usize) {
        self.proving_range = Self::analyse_times(&self.proving_ranges, elements);
        self.four_squares = Self::analyse_times(&self.computing_four_squares, elements);
    }
}

fn run_prover<C: AbstractChannel>(
    channel: &mut TrackChannel<C>,
    inp_vec: &Vec<F61p>,
    lower_bounds: &Vec<u64>,
    upper_bounds: &Vec<u64>,
    lpn_setup: LpnParams,
    lpn_extend: LpnParams,
) -> (Duration, Duration) {
    let mut rng = AesRng::new();
    let mut prover = RangeProver::init(channel, &mut rng, lpn_setup, lpn_extend);

    let xs_m = prover
        .fcom
        .input_with_macprover(channel, &mut rng, &inp_vec)
        .unwrap();

    channel.clear();

    let start = Instant::now();

    let computing_squares = prover.prove_range::<_, _, N, PROD, C_VEC>(
        channel,
        &mut rng,
        xs_m,
        lower_bounds,
        upper_bounds,
        STAT_SEC,
        compute_bound(),
    );

    let proving_range_time = start.elapsed();
    (proving_range_time, computing_squares)
}

fn run_verifier<C: AbstractChannel>(
    channel: &mut C,
    lower_bounds: &Vec<u64>,
    upper_bounds: &Vec<u64>,
    lpn_setup: LpnParams,
    lpn_extend: LpnParams,
) -> Duration {
    let mut rng = AesRng::new();
    let mut verifier = RangeVerifier::init(channel, &mut rng, lpn_setup, lpn_extend);

    let xs_m = verifier.fcom.input(channel, &mut rng, N).unwrap();

    let start = Instant::now();

    verifier.prove_range(
        channel,
        &mut rng,
        xs_m,
        &lower_bounds,
        &upper_bounds,
        STAT_SEC,
        compute_bound(),
    );

    start.elapsed()
}

fn run(iterations: usize) {
    let (bound, mask, mask_bit_size) = compute_bound();
    println!(
        "BOUND: {:?}\nMASK: {:?},\nMASK_BIT_SIZE: {:?}",
        bound, mask, mask_bit_size
    );

    let (mut sender, mut receiver) = track_unix_channel_pair();
    let mut sampling_rng = AesRng::new();
    let mut inp_vec: Vec<F61p> = Vec::with_capacity(N);
    let mut lower_bounds: Vec<u64> = Vec::with_capacity(N);
    let mut upper_bounds: Vec<u64> = Vec::with_capacity(N);

    let (lpn_setup, lpn_extend) = (LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM);
    let verifier_lpn_setup = lpn_setup.clone();
    let verifier_lpn_extend = lpn_extend.clone();

    for _ in 0..N {
        inp_vec.push(F61p(
            sampling_rng.sample::<u64, _>(Uniform::new(10000, 15000)),
        ));
        lower_bounds.push(sampling_rng.sample::<u64, _>(Uniform::new(5000, 6000)));
        upper_bounds.push(sampling_rng.sample::<u64, _>(Uniform::new(16000, 17000)));
    }

    let verifier_lower = lower_bounds.clone();
    let verifier_upper = upper_bounds.clone();

    let mut results_p = RuntimeStats::init();
    let prover_handle = std::thread::spawn(move || {
        for i in 0..iterations {
            let (prover_proving_time, computing_squares) = run_prover(
                &mut sender,
                &inp_vec,
                &lower_bounds,
                &upper_bounds,
                lpn_setup,
                lpn_extend,
            );
            results_p.proving_ranges.push(prover_proving_time);
            results_p.computing_four_squares.push(computing_squares);
            if i == 0 {
                results_p.kilobytes_sent = sender.kilobytes_written();
                results_p.kilobytes_read = sender.kilobytes_read();
            }
        }
        results_p
    });

    let verifier_handle = std::thread::spawn(move || {
        for _ in 0..iterations {
            let verifier_proving_time = run_verifier(
                &mut receiver,
                &verifier_lower,
                &verifier_upper,
                verifier_lpn_setup,
                verifier_lpn_extend,
            );
        }
    });

    let mut results_p = prover_handle.join().unwrap();
    verifier_handle.join().unwrap();

    results_p.compute_stats(iterations, N);
    println!("results prover: {:?}", results_p);
}

pub fn main() {
    run(1);
}
