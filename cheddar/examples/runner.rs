use num_traits::pow;
use ocelot::svole::wykw::{
    LpnParams, LPN_EXTEND_MEDIUM, LPN_EXTEND_SMALL, LPN_SETUP_MEDIUM, LPN_SETUP_SMALL,
};
use rand::distributions::Uniform;
use rand::Rng;
use scuttlebutt::field::F61p;
use scuttlebutt::{
    track_unix_channel_pair, unix_channel_pair, AbstractChannel, AesRng, TrackChannel,
};
use square_decomp::four_squares::{FourSquaresRangeProver, FourSquaresRangeVerifier};
use std::time::{Duration, Instant};

use serde::Serialize;
use serde_json;

// todo: something is completely off with the numbers. Look at last run
// the ones that break are above the mask_bit_size, so they can't be masked.

const N: usize = 2000;
const STAT_SEC: usize = 20;
const PROD: usize = N * STAT_SEC;
const C_VEC: usize = N * 4;
const MODULUS: u64 = (1 << 61) - 1;

const VAL_BOUNDS: [u64; 2] = [2000, 4000];
const LOWER_BOUNDS: [u64; 2] = [49, 50];
const UPPER_BOUNDS: [u64; 2] = [5999, 6000];

fn compute_bound() -> (u64, u64, u64) {
    // bound = sqrt(p) / (39 * sqrt(4m) + 4))
    let infinity_bound = f64::ceil(f64::sqrt((MODULUS as f64 - 1f64) / (2f64 * 4f64)));
    let numerator: f64 = infinity_bound;
    // let denominator: f64 = 2f64 * (9.75 * f64::sqrt(4f64 * N as f64) + 2f64);
    let denominator: f64 = (9.75 * f64::sqrt(4f64 * N as f64) + 2f64);
    let old_denominator: f64 = (4f64 * N as f64) + 2f64;

    println!("num: {:?}, denom: {:?}", numerator, denominator);

    // We have to add a mask that essentially adds an extra bit, so we what we have to add is really
    // just the size of the infinity_bound / denominator. Now, the check changes to ensuring that
    // our number is really below twice that.

    let bound = f64::ceil(numerator / denominator) as u64;
    let old_bound = f64::ceil(numerator / old_denominator) as u64;

    println!("bound: {:?}, old_bound: {:?}", bound / 2, old_bound / 2);

    /* todo: since we are now adding much more than the bound, do we allow the value to potentially go very negative, since we push it above 0?
        does this allow the prover to cheat if a lot of the random values are -1? We probably need a sign bit :(
        giving a hidden bit means both proving it's actually a bit, committing to it as well as computing the multiplication

        We can't multiply with -1, what happens if we send 0/1 and then multiply this value with bound
        so that if the value is negative, P sends 1 and adds the result to res. Otherwise P sends 0, which
        adds 0 to res now. Can P use this to cheat? If the value is already positive, it will
        likely not be representable by the bits anymore.

        Now, if res is already too large, P can make it overflow, but then adding the resulting stuff
        to fill the bitlength, won't make it loop again, so now the number is very negative.
        Other avenue is to have the value be too negative. Now, adding the bound + the buffer, can actually
        make it go positive and thus work. Welp. We need to solve this problem :(

        Can we add 2*res to res, if b = 1? This is equivalent to multiplying with -1, but it'll end up costing
        one multiplication PER output (of which there are s, usually 40 or 50). So if we do
        res' = res + (-2 * res), we get res' = |res|. We only do this res < 0. This does not require
        two multiplications, since -2 is just a constant. What happens if P cheats here and let
        b = 1 but res > 0. Now we get r' = res + (-2 * res) so this becomes -res instead and can no longer
        be represented. What if it's too large though and P does this? It just becomes very negative.
        If it's negative, since we compute the absolute value, it could technically overflow when we add
        the buffer? Nope, since we have negatives after the positives, it becomes very negative instead?

    */
    // compute the buffer to the nearest power of two for cleaner proof (note, closest to 2 * bound)
    let mask_bit_size = f64::ceil((bound as f64 / 2.0).log2());

    let buffer = pow(2, mask_bit_size as usize) - bound;
    let infinity_bound_bit = f64::ceil(infinity_bound as f64).log2();
    let buffer = pow(2, infinity_bound_bit as usize) - infinity_bound as u64;
    let mask = bound + buffer;
    let buffer = 0;

    return (bound, buffer, mask_bit_size as u64);
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
    let mut prover = FourSquaresRangeProver::init(channel, &mut rng, lpn_setup, lpn_extend);

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
    let mut verifier = FourSquaresRangeVerifier::init(channel, &mut rng, lpn_setup, lpn_extend);

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
            sampling_rng.sample::<u64, _>(Uniform::new(VAL_BOUNDS[0], VAL_BOUNDS[1])),
        ));
        lower_bounds
            .push(sampling_rng.sample::<u64, _>(Uniform::new(LOWER_BOUNDS[0], LOWER_BOUNDS[1])));
        upper_bounds
            .push(sampling_rng.sample::<u64, _>(Uniform::new(UPPER_BOUNDS[0], UPPER_BOUNDS[1])));
    }

    //inp_vec[10] = F61p(upper_bounds[1] + 2);

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
