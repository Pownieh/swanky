use ocelot::edabits::MacProver;
use ocelot::svole::wykw::{
    LpnParams, LPN_EXTEND_MEDIUM, LPN_EXTEND_SMALL, LPN_SETUP_MEDIUM, LPN_SETUP_SMALL,
};
use rand::distributions::Uniform;
use rand::Rng;
use scuttlebutt::field::F61p;
use scuttlebutt::{unix_channel_pair, AesRng};
use square_decomp::{prover::RangeProver, verifier::RangeVerifier};

const N: usize = 2500;
const ITERATIONS: usize = 30;
const MODULUS: u64 = (1 << 61) - 1;

fn compute_bound() -> u64 {
    // bound = sqrt(p) / (39 * sqrt(4m) + 4))
    let numerator: f64 = f64::ceil(f64::sqrt((MODULUS as f64 - 1f64) / (2f64 * 4f64)));
    let denominator: f64 = 2f64 * 9.75 * f64::sqrt(4f64 * N as f64) + 2f64;
    return f64::ceil(numerator / denominator) as u64;
}

pub fn main() {
    let bound = compute_bound();
    println!("BOUND: {:?}", bound);
    let (mut sender, mut receiver) = unix_channel_pair();
    let mut sampling_rng = AesRng::new();
    let mut inp_vec: Vec<F61p> = Vec::with_capacity(N);
    let mut lower_bounds: Vec<u64> = Vec::with_capacity(N);
    let mut upper_bounds: Vec<u64> = Vec::with_capacity(N);

    let (lpn_setup, lpn_extend) = (LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM);
    let verifier_lpn_setup = lpn_setup.clone();
    let verifier_lpn_extend = lpn_extend.clone();

    for _ in 0..N {
        inp_vec.push(F61p(
            sampling_rng.sample::<u64, _>(Uniform::new(1000, 1500)),
        ));
        lower_bounds.push(sampling_rng.sample::<u64, _>(Uniform::new(500, 600)));
        upper_bounds.push(sampling_rng.sample::<u64, _>(Uniform::new(3000, 4000)));
    }

    let verifier_lower = lower_bounds.clone();
    let verifier_upper = upper_bounds.clone();

    let prover_handle = std::thread::spawn(move || {
        let mut rng = AesRng::new();
        let mut prover = RangeProver::init(&mut sender, &mut rng, lpn_setup, lpn_extend);

        let xs_m = prover
            .fcom
            .input_with_macprover(&mut sender, &mut rng, &inp_vec)
            .unwrap();

        prover.prove_range(
            &mut sender,
            &mut rng,
            xs_m,
            lower_bounds,
            upper_bounds,
            ITERATIONS,
        );
    });

    let mut rng = AesRng::new();
    let mut verifier = RangeVerifier::init(
        &mut receiver,
        &mut rng,
        verifier_lpn_setup,
        verifier_lpn_extend,
    );

    let xs_m = verifier.fcom.input(&mut receiver, &mut rng, N).unwrap();

    verifier.prove_range(
        &mut receiver,
        &mut rng,
        xs_m,
        verifier_lower,
        verifier_upper,
        ITERATIONS,
        bound,
    );

    prover_handle.join().unwrap();
}
