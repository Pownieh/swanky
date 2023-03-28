use num_traits::pow;
use ocelot::edabits::{FComProver, FComVerifier, MacVerifier};
use ocelot::svole::wykw::LpnParams;
use ocelot::Error;
use rand::distributions::Uniform;
use rand::{CryptoRng, Rng, SeedableRng};
use scuttlebutt::field::{F40b, F61p, FiniteField};
use scuttlebutt::{AbstractChannel, AesRng, Block};

const MODULUS: u64 = (1 << 61) - 1;

pub struct RangeVerifier {
    pub fcom: FComVerifier<F61p>,
    pub fcom_f2: FComVerifier<F40b>,
}
impl RangeVerifier {
    pub fn init<C: AbstractChannel, RNG: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut RNG,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
    ) -> Self {
        let verifier = FComVerifier::init(channel, rng, lpn_setup, lpn_extend).unwrap();
        let verifier_f2 = FComVerifier::init(channel, rng, lpn_setup, lpn_extend).unwrap();
        Self {
            fcom: verifier,
            fcom_f2: verifier_f2,
        }
    }

    fn check_results(
        &mut self,
        results: Vec<F61p>,
        bound: u64,
        mask_bit_size: u64,
    ) -> Result<(), Error> {
        let mask = pow(2, mask_bit_size as usize);

        for x in results {
            if !(x.compute_signed() < mask + bound as i64) {
                return Err(Error::Other("Prover fucked up".to_string()));
            }
        }
        Ok(())
    }

    pub fn prove_range<C: AbstractChannel, RNG: CryptoRng + Rng>(
        &mut self,
        channel: &mut C,
        rng: &mut RNG,
        xs: Vec<MacVerifier<F61p>>,
        lower_bounds: &Vec<u64>,
        upper_bounds: &Vec<u64>,
        iterations: usize,
        (bound, mask, mask_bit_size): (u64, u64, u64),
    ) {
        let seed = rng.gen::<Block>();

        let _ = channel.write_block(&seed).unwrap();

        let mut vector_rng = AesRng::from_seed(seed);

        let mut triples: Vec<(MacVerifier<F61p>, MacVerifier<F61p>, MacVerifier<F61p>)> =
            Vec::with_capacity(iterations * xs.len());

        let mut results: Vec<MacVerifier<F61p>> = Vec::with_capacity(iterations);

        let mut equality_checks: Vec<MacVerifier<F61p>> = Vec::with_capacity(xs.len() * iterations);

        let mut hiding_macs: Vec<Vec<MacVerifier<F61p>>> =
            vec![Vec::with_capacity(mask_bit_size as usize); iterations];

        for i in 0..iterations {
            hiding_macs[i] = self
                .fcom
                .input(channel, rng, mask_bit_size as usize)
                .unwrap();
        }

        for j in 0..iterations {
            let out_decomposed = self.fcom.input(channel, rng, xs.len() * 4).unwrap();
            let out_squares = self.fcom.input(channel, rng, xs.len() * 4).unwrap();
            let out_bounds = self.fcom.input(channel, rng, xs.len()).unwrap();

            let mut c_vec: Vec<i8> = Vec::with_capacity(4);
            for _ in 0..(xs.len() * 4) {
                c_vec.push(vector_rng.sample::<i8, _>(Uniform::new(-1, 2)));
            }

            for i in 0..xs.len() {
                // flip lower_bound to become a negative: F61p(lb[i]) - modulus
                let x_minus_l = self.fcom.affine_add_cst(-F61p(lower_bounds[i]), xs[i]);
                let u_minus_x = self
                    .fcom
                    .affine_add_cst(F61p(upper_bounds[i]), self.fcom.neg(xs[i]));
                triples.push((x_minus_l, u_minus_x, out_bounds[i]));

                let mut right_side = out_squares[i * 4];
                triples.push((
                    out_decomposed[i * 4],
                    out_decomposed[i * 4],
                    out_squares[i * 4],
                ));
                for j in 1..4 {
                    triples.push((
                        out_decomposed[i * 4 + j],
                        out_decomposed[i * 4 + j],
                        out_squares[i * 4 + j],
                    ));

                    right_side = self.fcom.add(right_side, out_squares[i * 4 + j]);
                }

                equality_checks.push(self.fcom.sub(out_bounds[i], right_side));
            }

            let mut res: MacVerifier<F61p>;
            if c_vec[0] == -1 {
                res = self.fcom.neg(out_decomposed[0]);
            } else if c_vec[0] == 1 {
                res = out_decomposed[0];
            } else {
                res = self.fcom.sub(out_decomposed[0], out_decomposed[0]);
            }

            for i in 1..out_decomposed.len() {
                if c_vec[i] == -1 {
                    res = self.fcom.sub(res, out_decomposed[i]);
                } else if c_vec[i] == 1 {
                    res = self.fcom.add(res, out_decomposed[i]);
                }
            }

            for k in 0..mask_bit_size {
                res = self.fcom.add(
                    res,
                    self.fcom
                        .affine_mult_cst(pow(F61p(2), k as usize), hiding_macs[j][k as usize]),
                )
            }

            results.push(res);
        }
        let mut out: Vec<F61p> = Vec::with_capacity(iterations);
        self.fcom.open(channel, &results, &mut out).unwrap();

        self.check_results(out, bound, mask_bit_size).unwrap();

        self.fcom
            .quicksilver_check_multiply(channel, rng, &triples)
            .unwrap();

        self.fcom
            .check_zero(channel, rng, &equality_checks)
            .unwrap();
    }
}
