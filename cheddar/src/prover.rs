use crate::square_decomp::decomposition::decompose_four_squares;
use num_traits::pow;
use ocelot::edabits::{FComProver, MacProver};
use ocelot::svole::wykw::LpnParams;
use ocelot::Error;
use rand::distributions::Uniform;
use rand::{CryptoRng, Rng, SeedableRng};
use scuttlebutt::field::F61p;
use scuttlebutt::{AbstractChannel, AesRng};
use std::ptr::eq;

const MODULUS: u64 = (1 << 61) - 1;

pub struct RangeProver {
    pub fcom: FComProver<F61p>,
}

#[allow(non_snake_case)]
impl RangeProver {
    pub fn init<C: AbstractChannel, RNG: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut RNG,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
    ) -> Self {
        let prov: FComProver<F61p> = FComProver::init(channel, rng, lpn_setup, lpn_extend).unwrap();
        println!("PROVER> Init");
        Self { fcom: prov }
    }

    fn check_results(results: Vec<F61p>, bound: u64) -> Result<(), Error> {
        for x in results {
            if !(x.0 < bound) {
                return Err(Error::Other("Prover fucked up".to_string()));
            }
        }

        Ok(())
    }

    pub fn prove_range<C: AbstractChannel, RNG: CryptoRng + Rng>(
        &mut self,
        channel: &mut C,
        rng: &mut RNG,
        xs: Vec<MacProver<F61p>>,
        lower_bounds: Vec<u64>,
        upper_bounds: Vec<u64>,
        iterations: usize,
    ) {
        let seed = channel.read_block().unwrap();
        let mut vector_rng = AesRng::from_seed(seed);

        //println!("c_vec={:?}", c_vec);

        //let u = u64::from(x.0 .0);

        let mut doubles: Vec<(MacProver<F61p>, MacProver<F61p>)> =
            Vec::with_capacity(iterations * xs.len());

        let mut triples: Vec<(MacProver<F61p>, MacProver<F61p>, MacProver<F61p>)> =
            Vec::with_capacity(iterations * xs.len()); // todo: fix the size of this thing to also include the squares

        let mut equality_checks: Vec<MacProver<F61p>> = Vec::with_capacity(xs.len() * iterations);

        let mut results: Vec<MacProver<F61p>> = Vec::with_capacity(iterations);

        for _ in 0..iterations {
            let mut sum = 0u64;
            let mut decomposed_vals = Vec::with_capacity(xs.len() * 4);
            let mut decomposed_vals_squared: Vec<F61p> = Vec::with_capacity(xs.len() * 4);
            let mut combined_bounds: Vec<F61p> = Vec::with_capacity(xs.len());

            let mut c_vec: Vec<i8> = Vec::with_capacity(xs.len());
            for _ in 0..(xs.len() * 4) {
                c_vec.push(vector_rng.sample::<i8, _>(Uniform::new(-1, 2)));
            }
            for i in 0..xs.len() {
                let x_minus_l = self.fcom.affine_add_cst(-F61p(lower_bounds[i]), xs[i]);
                let u_minus_x = self
                    .fcom
                    .affine_add_cst(F61p(upper_bounds[i]), self.fcom.neg(xs[i]));
                doubles.push((x_minus_l, u_minus_x));
                let val_to_prove = x_minus_l.0 * u_minus_x.0;
                combined_bounds.push(val_to_prove);
                let (a, b, c, d) = decompose_four_squares(val_to_prove.0);
                assert_eq!(a * a + b * b + c * c + d * d, val_to_prove.0);
                sum += a + b + c + d;
                decomposed_vals.extend(vec![F61p(a), F61p(b), F61p(c), F61p(d)]);
                decomposed_vals_squared.extend(vec![
                    F61p(a) * F61p(a),
                    F61p(b) * F61p(b),
                    F61p(c) * F61p(c),
                    F61p(d) * F61p(d),
                ]);
            }

            // todo: prove the equality of the decomposed values and the x

            let out_decomposed = self
                .fcom
                .input_with_macprover(channel, rng, &decomposed_vals)
                .unwrap();

            let out_squares = self
                .fcom
                .input_with_macprover(channel, rng, &decomposed_vals_squared)
                .unwrap();

            let out_bounds = self
                .fcom
                .input_with_macprover(channel, rng, &combined_bounds)
                .unwrap();

            for i in 0..out_bounds.len() {
                triples.push((doubles[i].0, doubles[i].1, out_bounds[i]));
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

            let mut res: MacProver<F61p>;
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

            let mut actual_val: i64 = res.0 .0 as i64;
            if res.0 .0 > (MODULUS - 1) / 2 {
                actual_val = (MODULUS - res.0 .0) as i64 * -1;
            }
            println!("The final value: {:?}, sum: {:?}", actual_val, sum);
            results.push(res);
        }
        self.fcom.open(channel, &results).unwrap();

        self.fcom
            .quicksilver_check_multiply(channel, rng, &triples)
            .unwrap();

        self.fcom.check_zero(channel, &equality_checks).unwrap();
    }
}
