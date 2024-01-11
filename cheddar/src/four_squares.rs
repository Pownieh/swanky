use crate::square_decomp::decomposition::decompose_four_squares;
use num_traits::pow;
use ocelot::edabits::{FComProver, FComVerifier, MacProver, MacVerifier};
use ocelot::svole::wykw::LpnParams;
use ocelot::Error;
use rand::distributions::Uniform;
use rand::{CryptoRng, Rng, SeedableRng};
use scuttlebutt::field::{F61p, FiniteField};
use scuttlebutt::ring::FiniteRing;
use scuttlebutt::Block;
use scuttlebutt::{AbstractChannel, AesRng};
use std::ops::Sub;
use std::time::{Duration, Instant};

const MODULUS: u64 = (1 << 61) - 1;

pub struct FourSquaresRangeProver {
    pub fcom: FComProver<F61p>,
}

#[allow(non_snake_case)]
impl FourSquaresRangeProver {
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

    fn compute_bit_decomposition(&mut self, x: F61p, size: usize) -> Vec<F61p> {
        let mut val = x.0;

        // all uneven integers are 1 off lol, so a lot of the equality checks fail.. fix lol

        let mut res = Vec::with_capacity(size);
        for _ in 0..size {
            res.push(if val & 1 == 1 { F61p::ONE } else { F61p::ZERO });
            val = val >> 1;
        }
        // println!("x: {:?}, bit_decomp in bits: {:?}", x, res);
        return res;
    }

    pub fn prove_range<
        C: AbstractChannel,
        RNG: CryptoRng + Rng,
        const N: usize,
        const PROD: usize, // iterations * N
        const C_VEC: usize,
    >(
        &mut self,
        channel: &mut C,
        rng: &mut RNG,
        xs: Vec<MacProver<F61p>>,
        lower_bounds: &Vec<u64>,
        upper_bounds: &Vec<u64>,
        iterations: usize,
        (bound, buffer, mask_bit_size): (u64, u64, u64),
    ) -> Duration {
        let seed = channel.read_block().unwrap();
        let mut vector_rng = AesRng::from_seed(seed);

        let mut doubles: [(MacProver<F61p>, MacProver<F61p>); N] =
            [(MacProver::default(), MacProver::default()); N];

        /*let mut triples: [(MacProver<F61p>, MacProver<F61p>, MacProver<F61p>); MORE_PROD] = [(
            MacProver::default(),
            MacProver::default(),
            MacProver::default(),
        );
            MORE_PROD];*/
        let mut triples: Vec<(MacProver<F61p>, MacProver<F61p>, MacProver<F61p>)> =
            Vec::with_capacity(5 * xs.len());

        let mut equality_checks: Vec<MacProver<F61p>> = Vec::with_capacity(xs.len() * iterations);

        let mut results: Vec<MacProver<F61p>> = Vec::with_capacity(iterations);

        let mut computing_ranges = Duration::ZERO;
        let mut sum = 0u64;
        let mut decomposed_vals = Vec::with_capacity(xs.len() * 4);
        let mut decomposed_vals_squared: Vec<F61p> = Vec::with_capacity(xs.len() * 4);
        let mut combined_bounds: Vec<F61p> = Vec::with_capacity(xs.len());

        let low: u64 = 50;
        let upp: u64 = 6000;
        let n_x: u64 = (3000 - low) * (upp - 3000);

        let decomposed = decompose_four_squares(n_x);
        //println!("decomposed vals: {:?}", decomposed);

        for i in 0..xs.len() {
            let x_minus_l = self.fcom.affine_add_cst(-F61p(lower_bounds[i]), xs[i]);
            let u_minus_x = self
                .fcom
                .affine_add_cst(F61p(upper_bounds[i]), self.fcom.neg(xs[i]));
            doubles[i] = (x_minus_l, u_minus_x);
            let val_to_prove = x_minus_l.0 * u_minus_x.0;
            combined_bounds.push(val_to_prove);
            let start = Instant::now();

            let (a, b, c, d) = decompose_four_squares(val_to_prove.0);
            println!(
                "decomposed: {:?} into vals: ({:?},{:?},{:?},{:?})",
                val_to_prove.0, a, b, c, d
            );

            computing_ranges += start.elapsed();
            //assert_eq!(a * a + b * b + c * c + d * d, val_to_prove.0);
            //println!("{:?}, {:?}, {:?}, {:?}", a, b, c, d);
            sum += a + b + c + d;
            decomposed_vals.extend(vec![F61p(a), F61p(b), F61p(c), F61p(d)]);
            decomposed_vals_squared.extend(vec![
                F61p(a) * F61p(a),
                F61p(b) * F61p(b),
                F61p(c) * F61p(c),
                F61p(d) * F61p(d),
            ]);
        }

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

        let mut bit_decompositions: Vec<F61p> = Vec::with_capacity(iterations * 6);

        // Compute the inner product thing
        for j in 0..iterations {
            let mut c_vec: [i8; C_VEC] = [0i8; C_VEC];
            for i in 0..C_VEC {
                c_vec[i] = (vector_rng.sample::<i8, _>(Uniform::new(-1, 2)));
            }

            let mut res: MacProver<F61p>;
            //let zeroth = c_vec[j * out_decomposed.len()];
            let zeroth = c_vec[0];
            if zeroth == -1 {
                res = self.fcom.neg(out_decomposed[0]);
            } else if zeroth == 1 {
                res = out_decomposed[0];
            } else {
                res = self.fcom.sub(out_decomposed[0], out_decomposed[0]);
            }

            for i in 1..out_decomposed.len() {
                //let idx = c_vec[(j * out_decomposed.len()) + i];
                let idx = c_vec[i];
                if idx == -1 {
                    res = self.fcom.sub(res, out_decomposed[i]);
                } else if idx == 1 {
                    res = self.fcom.add(res, out_decomposed[i]);
                }
            }

            /*println!(
                "The final value: {:?}, sum: {:?}",
                res.0.compute_signed(),
                sum
            );*/

            // todo: Do this in a batch lol
            let mut sign_bit = Vec::with_capacity(1);
            if res.0.compute_signed() < 0 {
                sign_bit = self
                    .fcom
                    .input_with_macprover(channel, rng, &[F61p::ONE])
                    .unwrap();
            } else {
                sign_bit = self
                    .fcom
                    .input_with_macprover(channel, rng, &[F61p::ZERO])
                    .unwrap();
            }

            let x = self.fcom.affine_mult_cst(-F61p(2), res);
            let adder = sign_bit[0].0 * x.0;
            let adder_m = self
                .fcom
                .input_with_macprover(channel, rng, &[adder])
                .unwrap();

            triples.push((sign_bit[0], x, adder_m[0]));

            res = self.fcom.add(adder_m[0], res);
            res = self.fcom.affine_add_cst(F61p(buffer), res);

            bit_decompositions
                .extend(&self.compute_bit_decomposition(res.0, mask_bit_size as usize));

            results.push(res);
        }

        let bit_decomposition_macs = self
            .fcom
            .input_with_macprover(channel, rng, &bit_decompositions)
            .unwrap();

        for i in &bit_decomposition_macs {
            let x_minus_one = self.fcom.affine_add_cst(F61p(0).sub(F61p::ONE), *i);
            let out = x_minus_one.0 * i.0;
            let out_m = self
                .fcom
                .input_with_macprover(channel, rng, &[out])
                .unwrap();
            triples.push((x_minus_one, *i, out_m[0]));
            equality_checks.push(out_m[0]);
        }

        // todo: prove that this bit_decomposition stuff is actually bits
        let mut length_checks: Vec<MacProver<F61p>> = Vec::with_capacity(iterations);
        for i in 0..iterations {
            let mut res = Default::default();
            for k in 0..mask_bit_size {
                if k == 0 {
                    res = self.fcom.affine_mult_cst(
                        pow(F61p(2), k as usize),
                        bit_decomposition_macs[i * mask_bit_size as usize + k as usize],
                    );
                    continue;
                }
                res = self.fcom.add(
                    res,
                    self.fcom.affine_mult_cst(
                        pow(F61p(2), k as usize),
                        bit_decomposition_macs[i * mask_bit_size as usize + k as usize],
                    ),
                )
            }
            // println!("result: {:?}, bit_decomp: {:?}", results[i], res);
            length_checks.push(self.fcom.sub(results[i], res));
        }

        // todo: this open is not required by the protocol, it is debug
        self.fcom.open(channel, &results).unwrap();

        println!("triples len: {:?}", triples.len());
        self.fcom
            .quicksilver_check_multiply(channel, rng, &triples)
            .unwrap();

        equality_checks.append(&mut length_checks);
        self.fcom.check_zero(channel, &equality_checks).unwrap();

        computing_ranges
    }
}

pub struct FourSquaresRangeVerifier {
    pub fcom: FComVerifier<F61p>,
}

impl FourSquaresRangeVerifier {
    pub fn init<C: AbstractChannel, RNG: CryptoRng + Rng>(
        channel: &mut C,
        rng: &mut RNG,
        lpn_setup: LpnParams,
        lpn_extend: LpnParams,
    ) -> Self {
        let verifier = FComVerifier::init(channel, rng, lpn_setup, lpn_extend).unwrap();
        Self { fcom: verifier }
    }

    fn check_results(
        &mut self,
        results: Vec<F61p>,
        bound: u64,
        mask_bit_size: u64,
    ) -> Result<(), Error> {
        let mask = pow(2, mask_bit_size as usize);

        for x in results {
            println!(
                "x={:?}, x_signed={:?}, mask={:?}",
                x,
                x.compute_signed(),
                mask
            );
            if !(x.compute_signed() < mask as i64) {
                println!(
                    "x={:?}, x_signed={:?}, mask={:?}",
                    x,
                    x.compute_signed(),
                    mask
                );
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
        (bound, buffer, mask_bit_size): (u64, u64, u64),
    ) {
        let seed = rng.gen::<Block>();

        let _ = channel.write_block(&seed).unwrap();

        let mut vector_rng = AesRng::from_seed(seed);

        let mut triples: Vec<(MacVerifier<F61p>, MacVerifier<F61p>, MacVerifier<F61p>)> =
            Vec::with_capacity(iterations * xs.len());

        let mut results: Vec<MacVerifier<F61p>> = Vec::with_capacity(iterations);

        let mut equality_checks: Vec<MacVerifier<F61p>> = Vec::with_capacity(xs.len() * iterations);

        let out_decomposed = self.fcom.input(channel, rng, xs.len() * 4).unwrap();
        let out_squares = self.fcom.input(channel, rng, xs.len() * 4).unwrap();
        let out_bounds = self.fcom.input(channel, rng, xs.len()).unwrap();

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

        for j in 0..iterations {
            let mut c_vec: Vec<i8> = Vec::with_capacity(xs.len() * 4);
            for _ in 0..(xs.len() * 4) {
                c_vec.push(vector_rng.sample::<i8, _>(Uniform::new(-1, 2)));
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

            let sign_bit = self.fcom.input(channel, rng, 1).unwrap()[0];
            let x = self.fcom.affine_mult_cst(-F61p(2), res);
            let adder_m = self.fcom.input(channel, rng, 1).unwrap()[0];
            triples.push((sign_bit, x, adder_m));

            res = self.fcom.add(adder_m, res);
            res = self.fcom.affine_add_cst(F61p(buffer), res);

            results.push(res);
        }

        let bit_decomposition_macs = self
            .fcom
            .input(channel, rng, iterations * mask_bit_size as usize)
            .unwrap();

        for i in &bit_decomposition_macs {
            let x_minus_one = self.fcom.affine_add_cst(F61p(0).sub(F61p::ONE), *i);
            let out_m = self.fcom.input(channel, rng, 1).unwrap();
            triples.push((x_minus_one, *i, out_m[0]));
            equality_checks.push(out_m[0]);
        }

        // prove that this bit_decomposition stuff is actually bits
        let mut length_checks: Vec<MacVerifier<F61p>> = Vec::with_capacity(iterations);
        for i in 0..iterations {
            let mut res = Default::default();
            for k in 0..mask_bit_size {
                if k == 0 {
                    res = self.fcom.affine_mult_cst(
                        pow(F61p(2), k as usize),
                        bit_decomposition_macs[i * mask_bit_size as usize + k as usize],
                    );
                    continue;
                }
                res = self.fcom.add(
                    res,
                    self.fcom.affine_mult_cst(
                        pow(F61p(2), k as usize),
                        bit_decomposition_macs[i * mask_bit_size as usize + k as usize],
                    ),
                )
            }
            length_checks.push(self.fcom.sub(results[i], res));
        }

        let mut out: Vec<F61p> = Vec::with_capacity(iterations);
        self.fcom.open(channel, &results, &mut out).unwrap();

        self.check_results(out, bound, mask_bit_size).unwrap();

        self.fcom
            .quicksilver_check_multiply(channel, rng, &triples)
            .unwrap();

        equality_checks.append(&mut length_checks);

        self.fcom
            .check_zero(channel, rng, &equality_checks)
            .unwrap();
    }
}
