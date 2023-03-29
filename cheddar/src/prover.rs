use crate::square_decomp::decomposition::decompose_four_squares;
use num_traits::{pow, zero};
use ocelot::edabits::{f2_to_fe, FComProver, MacProver};
use ocelot::svole::wykw::LpnParams;
use ocelot::Error;
use rand::distributions::Uniform;
use rand::{CryptoRng, Rng, SeedableRng};
use scuttlebutt::field::{F40b, F61p, FiniteField, F2};
use scuttlebutt::ring::FiniteRing;
use scuttlebutt::{AbstractChannel, AesRng};
use std::ptr::eq;
use std::time::{Duration, Instant};

const MODULUS: u64 = (1 << 61) - 1;

pub struct RangeProver {
    pub fcom: FComProver<F61p>,
    pub fcom_f2: FComProver<F40b>,
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
        let prov_f2 = FComProver::init(channel, rng, lpn_setup, lpn_extend).unwrap();

        println!("PROVER> Init");
        Self {
            fcom: prov,
            fcom_f2: prov_f2,
        }
    }

    fn check_results(results: Vec<F61p>, bound: u64) -> Result<(), Error> {
        for x in results {
            if !(x.0 < bound) {
                return Err(Error::Other("Prover fucked up".to_string()));
            }
        }

        Ok(())
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
        (bound, mask, mask_bit_size): (u64, u64, u64),
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

        let mut hiding_vals: Vec<Vec<F61p>> =
            vec![Vec::with_capacity(mask_bit_size as usize); iterations];
        let mut hiding_macs: Vec<Vec<MacProver<F61p>>> =
            vec![Vec::with_capacity(mask_bit_size as usize); iterations];

        for i in 0..iterations {
            for j in 0..mask_bit_size {
                let cm = f2_to_fe(F2::random(rng));
                hiding_vals[i].push(cm);
            }
        }

        // TODO: prove these are bits
        for i in 0..iterations {
            hiding_macs[i] = self
                .fcom
                .input_with_macprover(channel, rng, &hiding_vals[i])
                .unwrap();
        }

        let mut computing_ranges = Duration::ZERO;
        let mut sum = 0u64;
        let mut decomposed_vals = Vec::with_capacity(xs.len() * 4);
        let mut decomposed_vals_squared: Vec<F61p> = Vec::with_capacity(xs.len() * 4);
        let mut combined_bounds: Vec<F61p> = Vec::with_capacity(xs.len());

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

            for k in 0..mask_bit_size {
                res = self.fcom.add(
                    res,
                    self.fcom
                        .affine_mult_cst(pow(F61p(2), k as usize), hiding_macs[j][k as usize]),
                )
            }

            results.push(res);
        }
        // TODO: Remember to hide the bits before opening
        // TODO: Remember to do a bit decomposition to show that it fits, rather than just opening the value
        self.fcom.open(channel, &results).unwrap();

        println!("triples len: {:?}", triples.len());
        self.fcom
            .quicksilver_check_multiply(channel, rng, &triples)
            .unwrap();

        self.fcom.check_zero(channel, &equality_checks).unwrap();

        computing_ranges
    }
}
