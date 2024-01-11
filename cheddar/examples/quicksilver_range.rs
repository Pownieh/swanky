/*  Protocol:
    Given x of supposed k bits:
        Split x into k/m blocks of m bits
        Compute polynomials Prod_i X - i
        Run QS on the polynomials
        If result is 0, it's correct
*/

use ocelot::edabits::{
    EdabitsProver, FComProver, FComVerifier, MacProver, MacVerifier, ProverConv, VerifierConv,
};
use ocelot::svole::wykw::{LPN_EXTEND_MEDIUM, LPN_SETUP_MEDIUM};
use rand::Rng;
use rand_chacha::rand_core::CryptoRng;
use scuttlebutt::field::polynomial::Polynomial;
use scuttlebutt::field::{F61p, FiniteField};
use scuttlebutt::ring::FiniteRing;
use scuttlebutt::{track_unix_channel_pair, AbstractChannel, AesRng, Block};
use smallvec::{smallvec, SmallVec};
use std::sync::mpsc::channel;
use std::time::Instant;

/*fn decompose<FE: FiniteField<PrimeField = FE>>(
    xs: Vec<EdabitsProver<FE>>,
    m: usize,
    k: usize,
) -> Vec<Vec<FE>> {
    let mut out: Vec<Vec<FE>> = Vec::with_capacity(xs.len());
    let mut decompose_size = m / k;
    for x in xs {
        let mut res: Vec<FE> = Vec::with_capacity(decompose_size);
        for i in 0..decompose_size {
            let mut twos: FE = FE::ONE;
            let mut result = FE::ZERO;
            for j in i * m..(i + 1) * m {
                result += x.bits[j].0 * twos;
                twos += twos;
            }
            res.push(result);
        }
        out.push(res);
    }
    out
}
 */

fn sample_val<FE: FiniteField, RNG: CryptoRng + Rng>(
    rng: &mut RNG,
    m: usize,
    k: usize,
) -> (FE, Vec<FE>) {
    let chunks = k / m;
    let mut decomposition = Vec::with_capacity(chunks);
    let mut final_val = FE::ZERO;
    let mut outer_twos = FE::ONE;
    for j in 0..chunks {
        let mut val = FE::ZERO;
        let mut twos = FE::ONE;
        for i in 0..m {
            let tester = rng.gen_bool(0.5);
            let mut x = FE::ZERO;
            if tester {
                x = FE::ONE;
            }
            val += twos * x;
            twos += twos;
            if j > 0 {
                outer_twos += outer_twos
            };
        }
        decomposition.push(val);
        final_val += (outer_twos * val);
    }

    (final_val, decomposition)
}

/// variant of quicksilver
pub fn prover_quicksilver_variant<C: AbstractChannel, RNG: CryptoRng + Rng>(
    channel: &mut C,
    rng: &mut RNG,
    xs_chunks: Vec<Vec<MacProver<F61p>>>,
    range: usize,
    fcom: &mut FComProver<F61p>,
) {
    let chunks: usize = xs_chunks[0].len();
    let elements = xs_chunks.len();
    // compute the polynomials:
    let mut poly_xs: Vec<Vec<MacProver<F61p>>> = Vec::with_capacity(chunks * elements);
    for x_chunks in xs_chunks {
        //println!("x_chunks: {:?}", x_chunks);
        for i in 0..chunks {
            let mut cst = F61p::ZERO;
            let mut coeffs: Vec<MacProver<F61p>> = Vec::with_capacity(range);
            for _ in 0..range {
                let tmp = fcom.affine_sub_cst(cst, x_chunks[i]);
                coeffs.push(tmp);
                cst += F61p::ONE;
            }
            //println!("coeffs: {:?}", coeffs);
            poly_xs.push(coeffs);
        }
    }

    let mut res_polys: Vec<Polynomial<F61p>> = Vec::with_capacity(chunks * elements);

    for chunk in poly_xs {
        let mut coeff: SmallVec<[F61p; 3]> = Default::default();
        coeff.push(chunk[0].0);
        let mut res = Polynomial {
            constant: chunk[0].1,
            coefficients: coeff,
        };
        for i in 1..range {
            let mut coeff: SmallVec<[F61p; 3]> = Default::default();
            coeff.push(chunk[i].0);
            let tmp = Polynomial {
                constant: chunk[i].1,
                coefficients: coeff,
            };
            res *= &tmp;
        }

        res_polys.push(res);
    }

    println!("length of coeffs: {:?}", res_polys[0].coefficients.len());

    let mut coeffs: SmallVec<[F61p; 3]> = Default::default();
    coeffs.push(F61p::ONE);
    let mut y_poly: Polynomial<F61p> = Polynomial {
        constant: F61p::ZERO,
        coefficients: coeffs,
    };

    // degree d-1 poly that we multiply with the commitment to 0
    let mut res_y_poly = y_poly.clone();
    for i in 0..range - 2 {
        res_y_poly *= &y_poly;
    }

    // receive a seed that we don't use for the sake of communication
    let seed = channel.read_block().unwrap();

    for j in 0..(N / QUICKSILVER_LOOP) {
        let mut to_send: Polynomial<F61p> = Polynomial::constant(F61p::ZERO);
        for i in 0..QUICKSILVER_LOOP {
            to_send += &res_polys[(j * QUICKSILVER_LOOP) + i];
        }

        channel
            .write_serializable::<F61p>(&to_send.constant)
            .unwrap();
        for coeff in to_send.coefficients {
            channel.write_serializable::<F61p>(&coeff).unwrap();
        }
    }

    /*for poly in res_polys {
        channel.write_serializable::<F61p>(&poly.constant).unwrap();
        for coeff in poly.coefficients {
            channel.write_serializable::<F61p>(&coeff).unwrap();
        }
    }*/

    // we send over the polynomial as a way of "opening" it, the verifie can
    // verify using the delta polynomial thing we define.

    // q is a polynomial which is based on a commitment to the zs which are
    // the results of the multiplications; i.e. they are all 0. Can I just reuse the same
    // commitment here?

    //for poly in res_polys {
    //    println!("poly: {:?}", poly);
    //}
}

/// variant of quicksilver used for range proofs
pub fn verifier_quicksilver_variant<C: AbstractChannel, RNG: CryptoRng + Rng>(
    channel: &mut C,
    rng: &mut RNG,
    xs_chunks: Vec<Vec<MacVerifier<F61p>>>,
    range: usize,
    fcom: &mut FComVerifier<F61p>,
) {
    let chunks: usize = xs_chunks[0].len();
    let elements = xs_chunks.len();
    // compute the polynomials:
    let mut poly_xs: Vec<Vec<MacVerifier<F61p>>> = Vec::with_capacity(chunks * elements);
    for x_chunks in xs_chunks {
        //println!("x_chunks: {:?}", x_chunks);
        for i in 0..chunks {
            let mut cst = F61p::ZERO;
            let mut coeffs: Vec<MacVerifier<F61p>> = Vec::with_capacity(range);
            for _ in 0..range {
                let tmp = fcom.affine_sub_cst(cst, x_chunks[i]);
                coeffs.push(tmp);
                cst += F61p::ONE;
            }
            //println!("coeffs: {:?}", coeffs);
            poly_xs.push(coeffs);
        }
    }

    // todo: does some more work here

    channel.write_block(&Block::default()).unwrap();

    for j in 0..(N / QUICKSILVER_LOOP) {
        let _ = channel.read_serializable::<F61p>().unwrap();
        for coeff in 0..range - 1 {
            let _ = channel.read_serializable::<F61p>().unwrap();
        }
    }

    /*for poly in 0..chunks * elements {
        let _ = channel.read_serializable::<F61p>().unwrap();
        for coeff in 0..range - 1 {
            channel.read_serializable::<F61p>().unwrap();
        }
    }*/

    // we send over the polynomial as a way of "opening" it, the verifier can
    // verify using the delta polynomial thing we define.

    // q is a polynomial which is based on a commitment to the zs which are
    // the results of the multiplications; i.e. they are all 0. Can I just reuse the same
    // commitment here?

    //for poly in res_polys {
    //    println!("poly: {:?}", poly);
    //}
}

const CHUNK_LENGTH: usize = 6;
const BIT_LENGTH: usize = 24;
const N: usize = 1_000;
const QUICKSILVER_LOOP: usize = 1000;

fn run<FE: FiniteField<PrimeField = FE>>() {
    let (mut sender, mut receiver) = track_unix_channel_pair();
    let mut range = 1 << CHUNK_LENGTH;
    let handle = std::thread::spawn(move || {
        let mut rng = AesRng::new();
        let mut fcom: FComProver<F61p> =
            FComProver::init(&mut sender, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();
        let mut chunks_vec = Vec::with_capacity(N);
        let _ = fcom
            .input_with_macprover(&mut sender, &mut rng, &[F61p::ZERO])
            .unwrap();
        for _ in 0..N {
            let (x, x_decomp) = sample_val(&mut rng, CHUNK_LENGTH, BIT_LENGTH);

            //let x_mac = fcom
            //    .input_with_macprover(&mut sender, &mut rng, &[x])
            //    .unwrap();

            let x_decomp_mac = fcom
                .input_with_macprover(&mut sender, &mut rng, &x_decomp)
                .unwrap();
            chunks_vec.push(x_decomp_mac);
        }
        let start = Instant::now();
        prover_quicksilver_variant(&mut sender, &mut rng, chunks_vec, range, &mut fcom);
        let end = start.elapsed();
        println!("time in micros: {:?}", end.as_micros());
        println!(
            "time in micros on average: {:?}",
            end.as_micros() as f64 / N as f64
        );
    });

    let mut rng = AesRng::new();
    let mut fcom_receiver: FComVerifier<F61p> =
        FComVerifier::init(&mut receiver, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();

    let mut chunks_vec = Vec::with_capacity(N);
    let _ = fcom_receiver.input(&mut receiver, &mut rng, 1).unwrap();
    receiver.clear();
    for _ in 0..N {
        //let x_mac = fcom_receiver.input(&mut receiver, &mut rng, 1).unwrap();
        let x_decomp_mac = fcom_receiver
            .input(&mut receiver, &mut rng, BIT_LENGTH / CHUNK_LENGTH)
            .unwrap();

        chunks_vec.push(x_decomp_mac);
    }
    let mut communication_read = receiver.kilobytes_read();
    let mut communication_written = receiver.kilobytes_written();
    println!(
        "communication received for inputs: {:?}",
        communication_read
    );
    println!("communication sent for inputs: {:?}", communication_written);
    receiver.clear();
    verifier_quicksilver_variant(
        &mut receiver,
        &mut rng,
        chunks_vec,
        range,
        &mut fcom_receiver,
    );

    communication_read += receiver.kilobytes_read();
    communication_written += receiver.kilobytes_written();
    println!("--------------------------------------------------");

    println!(
        "communication received for proof: {:?}",
        receiver.kilobytes_read()
    );
    println!(
        "communication sent for proof: {:?}",
        receiver.kilobytes_written()
    );

    println!("--------------------------------------------------");

    println!("communication received in total: {:?}", communication_read);
    println!("communication sent in total: {:?}", communication_written);

    let _ = handle.join().unwrap();
}

fn main() {
    run::<F61p>();
}
