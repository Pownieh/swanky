extern crate core;

use ocelot::edabits::{EdabitsProver, EdabitsVerifier, ProverConv, VerifierConv};
use ocelot::svole::wykw::{LPN_EXTEND_MEDIUM, LPN_EXTEND_SMALL, LPN_SETUP_MEDIUM, LPN_SETUP_SMALL};
use rand::Rng;
use rand_chacha::rand_core::CryptoRng;
use scuttlebutt::field::F2;
use scuttlebutt::ring::FiniteRing;
use scuttlebutt::{channel::track_unix_channel_pair, field::F61p, AesRng};
use serde_json::to_vec;
use std::time::Instant;

type Prover = ProverConv<F61p>;
type Verifier = VerifierConv<F61p>;

/*  Protocol:
       Prove truncation size of x_trunc
       Prove total size of x
       prove that y = x - (2^trunc_bits * x_trunc + result) = 0
*/

fn bit_decompose<RNG: CryptoRng + Rng>(
    rng: &mut RNG,
    num: usize,
    range: usize,
) -> (Vec<F61p>, Vec<Vec<F2>>) {
    let mut x_bits = vec![Vec::with_capacity(range); num];
    let mut xs = Vec::with_capacity(num);

    let mut acc = F61p::ONE;
    let mut two = F61p::ONE + F61p::ONE;
    let mut twos = vec![F61p::ZERO; range];
    for item in twos.iter_mut().take(range) {
        *item = acc;
        acc *= two;
    }

    for j in 0..num {
        let mut x: F61p = F61p::ZERO;
        for i in 0..range {
            let bit: F2 = rng.gen();
            if bit == F2::ONE {
                x += twos[i];
            }
            x_bits[j].push(bit)
        }
        xs.push(x);
    }

    (xs, x_bits)
}

fn run_range() {
    let (mut sender, mut receiver) = track_unix_channel_pair();
    let nb_bits: usize = 13;
    let nb_range: usize = 13;
    let n = 2000;
    let num_bucket = 3;
    let num_cut = num_bucket;
    let with_quicksilver = true;
    let handle = std::thread::spawn(move || {
        #[cfg(target_os = "linux")]
        {
            let mut cpu_set = nix::sched::CpuSet::new();
            cpu_set.set(1).unwrap();
            nix::sched::sched_setaffinity(nix::unistd::Pid::from_raw(0), &cpu_set).unwrap();
        }
        let mut rng = AesRng::new();
        let start = Instant::now();
        let mut fconv_sender =
            Prover::init(&mut sender, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();
        let _ = fconv_sender
            .fcom
            .input(&mut sender, &mut rng, &[F61p::ONE])
            .unwrap();
        let _ = fconv_sender
            .fcom_f2
            .input(&mut sender, &mut rng, &[F2::ONE])
            .unwrap();

        println!("Send time (init): {:?}", start.elapsed());
        let start = Instant::now();
        let (xs, x_bits) = bit_decompose(&mut rng, n, nb_range);

        let mut xs_bits_mac = Vec::with_capacity(n);
        let xs_mac = fconv_sender
            .fcom
            .input_with_macprover(&mut sender, &mut rng, &xs)
            .unwrap();
        for i in 0..n {
            xs_bits_mac.push(
                fconv_sender
                    .fcom_f2
                    .input_with_macprover(&mut sender, &mut rng, &x_bits[i])
                    .unwrap(),
            )
        }

        let mut edabits = Vec::with_capacity(n);
        for i in 0..n {
            edabits.push(EdabitsProver {
                value: xs_mac[i],
                bits: xs_bits_mac[i].clone(),
            });
        }

        println!("Send time (random edabits): {:?}", start.elapsed());

        let start = Instant::now();
        let _ = fconv_sender
            .range_check(&mut sender, &mut rng, &edabits, nb_range, num_bucket)
            .unwrap();

        let end = start.elapsed();
        println!("Send time (conv): {:?}", end);
        println!(
            "Send time (conv) average: {:?}ns",
            end.as_nanos() as f64 / n as f64
        );
    });
    #[cfg(target_os = "linux")]
    {
        let mut cpu_set = nix::sched::CpuSet::new();
        cpu_set.set(3).unwrap();
        nix::sched::sched_setaffinity(nix::unistd::Pid::from_raw(0), &cpu_set).unwrap();
    }
    let mut rng = AesRng::new();
    let start = Instant::now();
    let mut fconv_receiver =
        Verifier::init(&mut receiver, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();

    let _ = fconv_receiver
        .fcom
        .input(&mut receiver, &mut rng, 1)
        .unwrap();
    let _ = fconv_receiver
        .fcom_f2
        .input(&mut receiver, &mut rng, 1)
        .unwrap();

    println!("Receive time (init): {:?}", start.elapsed());
    println!(
        "Send communication (init): {:.2} KB",
        receiver.kilobytes_read()
    );
    println!(
        "Receive communication (init): {:.2} KB",
        receiver.kilobytes_written()
    );

    let mut xs_bits_mac = Vec::with_capacity(n);
    let xs_mac = fconv_receiver
        .fcom
        .input(&mut receiver, &mut rng, n)
        .unwrap();

    receiver.clear();
    let start = Instant::now();

    for _ in 0..n {
        xs_bits_mac.push(
            fconv_receiver
                .fcom_f2
                .input(&mut receiver, &mut rng, nb_range)
                .unwrap(),
        )
    }

    let mut edabits = Vec::with_capacity(n);
    for i in 0..n {
        edabits.push(EdabitsVerifier {
            value: xs_mac[i],
            bits: xs_bits_mac[i].clone(),
        });
    }

    fconv_receiver
        .range_check(&mut receiver, &mut rng, &edabits, nb_range, num_bucket)
        .unwrap();

    println!("Receive time (range): {:?}", start.elapsed());
    println!(
        "Send communication (range): {:.2} KB",
        receiver.kilobytes_read()
    );
    println!(
        "Receive communication (range): {:.4} KB",
        receiver.kilobytes_written()
    );
    handle.join().unwrap();
}

fn run() {
    let (mut sender, mut receiver) = track_unix_channel_pair();
    let nb_bits: usize = 14;
    let nb_trunc_bits: usize = 5;
    let n = 3_000;
    let num_bucket = 5;
    let num_cut = num_bucket;
    let with_quicksilver = true;
    let handle = std::thread::spawn(move || {
        #[cfg(target_os = "linux")]
        {
            let mut cpu_set = nix::sched::CpuSet::new();
            cpu_set.set(1).unwrap();
            nix::sched::sched_setaffinity(nix::unistd::Pid::from_raw(0), &cpu_set).unwrap();
        }
        let mut rng = AesRng::new();
        let start = Instant::now();
        let mut fconv_sender =
            Prover::init(&mut sender, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();
        println!("Send time (init): {:?}", start.elapsed());
        let start = Instant::now();
        let xs = fconv_sender
            .random_edabits(&mut sender, &mut rng, nb_bits, n)
            .unwrap();

        let mut trunc_away_xs = Vec::with_capacity(n);
        let mut trunc_xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        for i in 0..n {
            let (truncated, truncated_away, y) =
                fconv_sender.truncate(&mut sender, &mut rng, &xs[i], nb_trunc_bits);
            trunc_xs.push(truncated);
            trunc_away_xs.push(truncated_away);
            ys.push(y);
        }

        println!("Send time (random edabits): {:?}", start.elapsed());
        let start = Instant::now();
        let _ = fconv_sender
            .conv(
                &mut sender,
                &mut rng,
                num_bucket,
                num_cut,
                &trunc_xs,
                None,
                with_quicksilver,
            )
            .unwrap();

        let _ = fconv_sender
            .conv(
                &mut sender,
                &mut rng,
                num_bucket,
                num_cut,
                &trunc_away_xs,
                None,
                with_quicksilver,
            )
            .unwrap();

        let mut ys_minus_xs = Vec::with_capacity(n);
        for i in 0..n {
            ys_minus_xs.push(fconv_sender.fcom.sub(xs[i].value, ys[i]));
        }

        assert_eq!(ys.len(), xs.len());

        fconv_sender
            .fcom
            .check_zero(&mut sender, ys_minus_xs.as_slice())
            .unwrap();

        let end = start.elapsed();
        println!("Send time (conv): {:?}", end);
        println!(
            "Send time (conv) average: {:?}ns",
            end.as_nanos() as f64 / n as f64
        );
    });
    #[cfg(target_os = "linux")]
    {
        let mut cpu_set = nix::sched::CpuSet::new();
        cpu_set.set(3).unwrap();
        nix::sched::sched_setaffinity(nix::unistd::Pid::from_raw(0), &cpu_set).unwrap();
    }
    let mut rng = AesRng::new();
    let start = Instant::now();
    let mut fconv_receiver =
        Verifier::init(&mut receiver, &mut rng, LPN_SETUP_MEDIUM, LPN_EXTEND_MEDIUM).unwrap();
    println!("Receive time (init): {:?}", start.elapsed());
    println!(
        "Send communication (init): {:.2} Mb",
        receiver.kilobits_read() / 1000.0
    );
    println!(
        "Receive communication (init): {:.2} Mb",
        receiver.kilobits_written() / 1000.0
    );
    receiver.clear();
    let start = Instant::now();

    let xs_mac = fconv_receiver
        .random_edabits(&mut receiver, &mut rng, nb_bits, n)
        .unwrap();

    let mut xs_trunc_away_mac = Vec::with_capacity(n);
    let mut xs_trunc_mac = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    for _ in 0..n {
        let (truncated, truncated_away, y) = fconv_receiver.truncate(
            &mut receiver,
            &mut rng,
            nb_trunc_bits,
            nb_bits - nb_trunc_bits,
        );
        xs_trunc_mac.push(truncated);
        xs_trunc_away_mac.push(truncated_away);
        ys.push(y);
    }

    println!("Receive time (random edabits): {:?}", start.elapsed());
    println!(
        "Send communication (random edabits): {:.2} Mb",
        receiver.kilobits_read() / 1000.0
    );
    println!(
        "Receive communication (random edabits): {:.2} Mb",
        receiver.kilobits_written() / 1000.0
    );
    receiver.clear();
    let start = Instant::now();
    fconv_receiver
        .conv(
            &mut receiver,
            &mut rng,
            num_bucket,
            num_cut,
            &xs_trunc_mac,
            None,
            with_quicksilver,
        )
        .unwrap();

    fconv_receiver
        .conv(
            &mut receiver,
            &mut rng,
            num_bucket,
            num_cut,
            &xs_trunc_away_mac,
            None,
            with_quicksilver,
        )
        .unwrap();

    let mut ys_minus_xs = Vec::with_capacity(n);
    for i in 0..n {
        ys_minus_xs.push(fconv_receiver.fcom.sub(xs_mac[i].value, ys[i]));
    }

    fconv_receiver
        .fcom
        .check_zero(&mut receiver, &mut rng, ys_minus_xs.as_slice())
        .unwrap();

    println!("Receive time (conv): {:?}", start.elapsed());
    println!(
        "Send communication (conv): {:.2} Mb",
        receiver.kilobits_read() / 1000.0
    );
    println!(
        "Receive communication (conv): {:.4} Mb",
        receiver.kilobits_written() / 1000.0
    );
    handle.join().unwrap();
}

fn main() {
    println!("\nField: F61p \n");
    run_range()
}
