extern crate core;

use ocelot::edabits::{EdabitsProver, ProverConv, VerifierConv};
use ocelot::svole::wykw::{LPN_EXTEND_MEDIUM, LPN_EXTEND_SMALL, LPN_SETUP_MEDIUM, LPN_SETUP_SMALL};
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

fn run() {
    let (mut sender, mut receiver) = track_unix_channel_pair();
    let nb_bits: usize = 45;
    let nb_trunc_bits: usize = 15;
    let n = 2_000;
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
            "Send time (conv) average: {:?}",
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
    run()
}
