use rand::{CryptoRng, Rng};
use scuttlebutt::ring::R64;
use scuttlebutt::ring::rx::RX;
use crate::ot::mozzarella::cache::prover::CachedProver;
use crate::ot::mozzarella::cache::verifier::CachedVerifier;



pub struct GenCache {}

impl GenCache {
    pub fn new<R: CryptoRng + Rng, const K: usize, const T: usize>(mut rng: R, delta: RX) -> (CachedProver, CachedVerifier) {

        // only produce K currently
        let mut prover_cache_u: Vec<RX> = Vec::with_capacity(K + (2 * T));
        let mut prover_cache_w: Vec<RX> = Vec::with_capacity(K + (2 * T));
        let mut verifier_cache: Vec<RX> = Vec::with_capacity(K + (2 * T));
        for _ in 0..K {
            let a1 = RX::from(rng.gen::<u128>());
            let b1 =  RX::from(rng.gen::<u128>());
            let mut tmp = a1;
            tmp *= delta;
            let mut c1 = tmp;
            c1 += b1;
            prover_cache_u.push(a1);
            prover_cache_w.push(c1);

            verifier_cache.push(b1);
        }


        // generate base voles for spsvole
        for _ in 0..T {
            let a1 =  RX::from(rng.gen::<u128>());
            let b1 = RX::from(rng.gen::<u128>());
            let mut tmp = a1;
            tmp *= delta;
            let mut c1 = tmp;
            c1 += b1;

            let a2 =  RX::from(rng.gen::<u128>());
            let b2 =  RX::from(rng.gen::<u128>());
            let mut tmp = a2;
            tmp *= delta;
            let mut c2 = tmp;
            c2 += b2;

            verifier_cache.push(b1);
            verifier_cache.push(b2);
            prover_cache_u.push(a1);
            prover_cache_w.push(c1);
            prover_cache_u.push(a2);
            prover_cache_w.push(c2);
        }

        (CachedProver::init(prover_cache_u, prover_cache_w), CachedVerifier::init(verifier_cache))
    }
}