use scuttlebutt::ring::{R64, Ring};
use scuttlebutt::ring::rx::RX;

pub struct CachedVerifier {
    v: Vec<RX>, // cache
}

impl CachedVerifier {

    pub fn init(v: Vec<RX>) -> Self {
        Self {
            v,
        }
    }

    pub fn append<I1: Iterator<Item=RX>>(&mut self, v: I1) {
        self.v.extend(v);
    }

    pub fn pop(&mut self) -> RX {
        self.v.pop().unwrap()
    }

    pub fn get(&mut self, amount: usize) -> Vec<RX> {
        self.v[..amount].to_vec()

    }

    pub fn capacity(&self) -> usize {
        self.v.len()
    }
}