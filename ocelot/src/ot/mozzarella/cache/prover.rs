use scuttlebutt::ring::{R64, Ring};
use scuttlebutt::ring::rx::RX;

/// A collection of correlated OT outputs
pub struct CachedProver {
    u: Vec<RX>, // cache
    w: Vec<RX>, // cache
}

impl CachedProver {
    pub fn init(u: Vec<RX>, w: Vec<RX>) -> Self {
        Self {
            u,
            w,
        }
    }

    pub fn get(&mut self, amount: usize) -> (Vec<RX>, Vec<RX>) {
        (
            self.u[..amount].to_vec(),
            self.w[..amount].to_vec(),
        )
    }

    pub fn pop(&mut self) -> (RX, RX) {
        let u = self.u.pop();
        let w = self.w.pop();
        (u.unwrap(), w.unwrap())
    }

    pub fn capacity(&self) -> usize {
        self.u.len()
    }

    pub fn append<I1: Iterator<Item=RX>, I2: Iterator<Item=RX>>(&mut self, u: I1, w: I2) {
        self.u.extend(u);
        self.w.extend(w);
        debug_assert_eq!(self.u.len(), self.w.len());
    }
}
