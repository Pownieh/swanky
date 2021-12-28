use rand::{CryptoRng, Rng, SeedableRng};

use crate::ot::mozzarella::utils::gen_column;
use rayon::prelude::*;
use scuttlebutt::{ring::R64, AesRng, Block};
use scuttlebutt::ring::rx::RX;
use scuttlebutt::utils::K_BIT_STRING;

// Z64 Local Linear Code with parameter D
pub struct LLCode<const ROWS: usize, const COLS: usize, const D: usize> {
    indices: Vec<[(usize, RX); D]>,
}

// columns have length of rows
impl<const ROWS: usize, const COLS: usize, const D: usize> LLCode<ROWS, COLS, D> {
    pub fn from_seed(seed: Block) -> Self {
        let mut rng = AesRng::from_seed(seed);
        Self::gen(&mut rng)
    }

    pub fn gen<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        let max_val =  K_BIT_STRING as usize; // K_BIT_STRING is k+s 1s
        let mut code = LLCode {
            indices: Vec::with_capacity(COLS),
        };
        for _ in 0..COLS {
            code.indices.push(gen_column::<_, D>(rng, ROWS, max_val));
        }

        println!("COLS:\t {}", COLS);
        code.indices.sort(); // TODO: test this - sorting the rows, seems to improve cache locality
        code
    }

    // TODO: Can likely be made more efficient somehow?
    pub fn mul(&self, v: &[RX]) -> Vec<RX> {
        (self.indices.par_iter().map(|col| {
            let mut cord: RX = RX::default();
            for i in col.iter().copied() {
                cord += i.1 * v[i.0];
            }
            cord
        }))
        .collect()
    }

    // takes the indices of the code (A) and adds them to elements of a.
    pub fn mul_add(&self, v: &[RX], a: &[RX]) -> Vec<RX> {
        (self.indices.par_iter().enumerate().map(|(j, col)| {
            let mut cord: RX = RX::default();
            for i in col.iter().copied() {
                cord += i.1 * v[i.0];
            }
            cord + a[j]
        }))
        .collect()
    }
}
