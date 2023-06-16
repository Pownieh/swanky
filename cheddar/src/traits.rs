use scuttlebutt::AbstractChannel;
use scuttlebutt::field::FiniteField;
use ocelot::edabits::{MacProver, MacVerifier};


pub trait RangeProver<FE: FiniteField> {
    // TODO:
    // - handle statistical security parameter
    // - handle possible preprocessing
    fn prove_range<C: AbstractChannel>(&mut self, channel: &mut C,
        lower_bound: u64, upper_bound: u64,
        commitments: &[MacProver<FE>]);
}

pub trait RangeVerifier<FE: FiniteField> {
    // TODO:
    // - handle statistical security parameter
    // - handle possible preprocessing
    fn verify_range<C: AbstractChannel>(&mut self, channel: &mut C,
        lower_bound: u64, upper_bound: u64,
        commitments: &[MacVerifier<FE>]);
}
