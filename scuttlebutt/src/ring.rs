mod r64;
pub mod rx;

use std::num::FpCategory::Zero;
pub use r64::R64;

use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};
use serde::{Deserialize, Serialize};
use generic_array::{ArrayLength, GenericArray};
use crate::Block;

pub trait Ring:
    'static
    + Send
    + Clone
    + Copy
    + Eq
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + Mul<Self, Output=Self>
    + Add<Self, Output=Self>
    + Sub<Self, Output=Self>
    + Ord
    + PartialOrd
    + PartialEq
    + Default
    + std::fmt::Debug
    + std::iter::Sum
    + rand::distributions::Distribution<Self>

{
    /// The number of bytes in the byte representation for this ring element.
    type ByteReprLen: ArrayLength<u8>;
    /// The error that can result from trying to decode an invalid byte sequence.
    type FromBytesError: std::error::Error + Send + Sync + 'static;
    /// Deserialize a field element from a byte array.
    ///
    /// NOTE: for security purposes, this function will accept exactly one byte sequence for each
    /// field element.
    fn from_bytes(
        bytes: &GenericArray<u8, Self::ByteReprLen>,
    ) -> Self;
    /// Serialize a field element into a byte array.
    fn to_bytes(&self) -> GenericArray<u8, Self::ByteReprLen>;

    fn from_block(b: Block) -> Self;

    fn from_u128(u: u128) -> Self;

    fn from_u64(u: u64) -> Self;

    fn as_mut_ptr(&mut self) -> *mut u8;

    fn as_ptr(&self) -> *const u8;

    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }

    fn reduce_to_delta(u: u128) -> Self;

    //fn compute_sum(inp: &[Self]) -> Self;

    const ZERO: Self;
}
