//! Defines a block as a 128-bit value, and implements block-related functions.

#[cfg(feature = "curve25519-dalek")]
use crate::Aes256;
#[cfg(feature = "curve25519-dalek")]
use curve25519_dalek::ristretto::RistrettoPoint;
use std::{
    hash::{Hash, Hasher},
};

use std::arch::aarch64::*;

/// A 128-bit chunk.
#[derive(Clone, Copy)]
pub struct ArmBlock(pub uint64x2_t);

impl ArmBlock {

    /// Convert into a pointer
    #[inline]
    pub fn as_ptr(&self) -> *const u8 {self.as_ref().as_ptr() }

    /// Extract i'th element as u64
    #[inline]
    pub fn extract_i_u64(&self, i: i8) -> u64 {
        unsafe {
            vgetq_lane_u64::<i>(self.0) as u64
        }
    }

    /// Convert into a mutable pointer.
    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut u8 {
        self.as_mut().as_mut_ptr()
    }

    /// Carryless multiplication.
    ///
    /// This code is adapted from the EMP toolkit's implementation.
    #[inline]
    pub fn clmul(self, rhs: Self) -> (Self, Self) {
        unsafe {
            let x = self.0;
            let y = rhs.0;


            (ArmBlock(x), ArmBlock(y))
        }
    }

impl AsRef<[u8]> for ArmBlock {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        unsafe { &*(self as *const ArmBlock as *const [u8; 16]) }
    }
}