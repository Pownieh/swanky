use num_traits::int::PrimInt;
use std::fmt::Display;

#[inline(always)]
pub fn invert<T: PrimInt + Display>(a: T, m: T) -> T {
    // xgcd adapted from wikipedia
    // find s such that r = s * a + t * m
    // for r = gcd(a, m) and some integer t
    let mut s = T::zero();
    let mut r = m;
    let mut old_s = T::one();
    let mut old_r = a;

    let mut n = 0;
    while r != T::zero() {
        let quotient = old_r / r;
        (old_r, r) = if old_r > quotient * r {
            (r, old_r - quotient * r)
        } else {
            (r, quotient * r - old_r)
        };
        (old_s, s) = (s, old_s + quotient * s);
        n += 1;
    }

    if n % 2 == 1 {
        old_s = m - old_s;
    }

    old_s
}

#[inline(always)]
pub fn full_mul(lhs: u64, rhs: u64) -> u128 {
    let (low, high) = lhs.widening_mul(rhs);
    ((high as u128) << 64) | (low as u128)
}

#[inline(always)]
pub fn mulmod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    (full_mul(lhs, rhs) % modulus as u128) as u64
}

#[allow(dead_code)]
pub fn powmod(base: u64, exponent: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1;
    let mut base = base % modulus;
    let mut exponent = exponent;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mulmod(result, base, modulus);
        }
        exponent >>= 1;
        base = mulmod(base, base, modulus);
    }
    result
}

#[inline(always)]
pub fn to_montgomery_form(x: u64, modulus: u64) -> u64 {
    (((x as u128) << 64) % modulus as u128) as u64
}

#[inline(always)]
#[allow(dead_code)]
pub fn from_montgomery_form(y: u64, r_inv: u64, modulus: u64) -> u64 {
    mulmod(y, r_inv, modulus)
}

#[inline(always)]
pub fn redc(t: u128, modulus: u64, neg_inv_modulus: u64) -> u64 {
    assert!(t < (modulus as u128) << 64);
    let m = (t as u64).wrapping_mul(neg_inv_modulus);
    let m_times_modulus = full_mul(m, modulus);
    let sum_of_lower_parts =
        (t & 0xffffffffffffffffu128) + (m_times_modulus & 0xffffffffffffffffu128);
    assert_eq!(sum_of_lower_parts as u64, 0);
    let new_t = (t >> 64) + (m_times_modulus >> 64) + (sum_of_lower_parts >> 64);
    (if new_t >= modulus as u128 {
        new_t - modulus as u128
    } else {
        new_t
    }) as u64
}

#[inline(always)]
#[allow(dead_code)]
pub fn montgomery_mulmod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    // assume lhs and rhs in montgomery form with R = 2^64
    let neg_inv_modulus = (invert(modulus as u128, 1u128 << 64) as u64).wrapping_neg();
    let m_lhs = to_montgomery_form(lhs, modulus);
    let m_rhs = to_montgomery_form(rhs, modulus);
    let m_res = redc(full_mul(m_lhs, m_rhs), modulus, neg_inv_modulus);
    redc(m_res as u128, modulus, neg_inv_modulus)
}

pub fn montgomery_powmod(base: u64, exponent: u64, modulus: u64) -> u64 {
    assert_eq!(modulus % 2, 1);
    if modulus == 1 {
        return 0;
    }
    let neg_inv_modulus = (invert(modulus as u128, 1u128 << 64) as u64).wrapping_neg();
    let mut m_result = to_montgomery_form(1, modulus);
    let mut m_base = to_montgomery_form(base, modulus);
    let mut exponent = exponent;
    while exponent > 0 {
        if exponent & 1 == 1 {
            m_result = redc(full_mul(m_result, m_base), modulus, neg_inv_modulus);
        }
        exponent >>= 1;
        m_base = redc(full_mul(m_base, m_base), modulus, neg_inv_modulus);
    }
    redc(m_result as u128, modulus, neg_inv_modulus)
}

pub fn isqrt(x: u64) -> u64 {
    // TODO: find better implementation
    (x as f64).sqrt() as u64
}

#[cfg(test)]
mod tests {
    use super::{invert, montgomery_mulmod, montgomery_powmod, mulmod, powmod};

    #[test]
    fn test_invert() {
        let result = invert::<u64>(4321357400652853522, 12246519611726154379);
        assert_eq!(result, 6146556101494203857);
    }

    #[test]
    fn test_mulmod() {
        let result = mulmod(
            4321357400652853522,
            8776890861113620337,
            12246519611726154379,
        );
        assert_eq!(result, 9252311878089139471);
    }

    #[test]
    fn test_montgomery_mulmod() {
        let result = montgomery_mulmod(
            4321357400652853522,
            8776890861113620337,
            12246519611726154379,
        );
        assert_eq!(result, 9252311878089139471);
    }

    #[test]
    fn test_powmod() {
        let result = powmod(10484, 2, 486737);
        assert_eq!(result, 398431);
        let result = powmod(
            12387419097989408532,
            810153116369050041,
            1915435252845314032,
        );
        assert_eq!(result, 1763535093261198032);
        let result = powmod(
            4561009538669186476,
            6444163453014084353,
            13600670555087413879,
        );
        assert_eq!(result, 9851387365051549754);
    }

    #[test]
    fn test_montgomery_powmod() {
        let result = montgomery_powmod(1, 1, 47);
        assert_eq!(result, 1);
        let result = montgomery_powmod(2, 4, 47);
        assert_eq!(result, 16);
        let result = montgomery_powmod(2, 6, 47);
        assert_eq!(result, 17);
        let result = montgomery_powmod(2, 64, 47);
        assert_eq!(result, 25);
        let result = montgomery_powmod(
            4561009538669186476,
            6444163453014084353,
            13600670555087413879,
        );
        assert_eq!(result, 9851387365051549754);
    }
}
