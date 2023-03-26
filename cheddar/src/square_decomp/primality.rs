use crate::square_decomp::arithmetic::{montgomery_powmod as powmod, mulmod};
use crate::square_decomp::data::{BASES32, BASES64, COMPRESSED_PRIMES16};

#[inline(always)]
pub fn is_odd_prime_16(n: u16) -> bool {
    assert_eq!(n & 1, 1);
    let nn = n >> 1;
    COMPRESSED_PRIMES16[(nn >> 6) as usize] & (1 << (nn & 0x3f)) != 0
}

#[inline(always)]
fn is_sprp(n: u64, d: u64, s: u32, b: u64) -> bool {
    assert_eq!(d & 1, 1);
    assert_eq!((1 << s) * d, n - 1);
    let mut test_value = powmod(b, d, n);
    if test_value == 1 || test_value == n - 1 {
        return true;
    }
    for _r in 1..s {
        test_value = mulmod(test_value, test_value, n);
        if test_value == n - 1 {
            return true;
        }
    }
    return false;
}

fn is_small_prime(n: u64) -> bool {
    n == 2 || n == 3 || n == 5 || n == 7
}

fn has_small_divisors(n: u64) -> bool {
    n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0
}

#[allow(dead_code)]
pub fn is_prime_32(n: u32) -> bool {
    if is_small_prime(n as u64) {
        return true;
    }
    if has_small_divisors(n as u64) {
        return false;
    }
    if n == 13 || n == 19 || n == 73 || n == 193 || n == 407521 || n == 299210837 {
        return true;
    }
    const BASES: [u64; 3] = [2, 7, 61];
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    for &b in BASES.iter() {
        if !is_sprp(n as u64, d as u64, s, b) {
            return false;
        }
    }
    true
}

#[allow(dead_code)]
pub fn is_prime_64(n: u64) -> bool {
    if n == 1 {
        return false;
    }

    if is_small_prime(n) {
        return true;
    }
    if has_small_divisors(n) {
        return false;
    }
    if n == 13 || n == 19 || n == 73 || n == 193 || n == 407521 || n == 299210837 {
        return true;
    }
    const BASES: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    let s = (n - 1).trailing_zeros();
    // println!("Size of shift: {:?} of {:?}", s, n);
    let d = (n - 1) >> s;
    for &b in BASES.iter() {
        if !is_sprp(n, d, s, b) {
            return false;
        }
    }
    true
}
#[allow(dead_code)]
pub fn is_prime_32_faster(n: u32) -> bool {
    is_prime_64_faster(n as u64)
}

pub fn is_prime_64_faster(n: u64) -> bool {
    if is_small_prime(n) {
        return true;
    }
    if n <= 1 || has_small_divisors(n) {
        return false;
    }
    const EXTRA_BASES_64: [u64; 8] = [0x0f, 0x87, 0x0d, 0x3c, 0x0f, 0x75, 0x41, 0x1d];
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    if !is_sprp(n, d, s, 2) {
        return false;
    }
    if n >> 32 != 0 {
        let h = (n as u32).wrapping_mul(0xad625b89) >> 18;
        let b = BASES64[h as usize];
        if !is_sprp(n, d, s, b as u64) {
            return false;
        }
        if n >> 49 != 0 {
            let b = EXTRA_BASES_64[(b >> 13) as usize];
            return is_sprp(n, d, s, b);
        }
    } else if n >> 16 != 0 {
        let h = (n as u32).wrapping_mul(0xad625b89) >> 24;
        let b = BASES32[h as usize];
        return is_sprp(n, d, s, b as u64);
    } else {
        return is_odd_prime_16(n as u16);
    }
    true
}

#[cfg(test)]
mod tests {
    use super::{is_prime_64, is_prime_64_faster};

    const PRIMES: [u64; 50] = [
        2,
        3,
        13,
        67,
        109,
        163,
        179,
        223,
        233,
        251,
        463,
        3967,
        9479,
        22961,
        26879,
        29587,
        48731,
        50909,
        54773,
        63179,
        420909421,
        650294921,
        783886867,
        1188435037,
        1672333667,
        1754489677,
        2845677833,
        3503458909,
        3611293409,
        4132872371,
        6094539728159,
        44283161804779,
        46549401506911,
        73583651457809,
        91718746998401,
        195769334221187,
        339362028841057,
        373083875718871,
        401622274612231,
        490523568034949,
        1234185778462221283,
        2408408657804543393,
        6839202078207594061,
        9507465300880812281,
        12295524207854767931,
        12899643504564021853,
        13213144882258043179,
        13484225183427939229,
        13946817839235265601,
        17575983250388076607,
    ];
    const COMPOSITES: [u64; 53] = [
        32,
        93,
        104,
        108,
        111,
        143,
        168,
        216,
        242,
        254,
        6908,
        17267,
        17286,
        25314,
        28747,
        29505,
        34665,
        42758,
        61769,
        63883,
        486737,
        1082401,
        2284453,
        467930341,
        1247332693,
        1372384182,
        1671609988,
        2385196022,
        2632148983,
        2929627100,
        3724056633,
        3932377135,
        4095899563,
        176383209028701,
        251167824017200,
        266120077698470,
        275553389651277,
        317859958937336,
        390866660108755,
        402512180110923,
        457494191908602,
        458631227188517,
        552412564479929,
        419580706053627838,
        1776017772331684059,
        4332266425710294913,
        5926659751551585706,
        7983848487050519290,
        10201683839938667304,
        13658694667273378390,
        14350244060091734518,
        14469213262370774219,
        16207761221996441327,
    ];

    #[test]
    fn test_is_prime_64() {
        for p in PRIMES {
            assert_eq!(is_prime_64(p), true, "classified prime {p} incorrectly");
        }
        for n in COMPOSITES {
            assert_eq!(
                is_prime_64(n),
                false,
                "classified composite {n} incorrectly"
            );
        }
    }

    #[test]
    fn test_is_prime_64_faster() {
        for p in PRIMES {
            assert_eq!(
                is_prime_64_faster(p),
                true,
                "classified prime {p} incorrectly"
            );
        }
        for n in COMPOSITES {
            assert_eq!(
                is_prime_64_faster(n),
                false,
                "classified composite {n} incorrectly"
            );
        }
    }
}
