use crate::square_decomp::arithmetic::{isqrt, montgomery_powmod as powmod, mulmod};
use crate::square_decomp::primality::{is_prime_64, is_prime_64_faster};
use rand;
use rand::Rng;

fn legendre_symbol(a: u64, p: u64) -> u64 {
    assert!(is_prime_64_faster(p));
    powmod(a, (p - 1) / 2, p)
}

#[allow(non_snake_case)]
pub fn sqrt_m1(p: u64) -> u64 {
    assert!(is_prime_64_faster(p));
    assert_eq!(p % 4, 1);

    let mut rng = rand::thread_rng();

    let n = p - 1;

    // find non-QR
    let mut z = rng.gen_range(2..p - 2);
    while legendre_symbol(z, p) == 1 {
        z = rng.gen_range(2..p - 2);
    }

    let S = (p - 1).trailing_zeros();
    let Q = (p - 1) >> S;

    let mut M = S;
    let mut c = powmod(z, Q, p);
    let mut t = powmod(n, Q, p);
    let mut R = powmod(n, (Q + 1) / 2, p);

    let mut it = 0;
    while t != 1 {
        let mut tt = t;
        let mut i = 0;
        while tt != 1 {
            tt = powmod(t, 2, p);
            i += 1;
        }
        let b = powmod(c, 1 << (M - i - 1), p);
        M = i;
        c = powmod(b, 2, p);
        t = mulmod(t, c, p);
        R = mulmod(R, b, p);
        it += 1;
    }

    assert_eq!(powmod(R, 2, p), n);
    R
}

/// Don't ask me what this is doing ...
fn cornacchia(p: u64) -> (u64, u64) {
    assert!(is_prime_64_faster(p));
    assert_eq!(p % 4, 1);
    assert!(legendre_symbol(p - 1, p) == 1); // equivalen to p = 1 (mod 4)

    let x0 = sqrt_m1(p);
    let mut a = p;
    let mut b = x0;
    let l = isqrt(p);
    let mut it = 0;
    while b > l {
        (a, b) = (b, a % b);
        it = it + 1;
    }

    let c = p - b * b;

    // assert c is a square
    (b, isqrt(c))
}

#[allow(dead_code)]
pub fn decompose_three_squares_rs(n: u64) -> (u64, u64, u64) {
    assert_eq!(n % 4, 1);

    match n {
        1 => return (1, 0, 0),
        5 => return (2, 1, 0),
        10 => return (3, 1, 0),
        13 => return (3, 2, 0),
        34 => return (5, 3, 0),
        37 => return (6, 1, 0),
        58 => return (7, 3, 0),
        61 => return (6, 5, 0),
        85 => return (9, 2, 0),
        130 => return (11, 3, 0),
        214 => return (14, 3, 3),
        226 => return (15, 1, 0),
        370 => return (19, 3, 0),
        526 => return (21, 9, 2),
        706 => return (25, 9, 0),
        730 => return (27, 1, 0),
        829 => return (27, 10, 0),
        1414 => return (33, 18, 1),
        1549 => return (35, 18, 0),
        1906 => return (41, 15, 0),
        2986 => return (45, 31, 0),
        7549 => return (85, 18, 0),
        9634 => return (97, 15, 0),
        _ => (),
    }
    let mut rng = rand::thread_rng();
    let mut p = 0u64;
    let mut x = 0u64;
    let n_sqrt = isqrt(n);
    if n_sqrt * n_sqrt == n {
        return (n_sqrt, 0, 0);
    }
    while !is_prime_64_faster(p) {
        x = rng.gen_range(0..n_sqrt + 1);
        p = n - x * x;
    }
    let (y, z) = cornacchia(p);
    (x, y, z)
}

pub fn decompose_four_squares(n: u64) -> (u64, u64, u64, u64) {
    assert!(!(n >> 63 == 1 && n & 1 == 1));

    if n == 0 {
        return (0, 0, 0, 0);
    }
    // find t, k s.t. n = 2^t * (2k+1)
    let t = n.trailing_zeros();
    let k = ((n >> t) - 1) / 2;
    // TODO: It likely gets stuck in here somewhere
    if t == 1 {
        let (a, b, p) = {
            let mut rng = rand::thread_rng();
            let mut a = 0;
            let mut b = 0;
            let mut p = 0;
            while (a % 2 == b % 2) || !is_prime_64_faster(p) {
                a = rng.gen_range(0..isqrt(n) + 1);
                b = rng.gen_range(0..isqrt(n - a * a) + 1);
                p = n - a * a - b * b;
                if p == 1 && a <= 1 && b <= 1 {
                    break;
                }
            }
            (a, b, p)
        };
        assert_eq!(p % 4, 1);
        let (c, d) = if p == 1 { (1, 0) } else { cornacchia(p) };
        return (a, b, c, d);
    } else if t % 2 == 1 {
        let (a, b, c, d) = decompose_four_squares(2 * (2 * k + 1));
        let s = 2u64.pow((t - 1) / 2);
        return (a * s, b * s, c * s, d * s);
    } else {
        let (a, b, c, d) = decompose_four_squares(2 * (2 * k + 1)); // overflows if n is large, odd
        let (w, x, y, z) = if a % 2 == b % 2 {
            (a + b, a.abs_diff(b), c + d, c.abs_diff(d))
        } else if a % 2 == c % 2 {
            (a + c, a.abs_diff(c), b + d, b.abs_diff(d))
        } else if a % 2 == d % 2 {
            (a + d, a.abs_diff(d), b + c, b.abs_diff(c))
        } else {
            assert!(false);
            (0, 0, 0, 0)
        };
        if t == 0 {
            return (w / 2, x / 2, y / 2, z / 2);
        } else {
            let s = 2u64.pow(t / 2 - 1);
            return (w * s, x * s, y * s, z * s);
        }
    }
    // let mut p = 0u64;
    // let mut x = 0u64;
    // while !is_prime_64_faster(p) {
    //     x = rng.gen_range(0..(n as f64).sqrt().floor() as u64);
    //     p = n - x * x;
    // }
    // let (y, z) = cornacchia(p);
}

#[cfg(test)]
mod tests {
    use super::{cornacchia, decompose_four_squares, decompose_three_squares_rs, sqrt_m1};
    use crate::square_decomp::arithmetic::montgomery_powmod as powmod;

    #[test]
    fn test_sqrt_m1() {
        let p = 233;
        let r = sqrt_m1(p);
        assert_eq!(powmod(r, 2, p), p - 1);
        assert!(r == 89 || r == 144);

        let p = 15148159568689304753;
        let r = sqrt_m1(p);
        assert_eq!(powmod(r, 2, p), p - 1);
        assert!(r == 4611594693747350390 || r == 10536564874941954363);
    }

    #[test]
    fn test_sqrt_m1_other() {
        //let p = 486737;
        let p = 486739;
        let r = sqrt_m1(p);
        //assert_eq!(powmod(r, 2, p), p - 1);
        //assert_eq!(r, 134996);
    }

    #[test]
    fn test_cornacchia() {
        let p = 233;
        let (a, b) = cornacchia(p);
        assert_eq!(p, a * a + b * b);

        let p = 15148159568689304753;
        let (a, b) = cornacchia(p);
        assert_eq!(p, a * a + b * b);
    }

    #[test]
    fn test_decompose_three_squares_rs() {
        const NS: [u64; 2] = [233, 15148159568689304753];
        for n in NS {
            let (a, b, c) = decompose_three_squares_rs(n);
            assert_eq!(n, a * a + b * b + c * c, "decomposition failed for n = {n}");
        }
        for k in 0..(1 << 7) {
            let n = 4 * k + 1;
            let (a, b, c) = decompose_three_squares_rs(n);
            assert_eq!(n, a * a + b * b + c * c, "decomposition failed for n = {n}");
        }
    }

    #[test]
    fn test_decompose_four_squares() {
        const NS: [u64; 3] = [233, 15148159568689304752, 1585975039535028632];

        for n in NS {
            let (a, b, c, d) = decompose_four_squares(n);
            assert_eq!(
                n,
                a * a + b * b + c * c + d * d,
                "decomposition failed for n = {n}"
            );
        }
        for n in 0..(1 << 16) {
            let (a, b, c, d) = decompose_four_squares(n);
            assert_eq!(
                n,
                a * a + b * b + c * c + d * d,
                "decomposition failed for n = {n}"
            );
        }
    }
}
