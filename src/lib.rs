#![expect(clippy::suspicious_arithmetic_impl)]
//! # ⚠️ WARNING ⚠️
//! Do not use this code for anything real! This crate is for educational purposes only.
//!
//! ## Shamir's Secret Sharing over GF(2⁸)
//! Split a secret into some number of shares (we call this number `n`). Collecting a threshold
//! number of those shares (we call this number `k`) reconstructs the secret. Collecting fewer than
//! `k` reveals nothing at all.
//!
//! ```
//! use how_to_share_a_secret::{split_secret, reconstruct_secret};
//!
//! let secret = b"Hello, world!"; // My favorite secret
//! let n_shares: u8 = 12;         // Max shares is 255
//! let k_threshold: u8 = 7;       // Threshold must not exceed the number of shares
//!
//! // We get `n_shares` shares. Each one holds a single point from each byte's polynomial.
//! let shares = split_secret(secret, k_threshold, n_shares).unwrap();
//! let result = reconstruct_secret(&shares[..(k_threshold as usize)]).unwrap();
//! assert_eq!(result, secret);
//! ```

use core::{
    convert::From,
    ops::{Add, Div, Mul, Sub},
};

/// Implementation of polynomial multiplication over GF(2⁸). All code using the [`u8`]
/// representation of the field is kept inside this module, so that ordinary `u8` arithmetic can't
/// be used by accident.
///
/// DANGER: This implementation branches on data, so it could be used in a side-channel attack to
/// retrieve information about the secret.
///
/// The implementation is naive, and done in two steps. First, multiply the two polynomials out
/// without any modulus, which gives a result of degree up to 14 -- hence the `u16`. Then reduce
/// that back below degree 8 by repeatedly subtracting shifted copies of the irreducible
/// polynomial.
mod u8_repr {
    /// MODULUS used to reduce polynomial to degree 8: m(x) = x^8 + x^4 + x^3 + x + 1
    const MODULUS: u16 = 0b100011011;

    const DEGREE: usize = 8;
    /// Multiply two polynomials of degree <= 7 into a single polynomial of degree <= 14.
    fn mul_expanded(a: u8, b: u8) -> u16 {
        let mut product: u16 = 0;
        for i in 0..DEGREE {
            if (b >> i) & 1 != 0 {
                product ^= (a as u16) << i;
            }
        }
        product
    }

    // NB: this branches on data. A real-world implementation would be constant time.
    /// Reduce a polynomial of degree <= 14 down to degree <= 7. Used for polynomial
    /// multiplication.
    fn poly_modulus(poly: u16) -> u8 {
        let mut out = poly;
        for i in (8..15).rev() {
            // If bit i is set, the current value has an x^i term to cancel.
            if (out >> i) & 1 != 0 {
                // NB: this is not constant time
                // Shift the modulus up so its leading term (x^8) aligns with x^i,
                // then XOR it in — exactly like "x^1 * m(x)" in the walkthrough.
                out ^= MODULUS << (i - 8);
            }
        }

        out as u8
    }

    pub fn poly_mul(a: u8, b: u8) -> u8 {
        let mid = mul_expanded(a, b);
        poly_modulus(mid)
    }

    #[cfg(test)]
    mod test {
        use super::mul_expanded;

        macro_rules! chk_poly {
            ($a:expr, $b:expr, $expected:expr) => {
                let res1 = mul_expanded($a, $b);
                assert_eq!(res1, $expected);
                let res2 = mul_expanded($b, $a);
                assert_eq!(res2, $expected);
            };
        }
        #[test]
        fn foo() {
            chk_poly!(128, 128, 1 << (7 + 7));
            chk_poly!(128, 1, 128);
            chk_poly!(128, 0, 0);
        }
    }
}

/// An element of the Galois field GF(2⁸).
///
/// The byte holds the coefficients of a polynomial of degree 7 or less over GF(2), where bit `i`
/// is the coefficient of `x^i`. So `0b0000_0011` is the polynomial `x + 1`.
///
/// This is **not** arithmetic modulo 256. ℤ/256 is not a field: `2 * 128 == 0` there, so it has
/// zero divisors and most of its elements have no inverse. Here every element except `0` is
/// invertible, which is exactly what Lagrange interpolation needs.
///
/// # Arithmetic
///
/// We define our field the same way as the [Advanced Encryption Standard (AES)](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197-upd1.pdf). See the linked document for a full explanation.
///
/// Addition are the same and they're both XOR.
///
/// Multiplication multiplies the two polynomials out, then takes the remainder modulo the
/// irreducible polynomial `x⁸ + x⁴ + x³ + x + 1` -- the same as AES. The following example is taken from [FIPS-197
/// (AES)](https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197-upd1.pdf).
///
/// Division is just is defined using multiplicative inverses so: `a / b == a * b.inv()`. Note that
/// division by zero returns zero but should be considered undefined behavior.
///
/// ```
/// use how_to_share_a_secret::GF28Element;
///
/// let a = GF28Element::from(42);
/// let b = GF28Element::from(11);
/// let additive_identity = GF28Element::from(0);
/// let mul_identity = GF28Element::from(1);
/// let expected = GF28Element::from(11 ^ 42);
///
/// assert_eq!(a + a, additive_identity);
/// assert_eq!(b + additive_identity, b);
/// assert_eq!(a + b, expected);
/// assert_eq!(a - b, expected);
/// assert_eq!(b - a, expected);
///
///
/// // Taken from FIPS-197
/// assert_eq!(
///     GF28Element::from(0x57) * GF28Element::from(0x83),
///     GF28Element::from(0xc1),
/// );
/// assert_eq!(a * mul_identity, a);
///
///
/// assert_eq!(a.inv(), GF28Element::from(152));
/// assert_eq!(a * a.inv(), GF28Element::from(1));
/// assert_eq!(a / a, GF28Element::from(1));
/// assert_eq!(a / mul_identity, a);
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GF28Element(pub(crate) u8);

impl GF28Element {
    /// The multiplicative inverse, computed as `self²⁵⁴` via Fermat's little theorem.
    ///
    /// Returns `0` for an input of `0`, which is undefined. See the type's documentation.
    pub fn inv(self) -> Self {
        poly_inv(self)
    }
    /// A uniformly random element, drawn from a cryptographically secure source.
    fn random() -> Self {
        Self(rand::random::<u8>())
    }
}
impl From<u8> for GF28Element {
    fn from(value: u8) -> Self {
        GF28Element(value)
    }
}
impl From<GF28Element> for u8 {
    fn from(value: GF28Element) -> Self {
        value.0
    }
}
impl Add for GF28Element {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        GF28Element(self.0 ^ rhs.0)
    }
}
impl Sub for GF28Element {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        GF28Element(self.0 ^ rhs.0)
    }
}
impl Mul for GF28Element {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        poly_mul(self, rhs)
    }
}
impl Div for GF28Element {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

// NB: we choose to use the same `x` value for every polynomial of our share. They could be
// different in theory, we just do this for conveinience.
/// A "share" is used to reconstruct a secret when combined with other shares.
pub struct Share {
    x: GF28Element,
    ys: Vec<GF28Element>,
}

impl Share {
    fn get_coord(&self, i: usize) -> Option<(GF28Element, GF28Element)> {
        self.ys.get(i).map(|y| (self.x, *y))
    }
}

#[derive(Debug, thiserror::Error)]
/// Error type for secret sharing
pub enum Error {
    #[error("Zero shares provided")]
    ZeroShares,
    #[error("Shares provided don't have matching number of coordinates")]
    MismatchedCoordinatesInShares,
    #[error("Threshold [{0}] is greater than the number of shares [{1}]")]
    ThresholdTooHigh(u8, u8),
    #[error("Threshold [{0}] is too small")]
    ThresholdTooSmall(u8),
}

// While using this I've accidentally swapped k & n. This will error but maybe these should be
// newtyped to prevent this.
/// Split provided `secret` into `n` shares so that it requires a threshold of `k` number of shares
/// to reconstruct the secret.
pub fn split_secret(secret: &[u8], k: u8, n: u8) -> Result<Vec<Share>, Error> {
    if k > n {
        return Err(Error::ThresholdTooHigh(k, n));
    }

    if k <= 1 {
        return Err(Error::ThresholdTooSmall(k));
    }

    let polynomial_degree = k as usize - 1;
    // NB: start at x = 1 because f(x = 0) is the secret value
    let mut shares_vec: Vec<Share> = (1..=n)
        .map(|i| Share {
            x: GF28Element(i),
            ys: Vec::with_capacity(secret.len()),
        })
        .collect();
    // create a random polynomial for each byte
    for sbyte in secret {
        let poly = create_polynomial((*sbyte).into(), polynomial_degree);
        // evaluate the polynomial at `n` random points for each share
        for share in shares_vec.iter_mut() {
            let y = evaluate_polynomial(share.x, &poly);
            share.ys.push(y);
        }
    }

    Ok(shares_vec)
}

fn poly_mul(a: GF28Element, b: GF28Element) -> GF28Element {
    u8_repr::poly_mul(a.into(), b.into()).into()
}

fn poly_pow(a: GF28Element, pow: usize) -> GF28Element {
    if pow == 0 {
        return 1.into();
    }
    let mut acc = a;
    for _ in 0..(pow - 1) {
        acc = poly_mul(acc, a);
    }
    acc
}

/// Recall that in a group of size N (here 255), raising any element to the power N gives the
/// identity (1). So if `a**255 == 1`, then `a**254 == a**-1`. This is Fermat's little theorem, and
/// we use it to calculate the inverse of an element.
fn poly_inv(a: GF28Element) -> GF28Element {
    let a2 = poly_pow(a, 2);
    let a3 = a2 * a;
    let a12 = poly_pow(a3, 4);
    let a14 = a12 * a2;
    let a15 = a12 * a3;
    let a240 = poly_pow(a15, 16);
    a240 * a14
}

/// Build a random polynomial of the given `degree` whose constant term is `secret`.
///
/// Coefficients are constant-first: `out[0]` is the secret, `out[degree]` is the
/// leading coefficient. A degree `d` polynomial has `d + 1` coefficients.
pub fn create_polynomial(secret: GF28Element, degree: usize) -> Vec<GF28Element> {
    (0..degree + 1)
        .map(|i| {
            if i == 0 {
                secret
            } else {
                GF28Element::random()
            }
        })
        .collect()
}

/// Use [Horner's rule](https://en.wikipedia.org/wiki/Horner%27s_method) to evaluate the polynomial.
pub fn evaluate_polynomial(x: GF28Element, poly: &[GF28Element]) -> GF28Element {
    poly.iter()
        .rev()
        .fold(GF28Element(0), |acc, &a| acc * x + a)
}

/// Interpolate a polynomial over the given points, returning the `y` value at `x = 0`.
pub fn lagrange_interpolation(coords: &[(GF28Element, GF28Element)]) -> GF28Element {
    let mut sum = GF28Element(0);
    for (i, (x_i, y_i)) in coords.iter().copied().enumerate() {
        let mut product = GF28Element(1);
        for (j, (x_j, _)) in coords.iter().copied().enumerate() {
            if i == j {
                continue;
            }
            product = product * (x_j / (x_j - x_i));
        }
        sum = sum + (y_i * product);
    }
    sum
}

/// Combine the given shares to recreate the secret.
///
/// The number of shares must be greater than or equal to the `k` threshold value. If not, the
/// function will still return `Ok`, but with an incorrect value.
pub fn reconstruct_secret(shares: &[Share]) -> Result<Vec<u8>, Error> {
    let Some(s) = shares.first() else {
        return Err(Error::ZeroShares);
    };

    let n_ys = s.ys.len();
    for s in shares {
        if n_ys != s.ys.len() {
            return Err(Error::MismatchedCoordinatesInShares);
        }
    }

    let mut out = vec![];

    for i in 0..n_ys {
        let mut coords = vec![];
        for s in shares {
            coords.push(
                s.get_coord(i)
                    .expect("length checked above where assert all ys are the same length"),
            );
        }
        let res = lagrange_interpolation(&coords);
        out.push(res.into());
    }
    Ok(out)
}

#[cfg(test)]
pub mod test {
    use crate::{
        create_polynomial, evaluate_polynomial, lagrange_interpolation, reconstruct_secret,
        split_secret,
    };

    use super::{GF28Element, poly_mul, poly_pow};

    macro_rules! chk_mul {
        ($a:expr, $b:expr, $expected:expr) => {
            let a = GF28Element($a);
            let b = GF28Element($b);
            let c = GF28Element($expected);

            // communtivit
            let res1 = poly_mul(a, b);
            assert_eq!(res1, c);
            let res2 = poly_mul(b, a);
            assert_eq!(res2, c);

            // associativity
            assert_eq!(a * (b * c), (a * b) * c);
            // distirbutivity
            assert_eq!((a * (b + c)), (a * b) + (a * c));
            // identity
            assert_eq!(a * 1.into(), a);
            assert_eq!(b * 1.into(), b);
            assert_eq!(c * 1.into(), c);

            // zero
            assert_eq!(a * 0.into(), 0.into());
            assert_eq!(b * 0.into(), 0.into());
            assert_eq!(c * 0.into(), 0.into());

            // frobenius
            assert_eq!((a + b) * (a + b), (a * a + b * b));
        };
    }

    macro_rules! chk_eval {
        ($poly:expr, $x:expr => $expected:expr) => {
            assert_eq!(
                evaluate_polynomial(GF28Element($x), &$poly),
                GF28Element($expected)
            );
        };
    }

    /// Deliberately naive evaluator: the literal sum of `a_i * x^i`.
    ///
    /// Exists only to cross-check [`evaluate_polynomial`], which uses Horner's
    /// rule. Two independent implementations agreeing is a stronger check than
    /// any hand-computed vector.
    fn naive_evaluate(x: GF28Element, poly: &[GF28Element]) -> GF28Element {
        poly.iter()
            .enumerate()
            .fold(GF28Element(0), |acc, (i, &a)| acc + a * poly_pow(x, i))
    }

    #[test]
    fn share_one_byte() {
        let secret = GF28Element(42);
        let threshold = 4;
        let poly_degree = threshold - 1;
        let polynomial = create_polynomial(secret, poly_degree);
        let mut shares = vec![];
        for x_i in 1..=(threshold + 1) {
            let x = GF28Element(x_i as u8);
            let y = evaluate_polynomial(x, &polynomial);
            shares.push((x, y));
        }
        let res = lagrange_interpolation(&shares[..threshold]);
        assert_eq!(res, secret);
    }

    #[test]
    fn split_and_reconstruct() {
        let secret = b"Hello, world!";
        let n = 20;
        let k = 10;
        let shares = split_secret(secret, k, n).unwrap();
        let result = reconstruct_secret(&shares[0..(k as usize)]).unwrap();
        assert_eq!(result, secret);

        // too few shares produces wrong result
        let wrong_result = reconstruct_secret(&shares[0..((k as usize) - 1)]).unwrap();
        assert_ne!(wrong_result, secret);
    }

    /// Vectors chosen so integer arithmetic on `u8` gives a *different* answer
    /// than GF(2⁸) does. Without that property a test still passes when `+` is
    /// an integer add and `*` an integer multiply.
    #[test]
    fn poly_eval_field_vectors() {
        // f(x) = 0x05*x + 0x09
        let deg1 = [0x09, 0x05].map(GF28Element);
        chk_eval!(deg1, 0x01 => 0x0C); // over the integers: 0x0E
        chk_eval!(deg1, 0x07 => 0x12); // over the integers: 0x2C
        // NB: 0x04 agrees with integer arithmetic -- 5*4 produces no carries --
        // so it discriminates nothing. Kept as an illustration.
        chk_eval!(deg1, 0x04 => 0x1D);

        // f(x) = 0x03*x^2 + 0x05*x + 0x09. Degree 2 exercises a second Horner
        // step, and x = 0x80 pushes intermediates past degree 7 so that the
        // reduction in poly_modulus actually runs.
        let deg2 = [0x09, 0x05, 0x03].map(GF28Element);
        chk_eval!(deg2, 0x02 => 0x0F);
        chk_eval!(deg2, 0x80 => 0x0A);
    }

    /// Horner must agree with the naive evaluator on every point of the field.
    #[test]
    fn horner_matches_naive() {
        for degree in 0..8 {
            let mut coeffs = vec![0u8; degree + 1];
            rand::fill(&mut coeffs[..]);
            let poly: Vec<GF28Element> = coeffs.into_iter().map(GF28Element).collect();

            for x in 0..=u8::MAX {
                let x = GF28Element(x);
                assert_eq!(
                    evaluate_polynomial(x, &poly),
                    naive_evaluate(x, &poly),
                    "degree {degree}, x = {x:?}, poly = {poly:?}"
                );
            }
        }
    }

    #[test]
    fn poly_eval_test() {
        assert_eq!(evaluate_polynomial(0.into(), &[0.into()]), GF28Element(0));
        assert_eq!(evaluate_polynomial(0.into(), &[4.into()]), GF28Element(4));
        assert_eq!(evaluate_polynomial(6.into(), &[4.into()]), GF28Element(4));
        assert_eq!(
            evaluate_polynomial(0.into(), &[0.into(), 1.into()]),
            GF28Element(0)
        );
        assert_eq!(
            evaluate_polynomial(1.into(), &[0.into(), 1.into()]),
            GF28Element(1)
        );
    }

    #[test]
    fn inv_test() {
        for i in 1..256 {
            let a = GF28Element::from(u8::try_from(i).unwrap());
            assert_eq!(a * a.inv(), GF28Element::from(1));
        }
    }

    #[test]
    fn mul() {
        chk_mul!(0x57, 0x01, 0x57);
        chk_mul!(0x57, 0x02, 0xae);
        chk_mul!(0x57, 0x04, 0x47);
        chk_mul!(0x57, 0x08, 0x8e);
        chk_mul!(0x57, 0x10, 0x07);
        chk_mul!(0x57, 0x13, 0xFE);
        chk_mul!(0x57, 0x83, 0xC1);
    }
}
