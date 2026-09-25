#![expect(clippy::suspicious_arithmetic_impl)]
//! Implementation of GF(2^8)
//! We need a new polynomial for every byte
//!
//! Split a secret into a given number of shares (we call this number `n`), so that by collecting some
//! threshold number of shares (we call this number `k`) you can reconstruct the secret.
//!
//! ```
//! //use crate::split_secret;
//!
//! //let s = b"Q"; // our favorite byte
//! //let n_shares = 12;
//! //let k_threshold = 7;
//!
//! //// we have `n_shares` shares, each share has a polynomial for each byte of the secret.
//! //let shares = split_secret(s, n_shares, k_threshold);
//! ```

use core::{
    convert::From,
    ops::{Add, Div, Mul, Sub},
};

/// Math with our polynomials in their [`u8`] representation. We hide it in a module so we don't use
/// them accidentally. All [`u8`] code lives here, no [`u8`] lives outside here.
mod u8_repr {
    /// MODULUS used to reduce polynomial to degree 8: m(x) = x^8 + x^4 + x^3 + x + 1
    const MODULUS: u16 = 0b100011011;

    const DEGREE: usize = 8;
    /// Multiply out two <= degree 7 polynomials into a <= degree 14 polynomial
    pub fn mul_expanded(a: u8, b: u8) -> u16 {
        let mut product: u16 = 0;
        for i in 0..DEGREE {
            if (b >> i) & 1 != 0 {
                product ^= (a as u16) << i;
            }
        }
        product
    }

    /// Reduce a degree 14 polynomial to degree 7
    pub fn poly_modulus(poly: u16) -> u8 {
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
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GF28Element(pub(crate) u8);

impl GF28Element {
    fn inv(self) -> Self {
        poly_inv(self)
    }
    /// cryptographically secure random value
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

#[derive(Debug, thiserror::Error)]
pub enum Error {}
// TODO Should we have a `Share` or `Poly` type? or others?
//pub struct Share(Vec<GF28Element>);

fn split_secret(
    secret: &[u8],
    k: u8,
    n: u8,
) -> Result<Vec<Vec<(GF28Element, GF28Element)>>, Error> {
    let mut shares_vec = vec![vec![]; n as usize];
    // create a random polynomial for each byte
    for sbyte in secret {
        let poly = create_polynomial((*sbyte).into(), (k - 1) as usize);
        // evaluate the polynomial at `n` random points for each share
        for share in &mut shares_vec {
            let x = GF28Element::random();
            let y = evaluate_polynomial(x, &poly);
            share.push((x, y));
        }
    }

    Ok(shares_vec)
}

pub fn poly_mul(a: GF28Element, b: GF28Element) -> GF28Element {
    let mid = u8_repr::mul_expanded(a.into(), b.into());
    u8_repr::poly_modulus(mid).into()
}

/// Recall that in a group of size N (here 255) multiplying an element by itself 255 times gives
/// the identity (1). So it stands that if : a**255 = 1 then a**254 = a**-1.
/// This is known as fermat's little theorem. We use it to calculate the inverse of elements.
pub fn poly_inv(a: GF28Element) -> GF28Element {
    //let a2: u8 = a.pow(2);
    let a2 = poly_pow(a, 2);
    let a3 = a2 * a;
    //let a12 = a3.pow(4);
    let a12 = poly_pow(a3, 4);
    let a14 = a12 * a2;
    let a15 = a12 * a3;
    //let a240 = a15.pow(16);
    let a240 = poly_pow(a15, 16);
    a240 * a14
}

pub fn poly_pow(a: GF28Element, pow: usize) -> GF28Element {
    if pow == 0 {
        return 1.into();
    }
    let mut acc = a;
    for _ in 0..(pow - 1) {
        acc = poly_mul(acc, a);
    }
    acc
}

pub fn poly_div(numerator: GF28Element, denomenator: GF28Element) -> GF28Element {
    poly_mul(numerator, poly_inv(denomenator))
}

/// Build a random polynomial of the given `degree` whose constant term is `secret`.
///
/// Coefficients are constant-first: `out[0]` is the secret, `out[degree]` is the
/// leading coefficient. A degree `d` polynomial has `d + 1` coefficients.
pub fn create_polynomial(secret: GF28Element, degree: usize) -> Vec<GF28Element> {
    let mut out = vec![8u8; degree + 1];
    out[0] = secret.into();
    rand::fill(&mut out[1..]);
    out.into_iter().map(GF28Element::from).collect()
}

/// Use [Horner's rule](https://en.wikipedia.org/wiki/Horner%27s_method) to evaluate the polynomial.
pub fn evaluate_polynomial(x: GF28Element, poly: &[GF28Element]) -> GF28Element {
    poly.iter()
        .rev()
        .fold(GF28Element(0), |acc, &a| acc * x + a)
}

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

#[cfg(test)]
pub mod test {
    use crate::evaluate_polynomial;

    use super::{GF28Element, poly_mul, poly_pow, u8_repr};
    macro_rules! chk_poly {
        ($a:expr, $b:expr, $expected:expr) => {
            let res1 = u8_repr::mul_expanded($a, $b);
            assert_eq!(res1, $expected);
            let res2 = u8_repr::mul_expanded($b, $a);
            assert_eq!(res2, $expected);
        };
    }

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
    #[test]
    fn foo() {
        chk_poly!(128, 128, 1 << (7 + 7));
        chk_poly!(128, 1, 128);
        chk_poly!(128, 0, 0);
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
