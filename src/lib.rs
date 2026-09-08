#![expect(clippy::suspicious_arithmetic_impl)]
//! Implementation of GF(2^8)
//! We need a new polynomial for every byte

use core::{
    convert::From,
    ops::{Add, Div, Mul, Sub},
};

/// MODULUS used to reduce polynomial to degree 8: m(x) = x^8 + x^4 + x^3 + x + 1
const MODULUS: u16 = 0b100011011;

const DEGREE: usize = 8;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GF28Element(pub(crate) u8);

impl GF28Element {
    fn inv(self) -> GF28Element {
        poly_inv(self).into()
    }
}

/// Multiply out two <= degree 7 polynomials into a <= degree 14 polynomial
pub fn u8_poly_mul_expanded(a: u8, b: u8) -> u16 {
    let mut product: u16 = 0;
    for i in 0..DEGREE {
        if (b >> i) & 1 != 0 {
            product ^= (a as u16) << i;
        }
    }
    product
}

/// Reduce a degree 14 polynomial to degree 7
pub fn u8_poly_modulus(poly: u16) -> u8 {
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

pub fn polynomial_mul(a: GF28Element, b: GF28Element) -> GF28Element {
    let mid = u8_poly_mul_expanded(a.into(), b.into());
    u8_poly_modulus(mid).into()
}
pub fn polynomial_mul2(a: u8, b: u8) -> u8 {
    let mid = u8_poly_mul_expanded(a, b);
    u8_poly_modulus(mid)
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
    for _ in 0..pow {
        acc = polynomial_mul(acc, a);
    }
    acc
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

    // a * b
    fn mul(self, rhs: Self) -> Self::Output {
        polynomial_mul(self.into(), rhs.into()).into()
    }
}
impl Div for GF28Element {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.inv()
    }
}

#[cfg(test)]
pub mod test {
    use crate::{GF28Element, polynomial_mul, u8_poly_mul_expanded};

    macro_rules! chk_poly {
        ($a:expr, $b:expr, $expected:expr) => {
            let res1 = u8_poly_mul_expanded($a, $b);
            assert_eq!(res1, $expected);
            let res2 = u8_poly_mul_expanded($b, $a);
            assert_eq!(res2, $expected);
        };
    }

    macro_rules! chk_mul {
        ($a:expr, $b:expr, $expected:expr) => {
            let a = GF28Element($a);
            let b = GF28Element($b);
            let c = GF28Element($expected);

            let res1 = polynomial_mul(a, b);
            assert_eq!(res1, c);
            let res2 = polynomial_mul(b, a);
            assert_eq!(res2, c);

            //
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
