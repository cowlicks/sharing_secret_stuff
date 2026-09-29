# share_secret_stuff

> ## ⚠️ WARNING ⚠️
>
> Do not use this code for anything real. This crate is for educational purposes only.

Shamir's Secret Sharing over GF(2⁸), written from scratch as a learning exercise.

Split a secret into some number of shares (we call this number `n`). Collecting a threshold number
of those shares (we call this number `k`) reconstructs the secret. Collecting fewer than `k` reveals
nothing at all — not "nothing practical", nothing. Any `k - 1` shares are consistent with every
possible secret, in equal measure.

## Example

```rust
use share_secret_stuff::{split_secret, reconstruct_secret};

let secret = b"Hello, world!"; // My favorite secret
let n_shares: u8 = 12;         // Max shares is 255
let k_threshold: u8 = 7;       // Threshold must not exceed the number of shares

// We get `n_shares` shares. Each one holds a single point from each byte's polynomial.
let shares = split_secret(secret, k_threshold, n_shares).unwrap();
let result = reconstruct_secret(&shares[..(k_threshold as usize)]).unwrap();
assert_eq!(result, secret);
```

## How it works

The secret is split one byte at a time. Each byte gets its own random polynomial of degree `k - 1`
whose constant term is that byte, so `f(0)` is the secret and every other point on the curve reveals
nothing on its own. Each shareholder is handed one point from each of those polynomials, all at the
same `x`.

All the arithmetic happens in GF(2⁸), the finite field of 256 elements. A byte holds the
coefficients of a polynomial of degree 7 or less over GF(2), so addition and subtraction are both
just XOR, and multiplication is polynomial multiplication reduced modulo AES's irreducible
polynomial `x⁸ + x⁴ + x³ + x + 1`.


## Tests

```
cargo test
```
