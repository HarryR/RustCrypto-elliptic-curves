//! Pure Rust implementation of group operations on secp384r1.
//!
//! Curve parameters can be found in [NIST SP 800-186] § G.1.3: Curve P-384.
//!
//! [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final

pub(crate) mod field;
#[cfg(feature = "hash2curve")]
mod hash2curve;
pub(crate) mod scalar;

use self::{field::FieldElement, scalar::Scalar};
use crate::NistP384;
use elliptic_curve::{CurveArithmetic, PrimeCurveArithmetic};
use primeorder::{point_arithmetic, PrimeCurveParams};

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
use primeorder::__risc0::FieldElement384;

/// Elliptic curve point in affine coordinates.
pub type AffinePoint = primeorder::AffinePoint<NistP384>;

/// Elliptic curve point in projective coordinates.
pub type ProjectivePoint = primeorder::ProjectivePoint<NistP384>;

impl CurveArithmetic for NistP384 {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

impl PrimeCurveArithmetic for NistP384 {
    type CurveGroup = ProjectivePoint;
}

/// Adapted from [NIST SP 800-186] § G.1.3: Curve P-384.
///
/// [NIST SP 800-186]: https://csrc.nist.gov/publications/detail/sp/800-186/final
#[cfg(not(all(target_os = "zkvm", target_arch = "riscv32")))]
impl PrimeCurveParams for NistP384 {
    type FieldElement = FieldElement;
    type PointArithmetic = point_arithmetic::EquationAIsMinusThree;

    const EQUATION_A: FieldElement = FieldElement::from_u64(3).neg();
    const EQUATION_B: FieldElement = FieldElement::from_hex("b3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef");

    const GENERATOR: (FieldElement, FieldElement) = (
        FieldElement::from_hex("aa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7"),
        FieldElement::from_hex("3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f"),
    );
}

/// zkVM-accelerated implementation of PrimeCurveParams for NistP384
#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
impl PrimeCurveParams for NistP384 {
    type FieldElement = FieldElement;
    type PointArithmetic = point_arithmetic::EquationAIsMinusThree;

    const EQUATION_A: FieldElement = FieldElement::from_u64(3).neg();
    const EQUATION_B: FieldElement = FieldElement::from_hex("b3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef");

    const GENERATOR: (FieldElement, FieldElement) = (
        FieldElement::from_hex("aa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7"),
        FieldElement::from_hex("3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f"),
    );

    const WIDTH_WORDS: usize = 12;

    fn to_affine_accelerated(p: &ProjectivePoint) -> AffinePoint {
        use primeorder::__risc0::{felt_to_u32_words_le_12, PrimeCurveParams384};

        if p.z.is_zero().into() {
            return AffinePoint::IDENTITY;
        }

        let z = felt_to_u32_words_le_12::<NistP384>(&p.z);
        let mut z_inv = [0u32; 12];
        risc0_bigint2::field::unchecked::modinv_384(
            &z,
            &<NistP384 as PrimeCurveParams384>::PRIME_LE_WORDS,
            &mut z_inv,
        );

        let mut buffer = [0u32; 12];
        let x_buffer = felt_to_u32_words_le_12::<NistP384>(&p.x);
        let y_buffer = felt_to_u32_words_le_12::<NistP384>(&p.y);

        risc0_bigint2::field::unchecked::modmul_384(
            &x_buffer,
            &z_inv,
            &<NistP384 as PrimeCurveParams384>::PRIME_LE_WORDS,
            &mut buffer,
        );
        let x = <NistP384 as PrimeCurveParams384>::from_u32_words_le(buffer);

        risc0_bigint2::field::unchecked::modmul_384(
            &y_buffer,
            &z_inv,
            &<NistP384 as PrimeCurveParams384>::PRIME_LE_WORDS,
            &mut buffer,
        );
        let y = <NistP384 as PrimeCurveParams384>::from_u32_words_le(buffer);

        AffinePoint { x, y, infinity: 0 }
    }

    fn mul_accelerated(p: &ProjectivePoint, k: &elliptic_curve::Scalar<NistP384>) -> ProjectivePoint {
        use primeorder::__risc0::{ec_impl_384, PrimeCurveParams384};
        ec_impl_384::mul::<NistP384>(p, k)
    }
}

#[cfg(all(target_os = "zkvm", target_arch = "riscv32"))]
impl primeorder::__risc0::PrimeCurveParams384 for NistP384 {
    const PRIME_LE_WORDS: [u32; 12] = crate::__risc0::SECP384R1_PRIME;
    const ORDER_LE_WORDS: [u32; 12] = crate::__risc0::SECP384R1_ORDER;
    const EQUATION_A_LE: FieldElement384<NistP384> =
        FieldElement384::new_unchecked(crate::__risc0::SECP384R1_EQUATION_A_LE);
    const EQUATION_B_LE: FieldElement384<NistP384> =
        FieldElement384::new_unchecked(crate::__risc0::SECP384R1_EQUATION_B_LE);

    fn from_u32_words_le(words: [u32; 12]) -> FieldElement {
        FieldElement::from_words_le(words)
    }
}
