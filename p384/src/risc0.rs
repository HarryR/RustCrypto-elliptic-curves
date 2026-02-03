use risc0_bigint2::ec::{Curve, WeierstrassCurve, EC_384_WIDTH_WORDS};

/// The secp384r1 curve's prime field characteristic (little-endian)
/// p = 2^{384} − 2^{128} − 2^{96} + 2^{32} − 1
/// = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff
pub(crate) const SECP384R1_PRIME: [u32; EC_384_WIDTH_WORDS] = [
    0xFFFFFFFF, 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];

/// The secp384r1 curve's order (little-endian)
/// n = ffffffff ffffffff ffffffff ffffffff ffffffff ffffffff c7634d81 f4372ddf 581a0db2 48b0a77a ecec196a ccc52973
pub(crate) const SECP384R1_ORDER: [u32; EC_384_WIDTH_WORDS] = [
    0xCCC52973, 0xECEC196A, 0x48B0A77A, 0x581A0DB2, 0xF4372DDF, 0xC7634D81,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];

/// Coefficient a = -3 (mod p) in little-endian
/// a = p - 3 = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000fffffffc
pub(crate) const SECP384R1_EQUATION_A_LE: [u32; EC_384_WIDTH_WORDS] = [
    0xFFFFFFFC, 0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFF,
    0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
];

/// Coefficient b in little-endian
/// b = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef
pub(crate) const SECP384R1_EQUATION_B_LE: [u32; EC_384_WIDTH_WORDS] = [
    0xD3EC2AEF, 0x2A85C8ED, 0x8A2ED19D, 0xC656398D, 0x5013875A, 0x0314088F,
    0xFE814112, 0x181D9C6E, 0xE3F82D19, 0x988E056B, 0xE23EE7E4, 0xB3312FA7,
];

const SECP384R1_CURVE: &WeierstrassCurve<EC_384_WIDTH_WORDS> =
    &WeierstrassCurve::<EC_384_WIDTH_WORDS>::new(
        SECP384R1_PRIME,
        SECP384R1_EQUATION_A_LE,
        SECP384R1_EQUATION_B_LE,
    );

impl Curve<EC_384_WIDTH_WORDS> for crate::NistP384 {
    const CURVE: &'static WeierstrassCurve<EC_384_WIDTH_WORDS> = SECP384R1_CURVE;
}
