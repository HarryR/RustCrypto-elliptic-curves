use crate::FieldBytes;
use crate::{affine::AffinePoint, projective::ProjectivePoint};
use core::marker::PhantomData;
use core::ops::Deref;
use elliptic_curve::consts::{U32, U48};
use elliptic_curve::generic_array::GenericArray;
use elliptic_curve::subtle::{Choice, ConditionallySelectable};
use elliptic_curve::{PrimeField, Scalar};
use risc0_bigint2::ec;

use crate::PrimeCurveParams;

// ============================================================================
// 256-bit Field Element
// ============================================================================

/// Representation of a 256-bit field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct FieldElement256<C> {
    pub data: [u32; 8],
    _phantom: PhantomData<C>,
}

impl<C> Deref for FieldElement256<C> {
    type Target = [u32; 8];

    fn deref(&self) -> &[u32; 8] {
        &self.data
    }
}

impl<C: elliptic_curve::Curve<FieldBytesSize = U32>> From<&FieldBytes<C>> for FieldElement256<C> {
    fn from(data: &FieldBytes<C>) -> Self {
        let mut words = [0u32; 8];

        // Process 4 bytes at a time to create little-endian u32 words
        for (i, chunk) in data.chunks(4).enumerate() {
            // Convert each big-endian chunk to a little-endian u32
            words[7 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
        }

        Self::new_unchecked(words)
    }
}

impl<C: PrimeCurveParams<FieldBytesSize = U32>> From<FieldElement256<C>> for GenericArray<u8, C::FieldBytesSize> {
    fn from(data: FieldElement256<C>) -> Self {
        let bytes_slice = bytemuck::cast_slice::<u32, u8>(&data.data);
        GenericArray::from_iter(bytes_slice.iter().copied().rev())
    }
}

impl<C: Copy> ConditionallySelectable for FieldElement256<C> {
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut output = *a;
        output.conditional_assign(b, choice);
        output
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for (a_i, b_i) in self.data.iter_mut().zip(other.data.iter()) {
            a_i.conditional_assign(b_i, choice)
        }
    }
}

impl<C> FieldElement256<C> {
    pub const fn new_unchecked(data: [u32; 8]) -> Self {
        Self {
            data,
            _phantom: PhantomData,
        }
    }
}

impl<C> FieldElement256<C>
where
    C: PrimeCurveParams256,
{
    #[inline]
    pub fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modmul_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn mul(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modadd_256(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    /// Calculate the square root of the field element, using the provided buffer as scratch, and
    /// writing the result to the `result` parameter.
    pub(crate) fn sqrt_unchecked(&self, scratch: &mut Self, result: &mut Self) {
        // New buffers to keep temporary values.
        let mut scratch_1 = Self::default();
        let mut scratch_2 = Self::default();

        // let t11 = self.mul(&self.square());
        self.square(scratch);
        self.mul_unchecked(scratch, result);

        // result = t11
        // let t1111 = t11.mul(&t11.sqn(2));
        result.sqn(2, (&mut scratch_1, &mut scratch_2), scratch);
        result.mul_unchecked(scratch, &mut scratch_1);

        // scratch_1 = t1111
        // let t11111111 = t1111.mul(&t1111.sqn(4));
        scratch_1.sqn(4, (&mut scratch_2, result), scratch);
        scratch_1.mul_unchecked(scratch, result);

        // result = t11111111
        // let x16 = t11111111.sqn(8).mul(&t11111111);
        result.sqn(8, (&mut scratch_1, &mut scratch_2), scratch);
        result.mul_unchecked(scratch, &mut scratch_1);

        // scratch_1 = x16
        // let sqrt = x16
        //     .sqn(16)
        scratch_1.sqn(16, (&mut scratch_2, result), scratch);
        //     .mul(&x16)
        scratch.mul_unchecked(&scratch_1, result);
        //     .sqn(32)
        result.sqn(32, (&mut scratch_1, &mut scratch_2), scratch);
        //     .mul(self)
        scratch.mul_unchecked(self, &mut scratch_1);
        //     .sqn(96)
        scratch_1.sqn(96, (&mut scratch_2, result), scratch);
        //     .mul(self)
        scratch.mul_unchecked(self, &mut scratch_2);
        //     .sqn(94);
        // Last result is written to the result buffer.
        scratch_2.sqn(94, (&mut scratch_1, scratch), result);
    }

    /// Returns self^(2^n) mod p.
    ///
    /// This implementation is designed to avoid any memcpy of the buffers for intermediate ops.
    fn sqn(&self, n: usize, scratch: (&mut Self, &mut Self), result: &mut Self) {
        let mut x = scratch.0;
        let mut buffer = scratch.1;

        if n == 1 {
            // self^(2^1) = self^2
            self.square(result);
            return;
        } else if n == 0 {
            // self^(1) = self
            *result = *self;
            return;
        }

        // write value to a scratch buffer.
        self.square(x);

        // Square n - 2 times.
        let mut i = 2;
        while i < n {
            x.square(buffer);
            i += 1;
            // Swap scratch buffers, to set x to the squared value.
            core::mem::swap(&mut x, &mut buffer);
        }

        // Write final square to result buffer.
        x.square(result);
    }

    /// Returns self^2 mod p
    pub(crate) fn square(&self, result: &mut Self) {
        self.mul_unchecked(self, result);
    }
}

// ============================================================================
// 384-bit Field Element
// ============================================================================

/// Representation of a 384-bit field element in raw bytes form. This is not in montgomery form.
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq)]
pub struct FieldElement384<C> {
    pub data: [u32; 12],
    _phantom: PhantomData<C>,
}

impl<C> Deref for FieldElement384<C> {
    type Target = [u32; 12];

    fn deref(&self) -> &[u32; 12] {
        &self.data
    }
}

impl<C: elliptic_curve::Curve<FieldBytesSize = U48>> From<&FieldBytes<C>> for FieldElement384<C> {
    fn from(data: &FieldBytes<C>) -> Self {
        let mut words = [0u32; 12];

        // Process 4 bytes at a time to create little-endian u32 words
        for (i, chunk) in data.chunks(4).enumerate() {
            // Convert each big-endian chunk to a little-endian u32
            words[11 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
        }

        Self::new_unchecked(words)
    }
}

impl<C: PrimeCurveParams<FieldBytesSize = U48>> From<FieldElement384<C>> for GenericArray<u8, C::FieldBytesSize> {
    fn from(data: FieldElement384<C>) -> Self {
        let bytes_slice = bytemuck::cast_slice::<u32, u8>(&data.data);
        GenericArray::from_iter(bytes_slice.iter().copied().rev())
    }
}

impl<C: Copy> ConditionallySelectable for FieldElement384<C> {
    #[inline]
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let mut output = *a;
        output.conditional_assign(b, choice);
        output
    }

    fn conditional_assign(&mut self, other: &Self, choice: Choice) {
        for (a_i, b_i) in self.data.iter_mut().zip(other.data.iter()) {
            a_i.conditional_assign(b_i, choice)
        }
    }
}

impl<C> FieldElement384<C> {
    pub const fn new_unchecked(data: [u32; 12]) -> Self {
        Self {
            data,
            _phantom: PhantomData,
        }
    }
}

impl<C> FieldElement384<C>
where
    C: PrimeCurveParams384,
{
    #[inline]
    pub fn mul_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modmul_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn mul(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::modmul_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    #[inline]
    pub fn add_unchecked(&self, rhs: &Self, result: &mut Self) {
        risc0_bigint2::field::unchecked::modadd_384(
            &self.data,
            &rhs.data,
            &C::PRIME_LE_WORDS,
            &mut result.data,
        );
    }

    /// Returns self^2 mod p
    pub(crate) fn square(&self, result: &mut Self) {
        self.mul_unchecked(self, result);
    }
}

// ============================================================================
// Helper functions
// ============================================================================

fn bytes_to_u32_words_le_8(bytes: &[u8]) -> [u32; 8] {
    let mut words = [0u32; 8];

    // Process 4 bytes at a time to create little-endian u32 words
    for (i, chunk) in bytes.chunks(4).enumerate() {
        // Convert each big-endian chunk to a little-endian u32
        words[7 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
    }

    words
}

fn bytes_to_u32_words_le_12(bytes: &[u8]) -> [u32; 12] {
    let mut words = [0u32; 12];

    // Process 4 bytes at a time to create little-endian u32 words
    for (i, chunk) in bytes.chunks(4).enumerate() {
        // Convert each big-endian chunk to a little-endian u32
        words[11 - i] = u32::from_be_bytes(chunk.try_into().unwrap());
    }

    words
}

pub fn felt_to_u32_words_le_8<C>(data: &C::FieldElement) -> [u32; 8]
where
    C: PrimeCurveParams256,
{
    bytes_to_u32_words_le_8(data.to_repr().as_slice())
}

pub fn felt_to_u32_words_le_12<C>(data: &C::FieldElement) -> [u32; 12]
where
    C: PrimeCurveParams384,
{
    bytes_to_u32_words_le_12(data.to_repr().as_slice())
}

// ============================================================================
// EC operations for 256-bit curves
// ============================================================================

#[inline]
fn affine_to_r0_affine_256<C>(affine: &AffinePoint<C>) -> ec::AffinePoint<8, C>
where
    C: PrimeCurveParams256,
{
    if bool::from(affine.is_identity()) {
        return ec::AffinePoint::IDENTITY;
    }

    let x = felt_to_u32_words_le_8::<C>(&affine.x);
    let y = felt_to_u32_words_le_8::<C>(&affine.y);
    ec::AffinePoint::new_unchecked(x, y)
}

pub(crate) fn projective_to_affine_256<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<8, C>
where
    C: PrimeCurveParams256,
{
    let aff = p.to_affine();
    affine_to_r0_affine_256(&aff)
}

pub(crate) fn affine_to_projective_256<C>(affine: &ec::AffinePoint<8, C>) -> ProjectivePoint<C>
where
    C: PrimeCurveParams256,
{
    if let Some(value) = affine.as_u32s() {
        // This should only not be within the modulus with a malicious host, panic in that case.
        let x = C::from_u32_words_le(value[0]);
        let y = C::from_u32_words_le(value[1]);

        let affine = AffinePoint { x, y, infinity: 0 };
        ProjectivePoint::from(affine)
    } else {
        ProjectivePoint::IDENTITY
    }
}

pub(crate) fn scalar_to_words_8<C>(s: &Scalar<C>) -> [u32; 8]
where
    C: PrimeCurveParams256,
{
    bytes_to_u32_words_le_8(s.to_repr().as_slice())
}

pub(crate) mod ec_impl_256 {
    use super::*;

    pub(crate) fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let scalar = scalar_to_words_8::<C>(rhs);
        let affine = projective_to_affine_256::<C>(lhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        affine.mul(&scalar, &mut result);
        return affine_to_projective_256(&result);
    }

    pub(crate) fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let lhs = projective_to_affine_256::<C>(lhs);
        let rhs = projective_to_affine_256::<C>(rhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_256(&result);
    }

    #[inline]
    pub(crate) fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let lhs = projective_to_affine_256::<C>(lhs);
        let rhs = affine_to_r0_affine_256(rhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_256(&result);
    }

    pub(crate) fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams256,
    {
        let point = projective_to_affine_256::<C>(point);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 8], [0u32; 8]);
        point.double(&mut result);
        return affine_to_projective_256(&result);
    }
}

// ============================================================================
// EC operations for 384-bit curves - Proper Projective Implementation
// ============================================================================

/// Projective point in Jacobian coordinates stored as raw words (standard form, not Montgomery).
/// This avoids Montgomery conversions during EC operations.
#[derive(Clone, Copy)]
pub struct Projective384 {
    pub x: [u32; 12],
    pub y: [u32; 12],
    pub z: [u32; 12],
}

impl Projective384 {
    /// Identity point (0, 1, 0)
    pub const IDENTITY: Self = Self {
        x: [0; 12],
        y: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        z: [0; 12],
    };

    /// Check if this is the identity point
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.z == [0u32; 12]
    }

    /// Point doubling using standard projective coordinates for a = -3
    /// Implements RCB 2015 Algorithm 6
    pub fn double(&self, prime: &[u32; 12], b: &[u32; 12]) -> Self {
        if self.is_identity() {
            return Self::IDENTITY;
        }

        let mut tmp = [0u32; 12];
        let mut tmp2 = [0u32; 12];

        // xx = X^2
        let mut xx = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.x, &self.x, prime, &mut xx);

        // yy = Y^2
        let mut yy = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.y, &self.y, prime, &mut yy);

        // zz = Z^2
        let mut zz = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.z, &self.z, prime, &mut zz);

        // xy2 = 2*X*Y
        let mut xy2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.x, &self.y, prime, &mut xy2);
        risc0_bigint2::field::unchecked::modadd_384(&xy2, &xy2, prime, &mut tmp);
        xy2 = tmp;

        // xz2 = 2*X*Z
        let mut xz2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.x, &self.z, prime, &mut xz2);
        risc0_bigint2::field::unchecked::modadd_384(&xz2, &xz2, prime, &mut tmp);
        xz2 = tmp;

        // bzz_part = b*zz - xz2
        let mut bzz = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(b, &zz, prime, &mut bzz);
        let mut bzz_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&bzz, &xz2, prime, &mut bzz_part);

        // bzz3_part = 3 * bzz_part
        let mut bzz3_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&bzz_part, &bzz_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &bzz_part, prime, &mut bzz3_part);

        // yy_m_bzz3 = yy - bzz3_part
        let mut yy_m_bzz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&yy, &bzz3_part, prime, &mut yy_m_bzz3);

        // yy_p_bzz3 = yy + bzz3_part
        let mut yy_p_bzz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&yy, &bzz3_part, prime, &mut yy_p_bzz3);

        // y_frag = yy_p_bzz3 * yy_m_bzz3
        let mut y_frag = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&yy_p_bzz3, &yy_m_bzz3, prime, &mut y_frag);

        // x_frag = yy_m_bzz3 * xy2
        let mut x_frag = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&yy_m_bzz3, &xy2, prime, &mut x_frag);

        // zz3 = 3 * zz
        let mut zz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&zz, &zz, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &zz, prime, &mut zz3);

        // bxz2_part = b*xz2 - zz3 - xx
        let mut bxz2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(b, &xz2, prime, &mut bxz2);
        risc0_bigint2::field::unchecked::modsub_384(&bxz2, &zz3, prime, &mut tmp);
        let mut bxz2_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&tmp, &xx, prime, &mut bxz2_part);

        // bxz6_part = 3 * bxz2_part
        let mut bxz6_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&bxz2_part, &bxz2_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &bxz2_part, prime, &mut bxz6_part);

        // xx3_m_zz3 = 3*xx - zz3
        let mut xx3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&xx, &xx, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &xx, prime, &mut xx3);
        let mut xx3_m_zz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&xx3, &zz3, prime, &mut xx3_m_zz3);

        // y3 = y_frag + xx3_m_zz3 * bxz6_part
        risc0_bigint2::field::unchecked::modmul_384(&xx3_m_zz3, &bxz6_part, prime, &mut tmp);
        let mut y3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&y_frag, &tmp, prime, &mut y3);

        // yz2 = 2*Y*Z
        let mut yz2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.y, &self.z, prime, &mut yz2);
        risc0_bigint2::field::unchecked::modadd_384(&yz2, &yz2, prime, &mut tmp);
        yz2 = tmp;

        // x3 = x_frag - bxz6_part * yz2
        risc0_bigint2::field::unchecked::modmul_384(&bxz6_part, &yz2, prime, &mut tmp);
        let mut x3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&x_frag, &tmp, prime, &mut x3);

        // z3 = 4 * yz2 * yy
        risc0_bigint2::field::unchecked::modmul_384(&yz2, &yy, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &tmp, prime, &mut tmp2);
        let mut z3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&tmp2, &tmp2, prime, &mut z3);

        Self { x: x3, y: y3, z: z3 }
    }

    /// Point addition using standard projective coordinates for a = -3
    /// Implements RCB 2015 Algorithm 4
    pub fn add(&self, other: &Self, prime: &[u32; 12], b: &[u32; 12]) -> Self {
        if self.is_identity() {
            return *other;
        }
        if other.is_identity() {
            return *self;
        }

        let mut tmp = [0u32; 12];
        let mut tmp2 = [0u32; 12];

        // xx = X1 * X2
        let mut xx = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.x, &other.x, prime, &mut xx);

        // yy = Y1 * Y2
        let mut yy = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.y, &other.y, prime, &mut yy);

        // zz = Z1 * Z2
        let mut zz = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.z, &other.z, prime, &mut zz);

        // xy_pairs = (X1+Y1)*(X2+Y2) - xx - yy
        let mut x1_p_y1 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&self.x, &self.y, prime, &mut x1_p_y1);
        let mut x2_p_y2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&other.x, &other.y, prime, &mut x2_p_y2);
        let mut xy_pairs = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&x1_p_y1, &x2_p_y2, prime, &mut xy_pairs);
        risc0_bigint2::field::unchecked::modsub_384(&xy_pairs, &xx, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modsub_384(&tmp, &yy, prime, &mut xy_pairs);

        // yz_pairs = (Y1+Z1)*(Y2+Z2) - yy - zz
        let mut y1_p_z1 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&self.y, &self.z, prime, &mut y1_p_z1);
        let mut y2_p_z2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&other.y, &other.z, prime, &mut y2_p_z2);
        let mut yz_pairs = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&y1_p_z1, &y2_p_z2, prime, &mut yz_pairs);
        risc0_bigint2::field::unchecked::modsub_384(&yz_pairs, &yy, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modsub_384(&tmp, &zz, prime, &mut yz_pairs);

        // xz_pairs = (X1+Z1)*(X2+Z2) - xx - zz
        let mut x1_p_z1 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&self.x, &self.z, prime, &mut x1_p_z1);
        let mut x2_p_z2 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&other.x, &other.z, prime, &mut x2_p_z2);
        let mut xz_pairs = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&x1_p_z1, &x2_p_z2, prime, &mut xz_pairs);
        risc0_bigint2::field::unchecked::modsub_384(&xz_pairs, &xx, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modsub_384(&tmp, &zz, prime, &mut xz_pairs);

        // bzz_part = xz_pairs - b*zz (uses b, not 3*b)
        let mut bzz = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(b, &zz, prime, &mut bzz);
        let mut bzz_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&xz_pairs, &bzz, prime, &mut bzz_part);

        // bzz3_part = 3 * bzz_part (double + add)
        let mut bzz3_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&bzz_part, &bzz_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &bzz_part, prime, &mut bzz3_part);

        // yy_m_bzz3 = yy - bzz3_part
        let mut yy_m_bzz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&yy, &bzz3_part, prime, &mut yy_m_bzz3);

        // yy_p_bzz3 = yy + bzz3_part
        let mut yy_p_bzz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&yy, &bzz3_part, prime, &mut yy_p_bzz3);

        // zz3 = 3 * zz
        let mut zz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&zz, &zz, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &zz, prime, &mut zz3);

        // bxz_part = b*xz_pairs - zz3 - xx (uses b, not 3*b)
        let mut bxz = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(b, &xz_pairs, prime, &mut bxz);
        risc0_bigint2::field::unchecked::modsub_384(&bxz, &zz3, prime, &mut tmp);
        let mut bxz_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&tmp, &xx, prime, &mut bxz_part);

        // bxz3_part = 3 * bxz_part
        let mut bxz3_part = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&bxz_part, &bxz_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &bxz_part, prime, &mut bxz3_part);

        // xx3_m_zz3 = 3*xx - zz3
        let mut xx3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modadd_384(&xx, &xx, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&tmp, &xx, prime, &mut xx3);
        let mut xx3_m_zz3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modsub_384(&xx3, &zz3, prime, &mut xx3_m_zz3);

        // x3 = yy_p_bzz3 * xy_pairs - yz_pairs * bxz3_part
        let mut x3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&yy_p_bzz3, &xy_pairs, prime, &mut x3);
        risc0_bigint2::field::unchecked::modmul_384(&yz_pairs, &bxz3_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modsub_384(&x3, &tmp, prime, &mut tmp2);
        x3 = tmp2;

        // y3 = yy_p_bzz3 * yy_m_bzz3 + xx3_m_zz3 * bxz3_part
        let mut y3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&yy_p_bzz3, &yy_m_bzz3, prime, &mut y3);
        risc0_bigint2::field::unchecked::modmul_384(&xx3_m_zz3, &bxz3_part, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&y3, &tmp, prime, &mut tmp2);
        y3 = tmp2;

        // z3 = yy_m_bzz3 * yz_pairs + xy_pairs * xx3_m_zz3
        let mut z3 = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&yy_m_bzz3, &yz_pairs, prime, &mut z3);
        risc0_bigint2::field::unchecked::modmul_384(&xy_pairs, &xx3_m_zz3, prime, &mut tmp);
        risc0_bigint2::field::unchecked::modadd_384(&z3, &tmp, prime, &mut tmp2);
        z3 = tmp2;

        Self { x: x3, y: y3, z: z3 }
    }

    /// Convert to affine coordinates (requires one inversion)
    /// Standard projective: x = X/Z, y = Y/Z
    pub fn to_affine(&self, prime: &[u32; 12]) -> ([u32; 12], [u32; 12], bool) {
        if self.is_identity() {
            return ([0; 12], [0; 12], true);
        }

        // z_inv = Z^(-1)
        let mut z_inv = [0u32; 12];
        risc0_bigint2::field::unchecked::modinv_384(&self.z, prime, &mut z_inv);

        // x = X * Z^(-1)
        let mut x = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.x, &z_inv, prime, &mut x);

        // y = Y * Z^(-1)
        let mut y = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&self.y, &z_inv, prime, &mut y);

        (x, y, false)
    }
}

/// Convert from primeorder ProjectivePoint to our Projective384
pub fn projective_to_proj384<C>(p: &ProjectivePoint<C>) -> Projective384
where
    C: PrimeCurveParams384,
{
    Projective384 {
        x: felt_to_u32_words_le_12::<C>(&p.x),
        y: felt_to_u32_words_le_12::<C>(&p.y),
        z: felt_to_u32_words_le_12::<C>(&p.z),
    }
}

/// Convert from Projective384 back to primeorder ProjectivePoint
pub fn proj384_to_projective<C>(p: &Projective384) -> ProjectivePoint<C>
where
    C: PrimeCurveParams384,
{
    if p.is_identity() {
        return ProjectivePoint::IDENTITY;
    }
    ProjectivePoint {
        x: C::from_u32_words_le(p.x),
        y: C::from_u32_words_le(p.y),
        z: C::from_u32_words_le(p.z),
    }
}

#[inline]
fn affine_to_r0_affine_384<C>(affine: &AffinePoint<C>) -> ec::AffinePoint<12, C>
where
    C: PrimeCurveParams384,
{
    if bool::from(affine.is_identity()) {
        return ec::AffinePoint::IDENTITY;
    }

    let x = felt_to_u32_words_le_12::<C>(&affine.x);
    let y = felt_to_u32_words_le_12::<C>(&affine.y);
    ec::AffinePoint::new_unchecked(x, y)
}

pub(crate) fn scalar_to_words_12<C>(s: &Scalar<C>) -> [u32; 12]
where
    C: PrimeCurveParams384,
{
    bytes_to_u32_words_le_12(s.to_repr().as_slice())
}

pub mod ec_impl_384 {
    use super::*;
    use elliptic_curve::Field;

    /// Scalar multiplication using risc0-bigint2's accelerated EC operations.
    /// This delegates to the AffinePoint::mul method which uses hardware-accelerated
    /// double and add operations via precompiled blobs.
    ///
    /// Optimization: converts directly from ProjectivePoint to risc0 AffinePoint
    /// without going through primeorder AffinePoint's Montgomery form.
    pub fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let scalar = scalar_to_words_12::<C>(rhs);
        let affine = projective_to_r0_affine_direct::<C>(lhs);

        if affine.is_identity() {
            return ProjectivePoint::IDENTITY;
        }

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        affine.mul(&scalar, &mut result);
        return r0_affine_to_projective_direct(&result);
    }

    /// z = 1 in standard form (little-endian)
    const Z_ONE: [u32; 12] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    /// Convert ProjectivePoint directly to risc0 AffinePoint.
    /// This avoids the intermediate primeorder AffinePoint and its Montgomery form conversions.
    /// Optimized: skips modinv/modmul when z=1 (common case for points from affine).
    fn projective_to_r0_affine_direct<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<12, C>
    where
        C: PrimeCurveParams384,
    {
        // Convert projective coordinates from Montgomery to standard form
        let x_std = felt_to_u32_words_le_12::<C>(&p.x);
        let y_std = felt_to_u32_words_le_12::<C>(&p.y);
        let z_std = felt_to_u32_words_le_12::<C>(&p.z);

        // Check if this is the identity (z == 0)
        if z_std == [0u32; 12] {
            return ec::AffinePoint::IDENTITY;
        }

        // Fast path: if z == 1, no division needed
        if z_std == Z_ONE {
            return ec::AffinePoint::new_unchecked(x_std, y_std);
        }

        // Compute affine coordinates: x = X/Z, y = Y/Z
        let prime = &C::PRIME_LE_WORDS;
        let mut z_inv = [0u32; 12];
        risc0_bigint2::field::unchecked::modinv_384(&z_std, prime, &mut z_inv);

        let mut x_aff = [0u32; 12];
        let mut y_aff = [0u32; 12];
        risc0_bigint2::field::unchecked::modmul_384(&x_std, &z_inv, prime, &mut x_aff);
        risc0_bigint2::field::unchecked::modmul_384(&y_std, &z_inv, prime, &mut y_aff);

        ec::AffinePoint::new_unchecked(x_aff, y_aff)
    }

    /// Convert risc0 AffinePoint directly back to ProjectivePoint.
    /// This avoids the intermediate primeorder AffinePoint.
    fn r0_affine_to_projective_direct<C>(affine: &ec::AffinePoint<12, C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        if let Some(coords) = affine.as_u32s() {
            // Convert from standard form to Montgomery form
            let x = C::from_u32_words_le(coords[0]);
            let y = C::from_u32_words_le(coords[1]);

            // Create ProjectivePoint with z = 1 (in Montgomery form)
            ProjectivePoint {
                x,
                y,
                z: C::FieldElement::ONE,
            }
        } else {
            ProjectivePoint::IDENTITY
        }
    }

    pub fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_r0_affine_direct::<C>(lhs);
        let rhs = projective_to_r0_affine_direct::<C>(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return r0_affine_to_projective_direct(&result);
    }

    #[inline]
    pub fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_r0_affine_direct::<C>(lhs);
        let rhs = affine_to_r0_affine_384(rhs);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return r0_affine_to_projective_direct(&result);
    }

    pub fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let point = projective_to_r0_affine_direct::<C>(point);

        let mut result = ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        point.double(&mut result);
        return r0_affine_to_projective_direct(&result);
    }
}

// ============================================================================
// Trait aliases for 256-bit and 384-bit curves
// ============================================================================

/// Alias trait for 256-bit PrimeCurveParams
pub trait PrimeCurveParams256: PrimeCurveParams<FieldBytesSize = U32> + ec::Curve<{ec::EC_256_WIDTH_WORDS}> {
    const PRIME_LE_WORDS: [u32; 8];
    const ORDER_LE_WORDS: [u32; 8];
    const EQUATION_A_LE: FieldElement256<Self>;
    const EQUATION_B_LE: FieldElement256<Self>;
    fn from_u32_words_le(words: [u32; 8]) -> Self::FieldElement;
}

/// Alias trait for 384-bit PrimeCurveParams
pub trait PrimeCurveParams384: PrimeCurveParams<FieldBytesSize = U48> + ec::Curve<{ec::EC_384_WIDTH_WORDS}> {
    const PRIME_LE_WORDS: [u32; 12];
    const ORDER_LE_WORDS: [u32; 12];
    const EQUATION_A_LE: FieldElement384<Self>;
    const EQUATION_B_LE: FieldElement384<Self>;
    fn from_u32_words_le(words: [u32; 12]) -> Self::FieldElement;
}
