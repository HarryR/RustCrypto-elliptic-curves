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
// EC operations for 384-bit curves
// ============================================================================

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

pub(crate) fn projective_to_affine_384<C>(p: &ProjectivePoint<C>) -> ec::AffinePoint<12, C>
where
    C: PrimeCurveParams384,
{
    let aff = p.to_affine();
    affine_to_r0_affine_384(&aff)
}

pub(crate) fn affine_to_projective_384<C>(affine: &ec::AffinePoint<12, C>) -> ProjectivePoint<C>
where
    C: PrimeCurveParams384,
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

pub(crate) fn scalar_to_words_12<C>(s: &Scalar<C>) -> [u32; 12]
where
    C: PrimeCurveParams384,
{
    bytes_to_u32_words_le_12(s.to_repr().as_slice())
}

pub mod ec_impl_384 {
    use super::*;

    pub fn mul<C>(lhs: &ProjectivePoint<C>, rhs: &Scalar<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let scalar = scalar_to_words_12::<C>(rhs);
        let affine = projective_to_affine_384::<C>(lhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        affine.mul(&scalar, &mut result);
        return affine_to_projective_384(&result);
    }

    pub fn add<C>(lhs: &ProjectivePoint<C>, rhs: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_affine_384::<C>(lhs);
        let rhs = projective_to_affine_384::<C>(rhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_384(&result);
    }

    #[inline]
    pub fn add_mixed<C>(lhs: &ProjectivePoint<C>, rhs: &AffinePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let lhs = projective_to_affine_384::<C>(lhs);
        let rhs = affine_to_r0_affine_384(rhs);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        lhs.add(&rhs, &mut result);
        return affine_to_projective_384(&result);
    }

    pub fn double<C>(point: &ProjectivePoint<C>) -> ProjectivePoint<C>
    where
        C: PrimeCurveParams384,
    {
        let point = projective_to_affine_384::<C>(point);

        let mut result = risc0_bigint2::ec::AffinePoint::new_unchecked([0u32; 12], [0u32; 12]);
        point.double(&mut result);
        return affine_to_projective_384(&result);
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
