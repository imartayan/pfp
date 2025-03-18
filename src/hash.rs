//! Rolling hash over windows of a sequence.

use core::iter::{Copied, Zip};
use core::slice::Iter;

/// Integer type for rolling hashes.
pub type HT = u64;

/// Multiplicative constant used for hashing (taken from [`FxHash`](https://docs.rs/rustc-hash/latest/rustc_hash/)).
pub const MUL_HASH: HT = 0xf1357aea2e62a9c5;

/// Hash a single integer.
#[inline(always)]
pub(crate) const fn hash_one(x: HT) -> HT {
    x.wrapping_mul(MUL_HASH)
}

/// Merges the hashes of two sequences.
///
/// `merge_hashes(hash("AB"), hash("CDE"), 3) == hash("ABCDE")`
#[inline(always)]
pub const fn merge_hashes(h1: HT, h2: HT, h2_length: u32) -> HT {
    h1.rotate_left(h2_length) ^ h2
}

/// An iterator of rolling hashes for each window of a sequence.
pub struct RollingHashIterator<'a> {
    in_out: Zip<Copied<Iter<'a, u8>>, Copied<Iter<'a, u8>>>,
    hash: HT,
    rot: u32,
}

impl<'a> RollingHashIterator<'a> {
    /// Builds a new [`RollingHashIterator`] over windows of size `window_size` in `seq`.
    #[inline(always)]
    pub fn new(seq: &'a [u8], window_size: usize) -> Self {
        assert!(window_size > 0, "window_size must be > 0");
        let mut hash = 0u64;
        let mut in_ = seq.iter().copied();
        let out_ = seq.iter().copied();
        in_.by_ref().take(window_size - 1).for_each(|c_in| {
            hash = hash.rotate_left(1) ^ hash_one(c_in as HT);
        });
        Self {
            in_out: in_.zip(out_),
            hash,
            rot: (window_size - 1) as u32,
        }
    }
}

impl Iterator for RollingHashIterator<'_> {
    type Item = HT;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.in_out.next().map(|(c_in, c_out)| {
            let res = self.hash.rotate_left(1) ^ hash_one(c_in as HT);
            self.hash = res ^ hash_one(c_out as HT).rotate_left(self.rot);
            res
        })
    }

    // Having an exact size hint helps collecting the iterator faster.
    #[inline(always)]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.in_out.size_hint()
    }

    #[inline(always)]
    fn count(self) -> usize {
        self.in_out.count()
    }
}

impl ExactSizeIterator for RollingHashIterator<'_> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_repeated_seq() {
        const LEN: usize = 1000;
        let mut seq = [0u8; LEN];
        fastrand::fill(&mut seq);
        let mut repeated = Vec::from(seq);
        repeated.extend_from_slice(&seq);
        repeated.extend_from_slice(&seq);

        for window_size in 1..=50 {
            let hashes: Vec<_> = RollingHashIterator::new(&repeated, window_size).collect();
            assert_eq!(hashes[0..LEN], hashes[LEN..(2 * LEN)]);
        }
    }
}
