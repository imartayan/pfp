//! Prefix-free parsing.

use crate::hash::{hash_one, merge_hashes, RollingHashIterator, HT};
use core::ops::Range;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Integer type for phrase length.
pub type LT = u32;

const SINGLE_THREAD_THRESHOLD: usize = 1 << 16;

/// Result of the prefix-free parsing.
/// - `prefix` contains the hash of the incomplete starting phrase and its length
/// - `suffix` contains the hash of the the incomplete ending phrase and its length
/// - `phrases` contains the hashes of the complete phrases, in their order of appearance
/// - `phrases_len` contains the lengths of the corresponding phrases
pub struct Parse {
    pub prefix: (HT, LT),
    pub suffix: (HT, LT),
    pub phrases: Vec<HT>,
    pub phrases_len: Vec<LT>,
}

impl Parse {
    /// Number of complete phrases in the parse.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.phrases.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.phrases.is_empty()
    }

    // TODO methods to generate/extend a dictionary
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters.
///
/// This is a wrapper for [`parse_seq_with_vec`] that allocates new vectors to store the parse.
#[inline(always)]
pub fn parse_seq(seq: &[u8], w: usize, p: usize) -> Parse {
    let mut phrases = Vec::new();
    let mut phrases_len = Vec::new();
    let (prefix, suffix) = parse_seq_with_vec(seq, w, p, &mut phrases, &mut phrases_len);

    Parse {
        prefix,
        suffix,
        phrases,
        phrases_len,
    }
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters, storing the resulting phrases in existing `Vec`.
pub fn parse_seq_with_vec(
    seq: &[u8],
    w: usize,
    p: usize,
    vec_hash: &mut Vec<HT>,
    vec_len: &mut Vec<LT>,
) -> ((HT, LT), (HT, LT)) {
    let expected_len = seq.len() / p * 11 / 10;
    let hash_bound = HT::MAX / p as HT + 1;
    vec_hash.reserve(expected_len);
    vec_len.reserve(expected_len);
    let mut window_hash_it = RollingHashIterator::new(seq, w);
    let mut prefix = (0, 0);
    let mut phrase_hash = 0u64;
    let mut phrase_len = 0;

    for window_hash in window_hash_it.by_ref() {
        let long_hash = hash_one(window_hash);
        phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
        phrase_len += 1;
        if window_hash < hash_bound {
            prefix = (phrase_hash, phrase_len);
            phrase_hash = long_hash;
            phrase_len = w as LT;
            break;
        }
    }
    for window_hash in window_hash_it {
        let long_hash = hash_one(window_hash);
        phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
        phrase_len += 1;
        if window_hash < hash_bound {
            vec_hash.push(phrase_hash);
            vec_len.push(phrase_len);
            phrase_hash = long_hash;
            phrase_len = w as LT;
        }
    }
    let suffix = (phrase_hash, phrase_len);
    (prefix, suffix)
}

/// Splits a range into `num_ranges` evenly-sized ranges overlapping by `overlap`.
///
/// When the length is not a multiple of `num_ranges`, the last range will contain slightly more values than the others.
#[inline(always)]
pub fn overlapping_ranges(
    range: Range<usize>,
    num_ranges: usize,
    overlap: usize,
) -> Vec<Range<usize>> {
    let start = range.start;
    let end = range.end;
    let n = (end - start).saturating_sub(overlap).div_ceil(num_ranges);
    let mut res: Vec<Range<usize>> = (0..num_ranges)
        .map(|i| (i * n)..((i + 1) * n + overlap))
        .collect();
    if let Some(r) = res.last_mut() {
        *r = r.start..end
    }
    res
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters, using multiple threads.
///
/// Note that this will automatically switch to [`parse_seq`] if the sequence is short.
///
/// This is a wrapper for [`parse_seq_par_with_vecs`] that allocates new vectors to store the parse.
#[inline(always)]
pub fn parse_seq_par(seq: &[u8], w: usize, p: usize, threads: usize) -> Parse {
    if threads <= 1 || seq.len() <= SINGLE_THREAD_THRESHOLD {
        return parse_seq(seq, w, p);
    }

    let mut vecs_hash = vec![vec![]; threads];
    let mut vecs_len = vec![vec![]; threads];
    parse_seq_par_with_vecs(seq, w, p, threads, &mut vecs_hash, &mut vecs_len)
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters, using multiple threads and storing the resulting phrases in existing `Vec`.
///
/// `vecs_hash` and `vecs_len` should contain at least `threads` vectors.
///
/// Note that this will automatically switch to [`parse_seq`] if the sequence is short.
pub fn parse_seq_par_with_vecs(
    seq: &[u8],
    w: usize,
    p: usize,
    threads: usize,
    vecs_hash: &mut [Vec<HT>],
    vecs_len: &mut [Vec<LT>],
) -> Parse {
    if threads <= 1 || seq.len() <= SINGLE_THREAD_THRESHOLD {
        return parse_seq(seq, w, p); // TODO avoid allocation?
    }

    let borders: Vec<_> = (
        overlapping_ranges(0..seq.len(), threads, w - 1),
        &mut vecs_hash[..threads],
        &mut vecs_len[..threads],
    )
        .into_par_iter()
        .map(|(range, vec_hash, vec_len)| parse_seq_with_vec(&seq[range], w, p, vec_hash, vec_len))
        .collect();

    let num_phrases = vecs_hash[..threads].iter().map(|v| v.len()).sum::<usize>() + threads - 1;
    let mut phrases = Vec::with_capacity(num_phrases); // TODO avoid allocation?
    let mut phrases_len = Vec::with_capacity(num_phrases);
    phrases.extend(&vecs_hash[0]);
    phrases_len.extend(&vecs_len[0]);
    for i in 1..threads {
        let (suffix_hash, suffix_len) = borders[i - 1].1;
        let (prefix_hash, prefix_len) = borders[i].1;
        #[allow(clippy::unnecessary_cast)]
        let border_hash = merge_hashes(suffix_hash, prefix_hash, prefix_len as u32);
        let border_len = suffix_len + prefix_len - w as LT + 1;
        phrases.push(border_hash);
        phrases_len.push(border_len);
        phrases.extend(&vecs_hash[i]);
        phrases_len.extend(&vecs_len[i]);
    }

    Parse {
        prefix: borders[0].0,
        suffix: borders[threads - 1].1,
        phrases,
        phrases_len,
    }
}

/// An iterator yielding hashes of complete phrases, along with their length.
///
/// When collecting all the phrases, using [`parse_seq`] or [`parse_seq_par`] should be faster.
pub struct ParseIterator<'a> {
    window_hash_it: RollingHashIterator<'a>,
    phrase_hash: HT,
    phrase_len: LT,
    hash_bound: HT,
    w: usize,
    p: usize,
    suffix: &'a mut (HT, LT),
}

impl<'a> ParseIterator<'a> {
    #[inline]
    pub fn new(
        seq: &'a [u8],
        w: usize,
        p: usize,
        prefix: &'a mut (HT, LT),
        suffix: &'a mut (HT, LT),
    ) -> Self {
        let mut window_hash_it = RollingHashIterator::new(seq, w);
        let mut phrase_hash = 0u64;
        let mut phrase_len = 0;
        let hash_bound = HT::MAX / p as HT + 1;
        for window_hash in window_hash_it.by_ref() {
            let long_hash = hash_one(window_hash);
            phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
            phrase_len += 1;
            if window_hash < hash_bound {
                *prefix = (phrase_hash, phrase_len);
                phrase_hash = long_hash;
                phrase_len = w as LT;
                break;
            }
        }
        Self {
            window_hash_it,
            phrase_hash,
            phrase_len,
            hash_bound,
            w,
            p,
            suffix,
        }
    }
}

impl Iterator for ParseIterator<'_> {
    type Item = (HT, LT);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        for window_hash in self.window_hash_it.by_ref() {
            let long_hash = hash_one(window_hash);
            self.phrase_hash = self.phrase_hash.rotate_left(1) ^ long_hash;
            self.phrase_len += 1;
            if window_hash < self.hash_bound {
                let res = (self.phrase_hash, self.phrase_len);
                self.phrase_hash = long_hash;
                self.phrase_len = self.w as LT;
                return Some(res);
            }
        }
        *self.suffix = (self.phrase_hash, self.phrase_len);
        None
    }

    #[inline(always)]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let (lo, hi) = self.window_hash_it.size_hint();
        (lo / self.p * 9 / 10, hi.map(|hi| hi / self.p * 11 / 10))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngCore;
    use rand::SeedableRng;

    #[test]
    fn test_parse_repeated_seq() {
        const LEN: usize = 10000;
        let mut seq = [0u8; LEN];
        rand_xoshiro::Xoshiro512StarStar::from_os_rng().fill_bytes(&mut seq);
        let mut repeated = Vec::from(seq);
        repeated.extend_from_slice(&seq);
        repeated.extend_from_slice(&seq);

        for window_size in 1..=50 {
            let parse = parse_seq(&repeated, window_size, 100);
            let mut idx = 0;
            let mut start = 0;
            while start < LEN {
                start += parse.phrases_len[idx] as usize - window_size;
                idx += 1;
            }
            assert_eq!(parse.phrases[0..idx], parse.phrases[idx..(2 * idx)]);
        }
    }
}
