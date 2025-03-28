//! Prefix-free parsing.

use crate::hash::{HT, RollingHashIterator, hash_one, merge_hashes};
use core::ops::Range;
use rayon::prelude::*;

#[cfg(feature = "mphf")]
use cacheline_ef::CachelineEfVec;
#[cfg(feature = "epserde")]
use epserde::prelude::*;
#[cfg(feature = "mphf")]
use ptr_hash::{PtrHash, PtrHashParams, bucket_fn::Linear, hash::FxHash};
#[cfg(feature = "radix")]
use rdst::RadixSort;

/// Integer type for phrase length.
pub type LT = u32;

#[cfg(feature = "mphf")]
pub type MPHF = PtrHash<u64, Linear, CachelineEfVec, FxHash>;

const SINGLE_THREAD_THRESHOLD: usize = 1 << 16;

#[cfg(feature = "mphf")]
const LEN_BITS: usize = 20;
#[cfg(feature = "mphf")]
const LEN_MASK: u64 = (1 << LEN_BITS) - 1;

/// Partial phrase, represented by its hash and length.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "epserde", derive(Epserde))]
#[cfg_attr(feature = "epserde", repr(C))]
#[cfg_attr(feature = "epserde", zero_copy)]
pub struct Phrase {
    hash: HT,
    len: LT,
}

impl Phrase {
    /// Creates a new [`Phrase`] with the given hash and length.
    #[inline(always)]
    pub const fn new(hash: HT, len: LT) -> Self {
        Self { hash, len }
    }

    /// Hash of the phrase.
    #[inline(always)]
    pub const fn hash(&self) -> HT {
        self.hash
    }

    /// Length of the phrase.
    #[inline(always)]
    pub const fn len(&self) -> LT {
        self.len
    }

    /// Returns true if the phrase is empty.
    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Merges two partial phrases overlapping by `overlap` characters.
    #[inline(always)]
    #[allow(clippy::unnecessary_cast)]
    pub const fn merge(&self, other: &Self, overlap: usize) -> Self {
        Self {
            hash: merge_hashes(self.hash, other.hash, other.len as u32 - overlap as u32),
            len: self.len + other.len - overlap as LT,
        }
    }
}

#[cfg(feature = "radix")]
impl rdst::RadixKey for Phrase {
    const LEVELS: usize = HT::BITS as usize / 8;

    #[inline]
    fn get_level(&self, level: usize) -> u8 {
        (self.hash() >> (level * 8)) as u8
    }
}

/// Prefix-free parse.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "epserde", derive(Epserde))]
pub struct Parse {
    /// Window size.
    pub w: usize,
    /// Sampling factor of the windows.
    pub p: usize,
    /// Incomplete starting phrase (with its hash and length).
    pub prefix: Phrase,
    /// Incomplete ending phrase (with its hash and length).
    pub suffix: Phrase,
    /// Hashes of the complete phrases, in their order of appearance.
    pub phrases: Vec<HT>,
    /// Lengths of the complete phrases, in their order of appearance.
    pub phrases_len: Vec<LT>,
}

impl Parse {
    /// Creates an empty [`Parse`] with parameters `w` and `p`.
    #[inline(always)]
    pub fn new(w: usize, p: usize) -> Self {
        Self {
            w,
            p,
            prefix: Phrase::default(),
            suffix: Phrase::default(),
            phrases: Vec::new(),
            phrases_len: Vec::new(),
        }
    }

    /// Number of complete phrases in the parse.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.phrases.len()
    }

    /// Returns true if the parse has no complete phrases.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.phrases.is_empty()
    }

    /// Clears the phrases without reducing the capacity.
    #[inline(always)]
    pub fn clear(&mut self) {
        self.prefix = Phrase::default();
        self.suffix = Phrase::default();
        self.phrases.clear();
        self.phrases_len.clear();
    }

    /// Reserves capacity for at least `additional` more phrases.
    #[inline(always)]
    pub fn reserve(&mut self, additional: usize) {
        self.phrases.reserve(additional);
        self.phrases_len.reserve(additional);
    }

    /// Extends with an existing [`Parse`], merging the current suffix with the next prefix in a new phrase.
    #[inline(always)]
    pub fn extend(&mut self, other: &Self) {
        if self.is_empty() {
            self.prefix = other.prefix;
            self.suffix = other.suffix;
            self.phrases.extend_from_slice(&other.phrases);
            self.phrases_len.extend_from_slice(&other.phrases_len);
        } else {
            let border = self.suffix.merge(&other.prefix, self.w - 1);
            self.suffix = if other.is_empty() {
                border
            } else {
                self.phrases.push(border.hash);
                self.phrases_len.push(border.len);
                self.phrases.extend_from_slice(&other.phrases);
                self.phrases_len.extend_from_slice(&other.phrases_len);
                other.suffix
            }
        }
    }

    /// Extends with the [`Parse`] computed from `seq`, merging the current suffix with the next prefix in a new phrase.
    pub fn extend_from_seq(&mut self, seq: &[u8]) {
        let num_phrases_upper_bound = seq.len() / self.p * 11 / 10;
        let hash_bound = HT::MAX / self.p as HT + 1;
        self.reserve(num_phrases_upper_bound);
        let mut window_hash_it = RollingHashIterator::new(seq, self.w);
        let mut phrase_hash = 0u64;
        // we start at w - 1 because the very first iteration of the loop
        // will always increment `phrase_len` even before we check if the
        // hash is a trigger, and after that increment, we want the current
        // phrase length to be `w`.
        let mut phrase_len = self.w as LT - 1;

        for window_hash in window_hash_it.by_ref() {
            let long_hash = hash_one(window_hash);
            phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
            phrase_len += 1;
            if window_hash < hash_bound {
                let prefix = Phrase::new(phrase_hash, phrase_len);
                if self.is_empty() {
                    self.prefix = prefix;
                } else {
                    let border = self.suffix.merge(&prefix, self.w - 1);
                    self.phrases.push(border.hash);
                    self.phrases_len.push(border.len);
                }
                phrase_hash = long_hash;
                phrase_len = self.w as LT;
                break;
            }
        }
        for window_hash in window_hash_it {
            let long_hash = hash_one(window_hash);
            phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
            phrase_len += 1;
            if window_hash < hash_bound {
                self.phrases.push(phrase_hash);
                self.phrases_len.push(phrase_len);
                phrase_hash = long_hash;
                phrase_len = self.w as LT;
            }
        }
        self.suffix = Phrase::new(phrase_hash, phrase_len);
    }

    /// Iterator over the complete phrases of the parse represented as [`Phrase`].
    #[inline(always)]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Phrase> {
        self.phrases
            .iter()
            .zip(self.phrases_len.iter())
            .map(|(&hash, &len)| Phrase::new(hash, len))
    }

    /// Parallel iterator over the complete phrases of the parse represented as [`Phrase`].
    #[inline(always)]
    pub fn par_iter(&self) -> impl IntoParallelIterator<Item = Phrase> {
        self.phrases
            .par_iter()
            .zip(self.phrases_len.par_iter())
            .map(|(&hash, &len)| Phrase::new(hash, len))
    }

    #[inline(always)]
    fn _unique_hashes(&self) -> Vec<HT> {
        let mut hashes = self.phrases.clone();
        #[cfg(feature = "radix")]
        hashes.radix_sort_unstable();
        #[cfg(not(feature = "radix"))]
        hashes.par_sort_unstable();
        hashes.dedup();
        hashes
    }

    #[cfg(feature = "mphf")]
    pub fn build_compact(&self) -> CompactParse {
        let unique_hashes = self._unique_hashes();
        let mphf = PtrHash::new(&unique_hashes, PtrHashParams::default_fast());
        let mut phrases = Vec::with_capacity(self.len());
        let mut hashes = vec![HT::MAX; unique_hashes.len()];
        let mut locations = vec![u64::MAX; unique_hashes.len()];
        let mut iter = self.iter();
        let idx_stream = mphf.index_stream::<32, true, _>(self.phrases.iter());
        idx_stream.fold(0, |start, idx| {
            let phrase = iter.next().unwrap();
            phrases.push(idx as u32);
            hashes[idx] = phrase.hash();
            debug_assert!(phrase.len() as u64 <= LEN_MASK);
            locations[idx] = (start << LEN_BITS) | phrase.len() as u64;
            start + phrase.len() as u64
        });
        CompactParse {
            prefix: self.prefix,
            suffix: self.suffix,
            mphf,
            phrases,
            hashes,
            locations,
        }
    }
}

#[cfg_attr(feature = "mphf", derive(Epserde))]
#[cfg(feature = "mphf")]
pub struct CompactParse {
    pub prefix: Phrase,
    pub suffix: Phrase,
    pub mphf: MPHF,
    pub phrases: Vec<u32>,
    pub hashes: Vec<HT>,
    pub locations: Vec<u64>,
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters.
///
/// This is a wrapper for [`Parse::extend_from_seq`] that creates a new [`Parse`].
#[inline(always)]
pub fn parse_seq(seq: &[u8], w: usize, p: usize) -> Parse {
    let mut parse = Parse::new(w, p);
    parse.extend_from_seq(seq);
    parse
}

/// Splits a range into `num_ranges` evenly-sized ranges overlapping by `overlap`.
///
/// When the length is not a multiple of `num_ranges`, the last range will contain slightly less values than the others.
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
/// This is a wrapper for [`parse_seq_par_with`] that creates a temporary [`Parse`] for each thread.
#[inline(always)]
pub fn parse_seq_par(seq: &[u8], w: usize, p: usize, threads: usize) -> Parse {
    let mut parse = Parse::new(w, p);
    let mut temp_parses = Vec::new();
    parse_seq_par_with(seq, threads, &mut parse, &mut temp_parses);
    parse
}

/// Computes the parse a sequence with windows of size `w` and a fraction 1/`p` of windows acting as delimiters, using multiple threads, storing temporary parses in `temp_parses` and extending `parse` with the final parse.
///
/// Note that this will automatically switch to [`parse_seq`] if the sequence is short.
pub fn parse_seq_par_with(
    seq: &[u8],
    threads: usize,
    parse: &mut Parse,
    temp_parses: &mut Vec<Parse>,
) {
    if threads <= 1 || seq.len() <= SINGLE_THREAD_THRESHOLD {
        parse.extend_from_seq(seq);
    }

    if temp_parses.len() < threads {
        temp_parses.resize_with(threads, || Parse::new(parse.w, parse.p));
    }

    (
        overlapping_ranges(0..seq.len(), threads, parse.w - 1),
        &mut temp_parses[..threads],
    )
        .into_par_iter()
        .for_each(|(range, temp_parse)| {
            temp_parse.clear();
            temp_parse.extend_from_seq(&seq[range]);
        });

    let num_phrases_upper_bound = temp_parses[..threads]
        .iter()
        .map(|temp_parse| temp_parse.len())
        .sum::<usize>()
        + threads
        - 1;
    parse.reserve(num_phrases_upper_bound);
    for temp_parse in temp_parses.iter() {
        parse.extend(temp_parse);
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
    suffix: &'a mut Phrase,
}

impl<'a> ParseIterator<'a> {
    /// Creates a new [`ParseIterator`], yielding complete phrases as [`Phrase`].
    /// The incomplete starting phrase is written to `prefix` at the creation of the iterator.
    /// The incomplete ending phrase is written to `suffix` when the iterator is **exhausted**.
    #[inline]
    pub fn new(
        seq: &'a [u8],
        w: usize,
        p: usize,
        prefix: &'a mut Phrase,
        suffix: &'a mut Phrase,
    ) -> Self {
        let mut window_hash_it = RollingHashIterator::new(seq, w);
        let mut phrase_hash = 0u64;
        let mut phrase_len = w as LT - 1;
        let hash_bound = HT::MAX / p as HT + 1;
        for window_hash in window_hash_it.by_ref() {
            let long_hash = hash_one(window_hash);
            phrase_hash = phrase_hash.rotate_left(1) ^ long_hash;
            phrase_len += 1;
            if window_hash < hash_bound {
                *prefix = Phrase::new(phrase_hash, phrase_len);
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
    type Item = Phrase;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        for window_hash in self.window_hash_it.by_ref() {
            let long_hash = hash_one(window_hash);
            self.phrase_hash = self.phrase_hash.rotate_left(1) ^ long_hash;
            self.phrase_len += 1;
            if window_hash < self.hash_bound {
                let res = Phrase::new(self.phrase_hash, self.phrase_len);
                self.phrase_hash = long_hash;
                self.phrase_len = self.w as LT;
                return Some(res);
            }
        }
        *self.suffix = Phrase::new(self.phrase_hash, self.phrase_len);
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

    #[test]
    fn test_parse_repeated_seq() {
        const LEN: usize = 10000;
        let mut seq = [0u8; LEN];
        fastrand::fill(&mut seq);
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

    #[test]
    fn test_parse_par() {
        const LEN: usize = 1_000_000;
        let mut seq = vec![0u8; LEN];
        fastrand::fill(&mut seq);

        for window_size in 1..=50 {
            let parse = parse_seq(&seq, window_size, 100);
            let parse_par = parse_seq_par(&seq, window_size, 100, 4);
            assert_eq!(parse.prefix, parse_par.prefix);
            assert_eq!(parse.suffix, parse_par.suffix);
            assert_eq!(parse.phrases, parse_par.phrases);
            assert_eq!(parse.phrases_len, parse_par.phrases_len);
        }
    }

    #[test]
    fn test_parse_iter() {
        const LEN: usize = 1_000_000;
        let mut seq = vec![0u8; LEN];
        fastrand::fill(&mut seq);

        for window_size in 1..=50 {
            let parse = parse_seq(&seq, window_size, 100);
            let prefix = &mut Phrase::default();
            let suffix = &mut Phrase::default();
            let mut parse_iter = ParseIterator::new(&seq, window_size, 100, prefix, suffix);
            for phrase in parse.iter() {
                assert_eq!(parse_iter.next(), Some(phrase));
            }
            assert_eq!(parse_iter.next(), None);
            assert_eq!(parse.prefix, *prefix);
            assert_eq!(parse.suffix, *suffix);
        }
    }
}
