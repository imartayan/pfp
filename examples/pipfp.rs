use clap::Parser;
use niffler::send::from_path;
use pfp::parse::{Parse, Phrase, parse_seq_par_with};
use rayon::prelude::*;
use rayon::{ThreadPoolBuilder, current_num_threads};
#[cfg(feature = "radix")]
use rdst::RadixSort;
use seq_io::fasta::{self, Record};
use std::path::Path;
use std::time::Instant;

// type Dict = std::collections::HashMap<HT, LT, rustc_hash::FxBuildHasher>;
type Dict = Vec<Phrase>;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed)
    #[arg(short, long)]
    input: String,
    /// Window size
    #[arg(short, default_value_t = 10)]
    w: usize,
    /// Sampling factor of the windows
    #[arg(short, default_value_t = 100)]
    p: usize,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

fn parse_file<P: AsRef<Path>>(w: usize, p: usize, path: P, threads: usize) -> (usize, Dict) {
    let mut num_phrases = 0;
    let mut dict = Dict::default();
    let mut parse = Parse::new(w, p);
    let mut temp_parses = Vec::new();
    let (reader, _) = from_path(path).expect("Failed to open input file");
    let mut reader = fasta::Reader::new(reader);
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let mut seq = Vec::with_capacity(record.seq().len());
        for line in record.seq_lines() {
            seq.extend_from_slice(line);
        }
        parse.clear();
        parse_seq_par_with(&seq, threads, &mut parse, &mut temp_parses);
        num_phrases += parse.phrases.len();
        // dict.extend(parse.phrases.iter().zip(parse.phrases_len.iter()));
        dict.par_extend(parse.par_iter());
    }
    #[cfg(feature = "radix")]
    dict.radix_sort_unstable();
    #[cfg(not(feature = "radix"))]
    dict.par_sort_unstable();
    dict.dedup_by_key(|phrase| phrase.hash());
    (num_phrases, dict)
}

fn main() {
    let args = Args::parse();
    let w = args.w;
    let p = args.p;
    let path = args.input;
    let threads = if let Some(t) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        current_num_threads()
    };
    eprintln!(
        "Start prefix-free parsing using {threads} thread{}",
        if threads > 1 { "s" } else { "" }
    );
    let start_parse = Instant::now();
    let (num_phrases, dict) = parse_file(w, p, path, threads);
    let elapsed = start_parse.elapsed().as_secs_f64();
    eprintln!("Parsed in {:.02} s", elapsed);
    let num_distinct_phrases = dict.len();
    // let len_distinct_phrases = dict.into_values().map(|len| len as usize).sum::<usize>();
    let len_distinct_phrases = dict
        .iter()
        .map(|phrase| phrase.len() as usize)
        .sum::<usize>();
    eprintln!("number of phrases (w/ rep) = {num_phrases}");
    eprintln!("number of distinct phrases = {num_distinct_phrases}");
    eprintln!("length of distinct phrases = {len_distinct_phrases}");
    eprintln!("pi = {}", num_phrases + len_distinct_phrases);
}
