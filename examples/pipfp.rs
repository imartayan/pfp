use clap::Parser;
use niffler::send::from_path;
use pfp::hash::HT;
use pfp::parse::{parse_seq, parse_seq_par, LT};
use rayon::{current_num_threads, ThreadPoolBuilder};
use seq_io::fasta::{self, Record};
use std::path::Path;
use std::time::Instant;

type Dict = std::collections::HashMap<HT, LT, rustc_hash::FxBuildHasher>;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed)
    #[arg(short, long)]
    input: String,
    /// window size
    #[arg(short, default_value_t = 10)]
    w: usize,
    /// sparsity?
    #[arg(short, default_value_t = 100)]
    p: usize,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

fn parse_file<P: AsRef<Path>>(w: usize, p: usize, path: P, threads: usize) -> Dict {
    let mut dict = Dict::default();
    let (reader, _) = from_path(path).expect("Failed to open input file");
    let mut reader = fasta::Reader::new(reader);
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let mut seq = Vec::with_capacity(record.seq().len());
        for line in record.seq_lines() {
            seq.extend_from_slice(line);
        }
        let parse = if threads > 1 {
            parse_seq_par(&seq, w, p, threads)
        } else {
            parse_seq(&seq, w, p)
        };
        dict.extend(parse.phrases.iter().copied().zip(parse.phrases_len));
    }
    dict
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
    let dict = parse_file(w, p, path, threads);
    let elapsed = start_parse.elapsed().as_secs_f64();
    eprintln!("Parsed in {:.02} s", elapsed);
    let num_phrases = dict.len();
    let len_phrases = dict.into_values().map(|len| len as usize).sum::<usize>();
    eprintln!("number of distinct phrases = {num_phrases}");
    eprintln!("length of distinct phrases = {len_phrases}");
    eprintln!("pi = {}", num_phrases + len_phrases);
}
