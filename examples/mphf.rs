use clap::Parser;
use epserde::prelude::*;
use niffler::send::from_path;
use pfp::parse::{Parse, parse_seq_par_with};
use rayon::{ThreadPoolBuilder, current_num_threads};
use seq_io::fasta::{self, Record};
use std::path::Path;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed)
    #[arg(short, long)]
    input: String,
    /// Output file to serialize the index
    #[arg(short, long)]
    output: Option<String>,
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

fn parse_file<P: AsRef<Path>>(w: usize, p: usize, path: P, threads: usize) -> Parse {
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
        parse_seq_par_with(&seq, threads, &mut parse, &mut temp_parses);
    }
    parse
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
    let parse = parse_file(w, p, path, threads);
    let elapsed = start_parse.elapsed().as_secs_f64();
    eprintln!("Parsed in {:.02} s", elapsed);
    eprintln!("Building MPHF using PtrHash");
    let start_mphf = Instant::now();
    let compact = parse.build_compact();
    let elapsed = start_mphf.elapsed().as_secs_f64();
    eprintln!("Built MPHF in {:.02} s", elapsed);
    if let Some(path) = args.output {
        let mut file = std::fs::File::create(path).unwrap();
        let start_ser = Instant::now();
        compact.serialize(&mut file).unwrap();
        let elapsed = start_ser.elapsed().as_secs_f64();
        eprintln!("Serialized MPHF in {:.02} s", elapsed);
    }
}
