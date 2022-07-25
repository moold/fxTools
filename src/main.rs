use byte_unit::Byte;
use clap::{crate_name, App, AppSettings, Arg};
use hashbrown::HashMap;
use kseq::{parse_path, record::Fastx};
use indoc::indoc;
use owo_colors::OwoColorize;
use regex::Regex;
use std::{
    cmp::max,
    fmt,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

mod utils;
use utils::stat::stat;

pub const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));

fn out_attr(
    record: Fastx,
    attr: &str,
    attr_lower: &str,
) {
    let mut is_first = true;
    let mut s = 0;

    for elm in attr_lower.split(':') {
        if is_first {
            is_first = false
        } else {
            print!("\t");
            s += 1;
        }
        match elm {
            "id" => print!("{}", record.head()),
            "len" => print!("{}", record.len()),
            "qs" => {
                let qual = record.qual();
                if !qual.is_empty() {
                     // \text{read Q} = -10\log_{10}\big[\tfrac{1}{N}\sum 10^{-q_i/10}\big]
                    let e_sum: f64 = qual.as_bytes().iter().map(|x| 10_f64.powf(- 0.1 * (x - 33) as f64)).sum();
                    print!("{}", -10.0 * (e_sum/record.len() as f64).log10());
                }else {
                    print!("NA");
                }
            }
            "" => (),
            _ => {
                let mut letter_counts: HashMap<u8, usize> = HashMap::with_capacity(10);
                record.seq().as_bytes().iter().for_each(|c| *letter_counts.entry(*c).or_insert(0) += 1);
                let x = &attr[s..s + elm.len()];
                print!(
                    "{}",
                    x.bytes()
                        .fold(0, |acc, x| acc + letter_counts.get(&x).unwrap_or(&0))
                );
            }
        }
        s += elm.len();
    }
    if !is_first {
        println!();
    }
}

fn complement_base(b: char) -> char {
    match b {
        'a' => 't',
        'A' => 'T',
        't' => 'a',
        'T' => 'A',
        'g' => 'c',
        'G' => 'C',
        'c' => 'g',
        'C' => 'G',
        _ => b,
    }
}

fn out_seq(seq: &str, qual:&str, reverse: bool, complement: bool) {
    if complement && reverse {
        seq.chars()
            .rev()
            .for_each(|b| print!("{}", complement_base(b)));

        if !qual.is_empty(){
            print!("\n+\n{}", qual.chars().rev().collect::<String>());
        }
    } else if complement {
        seq.chars().for_each(|b| print!("{}", complement_base(b)));
        if !qual.is_empty(){
            print!("\n+\n{}", qual);
        }
    } else if reverse {
        seq.chars().rev().for_each(|b| print!("{}", b));

        if !qual.is_empty(){
            print!("\n+\n{}", qual.chars().rev().collect::<String>());
        }
    } else {
        print!("{}", seq);
        if !qual.is_empty(){
            print!("\n+\n{}", qual);
        }
    }
    println!();
}

fn colour_out_base<T: Into<char>>(b: T, cs: &[i32; 6]) {
    let b = b.into();
    let c: i32 = cs.iter().sum();
    if c <= 1 { // only one base
        print!("{}", b)
    } else { // SNPs or INDELs
        match b {
            'a' | 'A' => print!("{}", b.bold().underline().green()),
            't' | 'T' => print!("{}", b.bold().underline().red()),
            'g' | 'G' => print!("{}", b.bold().underline().yellow()),
            'c' | 'C' => print!("{}", b.bold().underline().blue()),
            '*' => print!("{}", b.bold().underline().cyan()),
            '-' => print!("{}", b.bold().underline().magenta()),
            _ => unreachable!(),
        }
    }
}

fn main() {
    let args = App::new("fastx")
        // .version(crate_version!())
        .version(VERSION)
        .about("A toolset for processing sequences in FASTA/Q formats")
        .override_usage(format!("{} [SUBCOMMAND] [OPTIONS] <input>", crate_name!()).as_str())
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .global_setting(AppSettings::ColoredHelp)
        .global_setting(AppSettings::DisableHelpSubcommand)
        .global_setting(AppSettings::DeriveDisplayOrder)
        // .global_setting(AppSettings::ArgsNegateSubcommands)
        // .global_setting(AppSettings::UnifiedHelpMessage)
        .arg(
            Arg::new("input")
                .about("input file")
                .global(true)
        )
        .arg(
            Arg::new("attr")
                .short('a')
                .long("attr")
                .value_name("STR")
                .about(indoc!{"
                    get sequence attributes, id:len:x
                      len: sequences length
                      x:   count x, case sensitive, where x can be a single base
                           or multiple bases (gcGC means count g + c + G + C)
                      qs:  quality score"
                    })
                .takes_value(true),
        )
        .arg(
            Arg::new("reform")
                .short('r')
                .long("reform")
                .value_name("STR")
                .about(indoc!{"
                    reform/modify the sequence, accept values:
                      lower:    convert sequences into lowercase
                      upper:    convert sequences into uppercase
                      lineINT[c]:  
                                wrap sequences into INT characters per line, 
                                0 for no wrap, c is the alignment mode, 
                                color for SNVs and INDELs
                      linkINT:  link sequences with INT Ns
                      splitINT: split sequences at Ns with length >= INT
                      fq2fa:    converts FASTQ to FASTA
                      fa2fq:    converts FASTA to FASTQ"
                })
                .takes_value(true)
        )
        .arg(
            Arg::new("subsample")
                .short('s')
                .long("sample~")
                .value_name("int[G|M|K]|float")
                .about("subsample reads, int is the total bases, float is the fraction")
                .takes_value(true),
        )
        .subcommand(
            App::new("stat")
                .about("simple statistics of FASTA/Q files")
                .arg(
                    Arg::new("min_len")
                        .short('m')
                        .long("min_len")
                        .value_name("int[G|M|K]")
                        .default_value("0")
                        .about("minimum sequence length, shorter sequences are ignored")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("n_len")
                        .short('n')
                        .long("n_len")
                        .value_name("int[G|M|K]")
                        .default_value("0")
                        .about("minimum gap (N) length, 0 means do not stat contig length")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("genome_len")
                        .short('g')
                        .long("genome_len")
                        .value_name("int[G|M|K]")
                        .default_value("0")
                        // .requires("n_len")
                        .about("genome size, 0 means calculate genome size using the input file, co-used with -n")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("out_ctg")
                        .short('o')
                        .long("out_ctg")
                        .about("output contig sequences to <INPUT>.ctg.fa")
                )
        )
        .subcommand(
            App::new("findseq")
                .about("find subseq positions")
                .arg(
                    Arg::new("subseq")
                        .short('s')
                        .long("subseq")
                        .value_name("STR")
                        .about("sequence to be find")
                        .required(true)
                        .takes_value(true)
                )
                .arg(
                    Arg::new("ignore-case")
                        .short('i')
                        .long("ignore-case")
                        .about("ignore case")
                )
        )
        .subcommand(
            App::new("findgap")
            .about("find gap(Nn) regions")
            .arg(
                Arg::new("min_len")
                    .short('m')
                    .long("min_len")
                    .value_name("int[G|M|K]")
                    .default_value("0")
                    .about("minimum gap length, shorter gaps are ignored")
                    .takes_value(true),
            )
        )
        .subcommand(
            App::new("getseq")
                .about("get sequences or subsequences from a region or file")
                .arg(
                    Arg::new("region")
                        .short('r')
                        .long("region")
                        .value_name("STR|FILE")
                        .about("region (chr, chr:start-end, chr-start-end or chr1,chr2:start-end) or file (bed or ID list) \
                            to be extracted, format: 0-based, [start, end).")
                        .required(true)
                        .takes_value(true)
                )
                .arg(
                    Arg::new("reverse")
                        .short('v')
                        .long("reverse")
                        .about("reverse the sequence")
                )
                .arg(
                    Arg::new("complement")
                        .short('c')
                        .long("complement")
                        .about("complement the sequence")
                )
        )
        .subcommand(
            App::new("diff~")
                .about("compare sequences between files")
                .arg(
                    Arg::new("base")
                        .short('b')
                        .long("base")
                        .about("compare sequences with same ID base by base")
                )
        )
        .get_matches();

    let path: Option<String> = args.value_of("input").map(str::to_string);
    let mut records = parse_path(path.to_owned()).unwrap();
    if let Some(attr) = args.value_of("attr") {
        let attr = attr.trim_matches(':');
        let attr_lower = attr.to_ascii_lowercase();
        while let Some(record) = records.iter_record().unwrap() {
            out_attr(record, attr, &attr_lower);
        }
    } else if let Some(reform) = args.value_of("reform") {
        let reform = reform.to_ascii_lowercase();
        if reform == "lower" {
            while let Some(record) = records.iter_record().unwrap() {
                if record.sep().is_empty() {
                    println!(
                        ">{} {}\n{}",
                        record.head(),
                        record.des(),
                        record.seq().to_ascii_lowercase()
                    );
                } else {
                    println!(
                        "@{} {}\n{}\n{}\n{}",
                        record.head(),
                        record.des(),
                        record.seq().to_ascii_lowercase(),
                        record.sep(),
                        record.qual()
                    );
                }
            }
        } else if reform == "upper" {
            while let Some(record) = records.iter_record().unwrap() {
                if record.sep().is_empty() {
                    println!(
                        ">{} {}\n{}",
                        record.head(),
                        record.des(),
                        record.seq().to_ascii_uppercase()
                    );
                } else {
                    println!(
                        "@{} {}\n{}\n{}\n{}",
                        record.head(),
                        record.des(),
                        record.seq().to_ascii_uppercase(),
                        record.sep(),
                        record.qual()
                    );
                }
            }
        } else if reform == "fq2fa" {
            while let Some(record) = records.iter_record().unwrap() {
                println!(">{} {}\n{}", record.head(), record.des(), record.seq());
            }
        } else if reform == "fa2fq" {
            while let Some(record) = records.iter_record().unwrap() {
                println!(
                    "@{} {}\n{}\n+\n{:!<4$}",
                    record.head(),
                    record.des(),
                    record.seq(),
                    "",
                    record.len()
                );
            }
        } else if reform.starts_with("line") {
            let is_colour = reform.ends_with('c');
            let w: usize = reform
                .strip_suffix('c')
                .unwrap_or(&reform)
                .strip_prefix("line")
                .unwrap()
                .parse()
                .unwrap();
            if is_colour {
                let mut bases = Vec::new();
                while let Some(record) = records.iter_record().unwrap() {
                    if record.len() > bases.len() {
                        bases.resize(record.len(), [0; 6]);
                    }
                    for (p, c) in record.seq().char_indices() {
                        let cs = &mut bases[p];
                        match c {
                            'a' | 'A' => cs[0] = 1,
                            't' | 'T' => cs[1] = 1,
                            'g' | 'G' => cs[2] = 1,
                            'c' | 'C' => cs[3] = 1,
                            '*' => cs[4] = 1,
                            '-' => cs[5] = 1,
                            _ => unreachable!(),
                        }
                    }
                }

                let path: Option<String> = args.value_of("input").map(str::to_string);
                let mut records = parse_path(path).unwrap();
                while let Some(record) = records.iter_record().unwrap() {
                    println!(">{} {}", record.head(), record.des());
                    let (seq, len) = (record.seq().as_bytes(), record.len());
                    if w == 0 {
                        for i in 0..len {
                            colour_out_base(seq[i], &bases[i]);
                        }
                        println!();
                    } else {
                        for i in (0..len).step_by(w) {
                            if i + w < len {
                                for j in i..i + w {
                                    colour_out_base(seq[j], &bases[j]);
                                }
                            } else {
                                for j in i..len {
                                    colour_out_base(seq[j], &bases[j]);
                                }
                            }
                            println!();
                        }
                    }
                }
            } else {
                while let Some(record) = records.iter_record().unwrap() {
                    if record.sep().is_empty() {
                        println!(">{} {}", record.head(), record.des());
                        if w == 0 {
                            println!("{}", record.seq());
                        } else {
                            let (seq, len) = (record.seq(), record.len());
                            for i in (0..len).step_by(w) {
                                if i + w < len {
                                    println!("{}", &seq[i..i + w]);
                                } else {
                                    println!("{}", &seq[i..len]);
                                }
                            }
                        }
                    } else {
                        println!("@{} {}", record.head(), record.des());
                        if w == 0 {
                            println!("{}\n{}\n{}", record.seq(), record.sep(), record.qual());
                        } else {
                            let (seq, qual, len) = (record.seq(), record.qual(), record.len());
                            for i in (0..len).step_by(w) {
                                if i + w < len {
                                    println!("{}", &seq[i..i + w]);
                                } else {
                                    println!("{}", &seq[i..len]);
                                }
                            }
                            println!("{}", record.sep());
                            for i in (0..len).step_by(w) {
                                if i + w < len {
                                    println!("{}", &qual[i..i + w]);
                                } else {
                                    println!("{}", &qual[i..len]);
                                }
                            }
                        }
                    }
                }
            }
        } else if reform.starts_with("link") {
            let w: usize = reform.strip_prefix("link").unwrap().parse().unwrap();
            println!(">link_reads");
            let mut is_first = true;
            while let Some(record) = records.iter_record().unwrap() {
                if is_first {
                    is_first = false;
                    print!("{}", record.seq());
                } else {
                    print!("{:N<2$}{}", "", record.seq(), w);
                }
            }
            println!();
        } else if reform.starts_with("split") {
            let w: usize = reform.strip_prefix("split").unwrap().parse().unwrap();
            if w == 0 {
                return;
            };
            let re = Regex::new(&format!("(?i)N{{{w},}}", w = w)).unwrap();
            while let Some(record) = records.iter_record().unwrap() {
                let mut last_pos = 0;
                let len = record.len();
                let seq = record.seq();
                for mat in re.find_iter(seq) {
                    println!(
                        ">{}:{}_{}\n{}",
                        record.head(),
                        last_pos,
                        mat.start() - 1,
                        &seq[last_pos..mat.start()]
                    );
                    last_pos = mat.end();
                }
                if len > last_pos {
                    if last_pos == 0 {
                        println!(">{}\n{}", record.head(), &seq[last_pos..len]);
                    } else {
                        println!(
                            ">{}:{}_{}\n{}",
                            record.head(),
                            last_pos,
                            len - 1,
                            &seq[last_pos..len]
                        );
                    }
                }
            }
        } else {
            panic!("unknown values: {} for --reform", reform);
        }
    } else if let Some(subarg) = args.subcommand_matches("findseq") {
        let subseq = subarg.value_of("subseq").unwrap();
        let re = if subarg.is_present("ignore-case") {
            Regex::new(format!("(?i){}", subseq).as_str()).unwrap()
        } else {
            Regex::new(&subseq.to_string()).unwrap()
        };
        while let Some(record) = records.iter_record().unwrap() {
            for mat in re.find_iter(record.seq()) {
                println!("{}\t{}\t{}", record.head(), mat.start(), mat.end() - 1);
            }
        }
    } else if let Some(subarg) = args.subcommand_matches("findgap") {
        let min_len = Byte::from_str(subarg.value_of("min_len").unwrap())
            .unwrap()
            .get_bytes();
        let re = Regex::new(r"(?i)N+").unwrap();
        while let Some(record) = records.iter_record().unwrap() {
            for mat in re.find_iter(record.seq()) {
                if mat.end() - mat.start() >= min_len as usize {
                    println!("{}\t{}\t{}", record.head(), mat.start(), mat.end() - 1);
                }
            }
        }
    } else if let Some(subarg) = args.subcommand_matches("getseq") {
        let region = subarg.value_of("region").unwrap();
        let path = Path::new(region);
        let mut out_info: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

        if path.is_file() {
            let file = File::open(path).unwrap();
            for line in BufReader::new(file).lines().flatten() {
                if !line.starts_with('#') {
                    let v: Vec<&str> = line.split(char::is_whitespace).collect();
                    let out_info = out_info.entry(v[0].to_owned()).or_insert(vec![]);
                    if v.len() == 1 {
                        out_info.push((0, 0))
                    } else {
                        out_info.push((
                            v[1].parse::<u32>().unwrap_or(0),
                            v[2].parse::<u32>().unwrap_or(0),
                        ))
                    }
                }
            }
        } else {
            for region in region.split(&[',', ';'][..]) {
                let re = Regex::new(r"(?i)(\S+)[:-](start|\d+)[:-](end|\d+)$").unwrap();
                if let Some(caps) = re.captures(region) {
                    let chr = caps.get(1).unwrap().as_str();
                    let start = caps.get(2).unwrap().as_str().parse::<u32>().unwrap_or(0);
                    let end = caps.get(3).unwrap().as_str().parse::<u32>().unwrap_or(0);
                    out_info.entry(chr.to_owned()).or_insert(vec![]).push((start, end));
                } else {
                    out_info.entry(region.to_owned()).or_insert(vec![]).push((0, 0));
                }
            }
        }
        while let Some(record) = records.iter_record().unwrap() {
            let head = record.head();
            if let Some(regions) = out_info.get_mut(head) {
                let seqs = record.seq();
                let qual = record.qual();
                for (start, end) in regions.iter_mut() {
                    let mut sub_head: String;
                    if start == end {
                        *start = 0;
                        *end = record.len() as u32;
                        sub_head = head.to_owned();
                    } else {
                        if *end == 0 {
                            *end = record.len() as u32;
                        }
                        sub_head = format!("{}:{}-{}", head, start, *end - 1);
                    }
                    if subarg.is_present("reverse") {
                        sub_head += "_rev";
                    }
                    if subarg.is_present("complement") {
                        sub_head += "_com";
                    }
                    if qual.is_empty(){
                        println!(">{}", sub_head);
                    }else{
                        println!("@{}", sub_head);
                    }
                    out_seq(
                        &seqs[*start as usize..*end as usize],
                        if qual.is_empty(){
                            qual
                        }else{
                            &qual[*start as usize..*end as usize]
                        },
                        subarg.is_present("reverse"),
                        subarg.is_present("complement"),
                    );
                }
            }
            out_info.remove(head);
            if out_info.is_empty() {
                break;
            }
        }
        for (key, _) in out_info.iter() {
            eprintln!("Missing record: {}", key);
        }
    } else if let Some(subarg) = args.subcommand_matches("stat") {
        let min_len = Byte::from_str(subarg.value_of("min_len").unwrap())
            .unwrap()
            .get_bytes() as usize;
        let genome_len = Byte::from_str(subarg.value_of("genome_len").unwrap())
            .unwrap()
            .get_bytes() as usize;
        let n_len = Byte::from_str(subarg.value_of("n_len").unwrap())
            .unwrap()
            .get_bytes() as usize;
        stat(&[&path.unwrap()], min_len, genome_len, n_len, subarg.is_present("out_ctg"));
    }
}
