use byte_unit::Byte;
use clap::{crate_name, App, AppSettings, Arg};
use hashbrown::HashMap;
use kseq::parse_path;
// use lazy_static::lazy_static;
use indoc::indoc;
use owo_colors::OwoColorize;
use regex::Regex;
use std::path::Path;
use std::{
    cmp::max,
    fmt,
    fs::File,
    io::{BufRead, BufReader, Write},
};

pub const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));

#[derive(Default)]
struct NX {
    count: [usize; 10],
    len: [u32; 10],
}

impl NX {
    fn get_width(&self) -> (usize, usize) {
        let mut width1 = self.count[8].to_string().len();
        let mut width2 = self.len[0].to_string().len();
        if width1 < 9 {
            width1 = 9;
        }
        if width2 < 11 {
            width2 = 11;
        }
        (width1, width2)
    }

    fn fill_data(&mut self, lens: &[u32], total: usize) {
        let mut i = 0;
        let mut acc = 0;
        for pos in (0..lens.len()).rev() {
            acc += lens[pos] as usize;
            self.count[i] += 1;
            while acc as f64 > ((i + 1) * total) as f64 * 0.1 {
                self.len[i] = lens[pos];
                i += 1;
                if i >= 10 {
                    break;
                }
                self.count[i] += self.count[i - 1];
            }
        }
    }
}

impl fmt::Display for NX {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (width1, width2) = self.get_width();
        writeln!(
            f,
            "{:<5} {:^width1$} {:^width2$}",
            "Types",
            "Count (#)",
            "Length (bp)",
            width1 = width1,
            width2 = width2
        )?;
        for i in 0..9 {
            writeln!(
                f,
                "N{:<4} {:^width1$} {:^width2$}",
                (i + 1) * 10,
                self.count[i],
                self.len[i],
                width1 = width1,
                width2 = width2
            )?;
        }
        Ok(())
    }
}

fn out_stat(lens: &[u32], total: usize) {
    let total_count = lens.len();
    let bin = max(total_count / 200, 10);
    let bin_step = max(10, lens[total_count / 2] / 20) as usize;
    let mut bin_index: usize = max(lens[0] as usize / bin_step, 1);
    let mut bin_count: usize = 0;

    let mut stat: NX = Default::default();
    let mut stat_i = 0;
    let mut accumulate: usize = 0;

    println!("[length histogram ('*' =~ {:?} reads)]", bin);
    for (pos, len) in lens.iter().enumerate() {
        if (*len as usize) < bin_step * bin_index || total_count - 1 - pos < bin << 1 {
            bin_count += 1;
        } else {
            loop {
                if bin_count > 0 {
                    println!(
                        "{:>width1$} {:>width1$} {:>width2$} {:*>width3$}",
                        bin_step * (bin_index - 1),
                        bin_step * bin_index - 1,
                        bin_count,
                        "",
                        width1 = lens.last().unwrap_or(&0).to_string().len(),
                        width2 = total_count.to_string().len(),
                        width3 = bin_count / bin
                    );
                }
                bin_count = 0;
                bin_index += 1;
                if (*len as usize) < bin_step * bin_index {
                    //check <= or <
                    bin_count = 1;
                    break;
                }
            }
        }
        accumulate += lens[total_count - 1 - pos] as usize;
        stat.count[stat_i] += 1;
        while accumulate as f64 > ((stat_i + 1) * total) as f64 * 0.1 {
            stat.len[stat_i] = lens[total_count - 1 - pos];
            stat_i += 1;
            stat.count[stat_i] += stat.count[stat_i - 1];
        }
    }

    println!(
        "{:>width1$} {:>width1$} {:>width2$} {:*>width3$}",
        bin_step * (bin_index - 1),
        bin_step * bin_index - 1,
        bin_count,
        "",
        width1 = lens.last().unwrap_or(&0).to_string().len(),
        width2 = total_count.to_string().len(),
        width3 = bin_count / bin
    );

    println!("\n\n[length stat]\n{}", stat);
    let (stat_width1, stat_width2) = stat.get_width();
    println!(
        "{:<5} {:^width1$} {:^width2$}",
        "Min.",
        "-",
        lens[0],
        width1 = stat_width1,
        width2 = stat_width2
    );
    println!(
        "{:<5} {:^width1$} {:^width2$}",
        "Max.",
        "-",
        lens[total_count - 1],
        width1 = stat_width1,
        width2 = stat_width2
    );
    println!(
        "{:<5} {:^width1$} {:^width2$}",
        "Ave.",
        "-",
        total / total_count,
        width1 = stat_width1,
        width2 = stat_width2
    );
    println!(
        "{:<5} {:^width1$} {:^width2$}",
        "Total",
        total_count,
        total,
        width1 = stat_width1,
        width2 = stat_width2
    );
}

fn count_by_min(lens: &[u32], min: u32) -> (usize, usize) {
    let mut count = 0;
    let mut total = 0;
    for len in lens {
        if *len >= min {
            count += 1;
            total += *len;
        }
    }
    (total as usize, count)
}

fn out_stat_with_ctg(
    lens: &[u32],
    total: usize,
    ctg_lens: &[u32],
    ctg_total: usize,
    gap_lens: &[u32],
    gap_total: usize,
    genome_len: usize,
) {
    let mut stat: NX = Default::default();
    let mut ctg_stat: NX = Default::default();
    let mut gap_stat: NX = Default::default();
    stat.fill_data(lens, if genome_len > 0 { genome_len } else { total });
    ctg_stat.fill_data(
        ctg_lens,
        if genome_len > 0 {
            genome_len
        } else {
            ctg_total
        },
    );
    gap_stat.fill_data(gap_lens, gap_total);

    println!("{:=<7}{:=^26}{:=^26}{:=^25}", "", "", "", "");
    println!(
        "{:<7}{:^26}{:^26}{:^25}",
        "Types", "Scaffold", "Contig", "Gap"
    );
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "", "Length (bp)", "Count (#)", "Length (bp)", "Count (#)", "Length (bp)", "Count (#)"
    );
    println!("{:-<7}{:-^26}{:-^26}{:-^25}", "", "", "", "");
    for i in 0..9 {
        println!(
            "N{:<6}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
            (i + 1) * 10,
            stat.len[i],
            stat.count[i],
            ctg_stat.len[i],
            ctg_stat.count[i],
            gap_stat.len[i],
            gap_stat.count[i],
        );
    }
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "Longest",
        lens.last().unwrap_or(&0),
        1,
        ctg_lens.last().unwrap_or(&0),
        1,
        gap_lens.last().unwrap_or(&0),
        if gap_lens.is_empty() { 0 } else { 1 },
    );
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "Total",
        total,
        lens.len(),
        ctg_total,
        ctg_lens.len(),
        gap_total,
        gap_lens.len(),
    );
    let (len_total, len_count) = count_by_min(lens, 10000);
    let (ctg_total, ctg_count) = count_by_min(ctg_lens, 10000);
    let (gap_total, gap_count) = count_by_min(gap_lens, 10000);
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        ">=10kb", len_total, len_count, ctg_total, ctg_count, gap_total, gap_count
    );
    let (len_total, len_count) = count_by_min(lens, 100000);
    let (ctg_total, ctg_count) = count_by_min(ctg_lens, 100000);
    let (gap_total, gap_count) = count_by_min(gap_lens, 100000);
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        ">=100kb", len_total, len_count, ctg_total, ctg_count, gap_total, gap_count
    );
    let (len_total, len_count) = count_by_min(lens, 1000000);
    let (ctg_total, ctg_count) = count_by_min(ctg_lens, 1000000);
    let (gap_total, gap_count) = count_by_min(gap_lens, 1000000);
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        ">=1mb", len_total, len_count, ctg_total, ctg_count, gap_total, gap_count
    );
    println!("{:=<7}{:=^26}{:=^26}{:=^25}", "", "", "", "");
}

fn out_attr(
    head: &str,
    len: usize,
    attr: &str,
    attr_lower: &str,
    letter_counts: &HashMap<u8, usize>,
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
            //TODO change to use attr
            "id" => print!("{}", head),
            "len" => print!("{}", len),
            "" => (),
            _ => {
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
                           or multiple bases (gcGC means count g + c + G + C)."
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
    let mut records = parse_path(path).unwrap();
    if let Some(attr) = args.value_of("attr") {
        let attr = attr.trim_matches(':');
        let attr_lower = attr.to_ascii_lowercase();
        let is_total = !attr_lower.contains("id");
        let mut total_len = 0;
        let need_count = attr_lower.split(':').any(|x| x != "id" && x != "len");
        let mut letter_counts: HashMap<u8, usize> = HashMap::with_capacity(10);
        while let Some(record) = records.iter_record().unwrap() {
            if need_count {
                record
                    .seq()
                    .as_bytes()
                    .iter()
                    .for_each(|c| *letter_counts.entry(*c).or_insert(0) += 1);
            }

            if !is_total {
                out_attr(
                    record.head(),
                    record.len(),
                    attr,
                    &attr_lower,
                    &letter_counts,
                );
                letter_counts.clear();
            } else {
                total_len += record.len();
            }
        }

        if is_total {
            out_attr("Total", total_len, attr, &attr_lower, &letter_counts);
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
                    out_info.insert(chr.to_owned(), vec![(start, end)]);
                } else {
                    out_info.insert(region.to_owned(), vec![(0, 0)]);
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
            .get_bytes();
        let genome_len = Byte::from_str(subarg.value_of("genome_len").unwrap())
            .unwrap()
            .get_bytes();
        let n_len = Byte::from_str(subarg.value_of("n_len").unwrap())
            .unwrap()
            .get_bytes();
        let mut lens = Vec::with_capacity(1024);
        let mut total: usize = 0;

        let re = Regex::new(&format!("(?i)N{{{w},}}", w = n_len)).unwrap();
        let mut ctg_lens = Vec::with_capacity(1024);
        let mut ctg_total: usize = 0;

        let mut gap_lens = Vec::with_capacity(1024);
        let mut gap_total: usize = 0;

        let mut ctg_file = if subarg.is_present("out_ctg") {
            let file = args.value_of("input").unwrap().to_owned() + ".ctg.fa";
            Some(
                File::create(file.to_owned())
                    .unwrap_or_else(|_| panic!("failed create file: {}", file)),
            )
        } else {
            None
        };
        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len();
            if len < min_len as usize {
                continue;
            }
            lens.push(len as u32);
            total += len;

            if n_len > 0 {
                let mut last_pos = 0;
                let mut ctg_count = 1;
                let seq = record.seq();
                for mat in re.find_iter(seq) {
                    if mat.start() > last_pos {
                        ctg_lens.push((mat.start() - last_pos) as u32);
                        ctg_total += mat.start() - last_pos;
                        if let Some(ref mut file) = ctg_file {
                            writeln!(
                                file,
                                ">{}_ctg{}\n{}",
                                record.head(),
                                ctg_count,
                                &seq[last_pos..mat.start()]
                            )
                            .unwrap();
                            ctg_count += 1;
                        }
                    }
                    gap_lens.push((mat.end() - mat.start()) as u32);
                    gap_total += mat.end() - mat.start();
                    last_pos = mat.end();
                }
                if len > last_pos {
                    ctg_lens.push((len - last_pos) as u32);
                    ctg_total += len - last_pos;
                    if let Some(ref mut file) = ctg_file {
                        writeln!(
                            file,
                            ">{}_ctg{}\n{}",
                            record.head(),
                            ctg_count,
                            &seq[last_pos..len]
                        )
                        .unwrap();
                    }
                }
            }
        }
        lens.sort_unstable();
        if n_len > 0 {
            ctg_lens.sort_unstable();
            gap_lens.sort_unstable();
            out_stat_with_ctg(
                &lens,
                total,
                &ctg_lens,
                ctg_total,
                &gap_lens,
                gap_total,
                genome_len as usize,
            );
        } else {
            out_stat(&lens, total);
        }
    }
}
