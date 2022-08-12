use byte_unit::Byte;
use clap::{crate_name, App, AppSettings, Arg};
use indoc::indoc;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

mod utils;
use utils::{
    attr::attr,
    diff::diff,
    findgap::findgap,
    findseq::findseq,
    getseq::getseq,
    reform::reform,
    sample::sample,
    split::split,
    stat::{stat, sum_fx},
};

pub const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));

fn main() {
    let args = App::new("fxTools")
        .version(VERSION)
        .about("A toolset for processing sequences in FASTA/Q formats")
        .override_usage(format!("{} [SUBCOMMAND] [OPTIONS] <input>...", crate_name!()).as_str())
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .global_setting(AppSettings::ColoredHelp)
        .global_setting(AppSettings::DisableHelpSubcommand)
        .global_setting(AppSettings::DeriveDisplayOrder)
        .arg(
            Arg::new("input")
                .help("input file ...")
                .multiple_occurrences(true)
                .global(true)
        )
        .arg(
            Arg::new("attr")
                .short('a')
                .long("attr")
                .value_name("STR")
                .help(indoc!{"
                    get sequence attributes, id:len:x:qs
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
                .help(indoc!{"
                    reform/modify the sequence, accept values:
                      lower:    convert sequences into lowercase
                      upper:    convert sequences into uppercase
                      fq2fa:    converts FASTQ to FASTA
                      fa2fq:    converts FASTA to FASTQ
                      rev:      reverse the sequence
                      com:      complement the sequence
                      rc:       reverse and complement the sequence
                      lineINT[c]:  
                                wrap sequences into INT characters per line, 
                                0 for no wrap, c is the alignment mode, 
                                color for SNVs and INDELs
                      linkINT:  link sequences with INT Ns"
                      // splitINT: split sequences at Ns with length >= INT
                })
                .takes_value(true)
        )
        .arg(
            Arg::new("split")
                .short('p')
                .long("split")
                .value_name("INT")
                .help("split file with INT subfiles in total")
                .takes_value(true),
        )
        .arg(
            Arg::new("sample")
                .short('s')
                .long("sample")
                .value_name("int[G|M|K]|float")
                .help("subsample reads, int is the total bases, float is the fraction")
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
                        .help("minimum sequence length, shorter sequences are ignored")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("n_len")
                        .short('n')
                        .long("n_len")
                        .value_name("int[G|M|K]")
                        .default_value("0")
                        .help("minimum gap (N) length, 0 means do not stat contig length")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("genome_len")
                        .short('g')
                        .long("genome_len")
                        .value_name("int[G|M|K]")
                        .default_value("0")
                        // .requires("n_len")
                        .help("genome size, 0 means calculate genome size using the input file, co-used with -n")
                        .takes_value(true),
                )
                .arg(
                    Arg::new("out_ctg")
                        .short('o')
                        .long("out_ctg")
                        .help("output contig sequences to <INPUT>.ctg.fa")
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
                        .help("sequence to be find")
                        .required(true)
                        .takes_value(true)
                )
                .arg(
                    Arg::new("ignore-case")
                        .short('i')
                        .long("ignore-case")
                        .help("ignore case")
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
                    .default_value("1")
                    .help("minimum gap length, shorter gaps are ignored")
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
                        .help("region (chr, chr:start-end, chr-start-end or chr1,chr2:start-end) or file (bed or ID list) \
                            to be extracted, format: 0-based, [start, end), random access with a fai index file is support.")
                        .required(true)
                        .takes_value(true)
                )
                .arg(
                    Arg::new("reverse")
                        .short('v')
                        .long("reverse")
                        .help("reverse the sequence")
                )
                .arg(
                    Arg::new("complement")
                        .short('c')
                        .long("complement")
                        .help("complement the sequence")
                )
        )
        .subcommand(
            App::new("diff")
                .about("compare sequences between two files")
                .arg(
                    Arg::new("base")
                        .short('b')
                        .long("base")
                        .hide(true)
                        .help("compare sequences with same ID base by base")
                )
        )
        .get_matches();

    let paths = args
        .get_many::<String>("input")
        .map_or_else(|| vec!["-"], |v| v.map(|x| x.as_str()).collect());
    if let Some(v) = args.value_of("attr") {
        attr(&paths, v);
    } else if let Some(v) = args.value_of("reform") {
        reform(&paths, v);
    } else if let Some(v) = args.value_of("split") {
        let v = v.parse::<usize>().expect("not a valid split size");
        assert!(v > 0, "not a valid split size");
        split(&paths, v);
    } else if let Some(v) = args.value_of("sample") {
        let fra = Byte::from_str(&v).unwrap().get_bytes();
        let fra = if fra > 0 {
            let sum = sum_fx(&paths);
            fra as f64 / sum as f64
        } else {
            v.parse::<f64>().expect("not a valid sample size")
        };
        sample(&paths, fra as f64);
    } else if let Some(v) = args.subcommand_matches("findseq") {
        findseq(
            &paths,
            v.value_of("subseq").unwrap(),
            v.is_present("ignore-case"),
        );
    } else if let Some(subarg) = args.subcommand_matches("findgap") {
        let w = Byte::from_str(subarg.value_of("min_len").unwrap())
            .unwrap()
            .get_bytes();
        if w > 0 {
            findgap(&paths, w as usize);
        }
    } else if let Some(subarg) = args.subcommand_matches("getseq") {
        let region = subarg.value_of("region").unwrap();
        getseq(
            &paths,
            region,
            subarg.is_present("reverse"),
            subarg.is_present("complement"),
        );
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
        stat(
            &paths,
            min_len,
            genome_len,
            n_len,
            subarg.is_present("out_ctg"),
        );
    } else if let Some(_subarg) = args.subcommand_matches("diff") {
        diff(&paths);
    }
}
