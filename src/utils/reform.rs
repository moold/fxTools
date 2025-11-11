use super::common::{parse_fx, print_fx, print_seq};
use kseq::record::Fastx;
use owo_colors::OwoColorize;
use regex::Regex;

fn lower(r: Fastx) {
    if r.sep().is_empty() {
        println!(
            ">{} {}\n{}",
            r.head(),
            r.des(),
            r.seq().to_ascii_lowercase()
        );
    } else {
        println!(
            "@{} {}\n{}\n{}\n{}",
            r.head(),
            r.des(),
            r.seq().to_ascii_lowercase(),
            r.sep(),
            r.qual()
        );
    }
}

fn upper(r: Fastx) {
    if r.sep().is_empty() {
        println!(
            ">{} {}\n{}",
            r.head(),
            r.des(),
            r.seq().to_ascii_uppercase()
        );
    } else {
        println!(
            "@{} {}\n{}\n{}\n{}",
            r.head(),
            r.des(),
            r.seq().to_ascii_uppercase(),
            r.sep(),
            r.qual()
        );
    }
}

fn wrapc(path: &str) {
    fn out_cbase<T: Into<char>>(b: T, cs: &[i32; 6]) {
        let b = b.into();
        let c: i32 = cs.iter().sum();
        if c <= 1 {
            // only one base
            print!("{b}")
        } else {
            // SNPs or INDELs
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

    fn count_bases(path: &str) -> Vec<[i32; 6]> {
        let mut bases = Vec::new();
        let mut records = parse_fx(path);
        while let Ok(Some(record)) = records.iter_record() {
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
        bases
    }

    let bases = count_bases(path);
    let mut records = parse_fx(path);
    while let Ok(Some(record)) = records.iter_record() {
        println!(">{} {}", record.head(), record.des());
        let (seq, len) = (record.seq().as_bytes(), record.len());
        for i in 0..len {
            out_cbase(seq[i], &bases[i]);
        }
        println!();
    }
}

fn split(r: Fastx, re: &Regex) {
    let mut last_pos = 0;
    let (seq, len) = (r.seq(), r.len());
    for mat in re.find_iter(seq) {
        println!(
            ">{}:{}_{}\n{}",
            r.head(),
            last_pos,
            mat.start() - 1,
            &seq[last_pos..mat.start()]
        );
        last_pos = mat.end();
    }
    if len > last_pos {
        if last_pos == 0 {
            println!(">{}\n{}", r.head(), &seq[last_pos..len]);
        } else {
            println!(
                ">{}:{}_{}\n{}",
                r.head(),
                last_pos,
                len - 1,
                &seq[last_pos..len]
            );
        }
    }
}

pub fn reform(paths: &[&str], reform: &str) {
    for path in paths {
        let mut records = parse_fx(path);

        if reform == "lower" {
            while let Ok(Some(record)) = records.iter_record() {
                lower(record);
            }
        } else if reform == "upper" {
            while let Ok(Some(record)) = records.iter_record() {
                upper(record);
            }
        } else if reform == "fq2fa" {
            while let Ok(Some(record)) = records.iter_record() {
                println!(">{} {}\n{}", record.head(), record.des(), record.seq());
            }
        } else if reform == "fa2fq" {
            while let Ok(Some(record)) = records.iter_record() {
                println!(
                    "@{} {}\n{}\n+\n{}",
                    record.head(),
                    record.des(),
                    record.seq(),
                    "I".repeat(record.len())
                );
            }
        } else if reform.starts_with("line") {
            let is_align = reform.ends_with('c');
            let w: usize = reform
                .strip_suffix('c')
                .unwrap_or(reform)
                .strip_prefix("line")
                .unwrap()
                .parse()
                .unwrap();
            if is_align {
                wrapc(path);
            } else {
                while let Ok(Some(record)) = records.iter_record() {
                    print_fx(record, w);
                }
            }
        } else if reform.starts_with("link") {
            let w: usize = reform.strip_prefix("link").unwrap().parse().unwrap();
            println!(">link_reads");
            let mut is_head = true;
            while let Ok(Some(record)) = records.iter_record() {
                if is_head {
                    is_head = false;
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
            let re = Regex::new(&format!("(?i)N{{{w},}}")).unwrap();
            while let Ok(Some(record)) = records.iter_record() {
                split(record, &re);
            }
        } else if ["rev", "com", "rc"].contains(&reform) {
            let (rev, com) = if reform == "rev" {
                (true, false)
            } else if reform == "com" {
                (false, true)
            } else {
                (true, true)
            };
            while let Ok(Some(record)) = records.iter_record() {
                if record.sep().is_empty() {
                    println!(">{} {}", record.head(), record.des());
                } else {
                    println!("@{} {}", record.head(), record.des());
                }
                print_seq(record.seq(), rev, com);
                if !record.sep().is_empty() {
                    println!("{}", record.sep());
                    print_seq(record.qual(), rev, false);
                }
            }
        } else {
            panic!("unknown values: {reform} for --reform");
        }
    }
}
