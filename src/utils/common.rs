use kseq::{parse_path, parse_reader, record::Fastx, Paths};
use std::io::Cursor;
use std::path::Path;

const SEQ_COMP_TABLE: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, b'T', b'V', b'G', b'H', b'E', b'F',
    b'C', b'D', b'I', b'J', b'M', b'L', b'K', b'N', b'O', b'P', b'Q', b'Y', b'S', b'A', b'A', b'B',
    b'W', b'X', b'R', b'Z', 91, 92, 93, 94, 95, 64, b't', b'v', b'g', b'h', b'e', b'f', b'c', b'd',
    b'i', b'j', b'm', b'l', b'k', b'n', b'o', b'p', b'q', b'y', b's', b'a', b'a', b'b', b'w', b'x',
    b'r', b'z', 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
    139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
    158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176,
    177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195,
    196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214,
    215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
    234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252,
    253, 254, 255,
];

#[allow(dead_code)]
pub fn complement_base(b: char) -> char {
    SEQ_COMP_TABLE[b as usize] as char
}

#[allow(dead_code)]
pub fn complement_bases(seq: &str) -> String {
    let t = seq
        .chars()
        .map(|x| SEQ_COMP_TABLE[x as usize])
        .collect::<Vec<u8>>();
    unsafe { String::from_utf8_unchecked(t) }
}

#[allow(dead_code)]
pub fn reverse_complement_base(b: char) -> char {
    complement_base(b)
}

#[allow(dead_code)]
pub fn reverse_complement_bases(seq: &str) -> String {
    let t = seq
        .chars()
        .rev()
        .map(|x| SEQ_COMP_TABLE[x as usize])
        .collect::<Vec<u8>>();
    unsafe { String::from_utf8_unchecked(t) }
}

pub fn print_seq(seq: &str, reverse: bool, complement: bool) {
    let iter_seqs: Box<dyn Iterator<Item = char>> = if reverse {
        Box::new(seq.chars().rev())
    } else {
        Box::new(seq.chars())
    };

    if complement {
        iter_seqs.for_each(|x| print!("{}", SEQ_COMP_TABLE[x as usize] as char));
    } else {
        iter_seqs.for_each(|x| print!("{x}"));
    }
    println!();
}

pub fn print_fx(r: Fastx, w: usize) {
    if r.sep().is_empty() {
        println!(">{} {}", r.head(), r.des());
        if w == 0 {
            println!("{}", r.seq());
        } else {
            let (seq, len) = (r.seq(), r.len());
            for i in (0..len).step_by(w) {
                if i + w < len {
                    println!("{}", &seq[i..i + w]);
                } else {
                    println!("{}", &seq[i..len]);
                }
            }
        }
    } else {
        println!("@{} {}", r.head(), r.des());
        if w == 0 {
            println!("{}\n{}\n{}", r.seq(), r.sep(), r.qual());
        } else {
            let (seq, qual, len) = (r.seq(), r.qual(), r.len());
            for i in (0..len).step_by(w) {
                if i + w < len {
                    println!("{}", &seq[i..i + w]);
                } else {
                    println!("{}", &seq[i..len]);
                }
            }
            println!("{}", r.sep());
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

pub fn write_fx<W>(r: Fastx, w: usize, writer: &mut W)
where
    W: std::io::Write,
{
    if r.sep().is_empty() {
        writeln!(writer, ">{} {}", r.head(), r.des()).expect("failed to write result");
        if w == 0 {
            writeln!(writer, "{}", r.seq()).expect("failed to write result");
        } else {
            let (seq, len) = (r.seq(), r.len());
            for i in (0..len).step_by(w) {
                if i + w < len {
                    writeln!(writer, "{}", &seq[i..i + w]).expect("failed to write result");
                } else {
                    writeln!(writer, "{}", &seq[i..len]).expect("failed to write result");
                }
            }
        }
    } else {
        writeln!(writer, "@{} {}", r.head(), r.des()).expect("failed to write result");
        if w == 0 {
            writeln!(writer, "{}\n{}\n{}", r.seq(), r.sep(), r.qual())
                .expect("failed to write result");
        } else {
            let (seq, qual, len) = (r.seq(), r.qual(), r.len());
            for i in (0..len).step_by(w) {
                if i + w < len {
                    writeln!(writer, "{}", &seq[i..i + w]).expect("failed to write result");
                } else {
                    writeln!(writer, "{}", &seq[i..len]).expect("failed to write result");
                }
            }
            writeln!(writer, "{}", r.sep()).expect("failed to write result");
            for i in (0..len).step_by(w) {
                if i + w < len {
                    writeln!(writer, "{}", &qual[i..i + w]).expect("failed to write result");
                } else {
                    writeln!(writer, "{}", &qual[i..len]).expect("failed to write result");
                }
            }
        }
    }
}

pub fn is_fasta_record(r: &Fastx) -> bool {
    r.sep().is_empty()
}

pub fn is_fasta_file(path: &str) -> bool {
    let p = path.to_lowercase();
    if p.ends_with("fasta") || p.ends_with("fa") || p.ends_with("fa.gz") || p.ends_with("fasta.gz")
    {
        true
    } else if p.ends_with("fastq")
        || p.ends_with("fq")
        || p.ends_with("fq.gz")
        || p.ends_with("fastq.gz")
    {
        false
    } else {
        let mut records = parse_path(path).unwrap();
        if let Some(record) = records.iter_record().unwrap() {
            return is_fasta_record(&record);
        }
        false
    }
}

pub fn parse_fx(file: &str) -> Paths<'_> {
    if !Path::new(file).exists() && file.chars().all(|x| "ATGCNatgcn".contains(x)) {
        let file = format!(">unname\n{file}");
        parse_reader(Cursor::new(file)).unwrap()
    } else {
        parse_path(file).unwrap()
    }
}
