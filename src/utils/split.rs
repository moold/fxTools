use super::common::{is_fasta_file, parse_fx, write_fx};
use kseq::record::Fastx;
use std::{
    fs::{create_dir, File},
    io::BufWriter,
    path::Path,
};

struct Outfile {
    size: usize,
    handle: BufWriter<File>,
}

impl Outfile {
    fn new(name: &str) -> Self {
        let file = File::create(name).expect("failed create out file!");
        Outfile {
            size: 0,
            handle: BufWriter::with_capacity(1024000, file),
        }
    }

    fn write(&mut self, r: Fastx) {
        self.size += r.len();
        write_fx(r, 0, &mut self.handle);
    }
}

fn creat_outfiles(fname: &str, count: usize, is_fa: bool) -> Vec<Outfile> {
    let mut files = Vec::with_capacity(count);
    for i in 0..count {
        files.push(Outfile::new(&format!(
            "{fname}/subseq{:02}.{}",
            i,
            if is_fa { "fasta" } else { "fastq" }
        )));
    }
    files
}

pub fn splits(paths: &[&str], count: usize) {
    let outdir = &format!("{}.split{count}", paths[0]);
    if Path::new(outdir).exists() {
        panic!("failed create out directory {outdir:?}!");
    }
    create_dir(outdir).unwrap();

    let mut outfiles = creat_outfiles(&outdir, count, is_fasta_file(paths[0]));

    for path in paths {
        let mut records = parse_fx(path);
        while let Some(record) = records.iter_record().unwrap() {
            let w = outfiles.iter_mut().min_by_key(|v| v.size).unwrap();
            w.write(record);
        }
    }
}

pub fn splitr(paths: &[&str], count: usize) {
    let outdir = &format!("{}.split{count}", paths[0]);
    if Path::new(outdir).exists() {
        panic!("failed create out directory {outdir:?}!");
    }
    create_dir(outdir).unwrap();

    let is_fa = is_fasta_file(paths[0]);
    if count == 1 {
        for path in paths {
            let mut records = parse_fx(path);
            while let Some(record) = records.iter_record().unwrap() {
                let mut w = Outfile::new(&format!("{outdir}/{}.{}", record.head(), if is_fa { "fasta" } else { "fastq" }));
                w.write(record);
            }
        }
    }else{
        let mut i = 0;
        let mut curr_count = 0;
        let mut w: Option<Outfile> = None;
        for path in paths {
            let mut records = parse_fx(path);
            while let Some(record) = records.iter_record().unwrap() {
                if curr_count == 0 {
                    w = Some(Outfile::new(&format!("{outdir}/subseq{:02}.{}", i, if is_fa { "fasta" } else { "fastq" })));
                    i += 1;
                }
                w.as_mut().unwrap().write(record);
                curr_count += 1;
                if curr_count >= count {
                    curr_count = 0;
                }
            }
        }
    }
}
