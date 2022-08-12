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
    fn new(name: String) -> Self {
        let file = File::create(&name).expect("failed create out file!");
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
    if Path::new(&fname).exists() {
        panic!("failed create out directory {:?}!", fname);
    }
    create_dir(&fname).unwrap();
    let mut files = Vec::with_capacity(count);
    for i in 0..count {
        files.push(Outfile::new(format!(
            "{fname}/subseq{:02}.{}",
            i,
            if is_fa { "fasta" } else { "fastq" }
        )));
    }
    files
}

pub fn split(paths: &[&str], count: usize) {
    let outdir = format!("{}.split{count}", paths[0]);
    let mut outfiles = creat_outfiles(&outdir, count, is_fasta_file(paths[0]));

    for path in paths {
        let mut records = parse_fx(*path);
        while let Some(record) = records.iter_record().unwrap() {
            let w = outfiles.iter_mut().min_by_key(|v| v.size).unwrap();
            w.write(record);
        }
    }
}
