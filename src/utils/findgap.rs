use kseq::parse_path;
use regex::Regex;

pub fn findgap(paths: &[&str], w: usize) {
    let re = Regex::new(&format!("(?i)N{{{w},}}", w = w)).unwrap();
    for path in paths {
        let mut records = parse_path(*path).unwrap();
        while let Ok(Some(record)) = records.iter_record() {
            let seq = record.seq();
            for mat in re.find_iter(seq) {
                println!("{}\t{}\t{}", record.head(), mat.start(), mat.end() - 1);
            }
        }
    }
}