use super::common::parse_fx;
use regex::Regex;

pub fn findgap(paths: &[&str], w: usize) {
    let re = Regex::new(&format!("(?i)N{{{w},}}")).unwrap();
    for path in paths {
        let mut records = parse_fx(path);
        while let Ok(Some(record)) = records.iter_record() {
            let seq = record.seq();
            for mat in re.find_iter(seq) {
                println!("{}\t{}\t{}\t{}\t{}", record.head(), mat.start(), mat.end() - 1, mat.end() - mat.start(), record.len());
            }
        }
    }
}
