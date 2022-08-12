use super::common::parse_fx;
use regex::Regex;

pub fn findseq(paths: &[&str], subseq: &str, ignore_case: bool) {
    let re = if ignore_case {
        Regex::new(format!("(?i){}", subseq).as_str()).unwrap()
    } else {
        Regex::new(subseq).unwrap()
    };
    for path in paths {
        let mut records = parse_fx(*path);
        while let Ok(Some(record)) = records.iter_record() {
            let seq = record.seq();
            for mat in re.find_iter(seq) {
                println!("{}\t{}\t{}", record.head(), mat.start(), mat.end() - 1);
            }
        }
    }
}
