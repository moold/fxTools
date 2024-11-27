use super::common::parse_fx;
use hashbrown::HashMap;

fn read_to_dict(infile: &str) -> HashMap<String, String> {
    let mut seqs = HashMap::new();
    let mut records = parse_fx(infile);
    while let Ok(Some(record)) = records.iter_record() {
        seqs.insert(record.head().to_string(), record.seq().to_string());
    }
    seqs
}

pub fn diff(paths: &[&str]) {
    assert_eq!(paths.len(), 2, "require two input files only");
    let mut seqs = read_to_dict(paths[0]);
    let mut records = parse_fx(paths[1]);
    let (mut same, mut diff, mut file1, mut file2) = (0, 0, 0, 0);
    while let Ok(Some(record)) = records.iter_record() {
        if let Some(seq) = seqs.get(record.head()) {
            if seq.eq_ignore_ascii_case(record.seq()) {
                println!("{}\tsame", record.head());
                same += 1;
            } else {
                println!("{}\tdiff", record.head());
                diff += 1;
            }
            seqs.remove(record.head());
        } else {
            println!("{}\tfile2", record.head());
            file2 += 1;
        }
    }
    for (key, _) in seqs.iter() {
        println!("{key}\tfile1");
        file1 += 1;
    }

    eprintln!("same records between files: {same}\ndifferent records between files: {diff}\nuniq records in file1: {file1}\nuniq records in file2: {file2}");
}
