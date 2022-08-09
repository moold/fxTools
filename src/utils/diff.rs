use hashbrown::HashMap;
use kseq::parse_path;

fn read_to_dict(infile: &str) -> HashMap<String, String>{

	let mut seqs = HashMap::new();
	let mut records = parse_path(infile).unwrap();
	while let Ok(Some(record)) = records.iter_record() {
		seqs.insert(
			record.head().to_string(),
			record.seq().to_string(),
		);
	}
	seqs
}

pub fn diff(paths: &[&str]){
	assert_eq!(paths.len(), 2, "require two input files only");
	let mut seqs = read_to_dict(paths[0]);
	let mut records = parse_path(paths[1]).unwrap();
	while let Ok(Some(record)) = records.iter_record() {
		if let Some(seq) = seqs.get(record.head()){
			if seq == record.seq(){
				println!("{}\tsame", record.head());
			}else {
				println!("{}\tdiff", record.head());
			}
			seqs.remove(record.head());
		}else {
			println!("{}\tfile2", record.head());
		}
	}
	for (key, _) in seqs.iter() {
		println!("{}\tfile1", key);
	}
}
