use std::env::args;
use kseq::path::parse_path;

fn main(){
	let path: Option<String> = args().nth(1);
	let mut records = parse_path(path).unwrap();
	while let Some(record) = records.iter_record() {
		println!("head:{} des:{} seq:{} qual:{} len:{}",
			record.head(), record.des(), record.seq(), record.qual(), record.len());
	}
}