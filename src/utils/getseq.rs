use kseq::parse_path;
use hashbrown::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::Regex;
use super::common::out_seq;

fn get_out_info(path: &str) -> HashMap<String, Vec<(u32, u32)>> {
	let mut out_info: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

	if let Ok(file) = File::open(path){
		for line in BufReader::new(file).lines().flatten() {
			if !line.starts_with('#') {
				let v: Vec<&str> = line.split(char::is_whitespace).collect();
				let out_info = out_info.entry(v[0].to_owned()).or_insert(vec![]);
				if v.len() == 1 {
					out_info.push((0, 0))
				} else {
					out_info.push((
						v[1].parse::<u32>().unwrap_or(0),
						v[2].parse::<u32>().unwrap_or(0),
					))
				}
			}
		}
	} else {
		for region in path.split(&[',', ';'][..]) {
			let re = Regex::new(r"(?i)(\S+)[:-](start|\d+)[:-](end|\d+)$").unwrap();
			if let Some(caps) = re.captures(region) {
				let chr = caps.get(1).unwrap().as_str();
				let start = caps.get(2).unwrap().as_str().parse::<u32>().unwrap_or(0);
				let end = caps.get(3).unwrap().as_str().parse::<u32>().unwrap_or(0);
				out_info.entry(chr.to_owned()).or_insert(vec![]).push((start, end));
			} else {
				out_info.entry(region.to_owned()).or_insert(vec![]).push((0, 0));
			}
		}
	}
	out_info
}

pub fn getseq(paths: &[&str], region: &str, reverse: bool, complement: bool) {

	let mut infos = get_out_info(region);
	for path in paths {
		let mut records = parse_path(*path).unwrap();
		while let Ok(Some(record)) = records.iter_record() {
			let head = record.head();
			if let Some(regions) = infos.get_mut(head) {
				let des = record.des();
				let seqs = record.seq();
				let qual = record.qual();
				for (start, end) in regions.iter_mut() {
					let mut sub_head: String;
					if start == end {
						*start = 0;
						*end = record.len() as u32;
						sub_head = head.to_owned();
					} else {
						if *end == 0 {
							*end = record.len() as u32;
						}
						sub_head = format!("{}:{}-{}", head, start, *end - 1);
					}
					if reverse {
						sub_head += "_rev";
					}
					if complement {
						sub_head += "_com";
					}
					if qual.is_empty(){
						println!(">{} {}", sub_head, des);
					}else{
						println!("@{} {}", sub_head, des);
					}
					out_seq(
						&seqs[*start as usize..*end as usize],
						reverse,
						complement,
					);
					if !qual.is_empty(){
						println!("{}", record.sep());
						out_seq(
							&qual[*start as usize..*end as usize],
							reverse,
							false,
						);
					}    
				}
			}
			infos.remove(head);
			if infos.is_empty() {
				break;
			}
		}
		for (key, _) in infos.iter() {
			eprintln!("Missing record: {}", key);
		}
	}
}