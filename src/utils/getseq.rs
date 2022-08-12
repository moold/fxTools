use super::common::{is_fasta_record, parse_fx, print_seq};
use hashbrown::HashMap;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

fn get_out_info(path: &str) -> HashMap<String, Vec<(u32, u32)>> {
    let mut out_info: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    if let Ok(file) = File::open(path) {
        for line in BufReader::new(file).lines().flatten() {
            if !line.starts_with('#') {
                let v: Vec<&str> = line.split(char::is_whitespace).collect();
                let out_info = out_info.entry(v[0].to_owned()).or_insert(vec![]);
                if v.len() == 1 {
                    out_info.push((0, 0))
                } else {
                    let start = v[1].parse::<u32>().unwrap_or(0);
                    let end = v[2].parse::<u32>().unwrap_or(0);
                    if start != end || end == 0 {
                        out_info.push((start, end));
                    }
                }
            }
        }
    } else {
        let re = Regex::new(r"(?i)(\S+)[:-](start|\d+)[:-](end|\d+)$").unwrap();
        for region in path.split(&[',', ';'][..]) {
            if let Some(caps) = re.captures(region) {
                let chr = caps.get(1).unwrap().as_str();
                let start = caps.get(2).unwrap().as_str().parse::<u32>().unwrap_or(0);
                let end = caps.get(3).unwrap().as_str().parse::<u32>().unwrap_or(0);
                if start != end || end == 0 {
                    out_info
                        .entry(chr.to_owned())
                        .or_insert(vec![])
                        .push((start, end));
                }
            } else {
                out_info
                    .entry(region.to_owned())
                    .or_insert(vec![])
                    .push((0, 0));
            }
        }
    }
    out_info
}

fn read_fai(path: &str) -> Result<HashMap<String, (u32, u64)>, std::io::Error> {
    match File::open(path) {
        Ok(file) => {
            let mut out_info: HashMap<String, (u32, u64)> = HashMap::new();
            for line in BufReader::new(file).lines().flatten() {
                let v: Vec<&str> = line.split(char::is_whitespace).collect();
                out_info.insert(
                    v[0].to_owned(),
                    (v[1].parse::<u32>().unwrap(), v[2].parse::<u64>().unwrap()),
                );
            }
            Ok(out_info)
        }
        Err(e) => Err(e),
    }
}

fn read_next_seq<R: Read>(r: &mut BufReader<R>, mut seq: String, len: usize) -> String {
    seq.clear();
    while let Ok(n) = r.read_line(&mut seq) {
        if seq.ends_with('\n') {
            seq.pop();
        }
        if seq.len() >= len || n == 0 {
            break;
        }
    }
    seq
}

fn get_head_reset_region(
    start: &mut u32,
    end: &mut u32,
    len: u32,
    name: &str,
    reverse: bool,
    complement: bool,
) -> String {
    let mut sub_head;
    if *start == *end {
        *start = 0;
        *end = len;
        sub_head = name.to_owned();
    } else {
        if *end == 0 {
            *end = len;
        }
        sub_head = format!("{}:{}-{}", name, start, *end - 1);
    }
    if reverse {
        sub_head += "_rev";
    }
    if complement {
        sub_head += "_com";
    }
    sub_head
}

fn getseq_by_index(
    path: &str,
    fai: &HashMap<String, (u32, u64)>,
    infos: &mut HashMap<String, Vec<(u32, u32)>>,
    reverse: bool,
    complement: bool,
) {
    let mut seq = String::new();
    let mut file = BufReader::new(File::open(path).expect("failed read file"));
    let mut valid_seq = Vec::new();
    for (head, regions) in infos.iter() {
        if let Some((len, off)) = fai.get(head) {
            file.seek(SeekFrom::Start(*off))
                .expect("not a correct fai index file");
            seq = read_next_seq(&mut file, seq, *len as usize);
            for (mut start, mut end) in regions {
                let sub_head = get_head_reset_region(
                    &mut start,
                    &mut end,
                    *len as u32,
                    head,
                    reverse,
                    complement,
                );
                println!(">{}", sub_head);
                print_seq(&seq[start as usize..end as usize], reverse, complement);
            }
            valid_seq.push(head.to_owned());
        }
    }
    valid_seq.iter().for_each(|x| {
        infos.remove(x.as_str());
    });
}

pub fn getseq(paths: &[&str], region: &str, reverse: bool, complement: bool) {
    let mut infos = get_out_info(region);
    for path in paths {
        let fai = (*path).to_owned() + ".fai";
        if let Ok(fai_infos) = read_fai(&fai) {
            eprintln!("Note: detected fai file {}, enable random access mode", fai);
            getseq_by_index(path, &fai_infos, &mut infos, reverse, complement);
        } else {
            let mut records = parse_fx(*path);
            while let Ok(Some(record)) = records.iter_record() {
                let head = record.head();
                if let Some(regions) = infos.get_mut(head) {
                    let des = record.des();
                    let seqs = record.seq();
                    let qual = record.qual();
                    for (mut start, mut end) in regions {
                        let sub_head = get_head_reset_region(
                            &mut start,
                            &mut end,
                            record.len() as u32,
                            head,
                            reverse,
                            complement,
                        );
                        if is_fasta_record(&record) {
                            println!(">{} {}", sub_head, des);
                        } else {
                            println!("@{} {}", sub_head, des);
                        }
                        print_seq(&seqs[start as usize..end as usize], reverse, complement);
                        if !qual.is_empty() {
                            println!("{}", record.sep());
                            print_seq(&qual[start as usize..end as usize], reverse, false);
                        }
                    }
                }
                infos.remove(head);
            }
        }
        if infos.is_empty() {
            break;
        }
    }
    for (key, _) in infos.iter() {
        eprintln!("Missing record in the fastx/fai file: {}", key);
    }
}
