use kseq::{parse_path, record::Fastx};

fn out_attr(
    record: Fastx,
    attr: &str,
    attr_lower: &str,
) {
    let mut index = 0;
    let mut is_head = true;

    for elm in attr_lower.split(':') {
        if is_head {
            is_head = false
        } else {
            print!("\t");
            index += 1;
        }
        match elm {
            "id" => print!("{}", record.head()),
            "len" => print!("{}", record.len()),
            "qs" => {
                let qual = record.qual();
                if !qual.is_empty() {
                     // \text{read Q} = -10\log_{10}\big[\tfrac{1}{N}\sum 10^{-q_i/10}\big]
                    let e_sum: f64 = qual.as_bytes().iter().map(|x| 10_f64.powf(- 0.1 * (x - 33) as f64)).sum();
                    print!("{}", -10.0 * (e_sum/record.len() as f64).log10());
                }else {
                    print!("NA");
                }
            }
            "" => (),
            _ => {
                let mut t = [0; 256];
                record.seq().as_bytes().iter().for_each(|c| t[*c as usize] += 1);
                let x = &attr[index..index + elm.len()];
                print!("{}", x.bytes().fold(0, |acc, x| acc + t[x as usize]));
            }
        }
        index += elm.len();
    }
    if !is_head {
        println!();
    }
}

pub fn attr(paths: &[&str], attr: &str){
    for path in paths {
        let mut records = parse_path(*path).unwrap();
        let attr = attr.trim_matches(':');
        let attr_lower = attr.to_ascii_lowercase();
        while let Ok(Some(record)) = records.iter_record() {
            out_attr(record, attr, &attr_lower);
        }
    }
}
