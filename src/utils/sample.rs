use super::common::{parse_fx, print_fx};
use rand::seq::index::sample as rsample;

pub fn min_count(paths: &[&str]) -> usize {
    let mut n = 0;
    for path in paths {
        let mut records = parse_fx(path);
        while let Ok(Some(_)) = records.iter_record() {
            n += 1;
            if n >= 1000 {
                break;
            }
        }
    }
    n
}

pub fn sample_index(count: usize, fra: f64) -> Vec<u8> {
    let mut v = vec![0; count];
    let mut rng = rand::thread_rng();
    let total = (count as f64 * fra).ceil() as usize;
    rsample(&mut rng, count, total)
        .iter()
        .for_each(|i| v[i] = 1);
    v
}

pub fn sample(paths: &[&str], fra: f64) {
    let min_c = min_count(paths);
    let index = sample_index(min_c, fra);
    let mut i = 0;
    let mut total = 0;
    for path in paths {
        let mut records = parse_fx(path);
        while let Some(record) = records.iter_record().unwrap() {
            if index[i] > 0 {
                total += record.len();
                let w = if record.sep().is_empty() { 100 } else { 0 };
                print_fx(record, w);
            }
            i += 1;
            if i >= min_c {
                i = 0;
            }
        }
    }
    eprintln!("sample size: {total:?} bp");
}
