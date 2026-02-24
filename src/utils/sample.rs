use super::common::{parse_fx, print_fx};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

pub fn sample(paths: &[&str], fra: f64, seed_opt: Option<u64>) {
    let fraction = fra.clamp(0.0, 1.0); 
    let base_seed = seed_opt.unwrap_or_else(|| rand::random::<u64>());
    
    let mut total_bp = 0;
    let mut sampled_reads = 0;

    for path in paths {
        let mut rng = SmallRng::seed_from_u64(base_seed);
        
        let mut records = parse_fx(path);
        
        while let Ok(Some(record)) = records.iter_record() {
            if rng.gen_bool(fraction) {
                total_bp += record.len();
                sampled_reads += 1;
                
                let w = if record.sep().is_empty() { 100 } else { 0 };
                print_fx(record, w);
            }
        }
    }
    
    eprintln!("sample reads: {}, sample size: {:?} bp", sampled_reads, total_bp);
}