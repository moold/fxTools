use super::common::parse_fx;
use crossbeam_channel::{bounded, unbounded};
use crossbeam_utils::thread;
use rayon::prelude::*;
use regex::Regex;
use std::{cmp::max, fmt, fs::File, io::Write};

mod io;
mod path;
use io::Buffer;
use path::open_path;

const BUF_COUNT: usize = 2;

#[derive(Default)]
struct Nx {
    count: [usize; 10],
    len: [u32; 10],
}

impl Nx {
    fn new() -> Self {
        Default::default()
    }

    fn get_width(&self) -> (usize, usize) {
        let mut w1 = self.count[8].to_string().len();
        let mut w2 = self.len[0].to_string().len();
        if w1 < 9 {
            w1 = 9;
        }
        if w2 < 11 {
            w2 = 11;
        }
        (w1, w2)
    }

    fn fill(mut self, lens: &[u32], total: usize) -> Self {
        let mut i = 0;
        let mut acc = 0;
        for len in lens.iter().rev() {
            acc += *len as usize;
            self.count[i] += 1;
            while acc as f64 > ((i + 1) * total) as f64 * 0.1 {
                self.len[i] = *len;
                i += 1;
                if i >= 10 {
                    break;
                }
                self.count[i] += self.count[i - 1];
            }
        }
        self
    }
}

impl fmt::Display for Nx {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (w1, w2) = self.get_width();
        writeln!(
            f,
            "{:<5} {:^w1$} {:^w2$}",
            "Types", "Count (#)", "Length (bp)",
        )?;
        for i in 0..9 {
            writeln!(
                f,
                "N{:<4} {:^w1$} {:^w2$}",
                (i + 1) * 10,
                self.count[i],
                self.len[i],
            )?;
        }
        Ok(())
    }
}

const BIN_LEN: u32 = 30;
const UNIT_BIN: usize = 200;
struct His {
    pw: usize,   // pos width
    cw: usize,   // count width
    unit: usize, // count of a '*'
    step: u32,
    start: u32,
    min: u32,
    max: u32,
    count: [usize; BIN_LEN as usize],
}

impl His {
    fn new(s: u32, e: u32, c: usize) -> Self {
        let step = max((e - s) / (BIN_LEN - 2), 1);
        Self {
            pw: e.to_string().len(),
            cw: c.to_string().len(),
            unit: max(c / UNIT_BIN, 5),
            step,
            start: s / step * step,
            min: u32::MIN,
            max: u32::MAX,
            count: [0; BIN_LEN as usize],
        }
    }

    fn fill(mut self, lens: &[u32]) -> Self {
        if !lens.is_empty() {
            let mut idx: usize = 0;
            let mut max_len = self.start;
            for len in lens {
                while idx + 1 < BIN_LEN as usize && *len >= max_len {
                    idx += 1;
                    max_len += self.step;
                }
                self.count[idx] += 1;
            }
            self.min = lens[0];
            self.max = *lens.last().unwrap();
            self.pw = self.max.to_string().len();
        }
        self
    }
}

impl fmt::Display for His {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "[length histogram ('*' =~ {} reads)]", self.unit)?;
        for (p, v) in self.count.into_iter().enumerate() {
            let p = p as u32;
            if v > 0 {
                writeln!(
                    f,
                    "{:>pw$} {:>pw$} {:>cw$} {:*>cv$}",
                    if p == 0 {
                        self.min
                    } else {
                        self.start + self.step * (p - 1)
                    },
                    if p == BIN_LEN - 1 {
                        self.max
                    } else {
                        self.start + self.step * p - 1
                    },
                    v,
                    "",
                    pw = self.pw,
                    cw = self.cw,
                    cv = v / self.unit
                )?;
            }
        }
        Ok(())
    }
}

fn out_step(lens: &[u32], step_len: u32, genome_len: usize) {
    thread::scope(|work| {
        let (in_s, in_r) = bounded(1024);
        // input thread
        work.spawn(move |_| {
            let mut step = 0;
            for (i, len) in lens.iter().enumerate(){
                while *len >= step + step_len{
                    step += step_len;
                }

                if *len >= step {
                    in_s.send((step, &lens[i..])).unwrap();
                    step += step_len;
                }
            }
        });

        //work thread
        let mut nxs = work.spawn(move |_| {
            let work_thread = 2;
            thread::scope(|scoped| {
                let mut handles = Vec::with_capacity(work_thread);
                for _i in 0..work_thread {
                    let in_r = in_r.clone();
                    let handle = scoped.spawn(move |_| {
                        let mut nxs = Vec::new();
                        while let Ok((step, lens)) = in_r.recv() {
                            let total = lens.iter().sum::<u32>();
                           let nx = Nx::new().fill(lens, if genome_len > 0 { genome_len } else { total as usize });
                           nxs.push((step, lens.len(), total, nx.count[4], nx.len[4]));
                        }
                        nxs
                    });
                    handles.push(handle);
                }

                let mut nxs: Vec<(u32, usize, u32, usize, u32)> = Vec::new();
                for res in handles.into_iter().map(|h| h.join().unwrap()){
                    nxs.extend(res);
                }
                nxs
            }).unwrap()
        }).join().unwrap();

        if !nxs.is_empty(){
            // sort by step
            nxs.sort_unstable_by_key(|k| k.0);
            let w0 = max(nxs.last().unwrap().0.to_string().len(), 5) + 2; //safe unwrap
            let w1 = max(nxs[0].3.to_string().len(), 9);
            let w2 = max(nxs[0].4.to_string().len(), 10);
            println!("{:<w$} {:^w1$} {:^w2$} {:^w1$} {:^w2$}", "Types", "Count", "Length", "N50 Count", "N50 Length", w = w0 + 2);
            println!("{:<w$} {:^w1$} {:^w2$} {:^w1$} {:^w2$}", "", "(#)", "(bp)", "(#)", "(bp)", w = w0 + 2);
            for (step, n50_count, n50_base, total_count, toal_base) in nxs {
                println!(">={step:<w0$} {n50_count:^w1$} {n50_base:^w2$} {total_count:^w1$} {toal_base:^w2$}");
            }
        }
    }).unwrap();
}

fn out_stat(lens: &[u32], total: usize, genome_len: usize) {
    let total_count = lens.len();
    let (hist, nx) = thread::scope(|work| {
        let hist = work.spawn(move |_| {
            let dev = total_count / UNIT_BIN * 2;
            His::new(lens[dev], lens[total_count - max(1, dev)], total_count).fill(lens)
        });

        let nx = work
            .spawn(move |_| Nx::new().fill(lens, if genome_len > 0 { genome_len } else { total }));
        (
            hist.join().expect("Failed to generate histogram!"),
            nx.join().expect("Failed to generate Nx stats!"),
        )
    })
    .unwrap();

    println!("{hist}");
    println!("\n\n[length stat]\n{nx}");
    let (sw1, sw2) = nx.get_width();
    println!("{:<5} {:^sw1$} {:^sw2$}", "Min.", "-", lens[0],);
    println!(
        "{:<5} {:^sw1$} {:^sw2$}",
        "Max.",
        "-",
        lens[total_count - 1],
    );
    println!("{:<5} {:^sw1$} {:^sw2$}", "Ave.", "-", total / total_count,);
    println!("{:<5} {:^sw1$} {:^sw2$}", "Total", total_count, total,);
}

fn out_stats(
    lens: &[u32],
    total: usize,
    ctg_lens: &[u32],
    ctg_total: usize,
    gap_lens: &[u32],
    gap_total: usize,
    genome_len: usize,
) {
    let acc_min = |lens: &[u32], min| {
        lens.iter()
            .filter(|&&x| x >= min)
            .fold((0, 0), |acc, x| (acc.0 + 1, acc.1 + x))
    };
    let (nx, ctg_nx, gap_nx) = thread::scope(|work| {
        let nx = work
            .spawn(move |_| Nx::new().fill(lens, if genome_len > 0 { genome_len } else { total }));
        let ctg_nx = work.spawn(move |_| {
            Nx::new().fill(
                ctg_lens,
                if genome_len > 0 {
                    genome_len
                } else {
                    ctg_total
                },
            )
        });
        let gap_nx = work.spawn(move |_| Nx::new().fill(gap_lens, gap_total));
        (
            nx.join().expect("Failed to generate Nx stats!"),
            ctg_nx.join().expect("Failed to generate ctg_Nx stats!"),
            gap_nx.join().expect("Failed to generate gap_Nx stats!"),
        )
    })
    .unwrap();

    println!("{:=<7}{:=^26}{:=^26}{:=^25}", "", "", "", "");
    println!(
        "{:<7}{:^26}{:^26}{:^25}",
        "Types", "Scaffold", "Contig", "Gap"
    );
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "", "Length (bp)", "Count (#)", "Length (bp)", "Count (#)", "Length (bp)", "Count (#)"
    );
    println!("{:-<7}{:-^26}{:-^26}{:-^25}", "", "", "", "");
    for i in 0..9 {
        println!(
            "N{:<6}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
            (i + 1) * 10,
            nx.len[i],
            nx.count[i],
            ctg_nx.len[i],
            ctg_nx.count[i],
            gap_nx.len[i],
            gap_nx.count[i],
        );
    }
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "Longest",
        lens.last().unwrap_or(&0),
        1,
        ctg_lens.last().unwrap_or(&0),
        1,
        gap_lens.last().unwrap_or(&0),
        if gap_lens.is_empty() { 0 } else { 1 },
    );
    println!(
        "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
        "Total",
        total,
        lens.len(),
        ctg_total,
        ctg_lens.len(),
        gap_total,
        gap_lens.len(),
    );

    for min_len in [10000, 100000, 1000000] {
        let (len_total, len_count) = acc_min(lens, min_len);
        let (ctg_total, ctg_count) = acc_min(ctg_lens, min_len);
        let (gap_total, gap_count) = acc_min(gap_lens, min_len);
        println!(
            "{:<7}{:^16}{:^10}{:^16}{:^10}{:^16}{:^9}",
            if min_len == 1000000 {
                ">=1mb"
            } else if min_len == 100000 {
                ">=100kb"
            } else {
                ">=10kb"
            },
            len_total,
            len_count,
            ctg_total,
            ctg_count,
            gap_total,
            gap_count
        );
    }
    println!("{:=<7}{:=^26}{:=^26}{:=^25}", "", "", "", "");
}

fn stat_read(
    infiles: &[&str],
    min_len: usize,
    genome_len: usize,
    step_len: usize,
    out: bool,
) -> (usize, Vec<u32>) {
    // exit if any thread panics
    let orig_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |v| {
        orig_hook(v);
        std::process::exit(1);
    }));

    let (mut lens, total) = thread::scope(|work| {
        let (s1, r1) = unbounded();
        let (s2, r2) = unbounded();
        for _ in 0..BUF_COUNT {
            s1.send(Buffer::new()).unwrap();
        }

        // read thread
        work.spawn(move |_| {
            for infile in infiles {
                for mut reader in open_path(infile) {
                    loop {
                        let mut buf = r1.recv().unwrap();
                        match buf.fill(&mut reader) {
                            Ok(0) => {
                                let is_empty = buf.is_empty();
                                s2.send(Some(buf)).unwrap(); //send an empty buf to indicate this file reaches EOF
                                if is_empty {
                                    break;
                                }
                            }
                            Ok(_n) => s2.send(Some(buf)).unwrap(),
                            Err(e) => panic!("Failed to read file: {infile:?}, error: {e:?}"),
                        }
                    }
                }
            }
            s2.send(None).unwrap();
            while r1.len() != BUF_COUNT {} //wait the stat thread to finish
        });

        // statistics thread
        let stat = work.spawn(move |_| {
            let mut skip_lines = 0;
            let mut skip_bases = 0;
            let mut is_new_record = true;

            let mut len = 0;
            let mut lens = Vec::with_capacity(1024000);
            let mut total: usize = 0;
            while let Ok(Some(mut buf)) = r2.recv() {
                if buf.is_empty() {
                    // a file has reached EOF
                    if skip_bases != 0 || skip_lines != 0 {
                        panic!("truncate file");
                    } else if len > min_len {
                        // save the last fasta record
                        lens.push(len as u32);
                        total += len;
                    }
                    len = 0;
                    is_new_record = true;
                    s1.send(buf).unwrap();
                    continue;
                }

                if skip_lines > 0 {
                    let skip_line = buf.skip_lines(skip_lines);
                    skip_lines -= skip_line;
                    if skip_lines > 0 {
                        s1.send(buf).unwrap();
                        continue;
                    }
                }

                if skip_bases > 0 {
                    let skip_base = buf.skip_bases(skip_bases);
                    skip_bases -= skip_base;
                    if skip_bases > 0 {
                        s1.send(buf).unwrap();
                        continue;
                    }
                }

                'outer: while let Some(c) = buf.next_byte(true) {
                    if !is_new_record {
                        if c == b'>' || c == b'@' {
                            if len > min_len {
                                // save the previous fasta record
                                lens.push(len as u32);
                                total += len;
                            }
                            len = 0;
                            is_new_record = true;
                            continue;
                        } else if c == b'+' {
                            if len > min_len {
                                // save the previous fasta record
                                lens.push(len as u32);
                                total += len;
                            }
                            let skip_line = buf.skip_lines(1); //skip sep
                            if skip_line != 1 {
                                skip_lines = 1;
                                skip_bases = len;
                                break;
                            }

                            let l = buf.skip_bases(len); // skip qual
                            is_new_record = true;
                            if l != len {
                                skip_bases = len - l;
                                len = 0;
                                break;
                            } else {
                                len = 0;
                                continue;
                            }
                        }
                    } else if c == b'>' || c == b'@' {
                        is_new_record = false;
                        let skip_line = buf.skip_lines(1); //skip head
                        if skip_line != 1 {
                            skip_lines = 1;
                            break;
                        }
                    } else if is_new_record && c != b'>' && c != b'@' {
                        panic!("Not a correct fasta/fastq file");
                    }

                    while let Some((p, _a)) = buf.next_line_len() {
                        // iter seq
                        len += p;
                        if let Some(c) = buf.next_byte(false) {
                            if c == b'>' {
                                // fasta
                                break;
                            } else if c == b'+' {
                                // fastq
                                let skip_line = buf.skip_lines(1); //skip sep
                                if skip_line != 1 {
                                    skip_lines = 1;
                                    skip_bases = len;
                                    break;
                                }

                                let l = buf.skip_bases(len); // skip qual
                                is_new_record = true;
                                if l != len {
                                    skip_bases = len - l;
                                }
                                break;
                            }
                        } else {
                            //need continue to read from the next buf
                            break 'outer;
                        }
                    }
                    if len > 0 {
                        if len > min_len {
                            lens.push(len as u32);
                            total += len;
                        }
                        len = 0;
                        is_new_record = true;
                    }
                }
                s1.send(buf).unwrap();
            }
            (lens, total)
        });
        stat.join().expect("Failed to read from input file!")
    })
    .unwrap();

    if out && !lens.is_empty() {
        lens.par_sort_unstable();
        if step_len > 0 {
            out_step(&lens, step_len as u32, genome_len);
        } else {
            out_stat(&lens, total, genome_len);
        }
    };
    (total, lens)
}

pub fn stat(
    infiles: &[&str],
    min_len: usize,
    genome_len: usize,
    n_len: usize,
    step_len: usize,
    out_ctg: bool,
) {
    if n_len == 0 {
        stat_read(infiles, min_len, genome_len, step_len, true);
        return;
    }
    let mut lens = Vec::with_capacity(1024);
    let mut total: usize = 0;

    let re = Regex::new(&format!("(?i)N{{{n_len},}}")).unwrap();
    let mut ctg_lens = Vec::with_capacity(1024);
    let mut ctg_total: usize = 0;

    let mut gap_lens = Vec::with_capacity(1024);
    let mut gap_total: usize = 0;

    for infile in infiles {
        let mut records = parse_fx(infile);
        let out = out_ctg.then(|| {
            let out = infile.to_string() + ".ctg.fa";
            File::create(&out).unwrap_or_else(|_| panic!("failed create file: {out}"))
        });
        while let Some(record) = records.iter_record().unwrap() {
            let len = record.len();
            if len < min_len {
                continue;
            }
            lens.push(len as u32);
            total += len;

            let mut last_pos = 0;
            let mut ctg_count = 1;
            let seq = record.seq();
            for mat in re.find_iter(seq) {
                if mat.start() > last_pos {
                    ctg_lens.push((mat.start() - last_pos) as u32);
                    ctg_total += mat.start() - last_pos;
                    if out_ctg {
                        writeln!(
                            out.as_ref().unwrap(),
                            ">{}_ctg{}\n{}",
                            record.head(),
                            ctg_count,
                            &seq[last_pos..mat.start()]
                        )
                        .unwrap();
                        ctg_count += 1;
                    }
                }
                gap_lens.push((mat.end() - mat.start()) as u32);
                gap_total += mat.end() - mat.start();
                last_pos = mat.end();
            }
            if len > last_pos {
                ctg_lens.push((len - last_pos) as u32);
                ctg_total += len - last_pos;
                if out_ctg {
                    writeln!(
                        out.as_ref().unwrap(),
                        ">{}_ctg{}\n{}",
                        record.head(),
                        ctg_count,
                        &seq[last_pos..len]
                    )
                    .unwrap();
                }
            }
        }
        lens.par_sort_unstable();
        ctg_lens.par_sort_unstable();
        gap_lens.par_sort_unstable();
        out_stats(
            &lens, total, &ctg_lens, ctg_total, &gap_lens, gap_total, genome_len,
        );
    }
}

pub fn sum_fx(infiles: &[&str]) -> usize {
    stat_read(infiles, 0, 0, 0, false).0
}
