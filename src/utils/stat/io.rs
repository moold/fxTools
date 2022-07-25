use memchr::{
    memchr,
    // memchr3,
    memchr_iter,
};
use std::io::{self, ErrorKind};

const NEWLINE: u8 = b'\n';
// const GAPUPPER: u8 = b'N';
// const GAPLOWER: u8 = b'n';
const DEFAULT_CAPACITY: usize = 1024 * 100;

pub struct Buffer {
    buf: Vec<u8>,
    buf_len: usize, //valid buf data length
    pos: usize,     //pos of start to next serach
}

impl Buffer {
    pub fn new() -> Self {
        Self {
            buf: vec![0; DEFAULT_CAPACITY],
            buf_len: 0,
            pos: 0,
        }
    }

    pub fn fill<R: io::Read>(&mut self, r: &mut R) -> io::Result<usize> {
        self.pos = 0;
        self.buf_len = 0;
        let mut s = 0;
        loop {
            match r.read(&mut self.buf[s..]) {
                Ok(0) => return Ok(0),
                Ok(n) => {
                    s += n;
                    self.buf_len += n;
                    if s >= self.buf.capacity() {
                        return Ok(n);
                    }
                }
                Err(ref e) if e.kind() == ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.buf_len == 0
    }

    pub fn next_byte(&mut self, ignore_newline: bool) -> Option<u8> {
        if ignore_newline {
            self.skip_newlines();
        }

        if self.pos < self.buf_len {
            Some(self.buf[self.pos])
        } else {
            None
        }
    }

    pub fn next_line_len(&mut self) -> Option<(usize, bool)> {
        if self.pos < self.buf_len {
            match memchr(NEWLINE, &self.buf[self.pos..self.buf_len]) {
                Some(p) => {
                    self.pos += p + 1;
                    return Some((p, true));
                }
                None => {
                    let p = self.buf_len - self.pos;
                    self.pos = self.buf_len;
                    return Some((p, false));
                }
            }
        }

        None
    }

    pub fn skip_newlines(&mut self) -> usize {
        let mut skip_len = 0;
        while self.pos < self.buf_len && self.buf[self.pos] == NEWLINE {
            self.pos += 1;
            skip_len += 1;
        }
        skip_len
    }

    pub fn skip_lines(&mut self, len: usize) -> usize {
        let mut skip_line = 0;
        while let Some((_p, true)) = self.next_line_len() {
            skip_line += 1;
            if skip_line >= len {
                break;
            }
        }
        skip_line
    }

    pub fn skip_bases(&mut self, len: usize) -> usize {
        let mut skip_len = 0;
        if len < 10 {
            while let Some(b) = self.next_byte(false) {
                self.pos += 1;
                if b != NEWLINE {
                    skip_len += 1;
                    if skip_len >= len {
                        break;
                    }
                }
            }
        } else if self.pos + len - 1 < self.buf_len {
            let sep_count = memchr_iter(NEWLINE, &self.buf[self.pos..self.pos + len]).count();
            self.pos += len;
            skip_len = if sep_count > 0 {
                len - sep_count + self.skip_bases(sep_count)
            } else {
                len
            };
        } else {
            let sep_count = memchr_iter(NEWLINE, &self.buf[self.pos..self.buf_len]).count();
            skip_len = self.buf_len - self.pos - sep_count;
            self.pos = self.buf_len;
        }

        self.skip_newlines();
        skip_len
    }

    // pub fn skip_ns(&mut self) -> usize {
    //     let mut skip_n = 0;
    //     while self.buf[self.pos] == GAPUPPER || self.buf[self.pos] == GAPLOWER {
    //         self.pos += 1;
    //         skip_n += 1;
    //         if self.pos >= self.buf_len {
    //             match self.fill(true){
    //                 Ok(n) => if n == 0 {
    //                     break;
    //                 },
    //                 Err(_e) => break
    //             }
    //         }
    //     }
    //     skip_n
    // }
}
