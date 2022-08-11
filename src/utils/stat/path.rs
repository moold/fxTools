use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{stdin, BufRead, BufReader, Cursor, Read},
    path::Path,
};

//convert path to Read, accept auguments: file.txt/file.gz/file.fofn
pub fn open_path<T>(path: T) -> Vec<Box<dyn Read>>
where
    T: AsRef<str>,
{
    let mut readers = Vec::new();
    let path = path.as_ref();
    let mut reader: Box<dyn Read> = match path {
        "-" => {
            if atty::is(atty::Stream::Stdin) {
                panic!("Missing input from Stdin");
            }
            Box::new(stdin())
        }
        path => {
            Box::new(File::open(path).unwrap_or_else(|_| panic!("Failed open file {:?}", path)))
        }
    };

    let mut format_bytes = [0u8; 2];
    reader
        .read_exact(&mut format_bytes)
        .unwrap_or_else(|_| panic!("Failed read file {:?}", path));
    reader = Box::new(Cursor::new(format_bytes.to_vec()).chain(reader));
    if &format_bytes[..2] == b"\x1f\x8b" {
        // for gz foramt
        reader = Box::new(MultiGzDecoder::new(reader));
        reader.read_exact(&mut format_bytes).unwrap();
        reader = Box::new(Cursor::new(format_bytes.to_vec()).chain(reader));
    }
    match format_bytes[0] {
        b'@' | b'>' => {
            readers.push(reader);
        }
        _ => {
            // for a fofn file
            let reader = BufReader::with_capacity(65536, reader);
            let parent = Path::new(path).parent().unwrap_or_else(|| Path::new(""));
            for _line in reader.lines().map(|l| l.unwrap()) {
                let line = _line.trim();
                if line.starts_with('#') || line.is_empty() {
                    continue;
                }
                let _path = parent.join(line); // convert to a absolute path
                readers.extend(open_path(_path.to_string_lossy()));
            }
        }
    }
    readers
}
