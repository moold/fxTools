# fxTools
A toolset for processing sequences in FASTA/Q formats

## Installation

#### Dependencies

`fxTools` is written in rust, see [here](https://www.rust-lang.org/tools/install) to install `Rust` first.

#### Download and install

```sh
git clone git@github.com:moold/fxTools.git
cd fxTools && cargo build --release
# ./target/release/fxTools
```

## Usage and Parameters
```
fxTools
A toolset for processing sequences in FASTA/Q formats

USAGE:
    fxTools [SUBCOMMAND] [OPTIONS] <input>...

ARGS:
    <input>...    input file ...

OPTIONS:
    -a, --attr <STR>                   get sequence attributes, id:len:x:qs
                                         len: sequences length
                                         x:   count x, case sensitive, where x can be a single base
                                              or multiple bases (gcGC means count g + c + G + C)
                                         qs:  quality score, for FASTQ
    -r, --reform <STR>                 reform/modify the sequence, accept values:
                                         lower:    convert sequences into lowercase
                                         upper:    convert sequences into uppercase
                                         fq2fa:    converts FASTQ to FASTA
                                         fa2fq:    converts FASTA to FASTQ
                                         rev:      reverse the sequence
                                         com:      complement the sequence
                                         rc:       reverse and complement the sequence
                                         lineINT[c]:
                                                   wrap sequences into INT characters per line,
                                                   0 for no wrap, c is the alignment mode,
                                                   color for SNVs and INDELs
                                         linkINT:  link sequences with INT Ns
    -p, --split <INT>                  split file with INT subfiles in total
    -s, --sample <int[G|M|K]|float>    subsample reads, int is the total bases, float is the
                                       fraction
    -h, --help                         Print help information
    -V, --version                      Print version information

SUBCOMMANDS:
    stat       simple statistics of FASTA/Q files
    findseq    find subseq positions
    findgap    find gap(Nn) regions
    getseq     get sequences or subsequences from a region or file
    diff       compare sequences between two files
```