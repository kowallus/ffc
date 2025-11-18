# Fast FASTA Compressor (FFC)

[![GitHub downloads](https://img.shields.io/github/downloads/kowallus/ffc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/kowallus/ffc/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/fast-fasta-compressor.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/fast-fasta-compressor)

Fast FASTA Compressor (FFC) is a tool designed for compressing FASTA files.
It is focused on speed and for non-redundant (or mildly redundant) DNA sequences it basically packs fours of nucleotides into bytes.
If, however, data seem compressible, like for e.g. virus strains, it makes use of a fast zstd mode.
We point, however, that FFC accepts any input,
and for example its performance on a tar genomic archive is essentially not worse than on its raw (FASTA) contents
(and if the genomes within the tar are similar, and small enough to fit a block, may even be better due to exploited repetitions).
FFC compression and decompression speed is often around 1500-2000 MB/s on an SSD, and above 3000 MB/s in RAM memory.

Major features and properties:
* supports FASTA files as input and output,
* allows for ambiguous IUPAC codes, masked sequences, and has no limit on sequence length or the number of sequences,
* is fully lossless (even for input with irregular EOL locations),
* supports standard input and output during (de)compression, which makes it easy to integrate into pipelines,
* applies zstd backend compression for selected streams,
* allows to set the number of (de)compression threads; the output archive is the same no matter the number of chosen compression worker threads,
* works on moderately-sized blocks of data (4 MB per thread by default), keeping thus the memory usage low, regardless of the input/output size,
* supports higher memory usage in order to (hopefully) boost compression ratios.

FFC is available on Bioconda as [fast-fasta-compressor](https://anaconda.org/bioconda/fast-fasta-compressor).

### FFC format specification

To understand the FFC format, please consult its [specification](specification.pdf). We
also provide a Python [decompression script](example-scripts/ffc_decompressor.py), generated according to
the specification.

### Compile:
```shell
git clone https://github.com/kowallus/ffc.git
cd ffc
make -C zstd/lib ZSTD_LEGACY_SUPPORT=0 ZSTD_LIB_DEPRECATED=0 ZSTD_LIB_DICTBUILDER=0 libzstd.a
mkdir build
cd build
cmake ..
make
```

### Usage

FFC can compress or decompress.
The compression mode is the default one, so to run it type, e.g.:
```
ffc -i human.fa -o human.ffc
```
or simply:
```
ffc -i human.fa
```
(which will create human.fa.ffc, assuming that no such file exists).

To decompress, use `-d` or `--decompress`.
For example:
```
ffc -d -i human.fa.ffc -o human.fa
```
(assuming that human.fa does not exist).

To overwrite an existing output file, use `-f` or `--force`
So the former example can be modified to:
```
ffc -d -i human.fa.ffc -o human.fa -f
```

The order of parameters is arbitrary, so you can write, e.g.
```
ffc -d -f -o human.fa -i human.fa.ffc
```

The parameters for providing input and output also have their longer forms:
`-i` or `--input`, and
`-o` or `--output`

FFC handles also standard input/output.
Following POSIX convention, a single hyphen character can be used to specify input from or output to the standard input and output streams.
To read a file from the standard input, use e.g.:
```
ffc -i - -o comp.ffc
```
To decompress to the standard output:
```
ffc -d -i comp.ffc -o -
```


The backend compression is based on zstd and by default it detects if zstd compression is needed for the actual data (we call this mode "adaptive").
Alternatively, you can specify the compression level as
`-l` or `--level`
followed by an integer from the range \[0, 22\], where 0 is no zstd compression and 22 is max zstd compression (slow!).
We do not recommend a value greater than 4 in most cases.

FFC works on equal-length blocks (the last block of the input may be shorter). The default block size is $2^{22}$ bytes = 4 MiBytes. This parameter can be changed with
`-b` or `--block`
followed by an integer from the range \[20, 30\] (we point out that, for technical reasons, -b 30 corresponds to $2^{30} - 1$, not $2^{30}$ bytes!).
Using large blocks is beneficial for the compression ratio on redundant data (e.g., bacterial genome collections), but requires more memory and may be slower as it might limit the number of threads working in parallel.

By default, FFC uses 12 threads for compression and 4 threads for decompression. This parameter can be set with:
`-t` or `--threads`

So, for example, to compress a file with 12 threads, you can write:
```
ffc -i influenza.fa -f -t 12
```
(which will create influenza.fa.ffc, even if it previously existed).

If you are curious which FFC version you have, use:
`-v` or `--version`

Finally, the list of commands and parameters is printed with
`-h` or `--help`

Here's a rundown of the options:
```
  -h,          --help                         Print this help message and exit 
  -d,          --decompress                   Decompress mode 
  -f,          --force                        Overwrite the output file if exists 
  -i,          --input FNAME                  Input file 
  -o,          --output FNAME                 Output file 
  -l,          --level [0 - 22]               Backend compr. level, default: adaptive 
  -b,          --block [20 - 30]              Block size order, default: 22 
  -t,          --threads                      Number of threads, default: 12c / 4d 
  -v,          --version                      Show version information 
```
