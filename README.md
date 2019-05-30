# RNAfiltering
RNAfilter is a Python tool that was designed to preprocess RNA sequences by trimming, filtering, and adapter removal.
The user has the option of specifying the criterial of RNA-seq analysis.

### How to use
RNAfilter was designed to run from the command line. Use the following code to process RNA sequences file using the standard criteria:
```bash
python RNAfilter.py -F input.fastq 
```

### Manual
RNAfilter tool offeres a set of parameters that can be used to tailer the processing step.

- `-F <path>`: takes the input file in fastq format. Note that it's the only required argument and the rest are optional.
- `-A <str>`: takes the adapter sequence to be trimmed from the beginning of the read. If no sequence is provided, no thing will be trimmed.
- `--read_qual <int>`: takes the minimum average read quality score to filter out the read. The default value is 20.
- `--base_qual <int>`: takes the minimum base quality score to trim the base from both ends. The default value is 15.
- `--read_len <float>`: takes the minimum accepted fraction of the trimmed read to be considered in the output. The default is that the read is filtered out if its length after trimming is less that 65% of its initial length.
- `-h|--help`: show the help message and exits. Don't use this argument while running the tool.

### Example usage
```bash
 python RNAfilter.py -F input.fastq --read_qual 20 --base_qual 15 --read_len 0.65 -A ACGT
```
