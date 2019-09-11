# transcriptomics-useful
Useful scripts and code snippets for bioinformatics and transcriptomics


## Python Scripts:
### filterReadsBySequenceFasta.py

Python script to filter (in or out) reads in FASTA format containing complete or partial sequences

Options:
```
-p path to input FASTA file.
-s input sequence for matching
-o path to output directory
-f filter read "in" or "out" based on sequence
-m number of allowed mismatches with input sequence (use 0 for exact matches only)
-n maximum number of reads desired in destination file
```
Usage example:\
`python filterReadsBySequenceFasta.py -p input.fasta -s filter_sequence -o output_path -f in -m 0 -n 10000`

### makeRandomSequences.py

Python script to generate user-defined number of random DNA sequences of input length

Options:
```
-n number of desired sequences
-l length of desired sequences
```

Usage example:\
`python makeRandomSequences.py -n 100 -l 10`
