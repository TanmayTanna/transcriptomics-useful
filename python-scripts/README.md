# Python Scripts:

## SE_rmdup.py

Python script to remove duplicate alignments from SAM file for SE sequencing. Dupicate alignments are defined as having the same alignment start position and same length of aligning sequence.

Usage:\
`python SE_rmdup.py -i input.sam -o output.sam`

### filterReadsBySequenceFasta.py

Python script to filter (in or out) reads in FASTA format containing complete or partial sequences

Options:
```
--inFile	path to input FASTA file.
--seq 		input sequence for matching
--outPath 	path to output directory
--filter 	filter read "in" or "out" based on sequence
--mismatch			number of allowed mismatches with input sequence (use 0 for exact matches only)
--readNumber	maximum number of reads desired in destination file
```
Usage example:\
`python filterReadsBySequenceFasta.py --inFile input.fasta --seq filter_sequence --outPath output_path --filter in --mismatch 0 --readnumber 10000`

### makeRandomSequences.py

Python script to generate user-defined number of random DNA sequences of input length

Options:
```
-n number of desired sequences
-l length of desired sequences
```

Usage example:\
`python makeRandomSequences.py -n 100 -l 10`
