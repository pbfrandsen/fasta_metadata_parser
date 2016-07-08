# FASTA metadata parser

This parser takes a FASTA file as input and calculates contig statistics (N50, L50, etc.).

It requires NumPy and [scikit-bio](http://scikit-bio.org). Both of which can be easily obtained from the [Anaconda Python](https://www.continuum.io/downloads) distribution.

####Usage:

To run the parser, simply add it to the directory with your FASTA file and issue the command:

```
$ python fasta_meta_data_parser.py <fasta_file_name>
```

It will then write a file called ```genome_stats.txt```, which contains the relevant statistics.