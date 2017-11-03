# FASTA metadata parser

This parser takes a FASTA file as input and calculates both scaffold and contig statistics (N50, L50, etc.) from a scaffold FASTA file. It does this by breaking each scaffold wherever there is more than one `N` and then calculating statistics for both the scaffolds and contigs.

It requires NumPy and [scikit-bio](http://scikit-bio.org). Both of which can be easily obtained from the [Anaconda Python](https://www.continuum.io/downloads) distribution.

####Usage:

To run the parser, simply add it to the directory with your FASTA file and issue the command:

```
$ python fasta_meta_data_parser.py <fasta_file_name>
```

It will then write the stats to the screen. It will look something like this:

```
Contig statistics:
Total number of base pairs: 1106252391
Total number of contigs: 8104
N10: 2049786
N20: 1376365
N30: 1058529
N40: 829321
N50: 647974
L10: 43
L20: 112
L30: 204
L40: 321
L50: 472
GC content: 42.88%
Median contig size: 21973.5
Mean contig size: 136506.96
Longest contig is: 4781189.0
Shortest contig is: 1231.0

Scaffold statistics:
Total number of base pairs: 1106317551
Total number of scaffolds: 7918
N10: 2170181
N20: 1496579
N30: 1109545
N40: 872285
N50: 685849
L10: 41
L20: 105
L30: 192
L40: 304
L50: 447
GC content: 42.88%
Median scaffold size: 21818.5
Mean scaffold size: 139721.84
Longest scaffold is: 4781189.0
Shortest scaffold is: 1231.0
```