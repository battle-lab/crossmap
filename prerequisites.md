### Input requirements
The following items are required to compute cross-mappability genomewide.
- Genome sequences (one fasta file for each chromosome)
- Gene annotation files (Format: GTF)
- Mappability of each k-mer (Format: bedgraph or bigWig)
- Bowtie index for the genome

### Software requirements
- R
- bowtie v1
- Linux

#### Genome sequences

#### Gene annotation files

#### Mappability of each k-mer
A tab-delimitted [bedgraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) file containing the mappability of each k-mer is required. The bed file should have 4 columns: chr, start_pos, end_pos, mappability. The bed file should look like below:
```
chr1	0	10	0.25
chr1	10	27	1
chr1	27	43	0.5
```
Here, the first 10 k-mers have a mappability of 0.25, next 17 k-mers have a mappability of 1, and the next 16 k-mers have a mappability of 0.5. Note: a track definition line is not expected.

#### bigwig to bedgraph
If you have k-mer mappabilities in a [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) file, you may easily convert it into a bedgraph file using the following shell script:
```
bigwig_fn="mappability.bigwig"  # filename (with path) of your bigwig file
bedgraph_fn="mappability.bed"   # filename (with path) of the bed file
bigWigToBedGraph  "$bedgraph_fn" "$bedgraph_fn"
```
Further instructions to convert bigwig files are available [here](https://genome.ucsc.edu/goldenpath/help/bigWig.html). Mappabilities of k-mers in human genome hg19 are available [here](http://bit.ly/hg19_mappability) in bigwig files.


#### Bowtie index


