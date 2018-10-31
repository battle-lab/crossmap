# Prerequisites
The following items are required to compute cross-mappability genomewide.
#### Genome sequences (Format: FASTA)
The genome sequences should be stored in a directory. The genome sequence for each chromosome should be in a fasta file in a directory. The file name for each chromosome should be CHRNAME.fa. For example, file names for human genome would look like below:
```
chr10.fa  chr14.fa  chr18.fa  chr21.fa  chr4.fa  chr8.fa  chrY.fa
chr11.fa  chr15.fa  chr19.fa  chr22.fa  chr5.fa  chr9.fa
chr12.fa  chr16.fa  chr1.fa   chr2.fa   chr6.fa  chrM.fa
chr13.fa  chr17.fa  chr20.fa  chr3.fa   chr7.fa  chrX.fa
```
Each fasta file should look like below:
```
>chr17
AAGCTTCTCACCCTGTTCCTGCATAGATAATTGCATGACAATTGCCTTGT
CCCTGCTGAATGTGCTCTGGGGTCTCTGGGGTCTCACCCACGACCAACTC
CCTGGGCCTGGCACCAGGGAGCTTAACAAACATCTGTCCAGCGAATACCT
GCATCCCTAGAAGTGAAGCCACCGCCCAAAGACACGCCCATGTCCAGCTT
AACCTGCATCCCTAGAAGTGAAGGCACCGCCCAAAGACACGCCCATGTCC
AGCTTATTCTGCCCAGTTCCTCTCCAGAAAGGCTGCATGGTTGACACACA
GTGcctgcgacaaagctgaatgctatcatttaaaaactccttgctggttt
gagaggcagaaaatgatatctcatagttgctttactttgcatattttAAA
ATTGTGACTTTCATGGCATAAATAATACTGGTTTATTACAGAAGCACTAG
```

#### Gene annotation file (Format: GTF)
A tab-delimited [Gencode gtf](https://www.gencodegenes.org/pages/data_format.html) file containing gene annotations is required. You may download comprehensive gene annotations for human and mouse from [Gencode website](https://www.gencodegenes.org/). The gene annotation file should look like below:
```
##description: evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)
##provider: GENCODE
##contact: gencode@sanger.ac.uk
##format: gtf
##date: 2013-12-05
chr1    HAVANA  gene    11869   14412   .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
chr1    HAVANA  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
```

#### Mappability of each k-mer (Format: BEDGRAPH or BIGWIG)
A tab-delimitted [bedgraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) file containing the mappability of each k-mer is required. The bedgraph file should have 4 columns: chr, start_pos, end_pos, mappability. The bedgraph file should look like below:
```
chr1	0	10	0.25
chr1	10	27	1
chr1	27	43	0.5
```
Here, the first 10 k-mers have a mappability of 0.25, next 17 k-mers have a mappability of 1, and the next 16 k-mers have a mappability of 0.5. Note: a track definition line is not expected.

bigwig to bedgraph: If you have k-mer mappabilities in a [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) file, you may easily convert it into a bedgraph file using the following shell script:
```
bigwig_fn="mappability.bigwig"  # filename (with path) of your bigwig file
bedgraph_fn="mappability.bed"   # filename (with path) of the bed file
bigWigToBedGraph  "$bedgraph_fn" "$bedgraph_fn"
```
Further instructions to convert bigwig files are available [here](https://genome.ucsc.edu/goldenpath/help/bigWig.html). Mappabilities of k-mers in human genome hg19 are available [here](http://bit.ly/hg19_mappability) in bigwig files.


#### Bowtie index
Bowtie index files are required to align k-mers to the genome. You may either download prebuilt bowtie indexes from the [bowtie website](http://bowtie-bio.sourceforge.net/index.shtml) or create an index using for your genome following instructions from the bowtie wesbite.


#### Software requirements (Linux, R, bowtie)
- Linux
- R
  - Please make sure the path variable includes the location for `Rscript`.
  - Please make sure the following R packages are installed: data.table, stringr, intervals, argparser, stats, 
- bowtie v1
  - Please make sure the path variable includes the location for `bowtie`.

#### Notes: 
- Input parsing in the program is pretty basic, so please try to match the format with the given examples.
- The program has been tested on CentOS Linux (64 bit) using R v3.5.1 and bowtie v1.2.2.
