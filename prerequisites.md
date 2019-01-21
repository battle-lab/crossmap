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

You may download precomputed bedgraph (or bigwig) files for human genomes for certain settings from one the following links: 1) [hg19, k=2, mismatch=2](http://bit.ly/hg19_mappability), 2) [hg19, different k's and mismatches](https://figshare.com/articles/Cross_Mappability_hg19_gencode19/7315049).

If the k-mer mappability bed file for your reference genome with desired k and number of mismatches is not available in any of the above links, you may generate it by following instructions from [here](https://wiki.bits.vib.be/index.php/Create_a_mappability_track). For convenience, we summarized necessary instructions here. Before proceeding, please get the binary files of the [GEM Library](https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/) and [UCSC command-line tools (wigToBigWig and bigWigToBedGraph)](http://hgdownload.soe.ucsc.edu/admin/exe/). Please add the paths in the "$PATH" variable, and execute the following script with your settings.

```
export PATH=/your/GEM_library/bin:$PATH   # put the path of GEM Libary binaries
export PATH=/your/ucsc_binary:$PATH       # put the path of UCSC tools binaries


genome_fasta="your/reference/genome/ref_all_chr.fa" # put the path of the reference fasta file with all chromosomes.
k=75                                                # put the value of k
n_mismatch=2                                        # put the maximum number of mismatches allowed
n_threads=16                                        # put the number of threads
output_dir="your/output/directory"                  # put the output directory
gem_index_pref="ref_gem_index"                      # put the prefix of GEM index
gem_mappability_pref="mappability_75mer_2mismatch"  # put the prefix of mappability files

# index the refernce genome (execute once for one fasta file): -- approx time: ~40 min
gem-indexer -T $n_threads -c dna -i "$genome_fasta" -o "$output_dir/$gem_index_pref"   # ~40 min
# compute mappability -- approx time: ~6 hr
gem-mappability -m "$n_mismatch" -T $n_threads -I "$output_dir/$gem_index_pref.gem" -l $k -o "$output_dir/$gem_mappability_pref"
# convert mappability file to bed, step by step
gem-2-wig -I "$output_dir/$gem_index_pref.gem" -i "$output_dir/$gem_mappability_pref.mappability" -o "$output_dir/$gem_mappability_pref"
wigToBigWig "$output_dir/$gem_mappability_pref.wig" "$output_dir/$gem_mappability_pref.sizes" "$output_dir/$gem_mappability_pref.bigWig"
bigWigToBedGraph "$output_dir/$gem_mappability_pref.bigWig" "$output_dir/$gem_mappability_pref.bed"

```
The bedgraph file will be available in your output directory with .bed extension.


#### Bowtie index
Bowtie index files are required to align k-mers to the genome. You may either download prebuilt bowtie indexes from the [bowtie website](http://bowtie-bio.sourceforge.net/index.shtml) or create an index using for your genome following instructions from the bowtie wesbite.


#### Software requirements (Linux, R, bowtie)
- Linux
- R
  - Please make sure the path variable includes the location for `Rscript`.
  - Please make sure the following R packages are installed: data.table, intervals, argparser, stats, 
- bowtie v1
  - Please make sure the path variable includes the location for `bowtie`.

#### Notes: 
- Input parsing in the program is pretty basic, so please try to match the format with the given examples.
- The program has been tested on CentOS Linux (64 bit) using R v3.5.1 and bowtie v1.2.2.
