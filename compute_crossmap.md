# How to compute cross-mappability genome-wide?
There are 5 steps to compute cross-mappability genomewide.
1. Settings specification
2. Process annotation data
3. Generate gene-mappability from k-mer mappability
4. Generate ambiguous k-mers
5. Compute cross-mappability


### 1. Settings specification
We will run the following shell scripts in terminal. First, we need to specify inputs in variables, which we will use in later steps. Please look at the [prerequisite page](https://github.com/battle-lab/crossmap/blob/master/prerequisites.md) and specify appropriate settings in the following variables.

```shell
gene_annot_fn=""  # gene annotation file name with full path
genome_dir=""     # genome directory with full path
bowtie_index_prefix=""   # bowtie index prefix name with full path
exon_k=75  # 75-mers from exon, change if you use different k
utr_k=36   # 36-mers from UTR, change if you use different k
exon_kmer_mappability_fn=""  # bedgraph file with k-mer mappability to use for exons
utr_kmer_mappability_fn=""  # bedgraph file with k-mer mappability to use for UTRs
```

We need to perform a lot of file operations to compute cross-mappability. In this tutorial, we will put all files under a directory which we call 'computation directory'. Please do not edit/delete any file inside the computation direcotry unless you know the consequences.
```shell
comp_dir=""  # put your result directory here
```


Next, we set some variables that control the computation memory and time. 
```shell
n_chr_in_ram=7       # number of chromosome to load in memory at a time
n_genes_per_batch=1000  # number of genes to process in a batch
dir_name_len=7       # length of the sub-directory names. First dir_len letters from gene
```
Note: 
- If large n_chr_in_ram is chosen, total memory requirement will be high.
- Memory is cleaned after computing cross-mappability of every n_genes_per_batch genes.
- If there are thousands of file the same directory, reading/writing of a file get slow. To avoid such a situation, we store files in sub-directories where the sub-directory names come from the first few letters of gene ids. dir_name_len is the number of letters to use from gene ids for sub-directory names.


### 2. Process annotation data
First we need to crate separate annotation files for exons and UTRs from the gene annotation file, and put them in a subdirectory ("annot") inside our computation directory.

```shell
# specification
annot_dir="$comp_dir/annot"
exon_annot_fn="$annot_dir/annot.exon.gtf"
utr_annot_fn="$annot_dir/annot.utr.gtf"
exon_utr_annot_fn="$annot_dir/annot.exon_utr.gtf"

if [ ! -d $annot_dir ]; then mkdir -p $annot_dir; fi

# split the annotation file
awk '{FS = "\t" } $3=="exon" {print}' $gene_annot_fn > $exon_annot_fnn.gtf
awk '{FS = "\t" } $3=="UTR" {print}' $gene_annot_fn > $exon_annot_fn
awk '{FS = "\t" } $3=="UTR" || $3=='exon' {print}' $gene_annot_fn > $exon_annot_fn
```


### 3. Generate gene-mappability from k-mer mappability
Next, we need to compute mappability of genes from mappability of k-mers. We will need genes with mappability < 1 only in later steps to compute cross-mappability. We have an R script (compute_mappability.R) to compute gene mappabilities with the following arguments: 

 - -exon: exon annotation file
 - -utr: UTR annotation file
 - -k_exon: value of k for exon
 - -k_utr: value of k for UTR
 - -kmap_exon: bedgraph file containing k-mer mappabilities where k=k_exon
 - -kmap_utr: bedgraph file containing k-mer mappabilities where k=k_utr
 - -o : output file containing gene mappabilities

Here is the script to compute gene mappability:
```shell
# specification
mappability_dir="$comp_dir/gene_mappability"
mappability_fn="$mappability_dir/gene_mappability.txt"

if [ ! -d $mappability_dir ] ; then mkdir -p $mappability_dir ; fi

# call 
Rscript compute_mappability.R -exon $exon_annot_fn \
                              -utr $utr_annot_fn \
                              -k_exon $exon_k \
                              -k_utr $utr_k \
                              -kmap_exon $exon_kmer_mappability_fn \
                              -kmap_utr $utr_kmer_mappability_fn \
                              -o $mappability_fn
```

### 4. Generate ambiguous k-mers
At this stage, we need to generate a list of ambiguous k-mers that map to multiple regions of the genome for genes with mappability < 1. We can do so using one of our scripts generate_ambiguous_kmers.R with the following arguments:
 - -mappability: gene mappability file
 - -genome: the genome directory
 - -exon: exon annotation file
 - -utr: UTR annotation file
 - -k_exon: value of k for exon
 - -k_utr: value of k for UTR
 - -kmap_exon: bedgraph file containing k-mer mappabilities where k=k_exon
 - -kmap_utr: bedgraph file containing k-mer mappabilities where k=k_utr
 - -th1: starting gene mappability threshold. ambiguous k-mers created only from genes with mappability >= th1.
 - -th2: ending gene mappability threshold. ambiguous k-mers created only from genes with mappability < th2.
 - -dir_name_len: length of the sub-directory names. First dir_len letters from gene ids are used as sub-directory names.
 - -o : output folder where ambiguous k-mers and their frequencies will be stored

Note: You may choose to generate k-mers in multiple batches by running the above script with different mappability thresholds. But our experience suggests that it is OK to generate all k-mers in one batch with th1=0 and th2=1.

Here is the code to gene ambiguous k-mers.
```shell
ambiguous_kmer_dir="$comp_dir/ambiguous_kmers"
mappability_th1=0
mappability_th2=1

if [ ! -d $ambiguous_kmer_dir ] ; then mkdir -p $ambiguous_kmer_dir ; fi

Rscript generate_ambiguous_kmers.R  -mappability $mappability_fn \
                                    -genome $genome_dir \
                                    -exon $exon_annot_fn \
                                    -utr $utr_annot_fn \
                                    -k_exon $exon_k \
                                    -k_utr $utr_k \
                                    -kmap_exon $exon_kmer_mappability_fn \
                                    -kmap_utr $utr_kmer_mappability_fn \
                                    -th1 $mappability_th1 \
                                    -th2 $mappability_th2 \
                                    -dir_name_len $dir_name_len
                                    -o $ambiguous_kmer_dir
```
After running the above script, you would find two files for each gene with mappability<1 in a subdirectory under the directory for ambiguous kmers ($ambiguous_kmer_dir). One file (*.txt) contains the ambiguous k-mers, and another file (*.count.txt) contains the corresponding frequency of those k-mers.

We need ambiguous k-mers in a fasta file to align them to the genome using bowtie. Here is the code for this purpose.
```shell
awk '{i += 1 ; print ">r"i ; print}' < 75mers_unsorted.txt > 75mers_unsorted.fa
```
After running the above script, you will find fasta files (\*.fa) for each gene with mappability<1. The text files for ambiguous k-mers are removed.




### 5. Compute cross-mappability
Finally, we compute cross-mappability from every gene with mappability<1 to other genes by aligning the ambiguous k-mers to the genome using bowtie v1. We can do so using our script compute_cross_mappability.R with the following arguments:

 - -exon: exon annotation file
 - -utr: UTR annotation file
 - -mappability: gene mappability file
 - -kmer: ambigous k-mer directory
 - -align: directory where alignments will be temporarily stored
 - -index: bowtie index prefix
 - -th1: starting gene mappability threshold. cross-mappability is computed only from genes with mappability >= th1.
 - -th2: ending gene mappability threshold. cross-mappability is computed only from genes with mappability < th1.
 - -o: output directory where results and intermediate files will be stored

This part takes a long time. So, in our script, we kept option to compute cross-mappabilities in multiple batches, by spliting genes according to their mappabilities. 
```shell
alignment_dir="$comp_dir/ambiguous_kmers_alignment"
mappabiity_range=0.05

if [ ! -d $alignment_dir ] ; then mkdir -p $alignment_dir ; fi

for mappability_th1 in $(seq 0 $mappabiity_range 1); 
do
  mappability_th2=$(echo "$mappability_th1 + $mappabiity_range" | bc -l)
  if [[ $mappability_th1 = 1 ]]; then continue; fi
  Rscript compute_cross_mappability.R -exon $exon_annot_fn \
                                    -utr $utr_annot_fn \
                                    -mappability $mappability_fn \
                                    -kmer $ambiguous_kmer_dir \
                                    -align $alignment_dir \
                                    -index $bowtie_index_prefix \
                                    -th1 $mappability_th1 \
                                    -th2 $mappability_th2 \
                                    -p $n_threads \
                                    -o $cross_mappability_dir

done
```
Note: If possible (i.e., if there is enough memory), we would like to run cross-mappabilities for every batch (grouped by gene mappability) in parallel.

Once cross-mappabilities from every gene have been computed, we need to combine them.
```shell
cd $out_dir,'
rm mappability_conflicts.txt; 
for fn in ENSG*/ENSG*.txt
do 
  cat $fn >> mappability_conflicts.txt
done
```
