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
# gene annotation file name with full path
gene_annot_fn="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gtf" 
# genome directory with full path
genome_dir="/work-zfs/abattle4/lab_data/hg19"
# bowtie index prefix name with full path
bowtie_index_prefix="/work-zfs/abattle4/lab_data/hg19/prebuilt_bowtie_index/hg19"
# k for exon and utr, change if you use different k
exon_k=75
utr_k=36
# bedgraph files containing k-mer mappabilities with appropriate k
exon_kmer_mappability_fn="/work-zfs/abattle4/lab_data/annotation/kmer_alignability_hg19/wgEncodeCrgMapabilityAlign75mer.bed"
utr_kmer_mappability_fn="/work-zfs/abattle4/lab_data/annotation/kmer_alignability_hg19/wgEncodeCrgMapabilityAlign36mer.bed"
# maximum number of mismatches for an alignment
mismatch=2 
```

We need to call the scripts in our program (in `$prog_dir`) that perform a lot of file operations. In this pipeline, we will put all files under a directory which we call 'computation directory' (`$comp_dir`). Please do not edit/delete any file inside the computation direcotry unless you know the consequences.
```shell
# put the directory of this program here
prog_dir="/work-zfs/abattle4/ashis/prog/crossmap"
# put your computation directory here
comp_dir="/work-zfs/abattle4/ashis/progres/crossmap/test_hg19"
```


Next, we set some variables that control the computation memory, time, and status.
```shell
# maximum number of chromosomes to load in memory at a time
max_chr=7
# maximum number of genes to align before cleaning alignments
max_gene_alignment=200
# length of the sub-directory names. First dir_len letters from gene
dir_name_len=12
# verbose output? 1 for verbose, 0 for non-verbose
verbose=1
```
Note: 
- You may use a big number for `max_chr` if you can afford high memory.
- You may use a big number for `max_gene_alignment` if you have large disk space. Disk space is cleaned after computing cross-mappability from each of `max_gene_alignment` genes to other genes.
- If there are thousands of file the same directory, reading/writing of a file get slow. To avoid such a situation, we store files in sub-directories where the sub-directory names come from the first few letters of gene ids. dir_name_len is the number of letters to use from gene ids for sub-directory names.
- If `verbose` is 1, the program shows the current progress status.

Now, we create sub-directories in the computation directory to store scripts and logs.
```shell
# script and log directories
script_dir="$comp_dir/script"
log_dir="$comp_dir/log"
if [ ! -d $script_dir ]; then mkdir -p $script_dir; fi
if [ ! -d $log_dir ]; then mkdir -p $log_dir; fi
```

As the scripts in this program are likely to take relatively high memory and long time to run, you may find it easy to run them in computing clusters compared to run them in a stand-alone computer. Here, we will use [slurm](https://slurm.schedmd.com/overview.html).

```shell
# slurm utility
use_slurm=1  # use slurm (1: use, 0: don't use)

slurm_partition="shared"
source $prog_dir/slurm_util.sh
run_script()
{
  # if slurm is not used
  # this function expects only one argument: 
  # the script file.
  
  # if slurm is used,
  # this function expects 6 arguments: 
  # 1) script file, 2) partition, 3) no. of nodes, 
  # 4) no. of tasks, 5) time to run, 6) memory
  
  if [ $# -lt 1 ]; then 
    echo "no script to run"; 
    return 1; 
  fi
  script_fn=$1
  
  if [ $use_slurm -eq 1 ]; then 
    submit_slurm_job $@
  else
    sh $script_fn
  fi
  
  return 0
}
```

Note: 
- If you do not use slurm, then just change the value of `use_slurm` to 0: `use_slurm=0`, and ignore other settings.
- If you want to use a different type of computing cluster, please edit the `run_script()` function accordingly.
- If the slurm jobs cannot finish due to memory or time limit, please give appropriate resources when you call the `run_script()` function.


### 2. Process annotation data
First we need to process the gene annotation file. We will save the exon and UTR annotations in a tabular format in a subdirectory (`annot`) inside our computation directory.

```shell
# specification
annot_dir="$comp_dir/annot"
features="exon,UTR"
exon_utr_annot_fn="$annot_dir/annot.exon_utr.txt"
script_fn="${script_dir}/gtf_to_txt.sh"
log_fn="${log_dir}/gtf_to_txt.log"

if [ ! -d $annot_dir ]; then mkdir -p $annot_dir; fi

echo "Rscript \"$prog_dir/gtf_to_txt.R\" -gtf \"$gene_annot_fn\" \
                     -f \"$features\" \
                     -o \"$exon_utr_annot_fn\" \
                     2>&1 | tee \"$log_fn\"" > "$script_fn"
run_script "$script_fn" $slurm_partition 1 1 "0:10:0" "8GB"  
# Time: 3 min
```

Once the above script is completed, you should see a file in the `annot` subdirectory which should look like below:
```
gene_id chr     annotation_source       feature start_pos       end_pos strand  gene_name       gene_type
ENSG00000223972.4       chr1    HAVANA  exon    11869   12227   +       DDX11L1 pseudogene
ENSG00000223972.4       chr1    HAVANA  exon    12613   12721   +       DDX11L1 pseudogene
ENSG00000223972.4       chr1    HAVANA  exon    13221   14409   +       DDX11L1 pseudogene
```

### 3. Generate gene-mappability from k-mer mappability
Next, we need to compute mappability of genes from mappability of k-mers. We will need genes with mappability < 1 only in later steps to compute cross-mappability. We have an R script (`compute_mappability.R`) to compute gene mappabilities with the following arguments: 

 - -annot: exon and UTR annotation file (txt)
 - -k_exon: k-mer length for exon
 - -k_utr: k-mer length for UTR
 - -kmap_exon: bedgraph file containing k-mer mappabilities where k=k_exon
 - -kmap_utr: bedgraph file containing k-mer mappabilities where k=k_utr
 - -verbose: show computation status if verbose > 0
 - -o : output file containing gene mappabilities

Here is the script to compute gene mappability:
```shell
# specification
mappability_dir="$comp_dir/gene_mappability"
mappability_fn="$mappability_dir/gene_mappability.txt"
script_fn="${script_dir}/compute_mappability.sh"
log_fn="${log_dir}/compute_mappability.log"

if [ ! -d $mappability_dir ] ; then mkdir -p $mappability_dir ; fi

# create script and run it 
echo "Rscript \"$prog_dir/compute_mappability.R\" -annot \"$exon_utr_annot_fn\" \
                              -k_exon $exon_k \
                              -k_utr $utr_k \
                              -kmap_exon \"$exon_kmer_mappability_fn\" \
                              -kmap_utr \"$utr_kmer_mappability_fn\" \
                              -verbose $verbose \
                              -o \"$mappability_fn\" \
                              2>&1 | tee \"$log_fn\"" > "$script_fn"
run_script "$script_fn" $slurm_partition 1 1 "0:40:0" "20GB"  
# Memory: 20GB, Time: 26 min
```

### 4. Generate ambiguous k-mers
At this stage, we need to generate a list of ambiguous k-mers that map to multiple regions of the genome for genes with mappability < 1. We can do so using one of our scripts `generate_ambiguous_kmers.R` with the following arguments:
 - -mappability: gene mappability file
 - -genome: the genome directory
 - -annot: exon and UTR annotation file (txt)
 - -k_exon: k-mer length for exon
 - -k_utr: k-mer length for UTR
 - -kmap_exon: bedgraph file containing k-mer mappabilities where k=k_exon
 - -kmap_utr: bedgraph file containing k-mer mappabilities where k=k_utr
 - -th1: starting gene mappability threshold. ambiguous k-mers created only from genes with mappability >= th1.
 - -th2: ending gene mappability threshold. ambiguous k-mers created only from genes with mappability < th2.
 - -dir_name_len: length of the sub-directory names. First `dir_name_len` letters from gene form the subdirectory name.
 - -verbose: show computation status if verbose > 0
 - -o : output folder where ambiguous k-mers and their frequencies will be stored.

Note: You may choose to generate k-mers in multiple batches by running the above script with different mappability thresholds. But our experience suggests that it is OK to generate all k-mers in one batch with `th1=0` and `th2=1`.

Here is the code to gene ambiguous k-mers.
```
# specification
ambiguous_kmer_dir="$comp_dir/ambiguous_kmers"
mappability_th1=0
mappability_th2=1
script_fn="${script_dir}/generate_ambiguous_kmers.sh"
log_fn="${log_dir}/generate_ambiguous_kmers.log"

if [ ! -d $ambiguous_kmer_dir ] ; then mkdir -p $ambiguous_kmer_dir ; fi

# create script and run it
echo "Rscript \"$prog_dir/generate_ambiguous_kmers.R\"  -mappability \"$mappability_fn\" \
                                    -genome \"$genome_dir\" \
                                    -annot \"$exon_utr_annot_fn\" \
                                    -k_exon $exon_k \
                                    -k_utr $utr_k \
                                    -kmap_exon \"$exon_kmer_mappability_fn\" \
                                    -kmap_utr \"$utr_kmer_mappability_fn\" \
                                    -th1 $mappability_th1 \
                                    -th2 $mappability_th2 \
                                    -dir_name_len $dir_name_len \
                                    -verbose $verbose \
                                    -o \"$ambiguous_kmer_dir\" \
                                    2>&1 | tee \"$log_fn\"" > "$script_fn"
run_script "$script_fn" $slurm_partition 1 1 "1:0:0" "20GB" 
# Time: 40 min
```
After running the above script, you would find two files for each gene with mappability<1 in a subdirectory under the directory for ambiguous kmers ($ambiguous_kmer_dir). One file (\*.kmer.txt) contains the ambiguous k-mers, and another file (\*.count.txt) contains the corresponding frequency of those k-mers.

We need ambiguous k-mers in a fasta file to align them to the genome using bowtie. Here is the code for this purpose.
```shell
# create k-mer fasta file
for fn in "$ambiguous_kmer_dir"/*/*.kmer.txt
do
  fasta_fn=$(echo $fn | sed 's/.txt$/.fa/g')
  awk -v i=-1 '{i += 1 ; print ">"i ; print}' < $fn > $fasta_fn
done
# Time: 10 min
```
After running the above script, you will find fasta files (\*.kmer.fa) for each gene with mappability<1. 

### 5. Compute cross-mappability
Finally, we compute cross-mappability from every gene with mappability<1 to other genes by aligning the ambiguous k-mers to the genome using bowtie v1. We can do so using our script `compute_cross_mappability.R` with the following arguments:

 - -annot: exon and UTR annotation file (txt).
 - -mappability: gene mappability file.
 - -kmer: ambigous k-mer directory.
 - -align: directory where alignments will be temporarily stored.
 - -index: bowtie index prefix.
 - -n1: compute cross-mappability from n1-th to n2-th gene in the gene mappabiity file. negative number to use all genes.
 - -n2: compute cross-mappability from n1-th to n2-th gene in the mappabiity file.  negative number to use all genes.
 - -mismatch: maximum number of mismatch for alignment.
 - -max_chr: maximum number of chromosomes to load in memory.
 - -max_gene: maximum number of genes to align before cleaning alignments.
 - -initonly: initialize required resources only, without computing cross-mappability.
 - -dir_name_len: length of the sub-directory names in ambiguous kmer directory.
 - -verbose: show computation status if verbose > 0.
 - -o: output directory where results and intermediate files will be stored.
 

This part takes a long time. So, in our script, we kept option to compute cross-mappabilities in multiple batches, by spliting genes according to position in the gene mappability file (see argument `-n1` and `-n2`). However, we need to initialize necessary resources before actually computing cross-mappabilites. We can do so by setting `-initonly` to `TRUE`.

```shell
# specification
alignment_dir="$comp_dir/ambiguous_kmers_alignment"
cross_mappability_dir="$comp_dir/cross_mappability"
n_gene_per_crossmap_batch=2000

if [ ! -d "$alignment_dir" ] ; then mkdir -p "$alignment_dir" ; fi
if [ ! -d "$cross_mappability_dir" ] ; then mkdir -p "$cross_mappability_dir" ; fi

# initialize resources to compute cross-mappability (-initonly TRUE)
script_fn="${script_dir}/compute_cross_mappability_1_init.sh"
log_fn="${log_dir}/compute_cross_mappability_1_init.log"
echo "Rscript \"$prog_dir/compute_cross_mappability.R\" -annot \"$exon_utr_annot_fn\" \
                                    -mappability \"$mappability_fn\" \
                                    -kmer \"$ambiguous_kmer_dir\" \
                                    -align \"$alignment_dir\" \
                                    -index \"$bowtie_index_prefix\" \
                                    -n1 1 \
                                    -n2 $n_gene_per_crossmap_batch \
                                    -mismatch $mismatch \
                                    -max_chr $max_chr \
                                    -max_gene $max_gene_alignment \
                                    -initonly TRUE \
                                    -dir_name_len $dir_name_len \
                                    -verbose $verbose \
                                    -o \"$cross_mappability_dir\" \
                     2>&1 | tee \"$log_fn\"" > "$script_fn"
run_script "$script_fn" $slurm_partition 1 1 "1:30:0" "12GB"
# Memory: 6GB, Time: 45 min
```

Once the resources are initialized, you will see one file (\*.RData) for each chromosome in the `pos2gene` folder of `$cross_mappability_dir`. Then we run scripts to actually compute cross-mappabilities in the following way:
```shell
# actually compute cross-mappability (-initonly FALSE)
n_gene_in_mappability_file=$(wc -l $mappability_fn | sed 's/ .*//g')
for n1 in $(seq 1 $n_gene_per_crossmap_batch $n_gene_in_mappability_file)
do
  n2=$(($n1+$n_gene_per_crossmap_batch-1))
  script_fn="${script_dir}/compute_cross_mappability_2_${n1}_${n2}.sh"
  log_fn="${log_dir}/compute_cross_mappability_2_${n1}_${n2}.log"
  echo "Rscript \"$prog_dir/compute_cross_mappability.R\" -annot \"$exon_utr_annot_fn\" \
                                      -mappability \"$mappability_fn\" \
                                      -kmer \"$ambiguous_kmer_dir\" \
                                      -align \"$alignment_dir\" \
                                      -index \"$bowtie_index_prefix\" \
                                      -n1 $n1 \
                                      -n2 $n2 \
                                      -mismatch $mismatch \
                                      -max_chr $max_chr \
                                      -max_gene $max_gene_alignment \
                                      -initonly FALSE \
                                      -dir_name_len $dir_name_len \
                                      -verbose $verbose \
                                      -o \"$cross_mappability_dir\" \
                                      2>&1 | tee \"$log_fn\"" > "$script_fn"
  run_script "$script_fn" $slurm_partition 1 1 "48:0:0" "20GB"
  # Memory: 14GB, Time: 35 hours
done
```

Finally, once cross-mappabilities from every gene have been computed, we combine them using the following script.
```shell
# combine all cross-mappability results into one file
combined_cross_mappability_fn="$comp_dir/cross_mappability.txt"
if [ -f $combined_cross_mappability_fn ]; then rm "$combined_cross_mappability_fn";  fi
for fn in "$cross_mappability_dir"/*/*.crossmap.txt
do 
  cat $fn >> "$combined_cross_mappability_fn"
done
```

Genome-wide cross-mappabilities stored in the `cross_mappability.txt` file in the computation directory should look like below:
```
ENSG00000000003.10      ENSG00000008394.8       2
ENSG00000000003.10      ENSG00000011465.12      19
ENSG00000000003.10      ENSG00000053524.7       19
ENSG00000000003.10      ENSG00000059804.11      14
```
According to the above cross-mappability file, 2 k-mers from ENSG00000000003.10 map to ENSG00000008394.8, allowing for pre-configured number of mismatches. Note: Only genes with mappability < 1 (not NA) appear here. Cross-mappability to/from a gene with mappability of NA is indeterminable (NA). Cross-mappability between any gene pair not present in the above file is 0.
