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


### 2. Process annotation data
First we need to process the gene annotation file. We will save the exon and UTR annotations in a tabular format in a subdirectory ("annot") inside our computation directory.

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
Note: please wait for the above script to finish before moving to the next step.

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
