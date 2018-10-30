### this script computes cross-mappabilities genome-wide
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(intervals))

args <- arg_parser('program')
args <- add_argument(args, '-annot', 
                     help='exon and UTR annotation file (txt)',
                     default='/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/annot/annot.exon_utr.txt')
args <- add_argument(args, '-mappability',
                     help='mappability file',
                     default = '/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/gene_mappability/gene_mappability.txt')
args <- add_argument(args, '-kmer', 
                     help='ambiguous kmer directory',
                     default='/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/ambiguous_kmers')
args <- add_argument(args, '-align', 
                     help='kmer alignment directory',
                     default='/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/ambiguous_kmers_alignment')
args <- add_argument(args, '-index', 
                     help='bowtie index prefix',
                     default='/work-zfs/abattle4/lab_data/hg19/prebuilt_bowtie_index/hg19')
args <- add_argument(args, '-n1',
                     help='compute cross-mappability from n1-th to n2-th gene in the mappabiity file. negative number to use all genes.',
                     default=-1)
args <- add_argument(args, '-n2',
                     help='compute cross-mappability from n1-th to n2-th gene in the mappabiity file.  negative number to use all genes.',
                     default=-1)
args <- add_argument(args, '-mismatch',
                     help='maximum number of mismatch for alignment',
                     default=2)
args <- add_argument(args, '-max_chr', 
                     help='maximum number of chromosomes to load in memory',
                     default=7)
args <- add_argument(args, '-initonly', 
                     help='initialize required resources only, without computing cross-mappability',
                     default=FALSE)
args <- add_argument(args, '-dir_name_len',
                     help='length of the sub-directory names in ambiguous kmer directory.',
                     default=12)
args <- add_argument(args, '-verbose', 
                     help='show computation status if verbose > 0',
                     default=1)
args <- add_argument(args, '-o',
                     help='output directory',
                     default='results/cross_mappability')


argv <- parse_args(args)
annot_fn <- argv$annot
mappability_fn <- argv$mappability
kmer_dir <- argv$kmer
alignment_dir <- argv$align
bowtie_index_prefix <- argv$index
n1 <- argv$n1
n2 <- argv$n2
max_mismatch <- argv$mismatch
max_chr <- argv$max_chr
init_only_input <- argv$initonly
dir_name_len <- argv$dir_name_len
verbose = argv$verbose
out_dir <- argv$o

### check input
stopifnot(file.exists(mappability_fn))
stopifnot(file.exists(annot_fn))
stopifnot(dir.exists(out_dir))
stopifnot(dir.exists(kmer_dir))
stopifnot(dir.exists(alignment_dir))
stopifnot(dir.exists(dirname(bowtie_index_prefix)))
stopifnot(is.numeric(dir_name_len) & dir_name_len>=1)
stopifnot(is.numeric(max_chr) & max_chr>=1)
stopifnot(max_mismatch %in% 0:3)
stopifnot(is.numeric(n1))
stopifnot(is.numeric(n2))
stopifnot((n1<1 || n2 <1) || (n1>=1 && n2>=1 && n1<=n2))

max_chr = as.integer(max_chr)
resource_initialization_only = ifelse(is.na(as.logical(init_only_input)), FALSE, as.logical(init_only_input))
n1 = as.integer(n1)
n2 = as.integer(n2)

### function to verbose output
verbose_print <- function(msg){
  if(verbose > 0){
    print(sprintf("[%s] %s", format(Sys.time(), "%D %T"), msg))
  }
}

### function to get subdirectory for a gene
get_subdir <- function(gid, out_dir, dir_name_len){
  paste0(c(out_dir, '/', substr(gid, 1, dir_name_len)), collapse = "")
}

### function to get kmer files for a gene
get_kmer_file_names <- function(gid, kmer_dir, dir_name_len){
  out_subdir = get_subdir(gid, out_dir = kmer_dir, dir_name_len = dir_name_len)
  text_fn = paste0(out_subdir, '/',  gid, '.kmer.txt')
  fasta_fn = paste0(out_subdir, '/',  gid, '.kmer.fa')
  count_fn = paste0(out_subdir, '/',  gid, '.count.txt')
  return(list(text_fn = text_fn, fasta_fn = fasta_fn, count_fn = count_fn))
}

##### read inputs
verbose_print('reading input files ...')
mappability <- fread(input = mappability_fn, sep = '\t', header = F, stringsAsFactors = F, colClasses = c('character', 'numeric'), col.names = c('gene', 'mappability'), data.table = F)

annot_df = fread(input = annot_fn, sep='\t', header=T, stringsAsFactors = F, data.table = F, showProgress = verbose>0)
annotation_formatted <- annot_df[annot_df$feature == 'UTR' | annot_df$feature == 'exon', ]

##### filter genes: n1-th to n2-th gene in the mappability data
if(n1>=1 && n2>=1){
  target_mappability <- mappability[n1:n2, , drop = F]
  target_genes <- target_mappability[!is.na(target_mappability$mappability) &  target_mappability$mappability < 1, 'gene']
} else {
  target_genes <- mappability[!is.na(mappability$mappability)  &  mappability$mappability < 1, 'gene']
}

if(length(target_genes) == 0){
  warning('There is no gene to compute cross-mappability from.')
  quit(save = 'no', status=0)
}

### get imperfect genes (mappability < 1, and not NA)
imperfect_genes <- sort(mappability[!is.na(mappability$mappability) & mappability$mappability < 1, 'gene'])
imperfect_genes_2_idx <- new.env(hash = T)
tmp <- lapply(1:length(imperfect_genes), function(i) {imperfect_genes_2_idx[[imperfect_genes[i]]] <<- i; return()})

##### clear variables not required anymore
rm(annot_df)
gc(reset = T)

##### find genes in every position
### this function returns a list of vectors.
### i-th entry contains a vector of genes mapped to i-th position in cur_chr
get_pos2genes_mapping <- function(cur_chr){
  chr_annot_data <- annotation_formatted[annotation_formatted$chr == cur_chr, ]
  verbose_print('finding genome regions ...')
  genome_regions = by(chr_annot_data, as.factor(chr_annot_data$gene_id), simplify = F,
                      FUN=function(df){
                        idf <- Intervals(df[,c('start_pos', 'end_pos'), drop=F])
                        merged_intervals <- as.data.frame(interval_union(idf))
                        return(merged_intervals)
                      })
  verbose_print('initializing genes per position array ...')
  max_index <- max(chr_annot_data$end_pos)
  pos2genes <- lapply(1:max_index, function(item) NULL)
  verbose_print('populating genes per position array ...')
  tmp <- lapply(intersect(unique(chr_annot_data$gene_id), imperfect_genes), function(gid){
    #gid <- chr_annot_data$gene_id[1] # to debug
    gr <- genome_regions[[gid]]
    gid_index <- imperfect_genes_2_idx[[gid]]
    tmp2 <- apply(gr, MARGIN = 1, function(row){
      pos2genes[row[1]:row[2]] <<- mapply(append, pos2genes[row[1]:row[2]], gid_index, SIMPLIFY=F)
      return(NA)
    })
    return(NA)
  })
  return(pos2genes)
}

############### code: low memory version #########
geneids_pos2genes <- imperfect_genes
chromosomes <- unique(annotation_formatted$chr)
pos2genes_by_chr = list()

load_chromosomes <- function(chromosomes){
  pos2genes_by_chr <<- list()
  tmp <- lapply(chromosomes, function(cur_chr){
    gc()
    chr_pos2genes_data_dir = paste0(out_dir,'/pos2gene')
    if(!dir.exists(chr_pos2genes_data_dir)) 
      dir.create(chr_pos2genes_data_dir)
    chr_pos2genes_data_fn = paste0(chr_pos2genes_data_dir,'/pos2gene_',cur_chr,'.RData')
    create_and_save_pos2gene_mapping = T
    if(file.exists(chr_pos2genes_data_fn)){
      verbose_print(paste('loading', cur_chr))
      # loads two variable: 'pos2genes', 'geneids_pos2genes'
      load(chr_pos2genes_data_fn, envir = .GlobalEnv)
      # make sure geneids_pos2genes and imperfect_genes are exactly same. 
      # otherwise pos2genes points to wrong genes. so, recreate them.
      create_and_save_pos2gene_mapping = !all(geneids_pos2genes == imperfect_genes)
    }
    if(create_and_save_pos2gene_mapping==T) {
      verbose_print(paste('generating mapping of ', cur_chr))
      pos2genes <- get_pos2genes_mapping(cur_chr)
      save(list = c('pos2genes', 'geneids_pos2genes'), file = chr_pos2genes_data_fn)
    }
    pos2genes_by_chr[[cur_chr]] <<- pos2genes
    return(NA)
  })
}

##### align and compute cross-mappability from each gene
find_and_save_cross_mappability <- function(g, chromosomes, delete_alignment=F, append_conflict=F){
  verbose_print(paste('finding conflicts with gene ', g))
  ### align kmers of gene g
  kmer_subdir = get_subdir(gid = g, out_dir = kmer_dir, dir_name_len = dir_name_len)
  alignment_subdir = get_subdir(gid = g, out_dir = alignment_dir, dir_name_len = dir_name_len)
  align_fn <- paste0(alignment_subdir, '/', g, '.alignment.txt')
  kmer_fn <- paste0(kmer_subdir, '/', g, '.kmer.fa')
  if(!file.exists(align_fn)){
    if(!dir.exists(alignment_subdir))
      dir.create(alignment_subdir)
    align_cmd <- paste0('bowtie -v ', max_mismatch,' -B 1 --quiet -a ', bowtie_index_prefix,  ' -f ', kmer_fn, '  | cut -f 1,3-4  > ', align_fn)
    system(align_cmd)
  }
  
  align_dt <- fread(input = align_fn, sep = '\t', header = F, stringsAsFactors = T, colClasses = c('integer', 'character', 'integer'), col.names = c('kmer', 'chr', 'pos'), data.table = T)
  align_dt$kmer= align_dt$kmer+1  # convert to 1-based k-mer index
    
  # read k-mer occurences from file
  kmer_counts_fn = paste0(kmer_subdir, '/', g, '.count.txt')
  kmer_counts_df = fread(input=kmer_counts_fn, header = F, sep='\t', stringsAsFactors = F, colClasses = 'numeric', data.table = F)
  gene_kmer_occurences =  kmer_counts_df[,1]
  
  if(delete_alignment == T){
    file.remove(align_fn)
  }
  
  ### get conflicting genes
  get_chr_conflicts <- function(kmer, pos, chr){
    chr = as.character(chr)
    empty_return_list = list(G1=character(), G2=character(), S1=double())
    if(!(chr %in% chromosomes)){
      return(empty_return_list)
    }
    
    valid_pos_idx = pos <= length(pos2genes_by_chr[[chr]])
    pos <- pos[valid_pos_idx]
    kmer <- kmer[valid_pos_idx]
    
    gene_idx_list <- pos2genes_by_chr[[chr]][pos]
    gene_idx_expanded <- unlist(gene_idx_list)
    if(length(gene_idx_expanded) == 0)
      return(empty_return_list)
    
    kmer_expanded <- unlist(mapply(rep, kmer, sapply(gene_idx_list, length), SIMPLIFY=F))
    
    conflict_df = data.frame(gene_idx=gene_idx_expanded, kmer=kmer_expanded)
    unique_conflict_df = unique(conflict_df)
    # S1 is the total number of k-mers in gene1 that map to gene2 (duplicate counts multiple times)
    S1_counts = tapply(gene_kmer_occurences[unique_conflict_df$kmer], unique_conflict_df$gene_idx, sum)
    S1 = data.frame(gene_idx=as.numeric(names(S1_counts)), S1 = as.numeric(S1_counts), stringsAsFactors = F)
    S1 = S1[S1$gene_idx != imperfect_genes_2_idx[[g]], ] # exclude self-crossmapping
    rm(valid_pos_idx, pos, kmer, gene_idx_list, gene_idx_expanded, kmer_expanded, conflict_df, unique_conflict_df, S1_counts, empty_return_list)
    gc(reset = T)
    return(list(G1=rep(g, nrow(S1)), G2 = imperfect_genes[S1$gene_idx], S1=S1$S1))
  }
  
  conflicting_genes <- align_dt[,get_chr_conflicts(kmer, pos, chr), by=chr]
  
  ### save conflicting gene pairs
  out_subdir = get_subdir(gid = g, out_dir = out_dir, dir_name_len = dir_name_len)
  out_fn <- paste0(out_subdir, '/', g, '.crossmap.txt')
  
  if(append_conflict==FALSE && file.exists(out_fn))
    file.remove(out_fn)
  
  if(nrow(conflicting_genes) > 0 ){
    if(!dir.exists(out_subdir))
      dir.create(out_subdir)
    write.table(conflicting_genes[,c('G1','G2', 'S1')], file = out_fn, sep = '\t', quote = F, row.names = F, col.names = F, append = append_conflict)
  }
  
  # clean memory
  rm(align_dt, kmer_counts_df, gene_kmer_occurences, conflicting_genes)
  gc(reset = T)
  
  return()
}

### initialize resources only, without computing cross-mappability
if(resource_initialization_only){
  for(chr in chromosomes){
    load_chromosomes(chr)
    rm(pos2genes_by_chr)
    gc(reset = T)
    pos2genes_by_chr <- list()
  }
  quit(save = 'no', status=0)
}

n_genes_per_batch = 1000
for(batch in 1:ceiling(length(target_genes)/n_genes_per_batch)){
  batch_target_genes <- target_genes[((batch-1)*n_genes_per_batch+1):min((batch*n_genes_per_batch), length(target_genes))]
  for(chr_batch in 1:ceiling(length(chromosomes)/max_chr)){
    batch_chromosomes <- chromosomes[((chr_batch-1)*max_chr+1):min((chr_batch*max_chr), length(chromosomes))]
    load_chromosomes(batch_chromosomes)
    tmp <- lapply(batch_target_genes, FUN = find_and_save_cross_mappability, chromosomes=batch_chromosomes,  delete_alignment=T, append_conflict=chr_batch!=1)
    rm(pos2genes_by_chr)
    gc(reset = T)
    pos2genes_by_chr <- list()
  }
}
