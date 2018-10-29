##### this script finds imperfectly-mapped genes and saves the imperfect kmers in file for each gene.

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(intervals))
suppressMessages(library(seqinr))

args <- arg_parser('program')
args <- add_argument(args, '-mappability',
                     help='mappability file',
                     default = '/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/gene_mappability/gene_mappability.txt')
args <- add_argument(args, '-genome', 
                     help='genome sequence directory',
                     default='/work-zfs/abattle4/lab_data/hg19')
args <- add_argument(args, '-annot', 
                     help='exon and UTR annotation file (txt)',
                     default='/work-zfs/abattle4/ashis/progres/crossmap/test_hg19/annot/annot.exon_utr.txt')
args <- add_argument(args, '-k_exon', 
                     help='k-mer length for exon',
                     default=75)
args <- add_argument(args, '-k_utr', 
                     help='k-mer length for utr',
                     default=36)
args <- add_argument(args, '-kmap_exon', 
                     help='bedgraph file containing k-mer mappabilities where k=k_exon',
                     default='/work-zfs/abattle4/lab_data/annotation/kmer_alignability_hg19/wgEncodeCrgMapabilityAlign75mer.bed')
args <- add_argument(args, '-kmap_utr', 
                     help='bedgraph file containing k-mer mappabilities where k=k_utr',
                     default='/work-zfs/abattle4/lab_data/annotation/kmer_alignability_hg19/wgEncodeCrgMapabilityAlign36mer.bed')
args <- add_argument(args, '-th1',
                     help='starting appability threshold. genes w/ avg mappability >= th1 are included.',
                     default=0)
args <- add_argument(args, '-th2',
                     help='ending appability threshold. genes w/ avg mappability < th2 are included.',
                     default=1.0)
args <- add_argument(args, '-dir_name_len',
                     help='length of the sub-directory names. First dir_name_len letters from gene form the name.',
                     default=12)
args <- add_argument(args, '-verbose', 
                     help='show computation status if verbose > 0',
                     default=1)
args <- add_argument(args, '-o',
                     help='output directory',
                     default='results/multimapped_kmers')

argv <- parse_args(args)
mappability_fn <- argv$mappability
genome_dir <- argv$genome
annot_fn = argv$annot
k_exon = argv$k_exon
k_utr = argv$k_utr
bed75_fn = argv$kmap_exon
bed36_fn = argv$kmap_utr
mappability_threshold1 <- argv$th1
mappability_threshold2 <- argv$th2
dir_name_len <- argv$dir_name_len
verbose = argv$verbose
out_dir <- argv$o

### check input
stopifnot(file.exists(mappability_fn))
stopifnot(file.exists(annot_fn))
stopifnot(file.exists(bed75_fn))
stopifnot(file.exists(bed36_fn))
stopifnot(dir.exists(genome_dir))
stopifnot(dir.exists(out_dir))
stopifnot(is.numeric(dir_name_len) & dir_name_len>=1)


### function to verbose output
verbose_print <- function(msg){
  if(verbose > 0){
    print(sprintf("[%s] %s", format(Sys.time(), "%D %T"), msg))
  }
}


##### read inputs
verbose_print('reading input files ...')

mappability <- fread(input = mappability_fn, sep = '\t', header = F, stringsAsFactors = F, colClasses = c('character', 'numeric'), col.names = c('gene', 'mappability'), data.table = F)

annot_df = fread(input = annot_fn, sep='\t', header=T, stringsAsFactors = F, data.table = F, showProgress = verbose>0)
utr_data_formatted <- annot_df[annot_df$feature == 'UTR', ]
exon_data_formatted <- annot_df[annot_df$feature == 'exon', ]

bed75 = fread(input = bed75_fn, sep='\t', header=F, stringsAsFactors = F, data.table = F, colClasses = c('character', 'numeric', 'numeric', 'numeric'), showProgress = verbose>0 )
bed75$V2 = bed75$V2 + 1 # make 1 based index. also make end index including.

bed36 = fread(input = bed36_fn, sep='\t', header=F, stringsAsFactors = F, data.table = F, colClasses = c('character', 'numeric', 'numeric', 'numeric'), showProgress = verbose>0 )
bed36$V2 = bed36$V2 + 1 # make 1 based index. also make end index including.

##### find genes within mappability threshold and sort them by mappability
mappability_imperfect <- mappability[which(mappability$mappability >= mappability_threshold1 & mappability$mappability < mappability_threshold2), ]
mappability_imperfect_sorted <- mappability_imperfect[order(-mappability_imperfect$mappability), ]
target_genes <- mappability_imperfect_sorted$gene

rm(list=c('annot_df', 'mappability', 'mappability_imperfect', 'mappability_imperfect_sorted'))
gc(reset = T)

# filter annotations
exon_data_formatted <- exon_data_formatted[exon_data_formatted$gene_id %in% target_genes,] # annotation of target genes only
utr_data_formatted <- utr_data_formatted[utr_data_formatted$gene_id %in% target_genes,] # annotation of target genes only


#### sanity check: does any gene come from multiple chromosomes?
n_chr_per_exon_gene <- tapply(exon_data_formatted$chr, INDEX = exon_data_formatted$gene_id, FUN = function(x) length(unique(x)))
if (any(n_chr_per_exon_gene > 1))
  stop('some gene(s) comes from multiple chromosomes. following exon-interval merging will not work for it!')

n_chr_per_utr_gene <- tapply(utr_data_formatted$chr, INDEX = utr_data_formatted$gene_id, FUN = function(x) length(unique(x)))
if (any(n_chr_per_utr_gene > 1))
  stop('some gene(s) comes from multiple chromosomes. following utr-interval merging will not work for it!')

### function to get subdirectory for a gene
get_subdir <- function(gid, out_dir, dir_name_len){
  paste0(c(out_dir, '/', substr(gid, 1, dir_name_len)), collapse = "")
}

### function to get kmer files for a gene
get_kmer_file_names <- function(gid, out_dir, dir_name_len){
  out_subdir = get_subdir(gid, out_dir = out_dir, dir_name_len = dir_name_len)
  text_fn = paste0(out_subdir, '/',  gid, '.kmer.txt')
  fasta_fn = paste0(out_subdir, '/',  gid, '.kmer.fa')
  count_fn = paste0(out_subdir, '/',  gid, '.count.txt')
  return(list(text_fn = text_fn, fasta_fn = fasta_fn, count_fn = count_fn))
}

##### save target genes' kmers
save_kmers_of_genes_by_chr <- function(cur_chr, flag='UTR', append=F){
  verbose_print('subsetting chromosome data ...')

  if(flag=='UTR'){
    chr_annot_data <- utr_data_formatted[utr_data_formatted$chr == cur_chr, ]
    if (nrow(chr_annot_data) <= 0){
      verbose_print(paste0('no imperfect gene in ', cur_chr))
      return(NULL)
    }
    chr_bed_data <- bed36[bed36$V1 == cur_chr, ]
    K <- k_utr
  } else {
    chr_annot_data <- exon_data_formatted[exon_data_formatted$chr == cur_chr, ]
    if (nrow(chr_annot_data) <= 0){
      verbose_print(paste0('no imperfect gene in ', cur_chr))
      return(NULL)
    }
    chr_bed_data <- bed75[bed75$V1 == cur_chr, ]
    K <- k_exon
  }
  
  verbose_print('finding genome regions ...')
  genome_regions = by(chr_annot_data, as.factor(chr_annot_data$gene_id), simplify = F,
                      FUN=function(df){
                        idf <- Intervals(df[,c('start_pos', 'end_pos'), drop=F])
                        merged_intervals <- as.data.frame(interval_union(idf))
                        return(merged_intervals)
                      })
  
  ### hash to acccess mappability starting at position
  verbose_print('generating mappability per position array ...')
  max_index <- max(chr_bed_data$V3)
  pos2mappability <- array(1, dim = max_index)
  tmp <- apply(X = chr_bed_data[,2:4, drop=F], MARGIN = 1, FUN = function(row) {
    pos2mappability[row[1]:row[2]] <<- row[3]
    NA
  })
  
  ### reading genome
  verbose_print('reading genome ...')
  chr_genome_fn <- paste0(genome_dir, '/', cur_chr, '.fa')
  chr_dna <- read.fasta(file = chr_genome_fn, seqtype = 'DNA', as.string = T, forceDNAtolower = T, seqonly = T)[[1]]
  chr_dna <- tolower(chr_dna)  # lowercase sequence (read.fasta does not make lower case)
  chr_dna <- strsplit(chr_dna, '')[[1]]
  ### determine kmers with imperfect mappability of every gene
  verbose_print('determining kmers with imperfect mappability of every gene ...')
  n_imperfect_kmers_per_gene <- sapply(unique(chr_annot_data$gene_id), function(gid){
    #gid <- chr_annot_data$gene_id[1] # to debug
    gr <- genome_regions[[gid]]
    imperfect_kmer_positions <- apply(gr, MARGIN = 1, function(row){
      if((row[2] - row[1] + 1) < K){
        return(c())
      }
      indexes <- row[1]:(row[2]-K+1)
      mappabilities <- pos2mappability[indexes]
      kmer_positions <- indexes[which(mappabilities < 1)]
      return(kmer_positions)
    })
    imperfect_kmer_positions_vector = unlist(imperfect_kmer_positions)
    imperfect_kmers <- unlist(lapply(imperfect_kmer_positions_vector, function(st) paste0(chr_dna[st:(st+K-1)], collapse='') ))
    ### count imperfect k-mers
    imperfect_kmers_count = table(imperfect_kmers)
    if(length(imperfect_kmers_count) > 0){
      out_subdir = get_subdir(gid, out_dir = out_dir, dir_name_len = dir_name_len)
      if(!dir.exists(out_subdir)) 
        dir.create(out_subdir)
      out_files <- get_kmer_file_names(gid, out_dir = out_dir, dir_name_len = dir_name_len)
      write(names(imperfect_kmers_count), out_files$text_fn, sep="\n", append=append)
      write(as.numeric(imperfect_kmers_count), out_files$count_fn, sep="\n", append=append)
    }
    return(length(imperfect_kmers_count))
  })
  rm(chr_bed_data)
  rm(chr_annot_data)
  rm(chr_dna)
  gc(reset = T)
  return(n_imperfect_kmers_per_gene)
}


chromosomes <- unique(utr_data_formatted$chr)
all_n_imperfect_kmers <- lapply(chromosomes, function(cur_chr){
  # cur_chr = 'chr22' # to debug
  verbose_print( paste('########## saving UTR kmers of', cur_chr, '##########'))
  n_imperfect_kmers_utr <- save_kmers_of_genes_by_chr(cur_chr, 'UTR', append=F)
  verbose_print( paste('########## saving exon kmers of', cur_chr, '##########'))
  n_imperfect_kmers_exon <- save_kmers_of_genes_by_chr(cur_chr, 'exon', append=T)
  n_imperfect_kmers <- c(n_imperfect_kmers_exon, n_imperfect_kmers_utr)
  
  # ### create fasta file and if required (k_exon == k_utr), merge duplicate k-mers
  # chr_imperfect_genes <- union(names(n_imperfect_kmers_utr[n_imperfect_kmers_utr>0]), names(n_imperfect_kmers_exon[n_imperfect_kmers_exon>0])) 
  # tmp <- lapply(chr_imperfect_genes, function(igid){
  #   kmer_files <- get_kmer_file_names(gid = igid, out_dir = out_dir, dir_name_len = dir_name_len)
  #   kmer_dt <- fread(input = kmer_files$text_fn, sep = '\t', header = F, stringsAsFactors = F, showProgress = F, data.table = T, colClasses = 'character')
  #   
  #   # merge duplicate k-mers
  #   if(k_exon == k_utr){
  #     count_dt <- fread(input = kmer_files$count_fn, sep = '\t', header = F, stringsAsFactors = F, showProgress = F, data.table = T, colClasses = 'integer')
  #     stopifnot(nrow(kmer_dt)==nrow(count_dt))
  #     kmer_dt$count = count_dt[,V1]
  #     combined_kmer_dt = kmer_dt[,.(count.sum = sum(count)),by=V1]
  #     kmer_dt = combined_kmer_dt[,.(V1)]
  #     #write.table(combined_kmer_dt$V1, file = kmer_files$text_fn, col.names = F, row.names = F, quote = F)
  #     #write.table(combined_kmer_dt$count.sum, file = kmer_files$count_fn, col.names = F, row.names = F, quote = F)
  #     write(as.character(combined_kmer_dt$V1), file = kmer_files$text_fn, sep="\n")
  #     write(as.numeric(combined_kmer_dt$count.sum), file = kmer_files$count_fn, sep="\n")
  #   }
  #   
  #   # make fasta
  #   fasta_txt = paste0('>', 0:(nrow(kmer_dt)-1), '\n', kmer_dt$V1, collapse = '\n' )
  #   write(fasta_txt, kmer_files$fasta_fn, sep="\n")
  #   
  #   # clean memory
  #   rm(list=intersect(ls(), c('kmer_dt', 'count_dt', 'combined_kmer_dt', 'fasta_txt')))
  #   gc(reset = T)
  #   
  #   return()
  # })
  
  return(n_imperfect_kmers)
})
