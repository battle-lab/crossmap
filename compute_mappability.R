### this script computes average mappability of genes based on exons and UTRs.
### exon: mappability computed as average of all 75-mer's mappability (configurable k)
### UTR: mappability computed as average of all 36-mer's mappability (configurable k)
### took an weighted average of both mappabilities. weights are proportional to lengths.

suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(intervals))
suppressMessages(library(stats))
suppressMessages(library(argparser))

args <- arg_parser('program')
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
args <- add_argument(args, '-verbose', 
                     help='show computation status if verbose > 0',
                     default=1)
args <- add_argument(args, '-o', 
                     help='output file. extra files will also be generated with added suffixes.',
                     default='results/gene_mappability.txt')

argv <- parse_args(args)
annot_fn = argv$annot
k_exon = argv$k_exon
k_utr = argv$k_utr
bed75_fn = argv$kmap_exon
bed36_fn = argv$kmap_utr
verbose = argv$verbose
out_fn = argv$o

### function to verbose output
verbose_print <- function(msg){
  if(verbose > 0){
    print(sprintf("[%s] %s", format(Sys.time(), "%D %T"), msg))
  }
}

##### read input files
verbose_print('reading input files ...')
annot_df = fread(input = annot_fn, sep='\t', header=T, stringsAsFactors = F, data.table = F, showProgress = verbose>0)
utr_data_formatted <- annot_df[annot_df$feature == 'UTR', ]
exon_data_formatted <- annot_df[annot_df$feature == 'exon', ]

bed75 = fread(input = bed75_fn, sep='\t', header=F, stringsAsFactors = F, data.table = F, colClasses = c('character', 'numeric', 'numeric', 'numeric'), showProgress = verbose>0)
bed75$V2 = bed75$V2 + 1  # make 1 based index with inlcuded end index

bed36 = fread(input = bed36_fn, sep='\t', header=F, stringsAsFactors = F, data.table = F, colClasses = c('character', 'numeric', 'numeric', 'numeric'), showProgress = verbose>0 )
bed36$V2 = bed36$V2 + 1 # make 1 based index with inlcuded end index


#### sanity check: does any gene come from multiple chromosomes?
n_chr_per_utr_gene <- tapply(utr_data_formatted$chr, INDEX = utr_data_formatted$gene_id, FUN = function(x) length(unique(x)))
if (any(n_chr_per_utr_gene > 1))
  stop('some gene(s) comes from multiple chromosomes. following exon-interval merging will not work for it!')

n_chr_per_utr_gene <- tapply(exon_data_formatted$chr, INDEX = exon_data_formatted$gene_id, FUN = function(x) length(unique(x)))
if (any(n_chr_per_utr_gene > 1))
  stop('some gene(s) comes from multiple chromosomes. following utr-interval merging will not work for it!')


##### divide annotation data by chromosome and compute mappability (reason: processing time + memory)
get_mappability_by_chr <- function(cur_chr, flag='UTR'){
  verbose_print('subsetting chromosome data ...')
  if(flag == 'UTR'){
    chr_annot_data <- utr_data_formatted[utr_data_formatted$chr == cur_chr, ]
    chr_bed_data <- bed36[bed36$V1 == cur_chr, ]
    K <- k_utr
  } else {
    chr_annot_data <- exon_data_formatted[exon_data_formatted$chr == cur_chr, ]
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
  pos2mappability <- array(0, dim = max_index)
  tmp <- apply(X = chr_bed_data[,2:4], MARGIN = 1, FUN = function(row) {
    pos2mappability[row[1]:row[2]] <<- row[3]
    NA
  })
  
  ### calculate average mappability of every gene
  verbose_print('computing average mappability per gene ...')
  gene_mappability <- sapply(unique(chr_annot_data$gene_id), function(gid){
    #gid <- chr_annot_data$gene_id[1] # to debug
    gr <- genome_regions[[gid]]
    total_len <- sum(gr$V2 - gr$V1 + 1)
    mappabilities <- apply(gr, MARGIN = 1, function(row){
      if((row[2] - row[1]) + 1 < K){
        return(c())
      }      
      return(pos2mappability[row[1]:(row[2]-K+1)])
    })
    mappabilities_vector = unlist(mappabilities)
    avg_mappability = ifelse(length(mappabilities_vector)>0, mean(mappabilities_vector), NA)
    return(list(mappability=avg_mappability, length=total_len, n_kmer=length(mappabilities_vector)))
  })
  
  return(t(gene_mappability))
}


##### compute chromosome-wise mappability
chromosomes <- unique(exon_data_formatted$chr)
all_mappabilities <- lapply(chromosomes, function(cur_chr){
  # cur_chr = 'chr22' # to debug
  verbose_print( paste('########## computing utr-mappability of', cur_chr, '##########'))
  utr_mappability <- get_mappability_by_chr(cur_chr, flag='UTR')
  verbose_print( paste('########## computing exon-mappability of', cur_chr, '##########'))
  exon_mappability <- get_mappability_by_chr(cur_chr, flag='exon')
  
  # merge utr- and exon-mappability and compute weighted mean
  verbose_print( paste('########## computing weighted-mappability of', cur_chr, '##########'))
  utr_m = data.frame(gene=rownames(utr_mappability), mappability=unlist(utr_mappability[,'mappability']), length=unlist(utr_mappability[,'length']), kmer=unlist(utr_mappability[,'n_kmer']), stringsAsFactors = F)
  exon_m = data.frame(gene=rownames(exon_mappability), mappability=unlist(exon_mappability[,'mappability']), length=unlist(exon_mappability[,'length']), kmer=unlist(exon_mappability[,'n_kmer']), stringsAsFactors = F)
  both_m <- base::merge(exon_m, utr_m, by='gene', all.x=T, all.y=F, suffixes = c('.exon', '.utr') )
  weighted_mappability <- apply(both_m[c('mappability.exon','mappability.utr','length.exon','length.utr')], MARGIN = 1, FUN = function(row){
    weighted.mean(c(row['mappability.exon'],row['mappability.utr']), c(row['length.exon'],row['length.utr']), na.rm=T)
  })
  weighted_mappability_df <- data.frame(gene=both_m$gene, mappability=weighted_mappability, stringsAsFactors = F)
  both_m$weighted_mappability <- weighted_mappability
  
  return(list(weighted_mappability = weighted_mappability_df,
         info=both_m))
})


##### merge chromosome-wise mappabilities
verbose_print('merging all mappabilities ...')
if(length(all_mappabilities) == 0)
  stop('serious error - mappability could not be computed!')

final_df <- all_mappabilities[[1]][['weighted_mappability']]
for(i in 2:length(all_mappabilities)){
    final_df <- rbind(final_df, all_mappabilities[[i]][['weighted_mappability']])
}

info_df <- all_mappabilities[[1]][['info']]
for(i in 2:length(all_mappabilities)){
  info_df <- rbind(info_df, all_mappabilities[[i]][['info']])
}


##### save final data
write.table(x = final_df, file = out_fn, sep = '\t', quote = F, row.names = F, col.names = F)
write.table(x = info_df, file = paste0(out_fn, '.info'), sep = '\t', quote = F, row.names = F, col.names = T)

