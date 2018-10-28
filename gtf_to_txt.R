# this script converts gene and transcript gtf file to tab-delimitted txt file.
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

args <- arg_parser('program')
args <- add_argument(args, '-gtf', 
                     help='gtf file',
                     default='/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gtf')
args <- add_argument(args, '-f',
                     help='feautres to include, separated by comma',
                     default='exon,UTR')
args <- add_argument(args, '-o',
                     help='output text file',
                     default='results/annot.txt')


argv <- parse_args(args)
gtf_fn <- argv$gtf
features_str <- argv$f
out_fn <- argv$o

### input check
stopifnot(file.exists(gtf_fn))
stopifnot(dir.exists(dirname(out_fn)))

### parse features
parse_delimitted_param <- function(param_str, delim=',', rm.empty=T){
  stopifnot(class(param_str)=='character' && length(param_str)==1)
  parts = strsplit(param_str, split = delim)[[1]]
  if(rm.empty==T){
    is_valid_parts = sapply(parts, function(s) nchar(s)>0)
    parts = parts[is_valid_parts]
  }
  return(parts)
}

### define function to read a gtf file
read_gtf <- function(gtf_fn, sep = '\t'){
  ### determine how many header lines to skip
  header_line_count = 0
  con = file(gtf_fn, open="r")
  while(TRUE){
    line = readLines(con, n=1)
    if (substr(line, start = 1, stop = 1) == "#"){
      header_line_count = header_line_count + 1
    } else {
      break
    }
  }
  close(con)
  
  ### read gtf file without the header lines
  gtf_df = fread(input = gtf_fn, sep = sep, skip = header_line_count, header = F, stringsAsFactors = F, data.table = F, verbose = F, showProgress = F)
  return(gtf_df)
}

### define function to get a field value from the last column of a gtf file
get_field <- function(str, field_id){
  pattern = paste0(field_id, " \"[^\\\"]+\"")
  m = regexec(pattern, str)[[1]]
  if(m>0){
    field_val = substr(str, m+nchar(field_id)+2, m+attr(m, 'match.length')-2)
    field_val
  } else {
    field_val = NA
  }
  
  return(field_val)
}

### define function to convert gene gtf dataframe to a table
gene_gtf_to_table <- function(gtf_df){
  gene_id <- unlist(lapply(gtf_df[, 9], get_field, field_id='gene_id'))
  chr <- gtf_df[,1]
  annotation_source <- gtf_df[,2]
  feature <- gtf_df[,3]
  start_pos <- gtf_df[,4]
  end_pos <- gtf_df[,5]
  strand <- gtf_df[,7]
  gene_name <- unlist(lapply(gtf_df[, 9], get_field, field_id='gene_name'))
  gene_type <- unlist(lapply(gtf_df[, 9], get_field, field_id='gene_type'))
  
  gene_annotation_df <- data.frame(gene_id = gene_id,
                                   chr = chr,
                                   annotation_source = annotation_source,
                                   feature = feature,
                                   start_pos = start_pos,
                                   end_pos = end_pos,
                                   strand = strand,
                                   gene_name = gene_name,
                                   gene_type = gene_type,
                                   stringsAsFactors = F)
  
}

### read annotation files
gene_annotation = read_gtf(gtf_fn)

### inlcude given features only
features_to_include = parse_delimitted_param(features_str, delim = ',')
if(length(features_to_include)>0){
  gene_annotation = gene_annotation[gene_annotation$V3 %in% features_to_include, ]
}

### convert 
gene_annotation_df = gene_gtf_to_table(gtf_df = gene_annotation)

### save
write.table(gene_annotation_df, file = out_fn, sep = '\t', row.names = F, col.names = T, quote = F)
