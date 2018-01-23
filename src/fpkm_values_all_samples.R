library(dplyr)
library(tidyr)
library(ggplot2)

ROOT <- Sys.getenv('DANFORTH_HOME')

source(file.path(ROOT, 'src', 'counts_to_rpkm.R'))

COUNTS_DIR <- file.path(ROOT, 'work', 'differential_gene_expression', 'results', 'counts')
SAMPLE_INFO <- file.path(ROOT, 'sample_information', 'sample_info.txt')
TRANSCRIPT_LENGTHS_FILE <- file.path(ROOT, 'data', 'transcript_lengths', 'mm9.transcript_lengths')
OUT <- 'fpkm_values.csv'

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t', comment.char = '')
sample_info <- sample_info %>% dplyr::filter(experiment_type=='RNAseq')

# load in the expression values for each WT mouse
COUNTS_FILES <- list.files(COUNTS_DIR, full.names = T)

parse_file_name <- function(f) {
  regex <- '^(\\d+)\\.geneCounts.txt.gz$'
  seqcore_id <- gsub(regex, "\\1", basename(f))
  return(c('seqcore_id' = seqcore_id))
}

load_count_file <- function(f) {
  seqcore_id <- parse_file_name(f)['seqcore_id']
  
  tmp <- read.table(gzfile(f), head = F, as.is = T, sep = '\t')
  tmp$seqcore_id <- seqcore_id
  
  return(tmp)
}

counts <- lapply(COUNTS_FILES, load_count_file)
counts <- bind_rows(counts)
colnames(counts) <- c('gene', 'count', 'seqcore_id')
counts <- counts[counts$seqcore_id %in% sample_info$seqcore_id,]
counts <- counts[grep("^ENS", counts$gene),]
counts$gene <- gsub("\\..*$", "\\1", counts$gene)
counts <- tidyr::spread(counts, key = seqcore_id, value = count)

# convert to FPKM
transcript_lengths <- read.table(TRANSCRIPT_LENGTHS_FILE, head = F, as.is = T, sep = '\t')
colnames(transcript_lengths) <- c('gene', 'transcript', 'gene_name', 'length')
gene_lengths <- transcript_lengths %>% group_by(gene) %>%
  dplyr::summarize(length=mean(length))

counts <- left_join(counts, gene_lengths)

for(i in 2:(ncol(counts)-1)) {
  counts[,i] <- counts_to_rpkm(counts[,i], lengths = counts$length)
}

counts <- dplyr::select(counts, -length)
counts$avg <- apply(counts[,2:ncol(counts)], 1, mean)

counts <- counts[rev(order(counts$avg)),]
counts <- counts[!is.nan(counts$avg),] # throw out one gene whos only transcript has length = 0 (ENSMUSG00000092330 ENSMUST00000174222 Vmn2r-ps82)
counts <- dplyr::select(counts, -avg)

# change column names to encode genotype
for(i in colnames(counts)[2:ncol(counts)]) {
  genotype <- unique(sample_info$genotype[sample_info$seqcore_id==i])
  colnames(counts)[match(i, colnames(counts))] <- paste0(genotype, '_sample_', i, '_FPKM')
}

# add in human-readable names
human_readable_names <- unique(transcript_lengths[,c('gene', 'gene_name')])
counts <- left_join(counts, human_readable_names)

# reformat slightly
counts <- counts %>% dplyr::rename(gene_id=gene)
counts <- counts[,c('gene_id', 'gene_name', grep('sample_', colnames(counts), value = T))]

fpkms <- counts %>% tidyr::gather(key = sample, value = fpkm, WT_sample_67029_FPKM, WT_sample_67030_FPKM, WT_sample_67031_FPKM, Sd_sample_67032_FPKM, Sd_sample_67033_FPKM, Sd_sample_67034_FPKM)

#p <- ggplot(fpkms) + geom_density(aes(x = log10(fpkm), fill = sample), alpha = 0.2)
#p

write.table(counts, file = OUT, append = F, quote = F, row.names = F, col.names = T, sep = ',')
