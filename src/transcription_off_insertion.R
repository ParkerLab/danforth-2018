# plot the normalized RNA-seq coverage of the insertion

library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = T)
SIGNAL_DIR <- args[1]

coverage_files <- list.files(SIGNAL_DIR, pattern = 'signal.bed', full.names = T)

parse_file_name <- function(f) {
  regex <- "^(\\d+)\\.(fwd|rev).signal.bed$"
  seqcore_id <- gsub(regex, '\\1', basename(f))
  strand <- gsub(regex, '\\2', basename(f))
  
  return(c('seqcore_id' = seqcore_id, 'strand' = strand))
}

load_coverage_file <- function(f) {
  tmp <- read.table(f, head = F, as.is = T, sep = '\t')
  colnames(tmp) <- c('chrom', 'start', 'end', 'coverage')
  tmp$seqcore_id <- parse_file_name(f)['seqcore_id']
  tmp$strand <- parse_file_name(f)['strand']
    
  return(tmp)
}

coverage <- lapply(coverage_files, load_coverage_file)
coverage <- bind_rows(coverage)

coverage$position <- coverage$start
coverage$seqcore_id[coverage$seqcore_id=='67029'] <- 'WT_1'
coverage$seqcore_id[coverage$seqcore_id=='67030'] <- 'WT_2'
coverage$seqcore_id[coverage$seqcore_id=='67031'] <- 'WT_3'
coverage$seqcore_id[coverage$seqcore_id=='67032'] <- 'Sd_1'
coverage$seqcore_id[coverage$seqcore_id=='67033'] <- 'Sd_2'
coverage$seqcore_id[coverage$seqcore_id=='67034'] <- 'Sd_3'
coverage$replicate <- gsub(".*_", "", coverage$seqcore_id)
coverage$genotype <- gsub("_.*", "", coverage$seqcore_id)

coverage$coverage <- abs(coverage$coverage)
coverage$strand[coverage$strand=='rev'] <- "- strand"
coverage$strand[coverage$strand=='fwd'] <- "+ strand"

p <- ggplot(coverage) + geom_line(aes(x = position, y = coverage, color = genotype, alpha = replicate)) +
  theme_bw() +
  scale_alpha_manual(values = c('1' = 1, '2' = 0.75, '3' = 0.5)) +
  geom_vline(xintercept = 19355026, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 19363566, color = 'black', linetype = 'dashed') +
  ylab("RNA-seq signal (normalized for library size)") + xlab("Genomic position") +
  facet_wrap(~strand)

pdf("insertion_coverage.pdf", width = 9, height = 4)
p
dev.off()
