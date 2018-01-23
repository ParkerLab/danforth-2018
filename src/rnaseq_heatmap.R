library(dplyr)
library(ggplot2)
library(DESeq2)
library(cowplot)
library(ggrepel)

args <- commandArgs(trailingOnly = T)
RNASEQ <- args[1]
RNASEQ_COUNTS_DIR <- args[2]
OUT <- args[3]

shh_genes <- c('Ihh', 'Shh')

# First, prepare the RNA-seq data for the heatmap
rnaseq <- read.table(RNASEQ, head = T, sep = '\t', as.is = T)
rnaseq$in_hedgehog_pathway <- ifelse(rnaseq$gene_name %in% shh_genes, T, F)
rnaseq$in_hedgehog_pathway[rnaseq$gene_name %in% c("Smo", "Gli2")] <- T

# calculate the size factors
counts_files <- list.files(RNASEQ_COUNTS_DIR, pattern = 'geneCounts.txt.gz', full.names = T)
counts <- lapply(counts_files, function(f) {
  sample_id <- gsub("\\.geneCounts.txt.gz", "", basename(f))
  tmp <- read.table(f, head = F, as.is = T, sep = '\t')
  colnames(tmp) <- c('ensembl_gene_id', 'count')
  tmp$sample_id <- sample_id
  
  return(tmp)
})

counts <- bind_rows(counts)
counts <- counts[grep("ENS", counts$ensembl_gene_id),]
counts$ensembl_gene_id <- gsub("\\..*", "", counts$ensembl_gene_id)
counts.matrix <- counts %>% tidyr::spread(key = sample_id, value = count)
exp_info <- data.frame(sample = colnames(counts.matrix)[2:ncol(counts.matrix)], 
                       genotype = rep(c('WT', 'Sd'), each = 3))
counts.matrix <- DESeqDataSetFromMatrix(counts.matrix[,2:ncol(counts.matrix)], colData = exp_info, design = ~ genotype)
counts.matrix <- estimateSizeFactors(counts.matrix)
size_factors <- data.frame(size_factor = sizeFactors(counts.matrix),
                           sample_id = names(sizeFactors(counts.matrix)))

counts <- left_join(counts, size_factors) %>% mutate(count = count / size_factor) %>%
  dplyr::select(-size_factor)
counts$sample_id <- gsub(67029, "WT_1", counts$sample_id)
counts$sample_id <- gsub(67030, "WT_2", counts$sample_id)
counts$sample_id <- gsub(67031, "WT_3", counts$sample_id)
counts$sample_id <- gsub(67032, "Sd_1", counts$sample_id)
counts$sample_id <- gsub(67033, "Sd_2", counts$sample_id)
counts$sample_id <- gsub(67034, "Sd_3", counts$sample_id)

rnaseq <- left_join(rnaseq, counts)

other_genes_of_interest <- c('T', 'Noto', 'Cyp26a1', 'Wnt3a', 'Cdx2', 'Foxa2')
differential_counts <- filter(rnaseq, padj <= 0.05)
nondifferential_counts <- filter(rnaseq, padj > 0.05 & gene_name %in% other_genes_of_interest)

differential_counts$is_differential <- 'yes'
nondifferential_counts$is_differential <- 'no'

differential_counts <- rbind(differential_counts, nondifferential_counts)


# rescale the counts for each gene to be between 0 and 1?
differential_counts <- tidyr::spread(differential_counts, key = sample_id, value = count)
COUNT_COLS <- c("WT_1", "WT_2", "WT_3", "Sd_1", "Sd_2", "Sd_3")
differential_counts[,COUNT_COLS] <- differential_counts[,COUNT_COLS] / apply(differential_counts[,COUNT_COLS], 1, max)
differential_counts <- tidyr::gather(differential_counts, key = sample_id, value = count, Sd_1, Sd_2, Sd_3, WT_1, WT_2, WT_3)

# order by log2FC
differential_counts <- differential_counts[order(differential_counts$is_differential, differential_counts$lfcMLE),]
differential_counts$gene[is.na(differential_counts$gene)] <- 'Gm10664' # TODO: use the exact Ensembl gene ID...
differential_counts$gene <- factor(differential_counts$gene, ordered = T, levels = unique(differential_counts$gene))

differential_counts$sample_id <- gsub("_", "\n", differential_counts$sample_id)
differential_counts$sample_id <- factor(differential_counts$sample_id, levels = c("WT\n1", "WT\n2", "WT\n3", "Sd\n1", "Sd\n2", "Sd\n3"), ordered = T)

counts.heatmap <- ggplot(differential_counts %>% filter(!is.na(gene_name))) + geom_tile(aes(x = sample_id, y = gene, fill = count)) +
  theme_bw() +
  guides(fill = guide_colorbar(title = 'Expression\n(scaled 0-1)')) +
  scale_fill_gradient(high = 'red', low = 'white') +
  xlab("Sample") + ylab("Gene") +
  facet_grid(is_differential~., scales = 'free_y', space = 'free') +
  theme(axis.text.y = element_text(colour = 'black')) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(strip.background = element_blank()) + 
  theme(panel.spacing = unit(2, 'lines')) +
  theme(strip.text.y = element_blank()) +
  theme(legend.position = 'left')


tmp <- differential_counts %>% filter(!is.na(gene_name)) %>% dplyr::select(-count)
tmp$sample_id <- 'None\nNone'
tmp <- unique(tmp)

pvalue.heatmap <- ggplot(tmp) + geom_tile(aes(x = sample_id, y = gene, fill = neg_log_10_p)) +
  theme_bw() +
  guides(fill = guide_colorbar(title = '-log10(p)')) +
  scale_fill_gradient(high = 'blue', low = 'white') +
  xlab("Sample") + ylab('') +
  facet_grid(is_differential~., scales = 'free_y', space = 'free') +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_text(colour = 'white')) +
  theme(axis.title.x = element_text(colour = 'white')) +
  theme(axis.ticks = element_line(colour = 'white', )) +
  theme(strip.background = element_blank()) + 
  theme(panel.spacing = unit(2, 'lines')) +
  theme(strip.text.y = element_blank()) +
  theme(legend.position = 'right')


f <- OUT
png(f, width = 7, height = 7, res = 300, units = 'in')
cowplot::plot_grid(counts.heatmap, pvalue.heatmap, ncol = 2, rel_widths = c(1, 0.3))
dev.off()
