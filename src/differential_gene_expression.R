library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(optparse)
library(biomaRt) # for converting ENSEMBL gene names


option_list <- list(
  make_option(c("--sample_info"), action = 'store', type = 'character', help = '[Required] Path to the sample_info.txt file'),
  make_option(c("--counts_dir"), action = 'store', type = 'character', help = '[Required] Path to the directory containing the gene counts for each sample. Each file in the directory should be named {seqcore_id}.geneCounts.txt.gz'),
  make_option(c("--covariates"), action = 'store', type = 'character', help = '[Optional] Comma-separated list of covariates. May include "somites", "rin", "batch"'),
  make_option(c("--prefix"), action = 'store', type = 'character', help = '[Required] Prefix for the output file names')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

set.seed(627677)

prefix <- function(x, p = opts$prefix) {
  return(paste(p, x, sep = '.'))
}


# Load in the sample information and filter down to what we actually need
sample_info <- read.table(opts$sample_info, head = T, as.is = T, sep = '\t', comment.char = "")
sample_info <- dplyr::filter(sample_info, experiment_type == 'RNAseq') %>% 
  group_by(seqcore_id, genotype) %>%
  dplyr::summarize(rin = mean(rin), somites = mean(somites), batch = paste(unique(date), collapse = ':')) %>%
  as.data.frame()
sample_info$batch <- factor(sample_info$batch)




# load in the gene expression data
count_files <- list.files(opts$counts_dir, pattern = "geneCounts.txt.gz", full.names = T)

parse_file_name <- function(f) {
  regex <- '^(\\d+)\\.geneCounts\\.txt\\.gz$'
  seqcore_id <- gsub(regex, '\\1', basename(f))
  x <- c('seqcore_id' =  seqcore_id)
  return(x)
}

load_count_file <- function(f) {
  counts <- read.table(gzfile(f), header = F, as.is = T, sep = '\t')
  colnames(counts) <- c('gene', 'count')
  counts$seqcore_id <- parse_file_name(f)['seqcore_id']
  
  # There are a few items that do not correspond to genes. Kick these out:
  counts <- counts[grep('ENS', counts$gene),]
  
  return(counts)
}

counts <- bind_rows(lapply(count_files, load_count_file))
counts <- tidyr::spread(counts, key = seqcore_id, value = count, fill = NA)
rownames(counts) <- counts$gene
counts <- dplyr::select(counts, -gene)

# sanity check
if (sum(is.na(counts)) != 0) {
  warning('Missing data present in the counts data.frame; exiting')
  quit(save = 'no')
}

# filter out the genes that show 0 expression in all experiments
counts <- counts[apply(counts, 1, max) > 0,]


# prepare the sample information for DESeq2
exp_info <- sample_info
rownames(exp_info) <- exp_info$seqcore_id
exp_info <- dplyr::select(exp_info, -seqcore_id)
exp_info$genotype <- factor(exp_info$genotype, levels = c("WT", "Sd")) # set the levels such that positive log2FC means higher in Sd mice


# run DESeq2
des <- formula('~ genotype')
if (!is.null(opts$covariates)) {
  covariates <- strsplit(opts$covariates, ',')[[1]]
  stopifnot(all(covariates %in% colnames(exp_info)))
  des <- formula(paste('~', paste(covariates, collapse = " + "), '+', 'genotype'))
}


dds <- DESeqDataSetFromMatrix(countData = counts[,rownames(exp_info)], colData = exp_info, design = des)
featureData <- data.frame(gene=rownames(counts))
(mcols(dds) <- DataFrame(mcols(dds), featureData))

dds <- DESeq(dds, betaPrior = T)


# extract results and add some variables for plotting
res <- results(dds,  addMLE=T, alpha = 0.05)

res.df <- as.data.frame(res)
res.df$neg_log_10_p <- -1*log10(res.df$pvalue)
res.df$significance <- "n.s."
res.df$significance[!is.na(res.df$padj) & res.df$padj <= 0.05] <- "sign."
res.df$higher_in <- "neither"
res.df$higher_in[res.df$significance=="sign." & res.df$log2FoldChange>0] <- "Sd"
res.df$higher_in[res.df$significance=="sign." & res.df$log2FoldChange<0] <- "WT"

# Use biomart to convert ensembl gene ids to human-readable gene names
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
filter <- c("ensembl_gene_id")
attribute <- c("ensembl_gene_id", "external_gene_name")
values <- gsub("\\..*", "", rownames(res.df))
translations <- getBM(attributes = attribute, values = values, filters = filter, mart = ensembl)
translations <- dplyr::rename(translations, gene_name = external_gene_name)
res.df$ensembl_gene_id <- gsub("\\..*", "", rownames(res.df))
res.df <- left_join(res.df, translations)
# if no human-readable name found, keep the ensembl gene ID as the gene name
res.df$gene_name[is.na(res.df$gene_name)] <- res.df$ensembl_gene_id[is.na(res.df$gene_name)]


# make MA plot
png(prefix("MA_plot.png"), width = 4, height = 4, units = "in", res = 300)
plotMA(res, main="Shrunken LFC")
dev.off()


# look at the distribution of p-values
png(prefix("p_value_histogram.png"), width = 5, height = 5, units = "in", res = 600)
hist(res.df$pvalue, main = "P value histogram", xlab="p.value")
dev.off()

# create volcano plot
p <- ggplot(res.df, aes(x = log2FoldChange, y = neg_log_10_p, color = significance, alpha = significance)) + 
  theme_bw() + 
  geom_point(stroke = 0) + xlab("log2(Sd / WT)") + ylab("-log10(p)") + 
  xlim(c(-1,1) * max(abs(res.df$log2FoldChange))) +
  scale_color_manual(values = c("n.s." = "black", "sign." = "red"), guide = guide_legend(title = "Significance")) +
  scale_alpha_manual(values = c("n.s." = 0.1, "sign." = 1), guide = guide_legend(title = "Significance"))
png(prefix("volcano_plot.png"), width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()

# annotate volcano plot with the gene names of significantly DE genes
p <- p + geom_text_repel(aes(label = gene_name), data = subset(res.df, padj < 0.05 & !is.na(padj)), show.legend = F, size = 2, segment.size = 0.4, point.padding = unit(0.2, "lines"), nudge_y = 0.07)
png(prefix("volcano_plot.all_sig_annotated.png"), width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()


# plot the number of differential genes
res.counts <- res.df %>% group_by(higher_in) %>% dplyr::summarize(count=n())
res.counts <- res.counts[order(res.counts$higher_in),]
p <- ggplot(res.df) + 
  geom_bar(aes(x = higher_in, group = higher_in, fill = higher_in), stat = "count") +
  theme_bw() + 
  ylab("Number of genes") + 
  xlab("Genotype with greater signal") +
  geom_text(aes(x = higher_in, y = count, label = count), data = res.counts, vjust = 0) +
  scale_fill_manual(values = c("neither" = "grey", "Sd" = "red", "WT" = "black"), guide = guide_legend(title = "Higher in"))
png(prefix("counts.png"), width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()


# write out the results
write.table(x = res.df, file = prefix("deseq2_results.txt"), quote = F, sep = "\t", col.names = T, row.names = F)

sessionInfo()
