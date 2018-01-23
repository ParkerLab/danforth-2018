library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggrepel)

args <- commandArgs(trailingOnly = T)
GO <- args[1]  # path to the GO analysis results
ENTREZ_TO_ENSEMBL <- args[2]
OUT <- args[3]


kegg <- read.table(GO, head = T, sep = '\t', as.is = T, fill = T, quote = "")

entrez2ensembl_mouse <- read.table(gzfile(ENTREZ_TO_ENSEMBL), head = T, as.is = T, sep = '\t', comment.char = '')
entrez2ensembl_mouse <- entrez2ensembl_mouse %>%
  dplyr::rename(tax_id = X.tax_id) %>%
  filter(tax_id == 10090) %>%
  dplyr::select(GeneID, Ensembl_gene_identifier) %>%
  dplyr::rename(entrez=GeneID, ensembl = Ensembl_gene_identifier) %>%
  unique()

# translate the sig.genes to ensembl genes...

singleEntrez2Ensembl <- function(entrez_id) {
  ensembl_id <- NA
  if (entrez_id %in% entrez2ensembl_mouse$entrez) {
    ensembl_id <- entrez2ensembl_mouse$ensembl[entrez2ensembl_mouse$entrez==entrez_id]
  }
  
  if (length(ensembl_id) !=1 ) {
    return(NA)
  } else {
    return(ensembl_id)
  }
}

entrez2ensembl <- function(x) {
  # x is a list of entrez gene ids
  return(sapply(x, singleEntrez2Ensembl))
}


kegg$sig.genes.ensembl <- sapply(kegg$sig.genes, function(x) {
  y <- strsplit(x, ", ")[[1]]
  y <- entrez2ensembl(y)
  y <- paste(y, collapse = ', ')
  return(y)
})

kegg.for_plot <- kegg
label_colors <- rep('black', nrow(kegg.for_plot))
label_colors[kegg.for_plot$Concept.name=='Hedgehog signaling pathway'] <- 'red'
kegg.for_plot$Concept.name <- gsub(" \\(CAMs\\)", "", kegg.for_plot$Concept.name)
kegg.for_plot$Concept.name <- gsub(' ', "\n", kegg.for_plot$Concept.name)
kegg.for_plot$Concept.name <- factor(kegg.for_plot$Concept.name, levels = kegg.for_plot$Concept.name, ordered = T)

concepts_to_label <- as.character(kegg.for_plot$Concept.name[kegg.for_plot$FDR<=0.05])
concepts_to_label <- c(concepts_to_label, grep("wnt", kegg.for_plot$Concept.name, ignore.case = T, value = T))
concepts_to_label <- c(concepts_to_label, grep("retin", kegg.for_plot$Concept.name, ignore.case = T, value = T))
concepts_to_label <- c(concepts_to_label, grep("tgf", kegg.for_plot$Concept.name, ignore.case = T, value = T))
concepts_to_label <- c(concepts_to_label, grep("notch", kegg.for_plot$Concept.name, ignore.case = T, value = T))
concepts_to_label <- c(concepts_to_label, grep("signaling\npathway", kegg.for_plot$Concept.name, ignore.case = T, value = T))
concepts_to_label <- unique(concepts_to_label)
concepts_to_label <- concepts_to_label[grep("^T\ncell", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^B\ncell", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Chemokine", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Insulin", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Adipocytokine", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Calcium", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^p53", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^RIG-I", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^PPAR", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Fc\nepsilon", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^NOD-like", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Toll-like", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^GnRH", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Hepatitis", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Huntington's", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Leukocyte", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Neurotrophin", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^Jak-STAT", concepts_to_label, invert = T)]
concepts_to_label <- concepts_to_label[grep("^VEGF", concepts_to_label, invert = T)]

kegg.plot <- ggplot(kegg.for_plot %>% mutate(significant = ifelse(FDR <= 0.05, 'sign.', 'n.s.'))) +
  geom_point(aes(x = odds.ratio, y = -log10(p.value), color = significant, alpha = log2(n.genes))) +
  theme_bw() +
  xlab("Odds ratio") +
  ylab("-log10(p)") +
  scale_color_manual(values = c('sign.' = 'red', 'n.s.' = 'black'), labels = c('FDR >= 5%', 'FDR < 5%')) +
  geom_text_repel(aes(x = odds.ratio, y = -log10(p.value), label = Concept.name), data = kegg.for_plot %>% filter(Concept.name %in% concepts_to_label), size = 3, box.padding = unit(0.2, 'lines')) +
  guides(alpha = guide_legend(title = 'log2(number of genes)'), color = guide_legend(title = ''))


f <- OUT
png(f, width = 8, height = 6, res = 300, units = 'in')
print(kegg.plot)
dev.off()
