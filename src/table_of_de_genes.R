library(tidyr)
library(dplyr)

# compile the list of DE genes

ROOT <- Sys.getenv('DANFORTH_HOME')

deseq <- read.table(file.path(ROOT, 'work', 'differential_gene_expression', 'results', 'deseq2', 'differential_gene_expression.deseq2_results.txt'), head = T, as.is = T, sep = '\t')
stopifnot(is.numeric(deseq$padj))

output <- deseq[!is.na(deseq$pad) & deseq$padj<=0.05,c('gene_name', 'ensembl_gene_id', 'log2FoldChange', 'lfcMLE', "pvalue", "padj")]
colnames(output) <- c('Gene', 'Ensembl gene id', 'DESeq2 log2FoldChange', 'DESeq2 lfcMLE', 'unadjusted p-value', 'adjusted p-value')

write.table(output, file = 'DE_genes.csv', col.names = T, row.names = F, sep = ',', quote = F)
