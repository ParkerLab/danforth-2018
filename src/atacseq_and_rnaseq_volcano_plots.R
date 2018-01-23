library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)

args <- commandArgs(trailingOnly = T)
RNASEQ <- args[1]  # path to the RNA-seq DESeq2 results
ATACSEQ <- args[2]  # path to the ATAC-seq DESeq2 results
OUT <- args[3]

rnaseq <- read.table(RNASEQ, head = T, as.is = T, sep = '\t')
atacseq <- read.table(ATACSEQ, head = T, as.is = T, sep = '\t')

RNASEQ_LIMIT <- max(abs(rnaseq$log2FoldChange[is.finite(rnaseq$log2FoldChange)]))
ATACSEQ_LIMIT <- max(abs(atacseq$log2FoldChange[is.finite(atacseq$log2FoldChange)]))

# stack the plots, with legend in between
rnaseq.plot <- ggplot(rnaseq) +
  geom_point(aes(x = log2FoldChange, y = neg_log_10_p, alpha = significance, color = significance)) +
  scale_color_manual(values = c('n.s.' = 'black', 'sign.' = 'red')) +
  scale_alpha_manual(values = c('n.s.' = 0.2, 'sign.' = 1)) +
  geom_text_repel(aes(x = log2FoldChange, y = neg_log_10_p, label = gene_name),
                  data = filter(rnaseq, neg_log_10_p>=50 | gene_name %in% c('Ptf1a', "Gm13344")),
                  min.segment.length = unit(0.1, 'lines'),
                  point.padding = unit(0.5, 'lines'),
                  nudge_y = 2) +
  theme_bw() +
  xlab("Log2(Sd / WT), RNA-seq") +
  ylab("-log10(p)") + guides(color = F, alpha = F) +
  xlim(c(-1, 1) * RNASEQ_LIMIT)

atacseq.plot <- ggplot(atacseq) +
  geom_point(aes(x = log2FoldChange, y = neg_log_10_p, alpha = significance, color = significance)) +
  scale_color_manual(values = c('n.s.' = 'black', 'sign.' = 'red'), labels = c("FDR > 5%", "FDR < 5%"), guide = guide_legend(title=element_blank())) +
  scale_alpha_manual(values = c('n.s.' = 0.2, 'sign.' = 1), labels = c("FDR > 5%", "FDR < 5%"), guide = guide_legend(title=element_blank())) +
  geom_text_repel(aes(x = log2FoldChange, y = neg_log_10_p, label = gene),
                  data = filter(atacseq, neg_log_10_p>=5 | gene == 'Ptf1a'),
                  min.segment.length = unit(0.1, 'lines'),
                  point.padding = unit(0.5, 'lines'),
                  nudge_y = 1) +
  theme_bw() +
  xlab("Log2(Sd / WT), ATAC-seq") +
  ylab("-log10(p)") +
  xlim(c(-1, 1) * ATACSEQ_LIMIT) +
  theme(legend.position = 'right', legend.direction = 'vertical')

f <- OUT
png(f, height = 2.5, width = 7, units = 'in', res = 300)
plot_grid(atacseq.plot, rnaseq.plot, ncol = 2, rel_heights = c(1, 1), rel_widths = c(1.42, 1))
dev.off()
