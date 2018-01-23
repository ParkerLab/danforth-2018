library(ggplot2)
library(dplyr)

#ROOT <- '/lab/work/porchard/danforth.final'
ROOT <- Sys.getenv('DANFORTH_HOME')
h3k4me1_signal <- file.path(ROOT, "figures/h3k4me1_barplot/all_signals.bed")
ROADMAP_SAMPLE_INFO <- file.path(ROOT, 'control', 'work', 'roadmap_chipseq', 'setup', 'sample_info.txt')

sample_info <- read.table(ROADMAP_SAMPLE_INFO, head = T, as.is = T, sep = '\t')

h3k4 <- read.table(h3k4me1_signal, head = F, as.is = T, sep = '\t')
colnames(h3k4) <- c("gsm", "chrom", "start", "end", "signal")
h3k4 <- left_join(h3k4, sample_info)
h3k4$tissue[is.na(h3k4$tissue)] <- h3k4$gsm[is.na(h3k4$tissue)]
h3k4 <- dplyr::select(h3k4, tissue, signal)
h3k4 <- h3k4[order(h3k4$signal),]

tissue_counts <- table(h3k4$tissue)
tissue_counts <- tissue_counts[tissue_counts>1]
for(tissue in names(tissue_counts)) {
  h3k4$tissue[h3k4$tissue==tissue] <- paste(h3k4$tissue[h3k4$tissue==tissue], 1:nrow(h3k4[h3k4$tissue==tissue,]), sep = '_')
}

h3k4$tissue <- gsub(' ', '', h3k4$tissue)
h3k4$tissue <- factor(h3k4$tissue, levels = h3k4$tissue, ordered = T)

colors <- rep('grey', nrow(h3k4))
names(colors) <- unique(h3k4$tissue)
colors[grep('DevelopingPancreas', names(colors))] <- 'darkgreen'

p <- ggplot(h3k4) + geom_bar(aes(x = tissue, y = signal, fill = tissue), stat = 'identity') +
  scale_fill_manual(values = colors) +
  theme_bw() + coord_flip() +
  guides(fill = F) + ylab("H3K4me1 signal in region\northologous to Gm13344 promoter peak") +
  xlab("Tissue")

f <- file.path(dirname(h3k4me1_signal), 'h3k4me1_signals.png')
png(f, height = 5, width = 6, units = 'in', res = 300)
print(p)
dev.off()
