library(optparse)


option_list <- list(
		    make_option(c('--deseq2_results'), action = 'store', type = 'character', help = '[Required] Path to the differential gene expression analysis results'),
		    make_option(c('--rnaenrich'), action = 'store', type = 'character', help = '[Required] Path to the RNA-Enrich.r script'),
		    make_option(c('--entrez_to_ensembl'), action = 'store', type = 'character', help = '[Required] Path to the Ensembl-Entrez translation file'),
		    make_option(c('--prefix'), action = 'store', type = 'character', help = '[Required] Prefix for the output')
		    )

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

library(dplyr)

# need to install some databases for RNA-Enrich
# source("http://bioconductor.org/biocLite.R")
# biocLite("GO.db")
# biocLite("org.Mm.eg.db")
# biocLite("KEGG.db")

set.seed(873459)
RNAENRICH <- opts$rnaenrich
source(RNAENRICH)

prefix <- function(x, p = opts$prefix) {
  return(paste(p, x, sep = '.'))
}

# RNA-Enrich requires Entrez gene IDs rather than Ensembl
# Will use entrez2ensembl (downloaded from NCBI) to convert gene IDs
ensembl2entrez <- read.table(gzfile(opts$entrez_to_ensembl), head = T, as.is = T, sep = '\t', comment.char = '')
ensembl2entrez <- dplyr::rename(ensembl2entrez, tax_id = X.tax_id) %>%
  filter(tax_id == 10090) %>%
  dplyr::rename(entrez = GeneID,
         ensembl_gene_id = Ensembl_gene_identifier) %>%
  dplyr::select(entrez, ensembl_gene_id) %>%
  unique()


# Load differential gene expression results and add in the entrez IDs
deseq <- read.table(opts$deseq2_results, head = T, as.is = T, sep = '\t')
deseq <- left_join(deseq, ensembl2entrez)

# Run RNA-Enrich
RNAE.KEGG <- rna_enrich(sigvals = deseq$pvalue, geneids = deseq$entrez, species = 'mmu', min.g = 5, max.g = 500, avg_readcount = deseq$baseMean,
                        database = 'KEGG',
                        plot_file=prefix("rnaenrich_plot.KEGG.jpg"), plot_height=480, plot_width=480, results_file=prefix("rnaenrich_results.KEGG.txt"))

sessionInfo()
