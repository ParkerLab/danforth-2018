library("QoRTs")

# Get size factors from qorts output
# First argument should be the qorts path (i.e., to the main directory)
# Second argument should be the decoder file
# Third argument should be the name of the output file

args <- commandArgs(trailingOnly = T)

decoder <- read.table(args[2], head = T, as.is = T, sep = '\t')
decoder[,1] <- as.character(decoder[,1])

res <- read.qc.results.data(args[1], decoder = decoder, calc.DESeq2 = T)
get.size.factors(res, outfile = args[3])

sessionInfo()
