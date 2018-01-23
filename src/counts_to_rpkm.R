counts_to_rpkm <- function(counts_vector, lengths) {
  # lengths indicates the lengths (in bp) of each of the corresponding regions
  stopifnot(length(counts_vector) == length(lengths))
  
  counts_sum <- sum(counts_vector)
  
  # to reduce the probability of integer overflow,
  # enforce a certain order of operations
  rpkm <- counts_vector * (((10^9) / (lengths)) / counts_sum)
  
  return(rpkm)
}