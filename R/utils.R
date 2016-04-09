# Not exported
get_nuc_freqs <- function(seqs){
  #get nucleotide frequencies
  nucFreqs <- colSums(Biostrings::letterFrequency(seqs, c("A","C","G","T")))
  nucFreqs <- nucFreqs/sum(nucFreqs)
  return(nucFreqs)
}

# adjust_bg <- function(pwm){
#   lapply(pwm, function(x){
#     out = TFBSTools::as.matrix(x)
#     norm_mat = matrix(log(TFBSTools::bg(x)) - log(bg), nrow = 4, ncol = ncol(out))
#     return(out - norm_mat)
#   })
# }

