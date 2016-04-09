
#' match_pwms
#'
#' Find pwm matches
#' @import Biostrings
#' @import TFBSTools
#' @import Matrix
#'@export
setGeneric("match_pwms", function(pwms, subject,...) standardGeneric("match_pwms"))

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAStringSet"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions")){
            out = match.arg(out)
            motif_mats <- lapply(pwms, as.matrix)
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
            })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "character"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w =7){
            out = match.arg(out)

            motif_mats <- lapply(pwms, as.matrix)

            if (is.null(bg)){
              bg <- get_nuc_freqs(DNAStringSet(subject))
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }

            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAString"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){
            out = match.arg(out)
            motif_mats <- lapply(pwms, as.matrix)
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "GenomicRanges"),
          function(pwms, subject, genome, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){

            out = match.arg(out)
            seqs <- getSeq(genome, subject)
            motif_mats <- lapply(pwms, as.matrix)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })


### PFMatrixList ---------------------------------------------------------------


#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "DNAStringSet"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions")){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            seqs <- as.character(subject)

            if (out == "match"){
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "character"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(DNAStringSet(subject))
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "DNAString"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "GenomicRanges"),
          function(pwms, subject, genome, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){

            out = match.arg(out)
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "DNAStringSet"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions")){
            out = match.arg(out)
            motif_mats <- list(as.matrix(pwm))
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "character"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w =7){
            out = match.arg(out)

            motif_mats <- list(as.matrix(pwm))

            if (is.null(bg)){
              bg <- get_nuc_freqs(DNAStringSet(subject))
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }

            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,subject,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "DNAString"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){
            out = match.arg(out)
            motif_mats <- list(as.matrix(pwm))
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "GenomicRanges"),
          function(pwms, subject, genome, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){

            out = match.arg(out)
            seqs <- getSeq(genome, subject)
            motif_mats <- list(as.matrix(pwm))
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

# Single PFM -------------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "DNAStringSet"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions")){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "character"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(DNAStringSet(subject))
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,subject,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "DNAString"),
          function(pwms, subject, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "GenomicRanges"),
          function(pwms, subject, genome, p.cutoff = 0.00005, bg = NULL, out = c("match","scores","positions"), w = 7){

            out = match.arg(out)
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })


