% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{match_pwms}
\alias{match_pwms}
\alias{match_pwms,PFMatrix,DNAString-method}
\alias{match_pwms,PFMatrix,DNAStringSet-method}
\alias{match_pwms,PFMatrix,GenomicRanges-method}
\alias{match_pwms,PFMatrix,character-method}
\alias{match_pwms,PFMatrixList,DNAString-method}
\alias{match_pwms,PFMatrixList,DNAStringSet-method}
\alias{match_pwms,PFMatrixList,GenomicRanges-method}
\alias{match_pwms,PFMatrixList,character-method}
\alias{match_pwms,PWMatrix,DNAString-method}
\alias{match_pwms,PWMatrix,DNAStringSet-method}
\alias{match_pwms,PWMatrix,GenomicRanges-method}
\alias{match_pwms,PWMatrix,character-method}
\alias{match_pwms,PWMatrixList,DNAString-method}
\alias{match_pwms,PWMatrixList,DNAStringSet-method}
\alias{match_pwms,PWMatrixList,GenomicRanges-method}
\alias{match_pwms,PWMatrixList,character-method}
\title{match_pwms}
\usage{
match_pwms(pwms, subject, ...)

\S4method{match_pwms}{PWMatrixList,DNAStringSet}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrixList,character}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrixList,DNAString}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrixList,GenomicRanges}(pwms, subject,
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  bg = NULL, out = c("match", "scores", "positions"), p.cutoff = 5e-05,
  w = 7)

\S4method{match_pwms}{PFMatrixList,DNAStringSet}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrixList,character}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrixList,DNAString}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrixList,GenomicRanges}(pwms, subject,
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  bg = NULL, out = c("match", "scores", "positions"), p.cutoff = 5e-05,
  w = 7)

\S4method{match_pwms}{PWMatrix,DNAStringSet}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrix,character}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrix,DNAString}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PWMatrix,GenomicRanges}(pwms, subject,
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  bg = NULL, out = c("match", "scores", "positions"), p.cutoff = 5e-05,
  w = 7)

\S4method{match_pwms}{PFMatrix,DNAStringSet}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrix,character}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrix,DNAString}(pwms, subject, bg = NULL,
  out = c("match", "scores", "positions"), p.cutoff = 5e-05, w = 7)

\S4method{match_pwms}{PFMatrix,GenomicRanges}(pwms, subject,
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  bg = NULL, out = c("match", "scores", "positions"), p.cutoff = 5e-05,
  w = 7)
}
\arguments{
\item{pwms}{either \code{\link[TFBSTools]{PFMatrix}}, \code{\link[TFBSTools]{PFMatrixList}},
\code{\link[TFBSTools]{PWMatrix}}, \code{\link[TFBSTools]{PWMatrixList}}}

\item{subject}{either \code{\link[GenomicRanges]{GenomicRanges}}, \code{\link[Biostrings]{DNAStringSet}},
\code{\link[Biostrings]{DNAString}}, or character vector}

\item{bg}{background nucleotide frequencies. if not provided, computed from subject}

\item{out}{what to return? see details}

\item{p.cutoff}{p-value cutoff for returning motifs}

\item{w}{parameter controlling size of window for filtration; default is 7}

\item{genome}{BSgenome object, only used if subect is \code{\link[GenomicRanges]{GenomicRanges}}}
}
\description{
Find pwm matches
}
\details{
Can either return a sparse matrix with values set to 1 for a match (if return == "match"), a sparse matrix
with values set to the max motif score in each sequence (but zero for sequences with no score above minimum p value threshold),
or \code{\link[GenomicRanges]{GenomicRanges}} if positions
}
\section{Methods (by class)}{
\itemize{
\item \code{pwms = PWMatrixList,subject = DNAStringSet}: 

\item \code{pwms = PWMatrixList,subject = character}: 

\item \code{pwms = PWMatrixList,subject = DNAString}: 

\item \code{pwms = PWMatrixList,subject = GenomicRanges}: 

\item \code{pwms = PFMatrixList,subject = DNAStringSet}: 

\item \code{pwms = PFMatrixList,subject = character}: 

\item \code{pwms = PFMatrixList,subject = DNAString}: 

\item \code{pwms = PFMatrixList,subject = GenomicRanges}: 

\item \code{pwms = PWMatrix,subject = DNAStringSet}: 

\item \code{pwms = PWMatrix,subject = character}: 

\item \code{pwms = PWMatrix,subject = DNAString}: 

\item \code{pwms = PWMatrix,subject = GenomicRanges}: 

\item \code{pwms = PFMatrix,subject = DNAStringSet}: 

\item \code{pwms = PFMatrix,subject = character}: 

\item \code{pwms = PFMatrix,subject = DNAString}: 

\item \code{pwms = PFMatrix,subject = GenomicRanges}: 
}}

