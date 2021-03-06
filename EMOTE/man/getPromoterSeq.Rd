% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{getPromoterSeq}
\alias{getPromoterSeq}
\title{Retreive the promoter sequence around given genomic ranges}
\usage{
getPromoterSeq(fa, gr, upstream = 45, downstream = 5)
}
\arguments{
\item{fa}{a FaFile to search the sequence in}

\item{gr}{the genomic range used as reference to compute the promoter range}

\item{upstream}{number of required bases upstream}

\item{downstream}{number of required bases downstream}
}
\value{
a DNAStringSet with the promoters sequences
}
\description{
Retreive the promoter sequence around given genomic ranges
}
