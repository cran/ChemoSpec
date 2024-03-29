% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/averageReplicates.R
\name{averageReplicates}
\alias{averageReplicates}
\title{Average Replicates in a Spectra Object}
\usage{
averageReplicates(spectra, uniq)
}
\arguments{
\item{spectra}{An object of S3 class \code{\link[=Spectra]{Spectra()}}.}

\item{uniq}{Character.  A character vector containing strings representing unique sample
identifiers.  The sample names will be searched for these strings, and all samples
matching a given string will be averaged and put into a new Spectra object.  For
example, consider the case where samples are named
\code{S_1_01, S_1_02, ..., S_2_01, S_2_02, ...}
where \code{_01} and so forth signifies replicates of a particular sample.
With \code{uniq = c("S_1", "S_2")} all \code{S_1} replicates will be averaged and
all \code{S_2} replicates will be averaged.  N.B. the strings will be used as
regex pattern and grepped.}
}
\value{
An object of S3 class \code{\link{Spectra}}.
}
\description{
Average the replicates in a \code{\link{Spectra}} object and return a new Spectra object
with fewer samples.  One should probably not do this until each individual sample has been
visualized for quality control, in case it is a potential outlier.
}
\examples{
data(SrE.IR)
averaged <- averageReplicates(SrE.IR, uniq = c("EPO", "OO", "adSrE", "pSrE"))
sumSpectra(SrE.IR)
sumSpectra(averaged)

}
\author{
Bryan A. Hanson (DePauw University).
}
\keyword{utilities}
