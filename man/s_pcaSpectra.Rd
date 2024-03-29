% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s_pcaSpectra.R
\name{s_pcaSpectra}
\alias{s_pcaSpectra}
\title{Sparse PCA of Spectra Objects}
\usage{
s_pcaSpectra(spectra, choice = "noscale", K = 3, para = rep(0.5, K), ...)
}
\arguments{
\item{spectra}{An object of S3 class \code{\link[=Spectra]{Spectra()}}.}

\item{choice}{A character string indicating the choice of scaling.  One of
\code{c("noscale"}, \code{"autoscale"}, \code{"Pareto")}. \code{"autoscale"}
is called "standard normal variate" or "correlation matrix PCA" in some literature.}

\item{K}{Integer. The number of components desired.}

\item{para}{A vector of \code{length(K)} giving the tuning parameters.}

\item{...}{Other parameters to be passed to \code{\link[elasticnet]{arrayspc}}.}
}
\value{
An object of class \code{prcomp} and \code{converted_from_arrayspc},
which includes a list
element called \code{$method}, a character string describing the
pre-processing carried out and the type of PCA performed (used to annotate
plots).  A check is carried out to see if the computation was successful
and a warning issued if it failed.
}
\description{
A wrapper which carries out sparse PCA analysis on a
\code{\link{Spectra}} object.  The user can select various options for
scaling.  There is no normalization by rows - do this manually using
\code{\link{normSpectra}}. The data will be centered, as is required by PCA.
}
\details{
The scale choice \code{autoscale} scales the columns by their standard
deviation.  \code{Pareto} scales by the square root of the standard
deviation.
}
\examples{
\dontrun{
# This example assumes the graphics output is set to ggplot2 (see ?GraphicsOptions).
library("ggplot2")
data(SrE.NMR)
pca <- s_pcaSpectra(SrE.NMR)

p1 <- plotScree(pca)
p1

p2 <- plotScores(SrE.NMR, pca, pcs = c(1, 2), ellipse = "cls", tol = 0.05)
p2 <- p2 + ggtitle("Scores: SrE NMR Data")
p2

p3 <- plotLoadings(SrE.NMR, pca, loads = 1:2, ref = 1)
p3 <- p3 + ggtitle("Loadings: SrE NMR Data")
p3
}

}
\references{
H. Zou, T. Hastie and R. Tibshirani "Sparse Principal Components Analysis"
\emph{J. Comp. Stat. Graphics} vol. 15 no. 2 pgs. 265-286 (2006).
}
\seealso{
\code{\link[elasticnet]{arrayspc}} for the underlying function,
\code{\link{c_pcaSpectra}} for classical PCA calculations,
\code{\link{r_pcaSpectra}} for robust PCA calculations,
\code{\link{irlba_pcaSpectra}} for PCA via the IRLBA algorithm.
Additional documentation at \url{https://bryanhanson.github.io/ChemoSpec/}

For displaying the results, \code{\link[ChemoSpecUtils:plotScree]{ChemoSpecUtils::plotScree()}}, \code{\link[ChemoSpecUtils:plotScores]{ChemoSpecUtils::plotScores()}}, \code{\link[=plotLoadings]{plotLoadings()}}, \code{\link[=plot2Loadings]{plot2Loadings()}}, \code{\link[=sPlotSpectra]{sPlotSpectra()}}.
}
\author{
Bryan A. Hanson (DePauw University).
}
\keyword{multivariate}
