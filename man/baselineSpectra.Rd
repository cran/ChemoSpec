% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baselineSpectra.R
\name{baselineSpectra}
\alias{baselineSpectra}
\title{Baseline Correction of a Spectra Object}
\usage{
baselineSpectra(spectra, int = TRUE, retC = FALSE, show = 1, ...)
}
\arguments{
\item{spectra}{An object of S3 class \code{\link[=Spectra]{Spectra()}}.}

\item{int}{Logical; if \code{TRUE}, do the correction interactively using
widgets.  No results are saved.  Use this for inspection and exploration
only.  Automatically overriden to \code{FALSE} if \code{interactive()} returns
\code{FALSE}.  This is necessary so that plots appear in vignettes etc.}

\item{retC}{Logical: shall the baseline-corrected spectra be returned in the
\code{Spectra} object?}

\item{show}{Integer. A (single) sample number for which you wish
to see the results of the baseline correction.  By "sample number"
we mean the rows in the \code{spectra$data} matrix.  To find a specific
sample type \code{spectra$names} to see which row contains that sample.}

\item{\dots}{Other arguments passed downstream.  The relevant ones can be
found in \code{\link[baseline]{baseline}}.  Be sure to pay attention to
argument \code{method} as you will probably want to use it.  You can also
use \code{method = "linear"} for a simple linear fit, see Details.}
}
\value{
If \code{int = TRUE}, an interactive plot is created.  If \code{int
= FALSE} and \code{retC = FALSE}, an object of class \code{baseline} is
returned (see \code{\link[baseline]{baseline-class}}).  If \code{int =
FALSE} and \code{retC = TRUE}, a \code{Spectra} object containing the
corrected spectra is returned.  In these latter two cases plots are also
drawn.
}
\description{
This function mostly wraps functions in package \pkg{baseline} which
carries out a variety of baseline correction routines.  A simple linear
correction method is also available.
}
\details{
In plots using methods from the baseline package, the x axis ticks give the
data point index, not the original values from your data. Note that you
cannot zoom the non-interactive display of corrected spectra because the
underlying function hardwires the display.  Try the interactive version
instead (\code{int = TRUE}), or use \code{\link{plotSpectra}} on the
corrected data.
In addition to the methods provided by \code{baseline}, you can also use
\code{method = "linear"}.  This correction is handled locally, and is very
simple: a line is drawn from the first data point to the last, and this
becomes the new baseline.  This is most suitable for cases in which the
baseline rises or falls steadily, as is often seen in chromatograms.
}
\examples{
# You need to install package "baseline" for this example
if (requireNamespace("baseline", quietly = TRUE)) {
  data(SrE.IR)
  temp <- baselineSpectra(SrE.IR, int = FALSE, method = "modpolyfit")
}

}
\seealso{
Additional documentation at \url{https://bryanhanson.github.io/ChemoSpec/}
}
\author{
Bryan A. Hanson (DePauw University).
}
\keyword{hplot}
