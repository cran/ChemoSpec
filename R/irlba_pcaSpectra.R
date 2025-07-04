#' IRLBA PCA of Spectra Objects
#'
#' A wrapper which carries out IRLBA PCA analysis on a
#' \code{\link{Spectra}} object.  The user can select various options for
#' scaling.  There is no normalization by rows - do this manually using
#' \code{\link{normSpectra}}. The data can be supplied already centered
#' if desired.
#'
#' The scale choice \code{autoscale} scales the columns by their standard
#' deviation.  \code{Pareto} scales by the square root of the standard
#' deviation.
#'
#' @param spectra `r .writeDoc_Spectra1()`
#'
#' @param choice A character string indicating the choice of scaling.  One of
#' \code{c("noscale"}, \code{"autoscale"}, \code{"Pareto")}. \code{"autoscale"}
#' is called "standard normal variate" or "correlation matrix PCA" in some literature.
#'
#' @param n Integer. The number of components desired.
#'
#' @param center Logical.  Should the data be centered?  Data must be centered
#' for PCA, either before arriving here or via this argument.
#'
#' @param ... Other parameters to be passed to \code{\link[irlba]{irlba}}.
#'
#' @return A modified object of class \code{prcomp} and \code{computed_via_irlba},
#' which includes a list element called \code{$method}, a character string describing the
#' pre-processing carried out and the type of PCA performed (used to annotate
#' plots).
#'
#' @author `r .writeDoc_Authors("BH")`
#'
#' @seealso \code{\link[irlba]{prcomp_irlba}} for the underlying function,
#' \code{\link{c_pcaSpectra}} for classical PCA calculations,
#' \code{\link{r_pcaSpectra}} for robust PCA calculations,
#' \code{\link{s_pcaSpectra}} for sparse PCA calculations.
#' Additional documentation at \url{https://bryanhanson.github.io/ChemoSpec/}
#'
#' `r .writeDoc_LinksShowPCAResults()`
#'
#' @references J. Baglama and L. Reichel, "Augmented Implicitly Restarted Lanczos
#' Bidiagonalization Methods"  \emph{SIAM J. Sci. Comput.} (2005).
#'
#' @keywords multivariate
#'
#' @examples
#' \dontrun{
#' # This example assumes the graphics output is set to ggplot2 (see ?GraphicsOptions).
#' library("ggplot2")
#' data(SrE.NMR)
#' pca <- irlba_pcaSpectra(SrE.NMR)
#'
#' p1 <- plotScree(pca)
#' p1
#'
#' p2 <- plotScores(SrE.NMR, pca, pcs = c(1, 2), ellipse = "cls", tol = 0.05)
#' p2 <- p2 + ggtitle("Scores: SrE NMR Data")
#' p2
#'
#' p3 <- plotLoadings(SrE.NMR, pca, loads = 1:2, ref = 1)
#' p3 <- p3 + ggtitle("Loadings: SrE NMR Data")
#' p3
#' }
#'
#' @export irlba_pcaSpectra
#' @importFrom stats sd
#'
irlba_pcaSpectra <- function(spectra, choice = "noscale", n = 3, center = TRUE, ...) {
  .chkArgs(mode = 11L)
  chkSpectra(spectra)

  if (.chkReqPkgs("irlba")) {
    choices <- c("noscale", "autoscale", "Pareto") # trap for invalid scaling method
    check <- choice %in% choices
    if (!check) stop("The choice of scaling parameter was invalid")

    # irlba::prcomp_irlba does its own scaling, so we must plan for that

    if (choice == "noscale") {
      X <- scale(spectra$data, center = center, scale = FALSE)
    }

    if (choice == "Pareto") {
      col.sd <- apply(spectra$data, 2, sd)
      X <- scale(spectra$data, center = center, scale = col.sd^0.5)
    }

    if (choice == "autoscale") {
      col.sd <- apply(spectra$data, 2, sd)
      X <- scale(spectra$data, center = center, scale = col.sd)
    }

    pca <- irlba::prcomp_irlba(x = X, n = n, center = FALSE, scale. = FALSE, ...)

    # Modify the class
    pca$method <- paste("centered/", choice, "/", "irlba", sep = "")
    class(pca) <- c("computed_via_irlba", "prcomp")
    pca
  }
}
