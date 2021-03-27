## ----SetUp, echo = FALSE, eval = TRUE, results = "hide"-----------------------

# R options & configuration:
set.seed(9)
rm(list = ls())
suppressMessages(library("knitr"))
suppressMessages(library("ChemoSpec"))
suppressMessages(library("mclust"))
suppressMessages(library("RColorBrewer"))

desc <- packageDescription("ChemoSpec")

# Stuff specifically for knitr:

# Create a temp bib file w/citations of installed pkgs, on the fly
knitr::write_bib(c("knitr", "mclust", "baseline", "hyperSpec"), file = "manuals.bib", prefix = "R_")

# Hook for figure size control
knitr::opts_hooks$set(sq.fig = function(options) {

  if (isFALSE(options$sq.fig)) {
    # custom fig dimensions given, use w/o further delay
    if ((!is.null(options$fig.width)) & (!is.null(options$fig.height))) return(options)
    # otherwise set the default wide aspect ratio
    if ((is.null(options$fig.width)) & (is.null(options$fig.height))) {
      options$fig.width = 6
      options$fig.height = 3.5
    }
  }

  if (isTRUE(options$sq.fig)) {
    options$fig.width = 5
    options$fig.height = 5
  }
  options
})

# choices here are designed to work with the hook
knitr::opts_chunk$set(fig.align = "center", sq.fig = FALSE,
  fig.width = NULL, fig.height = NULL, out.width = "80%")


## ----eval = FALSE-------------------------------------------------------------
#  source("My_First_ChemoSpec.R")

## ----workflow, out.width = "80%", echo = FALSE, fig.cap = "A typical workflow.  For a given data set, some steps may be omitted and the order changed.  That is part of what is meant by exploratory data analysis!"----
knitr::include_graphics("MetabPreProcess.png")

## ----results = "hide", eval = FALSE-------------------------------------------
#  ssp <- files2SpectraObject(
#    gr.crit = c("sspA", "sspB"),
#    gr.cols = c("red", "blue"),
#    freq.unit = "ppm",
#    int.unit = "peak intensity",
#    descrip = "Subspecies Study",
#    out.file = "subsp")

## ----results = "hide", eval = FALSE-------------------------------------------
#  SubspeciesNMR <- loadObject("subsp.RData")

## -----------------------------------------------------------------------------
data(SrE.IR) # makes the data available
sumSpectra(SrE.IR)

## ----sample-plot, fig.cap = "Sample plot."------------------------------------
# We'll make a fancy title here
# and re-use in other plots
myt <- expression(
  bolditalic(Serenoa)~
  bolditalic(repens)~
  bold(Extract~IR~Spectra))
plotSpectra(SrE.IR,
  main = myt,
  which = c(1, 2, 14, 16),
  yrange = c(0, 1.6),
  offset = 0.4,
  lab.pos = 2200)

## ----subplot, fig.cap = "Detail of the carbonyl region."----------------------
plotSpectra(SrE.IR,
  main = myt,
  which = c(1, 2, 14, 16),
  yrange = c(0, 0.6),
  offset = 0.1,
  lab.pos = 1775,
  xlim = c(1650, 1800))

## -----------------------------------------------------------------------------
# if there are only a few spectra
# show all of the names
SrE.IR$names
# if there are a lot of spectra,
# grep for the desired names
grep("OO", SrE.IR$names)

## ----baseline, fig.cap = "Correcting baseline drift.", sq.fig = TRUE----------
SrE2.IR <- baselineSpectra(SrE.IR,
  int = FALSE,
  method = "modpolyfit",
  retC = TRUE)

## -----------------------------------------------------------------------------
tmp <- binSpectra(SrE.IR, bin.ratio = 4)
sumSpectra(tmp)

## -----------------------------------------------------------------------------
noTD <- removeSample(SrE2.IR,
  rem.sam = c("TD_adSrE"))
sumSpectra(noTD)
grep("TD_adSrE", noTD$names)

## -----------------------------------------------------------------------------
SrE <- grep("SrE", SrE2.IR$names)
# show the name(s) that contain "SrE"
SrE2.IR$names[SrE]
SrE # gives the corresponding indices

## ----survey-1, fig.cap = "Checking for regions of no interest."---------------
surveySpectra(SrE2.IR,
  method = "iqr",
  main = myt,
  by.gr = FALSE)

## ----survey-2, fig.cap = "Checking for regions of no interest."---------------
surveySpectra2(SrE2.IR,
  method = "iqr",
  main = myt)

## ----survey-3, fig.cap = "Detail of carbonyl region."-------------------------
surveySpectra(SrE2.IR,
  method = "iqr",
  main = "Detail of Carbonyl Region",
  by.gr = FALSE,
  xlim = c(1650, 1800))

## ----survey-4, fig.cap = "Detail of carbonyl region by group."----------------
surveySpectra(SrE2.IR,
  method = "iqr",
  main = "Detail of Carbonyl Region",
  by.gr = TRUE,
  xlim = c(1650, 1800))

## ----survey-5, fig.cap = "Inspection of an uninteresting spectral region."----
surveySpectra(SrE2.IR,
  method = "iqr",
  main = "Detail of Empty Region",
  by.gr = FALSE,
  xlim = c(1800, 2500),
  ylim = c(0.0, 0.05))

## -----------------------------------------------------------------------------
SrE3.IR <- removeFreq(SrE2.IR,
  rem.freq = SrE2.IR$freq > 1800 &
  SrE2.IR$freq < 2500)
sumSpectra(SrE3.IR)

## ----gaps, fig.cap = "Identifying gaps in a data set."------------------------
check4Gaps(SrE3.IR$freq, SrE3.IR$data[1,])

## ----hca-1, fig.cap = "Hierarchical cluster analysis.", sq.fig = TRUE---------
HCA <- hcaSpectra(SrE3.IR, main = myt)

## ----classPCA, fig.cap = "Classical PCA scores.", sq.fig = TRUE---------------
c_res <- c_pcaSpectra(SrE3.IR,
  choice = "noscale")
plotScores(SrE3.IR, c_res,
  main = myt,
  pcs = c(1,2),
  ellipse = "rob",
  tol = 0.01)

## ----robPCA, fig.cap = "Robust PCA scores.", sq.fig = TRUE--------------------
r_res <- r_pcaSpectra(SrE3.IR,
  choice = "noscale")
plotScores(SrE3.IR, r_res,
  main = myt,
  pcs = c(1,2),
  ellipse = "rob",
  tol = 0.01)

## ----OD, fig.cap = "Diagnostics: orthogonal distances.", sq.fig = TRUE--------
diagnostics <- pcaDiag(SrE3.IR, c_res,
  pcs = 2,
  plot = "OD")

## ----SD, fig.cap = "Diagnostics: score distances.", sq.fig = TRUE-------------
diagnostics <- pcaDiag(SrE3.IR, c_res,
  pcs = 2,
  plot = "SD")

## ----scree-1, fig.cap = "Scree plot."-----------------------------------------
plotScree(c_res, main = myt)

## ----scree-2, fig.cap = "Traditional style scree plot."-----------------------
plotScree(c_res, style = "trad", main = myt)

## ----boot, fig.cap = "Bootstrap analysis for no. of principal components.", sq.fig = TRUE----
out <- cv_pcaSpectra(SrE3.IR,
  pcs = 5)

## ----results = "hide", eval = FALSE-------------------------------------------
#  plotScoresRGL(SrE3.IR, c_res,
#    main = "S. repens IR Spectra",
#    leg.pos = "A",
#    t.pos = "B") # not run - it's interactive!

## ----scores3D, fig.cap = "Plotting scores in 3D using plotScores3D.", sq.fig = TRUE----
plotScores3D(SrE3.IR, c_res,
  main = myt,
  ellipse = FALSE)

## ----load1, fig.cap = "Loading plot.", sq.fig = TRUE--------------------------
plotLoadings(SrE3.IR, c_res,
  main = myt,
  loads = c(1, 2),
  ref = 1)

## ----load2, fig.cap = "Plotting one loading vs. another.", sq.fig = TRUE------
res <- plot2Loadings(SrE3.IR, c_res,
  main = myt,
  loads = c(1, 2),
  tol = 0.002)

## ----splot-1,  fig.cap = "s-Plot to identify influential frequencies.", sq.fig = TRUE----
spt <- sPlotSpectra(SrE3.IR, c_res,
  main = myt,
  pc = 1,
  tol = 0.001)

## ----splot-2,  fig.cap = "s-Plot detail.", sq.fig = TRUE----------------------
spt <- sPlotSpectra(SrE3.IR, c_res,
  main = "Detail of s-Plot",
  pc = 1,
  tol = 0.05,
  xlim = c(-0.04, -0.01),
  ylim = c(-1.05, -0.9))

## ----results = "hide", eval = FALSE-------------------------------------------
#  hcaScores(SrE3.IR,  c_res,
#    scores = c(1:5),
#    main = myt)

## ----aovPCA2, out.width = "80%", echo = FALSE, fig.cap = "aovPCA breaks the data into a series of submatrices."----
knitr::include_graphics("aovPCA2.png")

## ----aovPCA1, out.width = "80%", echo = FALSE, fig.cap = "Submatrices are composed of rows which are averages of each factor level."----
knitr::include_graphics("aovPCA1.png")

## ----mclust-1, fig.cap = "mclust chooses an optimal model.", results = "hide", sq.fig = TRUE----
model <- mclustSpectra(SrE3.IR, c_res,
  plot = "BIC",
  main = myt)

## ----mclust-2, fig.cap = "mclust's thoughts on the matter.", results = "hide", sq.fig = TRUE----
model <- mclustSpectra(SrE3.IR, c_res,
  plot = "proj",
  main = myt)

## ----mclust-3, fig.cap = "Comparing mclust results to the TRUTH.", results = "hide", sq.fig = TRUE----
model <- mclustSpectra(SrE3.IR, c_res,
  plot = "errors",
  main = myt,
  truth = SrE3.IR$groups)

## ----results = "hide", eval = FALSE-------------------------------------------
#  # not run - it's interactive!
#  mclust3dSpectra(SrE3.IR, c_res)

## ----colsym, echo = FALSE, fig.cap = "Color and symbol suggestions.", fig.width = 6, fig.height = 6----
data(Col7)
data(Col12)
data(Sym12)
data(Col8)
data(Sym8)
auto <- RColorBrewer::brewer.pal(8, "Set1")

sp <- 0.75 # space between major plot elements
tsp <- 0.15 # additional space between points and color swatches/descriptive text
h <- 0.25 # height of the swatch
y <- 0.0 # bottom of the plot, the reference point

# empty plot
plot(1:12, rep(0.0, 12),
  type = "n", yaxt = "n", xaxt = "n", bty = "n",
   xlab = "", ylab = "", ylim = c(0, 3.5)
 )
 text(6.5, y + h + tsp * 4 + sp * 3.5,
   labels = "Automatic Color & Symbol Options", cex = 1.25, font = 2
 )

# Col12
 for (i in 1:12) {
   rect(i - 0.5, y, i + 0.5, y + h, border = NA, col = Col12[i])
 }
 points(1:12, rep(y + h + tsp, 12), pch = Sym12)
 text(0.6, y + h + tsp * 2, adj = 0,
   labels = "gr.cols = 'Col12'     12 mostly paired distinct colors/symbols"
 )

# Col8
 for (i in 1:8) {
   rect(i - 0.5, y + sp, i + 0.5, y + sp + h, border = NA, col = Col8[i])
 }
 points(1:8, rep(y + h + tsp + sp, 8), pch = Sym8)
 text(0.6, y + h + tsp * 2 + sp, adj = 0,
   labels = "gr.cols = 'Col8'     8 distinct colors/symbols"
 )

# auto (original)
 for (i in 1:8) {
   rect(i - 0.5, y + sp * 2, i + 0.5, y + sp * 2 + h, border = NA, col = auto[i])
 }
 points(1:8, rep(y + h + tsp + sp * 2, 8), pch = Sym8)
 text(0.6, y + h + tsp * 2 + sp * 2, adj = 0,
   labels = "gr.cols = 'auto'     8 distinct colors/symbols"
 )

# colorblind-friendly
 for (i in 1:7) {
   rect(i - 0.5, y + sp * 3, i + 0.5, y + sp * 3 + h, border = NA, col = Col7[i])
 }
 points(1:7, rep(y + h + tsp + sp * 3, 7), pch = Sym8[1:7])
 text(0.6, y + h + tsp * 2 + sp * 3, adj = 0,
   labels = "gr.cols = 'Col7'     7 colorblind-friendly colors"
 )


