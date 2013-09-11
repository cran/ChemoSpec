
## ----SetUp, echo = FALSE, eval = TRUE, results = "hide"------------------

# You must knit this file with getwd() set to the directory it is in!

# R options & configuration:

rm(list = ls())
options(width =  50, show.signif.stars = FALSE)

suppressMessages(library("knitr"))
suppressMessages(library("ChemoSpec"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("mvbutils"))
suppressMessages(library("sna"))

desc <- packageDescription("ChemoSpec")
vers <- paste("(Package Version ", desc$Version, ")", sep = "")

# Stuff specifically for knitr:

opts_chunk$set(out.width = "0.8\\textwidth", fig.align = "center", fig.width = 7, fig.height = 7, cache = FALSE)

# Note: defaults are eval = TRUE, echo = TRUE


## ----Chunk1,  results = "hide", eval = FALSE-----------------------------
## source("My_First_ChemoSpec.R")


## ----Chunk2,  results = "hide", eval = FALSE-----------------------------
## getManyCsv(gr.crit = c("sspA", "sspB"), gr.cols = c("red3", "dodgerblue4"),
## freq.unit = "ppm", int.unit = "peak intensity", descrip = "Subspecies Study",
## out.file = "subspecies")


## ----Chunk3,  results = "hide", eval = FALSE-----------------------------
## SubspeciesNMR <- loadObject("subspecies.RData")


## ----Chunk4,  echo = FALSE, fig.width = 6, fig.height = 4, fig.cap = "\\label{colors}Recommended Color Sets in ChemoSpec"----
par(mfrow = c(2,1), mar = c(3, 1, 2, 1))
display.brewer.pal(8, "Set1")
title(main = "ChemoSpec Primary Scheme")
display.brewer.pal(8, "Set2")
title(main = "ChemoSpec Pastel Scheme")
par(mfrow = c(1,1))


## ----Chunk5, echo = TRUE-------------------------------------------------
data(SrE.IR) # makes the data available
sumSpectra(SrE.IR)


## ----Chunk8, fig.cap = "\\label{plot}Plotting Spectra"-------------------
# We'll make a fancy proper title just this once!
myt <- expression(paste(italic(S.), phantom(0), italic(repens), " Extract IR Spectra"))
plotSpectra(SrE.IR, title = myt,
which = c(1, 2, 14, 16),
yrange = c(0, 1.6), offset = 0.4, lab.pos = 2200)


## ----Chunk10, fig.cap = "\\label{subplot}Zooming in on a Spectral Region"----
plotSpectra(SrE.IR, title = "S. repens IR Spectra: Detail of Carbonyl Region",
which = c(1, 2, 14, 16), xlim = c(1650, 1800),
yrange = c(0, 0.6), offset = 0.1, lab.pos = 1775)


## ----Chunk9--------------------------------------------------------------
SrE.IR$names # suitable if there are not many spectra
grep("OO", SrE.IR$names) # use if there are more spectra


## ----Chunk11-------------------------------------------------------------
noTD <- removeSample(SrE.IR, rem.sam = c("TD_adSrE"))
sumSpectra(noTD)
grep("TD_adSrE", noTD$names)


## ----Chunk12-------------------------------------------------------------
SrE <- grep("SrE", SrE.IR$names)
SrE.IR$names[SrE] # gives the name(s) that contain "SrE"
SrE # gives the corresponding indicies


## ----Chunk10b, fig.cap = "\\label{baseline}Correcting baseline drift"----
SrE2.IR <- baselineSpec(SrE.IR, int = FALSE, method = "rfbaseline", retC = TRUE)


## ----Chunk14,  fig.width = 8, fig.height = 5, fig.cap = "\\label{surv}Checking for Regions of No Interest"----
specSurvey(SrE2.IR, method = "iqr", title = "S. repens Extract IR Spectra", by.gr = FALSE)


## ----Chunk14a, fig.cap = "\\label{survA}Detail of Carbonyl Region"-------
specSurvey(SrE2.IR, method = "iqr", title = "S. repens Detail of Carbonyl Region", by.gr = FALSE, xlim = c(1650, 1800))


## ----Chunk14b, fig.cap = "\\label{survB}Detail of Carbonyl Region by Group"----
specSurvey(SrE2.IR, method = "iqr", title = "S. repens Detail of Carbonyl Region", by.gr = TRUE, xlim = c(1650, 1800))


## ----Chunk14c, fig.cap = "\\label{survC}Inspection of an Uninteresting Spectral Region"----
specSurvey(SrE2.IR, method = "iqr", title = "S. repens Detail of Empty Region", by.gr = FALSE, xlim = c(1800, 2500))


## ----Chunk15-------------------------------------------------------------
SrE3.IR <- removeFreq(SrE2.IR, rem.freq = SrE2.IR$freq > 1800 & SrE2.IR$freq < 2500)
sumSpectra(SrE3.IR)


## ----Chunk7, fig.cap = "\\label{gaps}Procedure to Find Gaps in a Data Set" , tidy = FALSE----
check4Gaps(SrE3.IR$freq, SrE3.IR$data[1,], plot = TRUE)


## ----Chunk17, eval = FALSE, echo = TRUE----------------------------------
## SrE3.IR <- normSpectra(SrE3.IR)


## ----Chunk18-------------------------------------------------------------
tmp <- binBuck(SrE3.IR, bin.ratio = 4)
sumSpectra(tmp)


## ----Chunk19, fig.cap = "\\label{hca}Hierarchical Cluster Analysis"------
hcaSpectra(SrE3.IR, title = "S. repens IR Spectra")


## ----Chunk10a, fig.cap = "\\label{classPCA}Classical PCA"----------------
class <- classPCA(SrE3.IR, choice = "noscale")
plotScores(SrE3.IR, title = "S. repens IR Spectra", class,
pcs = c(1,2), ellipse = "rob", tol = 0.01)


## ----Chunk21, fig.cap = "\\label{robPCA}Robust PCA"----------------------
robust <- robPCA(SrE3.IR, choice = "noscale")
plotScores(SrE3.IR, title = "S. repens IR Spectra", robust,
pcs = c(1,2), ellipse = "rob", tol = 0.01)


## ----Chunk22, fig.cap = "\\label{OD}Diagnostics: Orthogonal Distances"----
diagnostics <- pcaDiag(SrE3.IR, class, pcs = 2, plot = "OD")


## ----Chunk23, fig.cap = "\\label{SD}Diagnostics: Score Distances"--------
diagnostics <- pcaDiag(SrE3.IR, class, pcs = 2, plot = "SD")


## ----Chunk24, fig.cap = "\\label{scree}Scree Plot"-----------------------
plotScree(class, title = "S. repens IR Spectra")


## ----Chunk24a, fig.cap = "\\label{scree2}Alternate Style Scree Plot"-----
plotScree2(class, title = "S. repens IR Spectra")


## ----Chunk25, fig.cap = "\\label{boot}Bootstrap Analysis for No. of PCs"----
out <- pcaBoot(SrE3.IR, pcs = 5, choice = "noscale")


## ----Chunk26,  results = "hide", eval = FALSE----------------------------
## plotScoresRGL(SrE3.IR, class, title = "S. repens IR Spectra",
## leg.pos = "A", t.pos = "B") # not run - it's interactive!


## ----Chunk27, fig.cap = "\\label{s3D}Plotting Scores in 3D using plotScores3D"----
plotScores3D(SrE3.IR, class, title = "S. repens IR Spectra", ellipse = FALSE)


## ----Chunk28,  results = "hide", eval = FALSE----------------------------
## plotScoresG(SrE3.IR, class) # not run - it's interactive!


## ----Chunk29, fig.cap = "\\label{load}Loading Plot"----------------------
plotLoadings(SrE3.IR, class, title = "S. repens IR Spectra",
loads = c(1, 2), ref = 1)


## ----Chunk30, fig.cap = "\\label{load2}Plotting One Loading vs. Another"----
plot2Loadings(SrE3.IR, class, title = "S. repens IR Spectra",
loads = c(1, 2), tol = 0.002)


## ----Chunk30a,  fig.cap = "\\label{splot}s-Plot to Identify Influential Frequencies"----
spt <- sPlotSpectra(SrE3.IR, class, title = "S. repens IR Spectra", pc = 1, tol = 0.001)


## ----Chunk30b,  fig.cap = "\\label{splot2}s-Plot Detail"-----------------
spt <- sPlotSpectra(SrE3.IR, class, title = "Detail of S. repens IR Spectra", pc = 1, tol = 0.05, xlim = c(-0.04, -0.01), ylim = c(-1.05, -0.9))


## ----Chunk31,  results = "hide", eval = FALSE----------------------------
## hcaScores(SrE3.IR,  class, scores = c(1:5), title = "S. repens IR Spectra")


## ----Chunk35, fig.cap = "\\label{mclust1}mclust Chooses an Optimal Model"----
model <- mclustSpectra(SrE3.IR, class, plot = "BIC",
title = "S. repens IR Spectra")


## ----Chunk36, fig.cap = "\\label{mclust2}mclust's Thoughts on the Matter"----
model <- mclustSpectra(SrE3.IR, class, plot = "proj",
title = "S. repens IR Spectra")


## ----Chunk37, fig.cap = "\\label{mclust3}Comparing mclust Results to the TRUTH"----
model <- mclustSpectra(SrE3.IR, class, plot = "errors",
title = "S. repens IR Spectra", truth = SrE3.IR$groups)


## ----Chunk33,  results = "hide", eval = FALSE----------------------------
## mclust3dSpectra(SrE3.IR, class) # not run - it's interactive!


## ----Chunk38,  echo = FALSE, fig.cap = "\\label{food}Map of Functions in ChemoSpec", fig.width = 10, fig.height = 10----
CS <- foodweb(where = "package:ChemoSpec", plotting = FALSE)
CS$funmat <- CS$funmat[-9, -9] # removes chkSpectra
set.seed(7) # 7, 9, 19 decent
tst <- gplot(CS$funmat, displaylabels = TRUE, label.cex = 0.6,
	arrowhead.cex = 0.5, label.pos = 6, pad = 1.5)


