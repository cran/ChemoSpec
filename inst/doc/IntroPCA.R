## ----SetUp, echo = FALSE, eval = TRUE, results = "hide"----
# R options & configuration:
set.seed(9)
rm(list = ls())
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("kableExtra"))
suppressPackageStartupMessages(library("chemometrics"))
options(width =  50) # for pinp specifically (narrow cols)

# Stuff specifically for knitr:
opts_chunk$set(eval = TRUE, echo = FALSE, results = "hide",
  fig.width = 3.5, fig.height = 3.5)

## ---- dataTaste, results = "asis"---------------
data(glass)
DF <- as.data.frame(glass[1:5, 1:8])
kable(DF, format = "latex", caption = "A portion of the archaeological glass data set. Values are percentages.\\label{tab:dataTaste}") %>% kable_styling(c("striped", "bordered"))

## ----glassScree, fig.cap = "Scree plot from PCA on the glass data set.\\label{fig:glassScree}"----
pca <- prcomp(glass)
plot(pca, main = "")

## ----glassScores, fig.cap = "Score plot from PCA on the glass data set.\\label{fig:glassScores}"----
plot(pca$x[,1], pca$x[,2], type = "p",
  xlab = "Principal Component 1", ylab = "Principal Component 2")

## ----glassLoadings, fig.cap = "Loadings plot for the first PC from PCA on the glass data set.\\label{fig:glassLoadings}"----
barplot(pca$rotation[,1], cex.names = 0.7, ylab = "Contribution")

## ----screeTable, results = "asis"---------------
eigensum <- sum(pca$sdev*pca$sdev)
variance <- 100*(pca$sdev*pca$sdev/eigensum)
cumvariance <- cumsum(variance)
labs <- paste("PC", 1:13, sep = " ")
DF <- data.frame(component = labs, variance = variance, cumulative = cumvariance)
kable(DF, format = "latex", digits = 0, caption = "Variance (signal) accounted for by PCs. Values in percent.\\label{tab:screeTable}")

## ----groups, fig.cap = "Score plot from PCA on the glass data set, with groups color-coded.\\label{fig:glassScores2}"----
gr <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L, 4L, 4L, 4L, 4L,
4L, 4L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Group.1", "Group.2",
"Group.3", "Group.4"), class = "factor")
plot(pca$x[,1], pca$x[,2], type = "p", col = gr, pch = 20,
  xlab = "Principal Component 1", ylab = "Principal Component 2")

## ----loadIR-------------------------------------
suppressPackageStartupMessages(library("ChemoSpec"))
data(SrE.IR)

## ----IRSpectrum, fig.cap = "Spectrum 1 from the IR data set.\\label{fig:IRSpecrum}"----
xl <- rev(range(SrE.IR$freq))
plot(SrE.IR$freq, SrE.IR$data[1,], type = "l", xlim = xl)

## ----IRScree, fig.cap = "Scree plot from PCA on the IR data set.\\label{fig:IRScree}"----
pca <- prcomp(SrE.IR$data)
plot(pca, main = "")

## ----IRScores, fig.cap = "Score plot from PCA on the IR data set.\\label{fig:IRScores}"----
plot(pca$x[,1], pca$x[,2], type = "p",
  xlab = "Principal Component 1", ylab = "Principal Component 2")

## ----IRLoadings, fig.cap = "Loadings plot for the first PC from PCA on the IR data set.\\label{fig:IRLoadings}"----
plot(SrE.IR$freq, pca$rotation[,1], type = "l", xlim = xl,
     xlab = "Wavelength", ylab = "Contribution")

## ----IRLoadings2, fig.cap = "Loadings plot for the first PC from PCA on the IR data set, carbonyl region. Reference spectrum shown in red.\\label{fig:IRLoadings2}"----
plot(SrE.IR$freq, pca$rotation[,1], type = "l", xlim = c(1800, 1650),
     xlab = "Wavelength", ylab = "Contribution", ylim = c(-0.3, 0.3))
lines(SrE.IR$freq, SrE.IR$data[1,], col = "red")
abline(v = c(1743, 1708), lty = 2, col = "gray50")

## ----IRLoadings3, fig.cap = "Loadings plot for the first PC from PCA on the IR data set, carbonyl region, shown as a bar plot.\\label{fig:IRLoadings3}"----
plot(SrE.IR$freq, pca$rotation[,1], type = "h", xlim = c(1800, 1650),
     xlab = "Wavelength", ylab = "Contribution")

