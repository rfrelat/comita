library(shiny)
library(comita)
require(ForeCA)
require(ade4)
require(MARSS)
require(RColorBrewer)
require(vegan)
require(freqdom)
require(fastICA)

met_choices <- list(
  "PCA" = "pca",
  "PCA+Varimax" = "pcav",
  "DPCA (long)" = "dpca",
  "Robust PCA" = "rpca",
  "Kernel PCA" = "kpca",
  "MAFA" = "mafa",
  "t-SNE" = "tsne",
  "isomap" = "isomap",
  "UMAP (long)" = "umap",
  "DFA (very long)" = "dfa",
  "ForeCA (very long)" = "fca"
)
