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

col_met <- rainbow(length(met_choices))

pal <- list(
 "1" = c("springgreen4", "chartreuse3", "yellow", "darkgoldenrod1", "red"),
 "2" = c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571"),
 "3" = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6")
)
