## Comparative tools for Integrated Trend Analysis

Integrated trend analyses (ITA) are commonly applied to summarise changes that have occurred in recent decades in ecosystems. Currently there are multiple methods used for ITA, raising question on which is the best method to use for ITA.
The aim of this application is to provide a tool to easily compare different ITA methods, and test the effect of data pre-processing and variables selection on the final results. 
Hopefully, this tool will raise awareness of ecologists to different dimension reduction techniques, as well as provide reason and confidence on the choice of ITA methods according to the characteristics of their dataset.
The application comes with a R-package \code{comita} that can reproduce and make the analysis in R-script.

## Dimension reduction techniques

It is impossible to list all the dimension reduction techniques that exist, but we implemented the most common ones that (1) are used by ecologists, and (2) have been implemented in \code{R}.

|  Abb   |                 Method                |                References               |
|:------:|:-------------------------------------:|:--------------------------------------:|
| PCA    | Principal Component Analysis          | [Pearson, 1901](https://doi.org/10.1080/14786440109462720) |
| PCAV   | PCA + Varimax                         | [Kaiser, 1958](https://doi.org/10.1007/BF02289233) |
| DPCA   | Dynamic principal components analysis | [Ku et al., 1995](https://www.sciencedirect.com/science/article/pii/0169743995000763), [De Ketelaere et al., 2015](https://doi.org/10.1080/00224065.2015.11918137) |
| RPCA   | Robust PCA (ROBPCA)                   | [Hubert et al., 2005](https://doi.org/10.1198/004017004000000563) |
| KPCA   | Kernel PCA                            | [Schölkopf et al., 1998](https://doi.org/10.1007/978-1-4471-1599-1_18) |
| MAFA   | Min/Max Autocorrelation Factor Analysis | [Solow, 1994](https://doi.org/10.2307/2533401) |
| t-SNE  | t-distributed Stochastic Neighbor Embedding | [van der Maaten and Hinton, 2008](https://www.jmlr.org/papers/v9/vandermaaten08a.html) |
| Isomap | Isometric Feature Mapping             | [Tenenbaum et al., 2000](https://doi.org/10.1126/science.290.5500.2319) |
| UMAP   | Uniform Manifold Approximation and Projection | [McInnes et al., 2018](https://arxiv.org/abs/1802.03426) |
| DFA    | Dynamic Factor Analysis               | [Zuur et al., 2003](https://doi.org/10.1002/env.611) |
| ForeCA | Forecastable Component Analysis       | [Goerg, 2013](https://arxiv.org/abs/1205.4591) |

