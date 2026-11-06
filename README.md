
# comita

<!-- badges: start -->
[![R-CMD-check](https://github.com/rfrelat/comita/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rfrelat/comita/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Integrated trend analyses (ITA) are commonly applied to summarise changes that have
occurred in recent decades in ecosystems. Currently there are multiple methods used
for ITA, raising question on which is the best method to use for ITA. The aim of this
application is to provide a tool to easily compare different ITA methods, and test
the effect of data pre-processing and variables selection on the final results.
Hopefully, this package will raise awareness of ecologists to different dimension
reduction techniques, as well as provide reason and confidence on the choice of ITA
methods according to the characteristics of their dataset.

## Installation

You can install the development version of comita from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("rfrelat/comita")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(comita)
## basic example code
data(baltic)

## Run a single multivariate analysis
mvar <- ita(baltic, npc=2, met = "pca")

# Visually interpret the results of the multivariate analysis
plot_eig(mvar) # show the eigen values (if any)
plot_var(mvar) # show the scores of the variables
plot_timeseries(mvar) # show the scores of the time series
plot_heatmap(mvar) # 2D representation of the dataset ordered by the first dimension

## Compare with other method
multi <- multita(baltic, npc=2, met=c("pca", "kpca"))
plot_compts(baltic, multi) # Compare scores of time series
plot_compvar(baltic, multi) # compare the score of the variable

## Compare to randomly shuffled dataset
mrand <- ita.nrand(baltic, metrand = "phase", nrep=10)
# Graphic comparing eigen values
plot_nrand(mvar, mrand)
```


## Related package

In parallel, the [`dimRed` package](https://github.com/gdkrmr/dimRed) offers another commun framework for multiple dimension reductions techniques, check out their [documentation](https://journal.r-project.org/articles/RJ-2018-039/index.html).

