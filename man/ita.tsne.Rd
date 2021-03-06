% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ita.R
\name{ita.tsne}
\alias{ita.tsne}
\title{Run t-distributed Stochastic Neighbor Embedding}
\usage{
ita.tsne(dat, npc, p = 2)
}
\arguments{
\item{dat}{matrix with}

\item{npc}{number of new dimensions}

\item{p}{perplexity parameter (not bigger than 3 * perplexity < nrow(X) - 1)}
}
\value{
Output of ITA with 
\itemize{
\item \code{ts} scores of the time series on the new dimensions 
\item \code{co} NULL - no scores of the variables on the new dimensions 
\item \code{eig} NULL 
\item \code{dat} original dataset
\item \code{npc} number of new dimensions
\item \code{mita} method used for integrated trend analysis 
}
}
\description{
Run t-distributed Stochastic Neighbor Embedding
}
\note{
Use function Rtsne from the Rtsne package.
}
\references{
Maaten, L. Van Der, 2014. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research, 15, p.3221-3245.
van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
}
